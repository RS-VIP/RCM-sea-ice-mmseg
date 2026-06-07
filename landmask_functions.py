import numpy as np
import rasterio

from rasterio.crs import CRS
from rasterio.transform import from_origin
from rasterio.control import GroundControlPoint
from rasterio.warp import reproject, Resampling, transform_bounds, transform as crs_transform
from rasterio.features import rasterize
from rasterio.env import Env

import geopandas as gpd
from shapely.geometry import box
from shapely.ops import unary_union


def build_land_masks_earth_geometry(shp_path: str, rcm_product: list[dict]) -> dict:
    """
      - The shapefile polygons represent OCEAN (ocean=1), so land is the inverse.
      - add boolean masks: True=land, False=other - to the dict  
    """

    # Read shapefile
    gdf_raw = gpd.read_file(shp_path)

    # Read SAR grid info from the rcm_product dict
    hh = rcm_product["hh"]
    transform = rcm_product["src_transform"]
    crs = rcm_product["src_crs"]
    bounds = rcm_product["src_bounds"]
    folder_name = rcm_product["folder_name"]
    shape = hh.shape

    # Reproject shapefile to this rcm_product CRS
    gdf = gdf_raw.to_crs(crs)

    # SAR bbox polygon
    sar_bbox = box(bounds.left, bounds.bottom, bounds.right, bounds.top)

    # Keep only shapefile features that intersect the SAR bbox + Clip the shapefile geometry to the SAR bbox 
    gdf = gdf[gdf.intersects(sar_bbox)].copy()
    gdf["geometry"] = gdf.intersection(sar_bbox)
    if len(gdf) == 0:
        print(f"warning: {folder_name}: shapefile does not intersect bbox")
        return None

    # Merge polygons
    geom = unary_union([g for g in gdf.geometry if g is not None and not g.is_empty])

    # Rasterize ocean polygons to mask
    ocean_mask = rasterize(
        [(geom, 1)],
        out_shape=shape,
        transform=transform,
        fill=0,
        dtype=np.uint8,
        all_touched=True,
    ).astype(bool)

    # Convert ocean-mask to land-mask (land=True, ocean=False)
    land_mask = ~ocean_mask

    return land_mask

def build_land_mask_sensor_geometry(
    rcm_product: dict,
    shp_path: str,
    mask_res_m: float = 200,
    simplify_m: float = 0,
    threshold: float = 0.5,
    chunk_rows: int = 512,
):
    """
    Build a high-accuracy land mask for an RCM image in sensor geometry.

    Pipeline:
        1) Geocode SAR image to EPSG:3413 at chosen resolution.
        2) Rasterize ocean shapefile onto that projected grid.
        3) Resample ocean mask back to sensor grid using bilinear interpolation.
        4) Invert to obtain land mask.

    Parameters
    ----------
    rcm_product : dict
        Output from load_rcm_product (must contain hh + tie arrays).
    shp_path : str
        Path to ocean shapefile (polygons represent ocean).
    mask_res_m : float
        Intermediate projected raster resolution in meters.
        Smaller = sharper coastline but larger intermediate raster.
    simplify_m : float
        Optional coastline simplification tolerance (meters).
        Use 0 for maximum accuracy.
    threshold : float
        Threshold applied after bilinear resampling (default 0.5).
    chunk_rows : int
        Number of sensor rows processed per chunk to control memory usage.

    Returns
    -------
    land_mask : np.ndarray (bool)
        Boolean mask in sensor geometry (True = land, False = ocean).
    """

    def geocode_hh_to_3413_in_memory(rcm_product, pad_m=2000, target_res_m=200):
        """
        Robust: build dst grid from tie-point lon/lat bounds (projected to EPSG:3413),
        then reproject using GCPs (TPS-like mapping handled by GDAL internally).
        """
        hh = rcm_product["hh"].astype(np.float32, copy=False)
        H, W = hh.shape

        tie_r  = np.asarray(rcm_product["tie_lines"],  dtype=np.float64)
        tie_c  = np.asarray(rcm_product["tie_pixels"], dtype=np.float64)
        tie_la = np.asarray(rcm_product["tie_lats"],   dtype=np.float64)
        tie_lo = np.asarray(rcm_product["tie_lons"],   dtype=np.float64)

        # GCPs: X=lon, Y=lat
        gcps = [
            GroundControlPoint(row=float(r), col=float(c), x=float(lon), y=float(lat))
            for r, c, lat, lon in zip(tie_r, tie_c, tie_la, tie_lo)
        ]

        src_crs = CRS.from_epsg(4326)
        dst_crs = CRS.from_epsg(3413)

        # Use tie-point lon/lat bounds (more robust than asking GDAL to infer from GCPs)
        lon_min, lon_max = float(tie_lo.min()), float(tie_lo.max())
        lat_min, lat_max = float(tie_la.min()), float(tie_la.max())

        # Project bounds to EPSG:3413
        with Env(OGR_CT_FORCE_TRADITIONAL_GIS_ORDER="YES"):
            left, bottom, right, top = transform_bounds(
                src_crs, dst_crs, lon_min, lat_min, lon_max, lat_max, densify_pts=21
            )

        # Pad a bit to be safe
        left -= pad_m; right += pad_m
        bottom -= pad_m; top += pad_m

        # Target resolution (meters)
        xres = target_res_m
        yres = target_res_m

        dst_width  = int(np.ceil((right - left) / xres))
        dst_height = int(np.ceil((top - bottom) / yres))

        # North-up affine (y pixel size is negative)
        dst_transform = from_origin(left, top, xres, yres)

        dst = np.zeros((dst_height, dst_width), dtype=np.float32)

        with Env(OGR_CT_FORCE_TRADITIONAL_GIS_ORDER="YES"):
            reproject(
                source=hh,
                destination=dst,
                src_crs=src_crs,
                gcps=gcps,
                dst_crs=dst_crs,
                dst_transform=dst_transform,
                resampling=Resampling.bilinear,
            )

        return dst, dst_transform, dst_crs

    def ocean_mask_on_geocoded_grid_fast(ocean_shp_path, out_shape, dst_transform, dst_crs, simplify_m=0):
        # Read + reproject
        gdf = gpd.read_file(ocean_shp_path).to_crs(dst_crs)

        # Raster bounds in dst CRS
        H, W = out_shape
        left, top = (dst_transform.c, dst_transform.f)
        right = left + dst_transform.a * W
        bottom = top + dst_transform.e * H  # dst_transform.e is negative for north-up
        rbbox = box(min(left, right), min(bottom, top), max(left, right), max(bottom, top))

        # Use spatial index to prefilter
        sidx = gdf.sindex
        hits = list(sidx.intersection(rbbox.bounds))
        gdf = gdf.iloc[hits].copy()

        # Precise clip to bounds
        gdf = gdf[gdf.intersects(rbbox)]
        if gdf.empty:
            return np.zeros(out_shape, dtype=bool)

        gdf["geometry"] = gdf.intersection(rbbox)

        # Optional: simplify coastline (meters, since CRS is EPSG:3413)
        if simplify_m and simplify_m > 0:
            gdf["geometry"] = gdf["geometry"].simplify(simplify_m, preserve_topology=True)

        # Dissolve/merge to reduce number of shapes
        geom = unary_union([g for g in gdf.geometry if g is not None and not g.is_empty])
        if geom.is_empty:
            return np.zeros(out_shape, dtype=bool)

        ocean = rasterize(
            [(geom, 1)],
            out_shape=out_shape,
            transform=dst_transform,
            fill=0,
            dtype=np.uint8,
            all_touched=False,
        ).astype(bool)

        return ocean

    def ocean_mask_back_to_sensor_bilinear(
        ocean_mask_3413: np.ndarray,         # bool or uint8 (H3413, W3413)
        transform_3413,                      # Affine
        rcm_product: dict,
        thr: float = 0.5,
        chunk_rows: int = 256,
    ):
        """
        Build full-res ocean mask on sensor grid by:
        - lon/lat per sensor pixel (from regular tie grid)
        - lon/lat -> EPSG:3413
        - rowcol -> fractional pixel coords in ocean_mask_3413
        - manual bilinear sampling
        - threshold

        Returns ocean_sensor bool (H,W).
        """
        Hs, Ws = rcm_product["hh"].shape

        tie_r  = np.asarray(rcm_product["tie_lines"],  dtype=np.float64)
        tie_c  = np.asarray(rcm_product["tie_pixels"], dtype=np.float64)
        tie_la = np.asarray(rcm_product["tie_lats"],   dtype=np.float64)
        tie_lo = np.asarray(rcm_product["tie_lons"],   dtype=np.float64)

        # regular grid axes
        r_axis = np.unique(tie_r)
        c_axis = np.unique(tie_c)
        Ny, Nx = len(r_axis), len(c_axis)

        ir = np.searchsorted(r_axis, tie_r)
        ic = np.searchsorted(c_axis, tie_c)
        lon_grid = np.empty((Ny, Nx), dtype=np.float64); lon_grid[ir, ic] = tie_lo
        lat_grid = np.empty((Ny, Nx), dtype=np.float64); lat_grid[ir, ic] = tie_la

        # ocean mask as float in [0,1]
        src = ocean_mask_3413.astype(np.float32, copy=False)
        Hg, Wg = src.shape

        def interp_lonlat(rows, cols):
            rr = np.interp(rows, r_axis, np.arange(Ny))
            cc = np.interp(cols, c_axis, np.arange(Nx))

            r0 = np.floor(rr).astype(int); r1 = np.clip(r0 + 1, 0, Ny - 1)
            c0 = np.floor(cc).astype(int); c1 = np.clip(c0 + 1, 0, Nx - 1)
            dr = (rr - r0)[:, None]
            dc = (cc - c0)[None, :]

            lon00 = lon_grid[r0[:, None], c0[None, :]]
            lon01 = lon_grid[r0[:, None], c1[None, :]]
            lon10 = lon_grid[r1[:, None], c0[None, :]]
            lon11 = lon_grid[r1[:, None], c1[None, :]]

            lat00 = lat_grid[r0[:, None], c0[None, :]]
            lat01 = lat_grid[r0[:, None], c1[None, :]]
            lat10 = lat_grid[r1[:, None], c0[None, :]]
            lat11 = lat_grid[r1[:, None], c1[None, :]]

            lon0 = lon00 * (1 - dc) + lon01 * dc
            lon1 = lon10 * (1 - dc) + lon11 * dc
            lon  = lon0 * (1 - dr) + lon1 * dr

            lat0 = lat00 * (1 - dc) + lat01 * dc
            lat1 = lat10 * (1 - dc) + lat11 * dc
            lat  = lat0 * (1 - dr) + lat1 * dr

            return lon, lat

        ocean_sensor = np.zeros((Hs, Ws), dtype=bool)
        cols = np.arange(Ws, dtype=np.float64)

        for r_start in range(0, Hs, chunk_rows):
            r_end = min(Hs, r_start + chunk_rows)
            rows = np.arange(r_start, r_end, dtype=np.float64)

            # (R,C) lon/lat for this chunk
            lon, lat = interp_lonlat(rows, cols)

            # lon/lat -> EPSG:3413 coords
            xs, ys = crs_transform("EPSG:4326", "EPSG:3413",
                                lon.ravel().tolist(), lat.ravel().tolist())
            xs = np.asarray(xs); ys = np.asarray(ys)

            # fractional (row, col) in geocoded raster
            rr, cc = rasterio.transform.rowcol(transform_3413, xs, ys, op=float)
            rr = np.asarray(rr); cc = np.asarray(cc)

            # bilinear sampling
            r0 = np.floor(rr).astype(np.int64)
            c0 = np.floor(cc).astype(np.int64)
            r1 = r0 + 1
            c1 = c0 + 1

            # valid where all neighbors are inside
            valid = (r0 >= 0) & (c0 >= 0) & (r1 < Hg) & (c1 < Wg)

            out = np.zeros(rr.shape, dtype=np.float32)
            if np.any(valid):
                dr = (rr[valid] - r0[valid]).astype(np.float32)
                dc = (cc[valid] - c0[valid]).astype(np.float32)

                v00 = src[r0[valid], c0[valid]]
                v01 = src[r0[valid], c1[valid]]
                v10 = src[r1[valid], c0[valid]]
                v11 = src[r1[valid], c1[valid]]

                out[valid] = (
                    v00 * (1 - dr) * (1 - dc) +
                    v01 * (1 - dr) * dc +
                    v10 * dr * (1 - dc) +
                    v11 * dr * dc
                )

            ocean_sensor[r_start:r_end, :] = (out.reshape(r_end - r_start, Ws) > thr)

        return ocean_sensor


    # ------------------------------------------------------------------
    # STEP 1 — Geocode SAR band to projected CRS (EPSG:3413)
    # ------------------------------------------------------------------
    # This creates an intermediate raster in map geometry
    # at the requested resolution (mask_res_m).
    # The result is used only to rasterize ocean polygons cleanly.
    # ------------------------------------------------------------------

    hh_3413, tr_3413, crs_3413 = geocode_hh_to_3413_in_memory(
        rcm_product,
        target_res_m=mask_res_m
    )

    # ------------------------------------------------------------------
    # STEP 2 — Rasterize ocean shapefile on projected grid
    # ------------------------------------------------------------------
    # - Reproject shapefile to EPSG:3413
    # - Clip to raster bounds
    # - Optionally simplify coastline
    # - Rasterize into binary ocean mask
    #
    # Result: ocean_3413 (True = ocean)
    # ------------------------------------------------------------------

    ocean_3413 = ocean_mask_on_geocoded_grid_fast(
        shp_path,
        out_shape=hh_3413.shape,
        dst_transform=tr_3413,
        dst_crs=crs_3413,
        simplify_m=simplify_m,
    )

    # ------------------------------------------------------------------
    # STEP 3 — Resample ocean mask back to sensor geometry
    # ------------------------------------------------------------------
    # For each sensor pixel:
    #   - Interpolate lon/lat from tie-point grid
    #   - Project to EPSG:3413
    #   - Sample ocean_3413 using bilinear interpolation
    #   - Apply threshold to obtain boolean mask
    #
    # chunk_rows controls memory usage.
    # ------------------------------------------------------------------

    ocean_sensor = ocean_mask_back_to_sensor_bilinear(
        ocean_3413,
        tr_3413,
        rcm_product,
        thr=threshold,
        chunk_rows=chunk_rows,
    )

    # ------------------------------------------------------------------
    # STEP 4 — Convert ocean mask to land mask
    # ------------------------------------------------------------------
    # Shapefile represents ocean polygons (ocean=True),
    # so land is simply the inverse.
    # ------------------------------------------------------------------

    land_mask = ~ocean_sensor

    return land_mask
