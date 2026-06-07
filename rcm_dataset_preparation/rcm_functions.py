import numpy as np

from pathlib import Path
import rasterio
from rasterio.warp import reproject, calculate_default_transform, Resampling
from rasterio.warp import transform as transform_coords
from rasterio.transform import rowcol
from rasterio.crs import CRS
from affine import Affine
from lxml import etree
from datetime import datetime

import matplotlib.pyplot as plt
import copy

# ============================================================
# Loader
# ============================================================

def load_rcm_product(data_dir):
    def _is_identity_transform(t, tol=1e-9) -> bool:
        return (
            abs(t.a - 1.0) < tol and abs(t.b) < tol and abs(t.c) < tol and
            abs(t.d) < tol and abs(t.e + 1.0) < tol and abs(t.f) < tol
        ) or (
            abs(t.a - 1.0) < tol and abs(t.b) < tol and abs(t.c) < tol and
            abs(t.d) < tol and abs(t.e - 1.0) < tol and abs(t.f) < tol
        )

    def _detect_geometry(src_crs, src_transform, gcps_count: int) -> str:
        if src_crs is not None:
            return "earth"
        return "sensor"

    def _parse_tie_points(xml_root, ns):
        lines = []
        pixels = []
        lats = []
        lons = []

        tps = xml_root.findall(".//rcm:geolocationGrid//rcm:imageTiePoint", ns)

        for tp in tps:
            img = tp.find(".//rcm:imageCoordinate", ns)
            geo = tp.find(".//rcm:geodeticCoordinate", ns)
            if img is None or geo is None:
                continue

            line = img.find("rcm:line", ns)
            pixel = img.find("rcm:pixel", ns)
            lat = geo.find("rcm:latitude", ns)
            lon = geo.find("rcm:longitude", ns)

            if None in (line, pixel, lat, lon):
                continue

            lines.append(float(line.text))
            pixels.append(float(pixel.text))
            lats.append(float(lat.text))
            lons.append(float(lon.text))

        return {
            "tie_lines": np.array(lines, dtype=np.float64),
            "tie_pixels": np.array(pixels, dtype=np.float64),
            "tie_lats": np.array(lats, dtype=np.float64),
            "tie_lons": np.array(lons, dtype=np.float64),
        }

    def _infer_axis_mapping_from_xml(pixel_spacing):
        return {
            "cols_are": "range",
            "rows_are": "azimuth",
            "range_spacing_m": pixel_spacing.get("range_m"),
            "azimuth_spacing_m": pixel_spacing.get("azimuth_m"),
            "confidence": "high (SAR convention: samples=pixels=cols, lines=rows)"
        }

    def extract_acquisition_time(xml_root, ns):
        """
        Extract acquisition time from RCM product XML.

        Priority:
            1. rawDataStartTime  (preferred for RCM)
            2. acquisitionStartTime
            3. sceneStartTime

        Returns
        -------
        str or None
            ISO timestamp string (e.g., '2021-03-15T14:23:08.123456Z')
            or None if not found.
        """

        # Preferred tag
        acq_time_elem = xml_root.find(".//rcm:rawDataStartTime", ns)

        # Fallbacks
        if acq_time_elem is None:
            acq_time_elem = xml_root.find(".//rcm:acquisitionStartTime", ns)

        if acq_time_elem is None:
            acq_time_elem = xml_root.find(".//rcm:sceneStartTime", ns)

        if acq_time_elem is None or acq_time_elem.text is None:
            return None

        return datetime.fromisoformat(acq_time_elem.text.strip().replace("Z", "+00:00"))

    img_file = list(Path(data_dir).glob("*.img"))
    if len(img_file) != 1:
        img_file = list(Path(data_dir).glob("*.tif"))
        if len(img_file) != 1:
            raise ValueError("expected one .img or .tif file")

    xml_file = list(Path(data_dir).glob("product.xml"))
    if len(xml_file) != 1:
        raise ValueError("expected one product.xml file")

    img_path = img_file[0]
    xml_path = xml_file[0]

    try:
        with rasterio.open(img_path) as src:
            if src.count != 2:
                raise ValueError(f"expected 2 bands, found {src.count}")

            hh = src.read(1).astype(np.float32, copy=False)
            hv = src.read(2).astype(np.float32, copy=False)

            nodata_hh = src.nodatavals[0]
            nodata_hv = src.nodatavals[1]

            hh[hh==nodata_hh] = np.nan
            hv[hv==nodata_hv] = np.nan

            src_transform = src.transform
            src_crs = src.crs
            src_bounds = src.bounds


            gcps, gcp_crs = src.gcps
            gcps_count = len(gcps) if gcps is not None else 0

        xml_root = etree.parse(str(xml_path)).getroot()
        
        ns = {"rcm": xml_root.tag.split("}")[0].strip("{")}

        acquisition_time = extract_acquisition_time(xml_root, ns)
        product_id_elem = xml_root.find(".//rcm:productId", ns)
        product_id = product_id_elem.text.strip() if product_id_elem is not None else None

        range_spacing_elem = xml_root.find(".//rcm:sampledPixelSpacing", ns)
        azimuth_spacing_elem = xml_root.find(".//rcm:sampledLineSpacing", ns)

        pixel_spacing = {
            "range_m": float(range_spacing_elem.text) if range_spacing_elem is not None else None,
            "azimuth_m": float(azimuth_spacing_elem.text) if azimuth_spacing_elem is not None else None,
        }

        geocoded_points = []
        geodetic_coords = xml_root.findall(".//rcm:geodeticCoordinate", ns)
        for coord in geodetic_coords:
            lat = coord.find("rcm:latitude", ns)
            lon = coord.find("rcm:longitude", ns)
            if lat is not None and lon is not None:
                geocoded_points.append({
                    "latitude": float(lat.text),
                    "longitude": float(lon.text)
                })

        tie_points = _parse_tie_points(xml_root, ns)
        geometry = _detect_geometry(src_crs, src_transform, gcps_count)
        axis_mapping = _infer_axis_mapping_from_xml(pixel_spacing)

        out = {
            "folder_name": str(data_dir),
            "product_id": product_id,
            "hh": hh,
            "hv": hv,
            "geocoded_points": geocoded_points,
            "xml": xml_root,
            "src_transform": src_transform,
            "src_crs": src_crs,
            "src_bounds": src_bounds,
            "nodata_hh": np.nan,
            "nodata_hv": np.nan,
            "geometry": geometry,
            "transform_is_identity": _is_identity_transform(src_transform),
            "gcps_count": gcps_count,
            "gcp_crs": gcp_crs.to_string() if gcp_crs is not None else None,
            "axis_mapping": axis_mapping,
        }

        out["range_pixel_spacing_m"] = pixel_spacing["range_m"]
        out["azimuth_pixel_spacing_m"] = pixel_spacing["azimuth_m"]
        out["acquisition_time"] = acquisition_time
        out.update(tie_points)

        return out

    except Exception as e:
        print(f"Skipping {data_dir}: {e}")
        return None

# ============================================================
# Scaling functions
# ============================================================

def sar_value_scale(arr, name="SAR band", verbose=False):
    """
    Classify likely SAR raster value scale.

    Returns:
        {
            "scale": "db", "linear", "dn", or "unknown",
            "stats": {...},
            "dtype": "...",
            "warning": "..." or None
        }
    """
    warning = None

    # ------------------------------------------------------------
    # 0. Check dtype (important!)
    # ------------------------------------------------------------
    if np.issubdtype(arr.dtype, np.integer):
        warning = (
            f"{name}: input dtype is integer ({arr.dtype}). "
            "Values may be quantized, scaled, or stored as DN. "
            "Metadata should be checked before assuming dB or linear scale."
        )
        if verbose:
            print(f"[WARNING] {warning}")

    # Keep only valid values
    vals = arr[np.isfinite(arr)]

    if vals.size == 0:
        return {
            "name": name,
            "scale": "unknown",
            "dtype": str(arr.dtype),
            "stats": None,
            "warning": warning,
            "reason": "no finite values",
        }

    # ------------------------------------------------------------
    # 1. Robust statistics (avoid outliers)
    # ------------------------------------------------------------
    p0, p1, p50, p99, p100 = np.percentile(vals, [0, 1, 50, 99, 100])

    stats = {
        "min": float(p0),
        "p1": float(p1),
        "median": float(p50),
        "p99": float(p99),
        "max": float(p100),
    }

    # ------------------------------------------------------------
    # 2. dB detection (log-scaled SAR backscatter)
    # ------------------------------------------------------------
    looks_db = (
        p50 < 0 and # Median is typically negative for SAR in dB (e.g., -25 to -5 dB)

        p1 < -5 and # Lower tail should be clearly negative (water/low backscatter)

        p99 < 50 and # Upper tail should not be extremely large; typical SAR rarely exceeds ~30 dB

        p0 < 0 # Minimum must be negative → ensures it's not linear/DN
    )

    # ------------------------------------------------------------
    # 3. Linear σ⁰ detection (calibrated, not log-transformed)
    # ------------------------------------------------------------
    looks_linear = (
        p0 >= 0 and # Linear σ⁰ cannot be negative

        p50 < 1.0 and # Median usually small (e.g., 0.01–0.2 for most scenes)

        p99 < 5.0 # Upper tail limited; even bright scatterers rarely exceed a few units
    )

    # ------------------------------------------------------------
    # 4. DN / uncalibrated detection
    # ------------------------------------------------------------
    looks_dn = (
        p0 >= 0 and
        # DN values are non-negative

        (
            p50 > 10 or # Median already large → not σ⁰

            p99 > 100 or # High values → typical of raw/scaled intensity

            p100 > 255 # Exceeds 8-bit range → likely raw or scaled sensor values
        )
    )

    # ------------------------------------------------------------
    # 5. Final decision
    # ------------------------------------------------------------
    if looks_db:
        scale = "db"
    elif looks_linear:
        scale = "linear"
    elif looks_dn:
        scale = "dn"
    else:
        scale = "unknown"

    return {
        "name": name,
        "scale": scale,
        "dtype": str(arr.dtype),
        "stats": stats,
        "warning": warning,
    }


def scale_hh_hv_earth_geometry(rcm_data, target_spacing_m=200):
    dst_transform, dst_width, dst_height = calculate_default_transform(
        rcm_data["src_crs"],
        rcm_data["src_crs"],
        rcm_data["hh"].shape[1],
        rcm_data["hh"].shape[0],
        *rcm_data["src_bounds"],
        resolution=target_spacing_m
    )

    hh_out = np.empty((dst_height, dst_width), dtype=np.float32)
    hv_out = np.empty((dst_height, dst_width), dtype=np.float32)

    reproject(
        source=rcm_data["hh"],
        destination=hh_out,
        src_transform=rcm_data["src_transform"],
        src_crs=rcm_data["src_crs"],
        dst_transform=dst_transform,
        dst_crs=rcm_data["src_crs"],
        resampling=Resampling.average,
        src_nodata=rcm_data.get("nodata_hh", None),
        dst_nodata=np.nan
    )

    reproject(
        source=rcm_data["hv"],
        destination=hv_out,
        src_transform=rcm_data["src_transform"],
        src_crs=rcm_data["src_crs"],
        dst_transform=dst_transform,
        dst_crs=rcm_data["src_crs"],
        resampling=Resampling.average,
        src_nodata=rcm_data.get("nodata_hv", None),
        dst_nodata=np.nan
    )

    out = {
        "hh": hh_out,
        "hv": hv_out,
        "src_transform": dst_transform,
        "src_crs": rcm_data["src_crs"],
        "src_bounds": rcm_data["src_bounds"],
        "folder_name": rcm_data["folder_name"],
        "geometry": rcm_data["geometry"],
    }

    # Recompute tie-point image coordinates in the resampled earth-geometry grid.
    # This is more accurate than simple row/column scaling because the new raster
    # is defined by dst_transform and src_crs.
    tie_lons = rcm_data["tie_lons"]
    tie_lats = rcm_data["tie_lats"]

    # Convert geographic lon/lat to the raster CRS, if needed.
    xs, ys = transform_coords(
        "EPSG:4326",
        rcm_data["src_crs"],
        tie_lons.tolist(),
        tie_lats.tolist(),
    )

    # Convert map coordinates to row/column in the new resampled raster.
    tie_rows, tie_cols = rowcol(
        dst_transform,
        xs,
        ys,
        op=float
    )

    out.update({
        "tie_lines":  np.array(tie_rows, dtype=np.float64),
        "tie_pixels": np.array(tie_cols, dtype=np.float64),
        "tie_lats":   tie_lats,
        "tie_lons":   tie_lons,
    })
    

    return out


def scale_hh_hv_sensor_geometry(
    rcm_data,
    target_spacing_m=200.0,
    range_spacing_m_key="range_pixel_spacing_m",
    azimuth_spacing_m_key="azimuth_pixel_spacing_m",
):
    hh = rcm_data["hh"]
    hv = rcm_data["hv"]

    range_spacing = float(rcm_data[range_spacing_m_key])
    az_spacing = float(rcm_data[azimuth_spacing_m_key])

    if range_spacing <= 0 or az_spacing <= 0:
        raise ValueError(f"Invalid spacings: range={range_spacing}, azimuth={az_spacing}")

    src_h, src_w = hh.shape

    sx = target_spacing_m / range_spacing
    sy = target_spacing_m / az_spacing

    dst_w = max(1, int(np.round(src_w / sx)))
    dst_h = max(1, int(np.round(src_h / sy)))

    src_transform = Affine(range_spacing, 0, 0, 0, -az_spacing, 0)
    dst_transform = Affine(target_spacing_m, 0, 0, 0, -target_spacing_m, 0)

    dummy_crs = CRS.from_epsg(3857)

    hh_out = np.empty((dst_h, dst_w), dtype=np.float32)
    hv_out = np.empty((dst_h, dst_w), dtype=np.float32)

    reproject(
        source=hh,
        destination=hh_out,
        src_transform=src_transform,
        src_crs=dummy_crs,
        dst_transform=dst_transform,
        dst_crs=dummy_crs,
        resampling=Resampling.average,
        src_nodata=rcm_data.get("nodata_hh", None),
        dst_nodata=np.nan,
    )

    reproject(
        source=hv,
        destination=hv_out,
        src_transform=src_transform,
        src_crs=dummy_crs,
        dst_transform=dst_transform,
        dst_crs=dummy_crs,
        resampling=Resampling.average,
        src_nodata=rcm_data.get("nodata_hv", None),
        dst_nodata=np.nan,
    )

    out = {
        "hh": hh_out,
        "hv": hv_out,
        "src_transform": dst_transform,
        "src_crs": None,
        "src_bounds": None,
        "sensor_transform": dst_transform,
        range_spacing_m_key: target_spacing_m,
        azimuth_spacing_m_key: target_spacing_m,
        "folder_name": rcm_data.get("folder_name", None),
        "geometry": rcm_data["geometry"],
    }

    # scale tie points (note: this is a simple scaling of pixel/line coordinates, 
    # which may not be perfectly accurate for large transformations or non-linear distortions, 
    # but can be a reasonable approximation for small resampling)
    rs = dst_h / src_h
    cs = dst_w / src_w
    out.update({
        "tie_lines":  rcm_data["tie_lines"] * rs,
        "tie_pixels": rcm_data["tie_pixels"] * cs,
        "tie_lats": rcm_data["tie_lats"],
        "tie_lons": rcm_data["tie_lons"],
        })

    return out


def scale_hh_hv(rcm_data, target_spacing_m):

    # ------------------------------------------------------------
    # Convert to linear scale if currently in dB, to ensure proper spatial averaging during resampling.
    # average( log(σ⁰) ) is NOT physically correct
    # average(σ⁰) SAR true averaging
    # log((x1 + x2)/2) != (log(x1) + log(x2))/2
    
    def db_to_linear(x, nodata=None):
        out = np.full(x.shape, np.nan, dtype=np.float32)

        valid = np.isfinite(x)

        if nodata is not None:
            valid &= x != nodata

        out[valid] = np.power(10.0, x[valid] / 10.0)

        return out

    def linear_to_db(x, eps=1e-12):
        return (10.0 * np.log10(np.maximum(x, eps))).astype(np.float32)


    rcm_data = copy.deepcopy(rcm_data)

    hh_scale = sar_value_scale(rcm_data["hh"], name="HH band")["scale"]
    hv_scale = sar_value_scale(rcm_data["hv"], name="HV band")["scale"]

    assert hh_scale == hv_scale, f"HH and HV bands appear to be on different scales (HH: {hh_scale}, HV: {hv_scale}). Please verify the value ranges and metadata for both bands."
    assert hh_scale != "dn" and hv_scale != "dn", "Data appears to be in DN or uncalibrated format. Proper calibration (sigma nought) and scaling to physical units (dB or linear) is required before spatial resampling. Please check metadata and apply appropriate scaling."
    assert hh_scale != "unknown" and hv_scale != "unknown", "Review SAR value ranges to determine correct scale (dB, linear, or DN)."

    if hh_scale == "db":
        rcm_data["hh"] = db_to_linear(rcm_data["hh"], rcm_data["nodata_hh"])
        rcm_data["hv"] = db_to_linear(rcm_data["hv"], rcm_data["nodata_hv"])
    
    # ------------------------------------------------------------

    # ------------------------------------------------------------
    # Scale
    if rcm_data["geometry"] == "earth":
        rcm_data =  scale_hh_hv_earth_geometry(
                    rcm_data,
                    target_spacing_m=target_spacing_m,
                    )
    else: # sensor geometry
        rcm_data = scale_hh_hv_sensor_geometry(
                    rcm_data,
                    target_spacing_m=target_spacing_m,
                    range_spacing_m_key="range_pixel_spacing_m",
                    azimuth_spacing_m_key="azimuth_pixel_spacing_m",
                    )
    
    # ------------------------------------------------------------
    # Convert back to dB if original data was in dB, to maintain consistency.
    if hh_scale == "db":
        rcm_data["hh"] = linear_to_db(rcm_data["hh"])
        rcm_data["hv"] = linear_to_db(rcm_data["hv"])

    return rcm_data
