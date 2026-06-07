"""
No@,
May 5th, 2026

Extract training patches from annotated RCM products.

Pipeline
--------
1. Load native RCM products.
2. Rescale HH/HV SAR channels to the requested pixel spacing.
3. Convert dB SAR values to linear scale before spatial averaging:
       
       dB -> linear -> averaging -> dB

4. Load and rescale annotation labels.
5. Build land masks on the target SAR grid.
6. Create sliding-window patch indexes and exclude fully invalid patches:
       
       invalid_mask = np.isnan(HH) | land_mask | (labels == 255)

7. Save one .pkl file per patch containing:
       
       HH, HV, label,
       product_id, product_name,
       indexes,
       pixel_spacing,
       month, day

8. Compute IRGS unsupervised segmentations for each patch and save:
       
       seg, n_t, bound

Notes
-----
- Spatial resampling uses area-weighted averaging based on geometric overlap.
- Patch indexes are stored on the rescaled image grid:
      
      [(x1, y1), (x2, y2)]

- Labels:
      
      0   -> WATER
      1   -> ICE
      255 -> invalid / unlabeled
"""

import os
import time
import shutil
import argparse
from datetime import datetime
import multiprocessing

import joblib
import numpy as np
from PIL import Image
from irgs import IRGS
from skimage import measure

from parallel_stuff import Parallel
from rcm_functions import load_rcm_product, scale_hh_hv
from landmask_functions import (
    build_land_masks_earth_geometry,
    build_land_mask_sensor_geometry,
)
import glob


# ============================================================
# Constants
# ============================================================

WATER_RGB = (0, 0, 128)
ICE_RGB = (128, 0, 0)
INVALID_LABEL = 255


# ============================================================
# Argument parsing
# ============================================================

def Arguments():
    parser = argparse.ArgumentParser(
        description="Extract patches from annotated RCM products."
    )

    parser.add_argument(
        "--data_root",
        type=str,
        default="D:/Temp/CIS - Images - ArcticScope/03_dB_noise_sub_OK/",
        help="Root folder containing RCM product folders"
    )

    parser.add_argument(
        "--labels_root",
        type=str,
        default="D:/Temp/ArcticScope-executable_200m/Custom_Annotation/",
        help="Root folder containing annotation subfolders"
    )

    parser.add_argument(
        "--annotated_products_txt",
        type=str,
        default="annotated_products.txt",
        help="Txt file listing annotated product folders"
    )

    parser.add_argument(
        "--output",
        type=str,
        default="./patches/",
        help="Output root folder"
    )

    parser.add_argument(
        "--shp_path",
        type=str,
        default="./StatCan_ocean/StatCan_ocean.shp",
        help="Path to ocean shapefile used to derive land mask"
    )

    parser.add_argument(
        "--pixel_spacing",
        type=float,
        default=200.0,
        help="Target SAR pixel spacing in meters"
    )

    parser.add_argument(
        "--labels_pixel_spacing",
        type=float,
        default=200.0,
        help="Pixel spacing at which annotation PNGs were created"
    )

    parser.add_argument(
        "--patch_size",
        type=int,
        default=384,
        help="Patch size in pixels after resampling"
    )

    parser.add_argument(
        "--overlap",
        type=float,
        default=0.0,
        help="Patch overlap fraction in [0, 1)"
    )

    parser.add_argument(
        "--n_cores",
        type=int,
        default=multiprocessing.cpu_count() - 2,
        help="Number of CPU cores for parallel processing"
    )

    parser.add_argument(
        "--mask_res_m",
        type=float,
        default=200.0,
        help="Intermediate map resolution for sensor-geometry land mask"
    )

    parser.add_argument(
        "--simplify_m",
        type=float,
        default=0.0,
        help="Optional coastline simplification tolerance in meters"
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=0.5,
        help="Threshold for projected ocean-mask back-projection"
    )

    parser.add_argument(
        "--chunk_rows",
        type=int,
        default=512,
        help="Chunk rows for projected ocean-mask back-projection"
    )

    parser.add_argument(
        "--validate_label_shape",
        action="store_true",
        help="Check whether label shape matches SAR shape at labels_pixel_spacing"
    )

    parser.add_argument(
        "--irgs_classes",
        type=int,
        default=70,
        help="Number of classes for IRGS segmentation"
    )

    return parser.parse_args()


# ============================================================
# Basic helpers
# ============================================================

def read_product_list(txt_path):
    with open(txt_path, "r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


def resolve_product_dirs(data_root, annotated_products_txt):
    product_names = read_product_list(annotated_products_txt)
    product_dirs = []

    for name in product_names:
        full_path = os.path.join(data_root, name)
        if os.path.isdir(full_path):
            product_dirs.append(full_path)
        else:
            print(f"[WARNING] Missing product folder: {full_path}")

    return product_dirs


def spacing_to_tag(pixel_spacing):
    return (
        f"{int(pixel_spacing)}m"
        if float(pixel_spacing).is_integer()
        else f"{pixel_spacing}m"
    )


def get_label_path(labels_root, product_dir):
    product_name = os.path.basename(product_dir.rstrip("/\\"))
    return os.path.join(labels_root, product_name, "custom_annotation.png")


def extract_month_day_from_acquisition_time(dt):
    """
    Returns:
        month: 1..12
        day: day of year in [0..365]
    """
    if dt is None:
        return None, None

    start_of_year = datetime(dt.year, 1, 1, tzinfo=dt.tzinfo)

    month = dt.month
    day_of_year = (dt - start_of_year).days
    return month, day_of_year


# ============================================================
# Label handling
# ============================================================

def load_binary_annotation(label_path):
    """
    Convert RGB annotation to:
        0   -> WATER
        1   -> ICE
        255 -> invalid / unlabeled
    """
    img = np.array(Image.open(label_path).convert("RGB"), dtype=np.uint8)

    label = np.full(img.shape[:2], INVALID_LABEL, dtype=np.uint8)

    water_mask = np.all(img == WATER_RGB, axis=-1)
    ice_mask = np.all(img == ICE_RGB, axis=-1)

    label[water_mask] = 0
    label[ice_mask] = 1

    return label


def resize_label_nearest(label, out_shape):
    """
    Resize label map to out_shape=(rows, cols) using nearest-neighbor.
    """
    pil_img = Image.fromarray(label)
    resized = pil_img.resize((out_shape[1], out_shape[0]), resample=Image.NEAREST)
    return np.array(resized, dtype=np.uint8)


def validate_label_against_reference_spacing(native_product, labels, labels_pixel_spacing):
    """
    Optional consistency check:
    compare label PNG shape against SAR shape rescaled to labels_pixel_spacing.
    """
    try:
        reference_product = scale_hh_hv(native_product, target_spacing_m=labels_pixel_spacing)
        reference_shape = reference_product["hh"].shape
        return labels.shape == reference_shape, reference_shape
    except Exception:
        return False, None


# ============================================================
# Patch padding
# ============================================================

def pad_2d(arr, patch_size, mode="symmetric", constant_values=0):
    """
    Pad a 2D array to (patch_size, patch_size).
    """
    h, w = arr.shape
    pad_h = max(0, patch_size - h)
    pad_w = max(0, patch_size - w)

    if pad_h == 0 and pad_w == 0:
        return arr

    if mode == "constant":
        return np.pad(
            arr,
            ((0, pad_h), (0, pad_w)),
            mode=mode,
            constant_values=constant_values
        )

    return np.pad(arr, ((0, pad_h), (0, pad_w)), mode=mode)


# ============================================================
# Patch index generator
# ============================================================

class SlidePatchesIndex:
    """
    Generate sliding patch indexes on the already-rescaled image grid.

    Patches that are entirely inside invalid_mask are excluded.
    """

    def __init__(self, h_img, w_img, patch_size, overlap_percent, invalid_mask):
        self.h_crop = min(patch_size, h_img)
        self.w_crop = min(patch_size, w_img)

        self.h_stride = (
            self.h_crop - round(self.h_crop * overlap_percent)
            if self.h_crop < h_img else h_img
        )
        self.w_stride = (
            self.w_crop - round(self.w_crop * overlap_percent)
            if self.w_crop < w_img else w_img
        )

        self.h_grids = max(h_img - self.h_crop + self.h_stride - 1, 0) // self.h_stride + 1
        self.w_grids = max(w_img - self.w_crop + self.w_stride - 1, 0) // self.w_stride + 1

        self.patches_list = []

        for h_idx in range(self.h_grids):
            for w_idx in range(self.w_grids):
                y1 = h_idx * self.h_stride
                x1 = w_idx * self.w_stride

                y2 = min(y1 + self.h_crop, h_img)
                x2 = min(x1 + self.w_crop, w_img)

                y1 = max(y2 - self.h_crop, 0)
                x1 = max(x2 - self.w_crop, 0)

                # exclude only if the full patch is invalid
                if not invalid_mask[y1:y2, x1:x2].all():
                    self.patches_list.append((y1, y2, x1, x2))

    def __getitem__(self, index):
        return self.patches_list[index]

    def __len__(self):
        return len(self.patches_list)


# ============================================================
# Scene preparation
# ============================================================

def build_land_mask_for_scaled_product(args, scaled_product):
    """
    Route land-mask computation depending on geometry.
    Returns boolean mask with True=land, False=other.
    """
    if scaled_product["geometry"] == "earth":
        land_mask = build_land_masks_earth_geometry(args.shp_path, scaled_product)
    else:
        land_mask = build_land_mask_sensor_geometry(
            scaled_product,
            args.shp_path,
            mask_res_m=args.mask_res_m,
            simplify_m=args.simplify_m,
            threshold=args.threshold,
            chunk_rows=args.chunk_rows,
        )

    if land_mask is None:
        # fail-safe: no land
        land_mask = np.zeros_like(scaled_product["hh"], dtype=bool)

    return land_mask


def prepare_product_for_patching(args, product_dir):
    """
    Load SAR + labels, rescale to target spacing, build land mask, and
    generate patch indexes.
    """
    native_product = load_rcm_product(product_dir)
    if native_product is None:
        raise RuntimeError(f"Could not load product: {product_dir}")

    product_name = os.path.basename(product_dir.rstrip("/\\"))
    product_id = native_product.get("product_id", None)
    acquisition_time = native_product.get("acquisition_time", None)

    # --- target SAR
    scaled_product = scale_hh_hv(native_product, target_spacing_m=args.pixel_spacing)
    hh = scaled_product["hh"]
    hv = scaled_product["hv"]

    # --- land mask on target grid
    land_mask = build_land_mask_for_scaled_product(args, scaled_product)

    # --- labels
    label_path = get_label_path(args.labels_root, product_dir)
    if not os.path.isfile(label_path):
        raise FileNotFoundError(f"Annotation not found: {label_path}")

    labels = load_binary_annotation(label_path)

    if args.validate_label_shape:
        is_match, reference_shape = validate_label_against_reference_spacing(
            native_product,
            labels,
            args.labels_pixel_spacing
        )
        if not is_match:
            print(
                f"[WARNING] Label shape mismatch for {product_name}: "
                f"label={labels.shape}, expected_at_{args.labels_pixel_spacing}m={reference_shape}"
            )

    # labels were created at a fixed pixel spacing; rescale to target image spacing
    labels = resize_label_nearest(labels, hh.shape)
    
    # make sure all invalid pixels have the same label value for easier patch exclusion
    invalid_mask = np.isnan(hh) | land_mask | (labels == INVALID_LABEL)
    labels[invalid_mask] = INVALID_LABEL

    patch_idx = SlidePatchesIndex(
        h_img=hh.shape[0],
        w_img=hh.shape[1],
        patch_size=args.patch_size,
        overlap_percent=args.overlap,
        invalid_mask=invalid_mask,
    )

    month, day = extract_month_day_from_acquisition_time(acquisition_time)

    return {
        "product_dir": product_dir,
        "product_name": product_name,
        "product_id": product_id,
        "hh": hh,
        "hv": hv,
        "labels": labels,
        "patch_idx": patch_idx,
        "pixel_spacing": float(args.pixel_spacing),
        "month": month,
        "day": day,
    }


def get_patch_index(args, product_dir):
    prepared = prepare_product_for_patching(args, product_dir)
    return prepared["patch_idx"]


# ============================================================
# Patch extraction
# ============================================================

def Extract_patches(args, item):
    """
    Extract patches for one product and save them.

    item:
        (product_dir, patch_idx)
    """
    product_dir, patch_idx = item
    prepared = prepare_product_for_patching(args, product_dir)

    print(f'from product: {prepared["product_name"]}')

    product_name = prepared["product_name"]
    product_id = prepared["product_id"]
    hh = prepared["hh"]
    hv = prepared["hv"]
    labels = prepared["labels"]
    pixel_spacing = prepared["pixel_spacing"]
    month = prepared["month"]
    day = prepared["day"]

    output_folder = os.path.join(
        args.output,
        f"pixel_spacing_{spacing_to_tag(pixel_spacing)}",
        product_name,
    )

    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    os.makedirs(output_folder, exist_ok=True)

    if len(patch_idx) == 0:
        print(f"Number of patches = 0 for scene {product_name}")
        return 0

    for i in range(len(patch_idx)):
        y1, y2, x1, x2 = patch_idx[i]

        hh_patch = hh[y1:y2, x1:x2]
        hv_patch = hv[y1:y2, x1:x2]
        label_patch = labels[y1:y2, x1:x2]

        hh_patch = pad_2d(hh_patch, args.patch_size, mode="symmetric")
        hv_patch = pad_2d(hv_patch, args.patch_size, mode="symmetric")
        label_patch = pad_2d(
            label_patch,
            args.patch_size,
            mode="constant",
            constant_values=INVALID_LABEL
        )
        
        data_patch = {
            "HH": hh_patch.astype(np.float32),
            "HV": hv_patch.astype(np.float32),
            "pixel_labels": label_patch.astype(np.uint8),
            "product_id": product_id,
            "product_name": product_name,
            "indexes": [(x1, y1), (x2, y2)],
            "pixel_spacing": float(pixel_spacing),
            "month": int(month) if month is not None else None,
            "day": int(day) if day is not None else None,
        }

        joblib.dump(data_patch, os.path.join(output_folder, f"{i:05d}.pkl"))

        # Image.fromarray(np.uint8((hh_patch - np.nanmin(hh_patch)) / (np.nanmax(hh_patch) - np.nanmin(hh_patch) + 1e-8) * 255)).save(f"{output_folder}/{i:05d}_hh.png" )
        # Image.fromarray(np.uint8((hv_patch - np.nanmin(hv_patch)) / (np.nanmax(hv_patch) - np.nanmin(hv_patch) + 1e-8) * 255)).save(f"{output_folder}/{i:05d}_hv.png" )
        # l = label_patch.copy(); l[l==INVALID_LABEL] = 255; l[l==0] = 0; l[l==1] = 128
        # Image.fromarray(np.uint8(l)).save(f"{output_folder}/{i:05d}_label.png" )

    return len(patch_idx)

def compute_irgs(args, file_path):

    data_patch = joblib.load(file_path)
    hh_patch = data_patch["HH"]
    hv_patch = data_patch["HV"]
    label_patch = data_patch["pixel_labels"]

    #============== Unsupervised segmentation using IRGS ==============
    # These variabes are used to train the ITT framework
    hh_normalized = (hh_patch - np.nanmin(hh_patch)) / (np.nanmax(hh_patch) - np.nanmin(hh_patch) + 1e-8)
    hv_normalized = (hv_patch - np.nanmin(hv_patch)) / (np.nanmax(hv_patch) - np.nanmin(hv_patch) + 1e-8)
    valid_mask = ~(label_patch == INVALID_LABEL)
    irgs_output, boundaries = IRGS(np.stack([hh_normalized, hv_normalized], axis=-1), 
                                    n_classes=args.irgs_classes, n_iter=200, 
                                    mask=~np.isnan(hh_normalized))  
    irgs_output[~valid_mask] = -1
    boundaries [~valid_mask] = -1
    
    # IRGS returns a segmentation map with arbitrary integer labels for each segment, 
    # so we need to re-map to assign a unique label to each segment and to make sure that all invalid pixels have the same label value (-1 in this case).
    components, n_tokens = measure.label(irgs_output, background=-1, return_num=True, connectivity=2)
    segments = components - 1
    # when the tokens are calculated on the fly (ITT framework)
    # the ids are changed so that they're exclusive among the sample images
    # within a batch. In that sense, -1 ids need to be changed first to infinite 
    # (the reason for this is restricted to the way the ITT code was implemented)
    segments = segments.astype('float')
    segments[segments==-1] = np.inf
    #============== ==================================== ==============
    
    data_patch.update({
        "segments": segments,
        "n_segments": n_tokens,
        "boundaries": boundaries,
    })
    
    joblib.dump(data_patch, file_path)

    # segments[segments==np.inf] = -1
    # Image.fromarray(np.uint8(segments*255/(n_tokens+1))).save(f"{output_folder}/{i:05d}_irgs.png")
    # Image.fromarray(255*np.uint8(boundaries==-1)).save(f"{output_folder}/{i:05d}_boundaries.png")        
    return file_path

# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    args = Arguments()

    product_dirs = resolve_product_dirs(args.data_root, args.annotated_products_txt)
    print(f"Number of annotated products: {len(product_dirs)}")

    if len(product_dirs) == 0:
        raise RuntimeError("No valid annotated products were found.")

    output_spacing_folder = os.path.join(
        args.output,
        f"pixel_spacing_{spacing_to_tag(args.pixel_spacing)}"
    )

    if os.path.isdir(output_spacing_folder):
        shutil.rmtree(output_spacing_folder)
    if not os.path.isdir(output_spacing_folder):
        os.makedirs(output_spacing_folder)

    # ---------------- GET INDEXES
    start_time = time.time()
    print("Calculating patch indexes ...")

    iterable = iter(product_dirs)
    if len(product_dirs) > 1:
        patches_idx = Parallel(get_patch_index, iterable, args, n_cores=args.n_cores)
    else:
        patches_idx = [get_patch_index(args, iter) for iter in iterable]

    print(f"Patches index generated! - time: {time.time() - start_time:.2f} s")
    print(
        f"Number of patches for pixel spacing {args.pixel_spacing} m: "
        f"{sum(len(l) for l in patches_idx)}"
    )

    # ---------------- EXTRACT PATCHES
    start_time = time.time()
    print("Saving patches ...")

    iterable = zip(product_dirs, patches_idx)
    if len(product_dirs) > 1:
        patch_counts = Parallel(Extract_patches, iterable, args, n_cores=args.n_cores)
    else:
        patch_counts = [Extract_patches(args, iter) for iter in iterable]

    print(f"Patches saved! - time: {time.time() - start_time:.2f} s")
    print(f"Total saved patches: {sum(patch_counts)}")

    # ---------------- COMPUTE IRGS
    start_time = time.time()
    print("Computing IRGS segments ...")
    file_paths = glob.glob(args.output + "/**/*.pkl", recursive=True)
    if len(file_paths) > 1:
        Parallel(compute_irgs, file_paths, args, n_cores=args.n_cores)
    else:
        for file_path in file_paths:
            compute_irgs(args, file_path)
    
    print(f"IRGS segments computed! - time: {time.time() - start_time:.2f} s")
