"""
No@
March 19, 2026

Create train/validation splits for annotated RCM SAR products using
multi-class pixel-level annotations, while balancing the split using only
WATER and ICE.

The script:
1. Reads annotated product folder names from a text file.
2. Loads annotations from labels_root/product_dir/custom_annotation.png.
3. Converts annotation colors into class labels:
       WATER   = (0, 0, 128)
       ICE     = (128, 0, 0)
       SHOALS  = (0, 255, 0)
       SHIP    = (255, 255, 0)
       ICEBERG = (255, 0, 255)
       UNKNOWN = (150, 150, 150)
4. Computes per-scene and global class statistics for all classes.
5. Selects a validation set that preserves WATER and ICE whenever possible.
6. Selects a training set to improve pixel-level balance between WATER and ICE.
7. Writes train_XX.txt, val_YY.txt, summary JSON, and class-distribution plots.
"""

import os
import json
import time
import argparse

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

from rcm_functions import load_rcm_product
from parallel_stuff import Parallel


# ============================================================
# Annotation utilities
# ============================================================

LABELS = {
    0: "WATER",
    1: "ICE",
    2: "SHOALS",
    3: "SHIP",
    4: "ICEBERG",
    5: "UNKNOWN",
}

WATER_RGB = (0, 0, 128)
ICE_RGB = (128, 0, 0)
SHOALS_RGB = (0, 255, 0)
SHIP_RGB = (255, 255, 0)
ICEBERG_RGB = (255, 0, 255)
UNKNOWN_RGB = (150, 150, 150)
INVALID_LABEL = 255

COLOR_TO_LABEL = {
    WATER_RGB: 0,
    ICE_RGB: 1,
    SHOALS_RGB: 2,
    SHIP_RGB: 3,
    ICEBERG_RGB: 4,
    UNKNOWN_RGB: 5,
}

N_CLASSES = len(LABELS)

BALANCE_CLASS_IDS = [0, 1]  # WATER, ICE


def read_product_list(txt_path):
    with open(txt_path, "r", encoding="utf-8") as f:
        product_dirs = [line.strip() for line in f if line.strip()]
    return product_dirs


def find_existing_products(data_root, product_dirs):
    out = []
    for p in product_dirs:
        full_path = os.path.join(data_root, p)
        if os.path.isdir(full_path):
            out.append(full_path)
        else:
            print(f"Warning: product folder not found: {full_path}")
    return out


def get_label_path(labels_root, product_dir):
    product_name = os.path.basename(product_dir.rstrip("/\\"))
    return os.path.join(labels_root, product_name, "custom_annotation.png")


def load_multiclass_annotation(label_path):
    """
    Convert RGB annotation to a class map:
        0 -> WATER
        1 -> ICE
        2 -> SHOALS
        3 -> SHIP
        4 -> ICEBERG
        5 -> UNKNOWN
        255 -> invalid / unlabeled
    """
    img = np.array(Image.open(label_path).convert("RGB"), dtype=np.uint8)

    class_map = np.full(img.shape[:2], INVALID_LABEL, dtype=np.uint8)

    for rgb, cls_id in COLOR_TO_LABEL.items():
        mask = np.all(img == rgb, axis=-1)
        class_map[mask] = cls_id

    return class_map


def compute_histograms_from_mask(mask):
    """
    Returns:
        hist_pixel: np.ndarray of shape (N_CLASSES,)
            Pixel-level counts for all classes.

        scene_presence: np.ndarray of shape (N_CLASSES,)
            Binary class presence for the current scene.
    """
    hist_pixel = np.zeros(N_CLASSES, dtype=np.int64)
    scene_presence = np.zeros(N_CLASSES, dtype=np.int64)

    valid = mask != INVALID_LABEL
    if not np.any(valid):
        return hist_pixel, scene_presence

    labels, counts = np.unique(mask[valid], return_counts=True)
    for lbl, cnt in zip(labels, counts):
        lbl = int(lbl)
        if 0 <= lbl < N_CLASSES:
            hist_pixel[lbl] += int(cnt)
            scene_presence[lbl] = 1

    return hist_pixel, scene_presence


# ============================================================
# Worker
# ============================================================

def Obtain_histograms(iterator):
    """
    iterator:
        (product_dir, idx, labels_root, validate_with_sar, n_scenes)

    Returns:
        hist_pixel: np.ndarray of shape (N_CLASSES,)
            Pixel-level class counts for the current scene.

        scene_count: np.ndarray of shape (N_CLASSES, n_scenes)
            Binary class presence matrix for the current scene, where only the
            column corresponding to this scene index is populated.

        scene_pixels_balance: np.ndarray of shape (2,)
            Per-scene pixel counts for WATER and ICE only:
            [num_water_pixels, num_ice_pixels].
    """
    product_dir, idx, labels_root, validate_with_sar, n_scenes = iterator

    print(f"Processing scene: {product_dir}")

    hist_pixel = np.zeros(N_CLASSES, dtype=np.int64)
    scene_count = np.zeros((N_CLASSES, n_scenes), dtype=np.int64)
    scene_pixels_balance = np.zeros(2, dtype=np.int64)

    label_path = get_label_path(labels_root, product_dir)
    if not os.path.isfile(label_path):
        print(f"Warning: annotation not found for {product_dir}")
        return hist_pixel, scene_count, scene_pixels_balance

    try:
        mask = load_multiclass_annotation(label_path)

        if validate_with_sar:
            rcm_product = load_rcm_product(product_dir)
            if rcm_product is None:
                print(f"Warning: could not load SAR product for {product_dir}")
                return hist_pixel, scene_count, scene_pixels_balance

            if mask.shape != rcm_product["hh"].shape:
                print(
                    f"Warning: label/SAR size mismatch for {product_dir} | "
                    f"label={mask.shape}, sar={rcm_product['hh'].shape}"
                )
                return hist_pixel, scene_count, scene_pixels_balance

        h, s = compute_histograms_from_mask(mask)
        hist_pixel += h
        scene_count[:, idx] = s

        scene_pixels_balance[0] = h[0]  # WATER
        scene_pixels_balance[1] = h[1]  # ICE

    except Exception as e:
        print(f"Error processing {product_dir}: {e}")

    return hist_pixel, scene_count, scene_pixels_balance


# ============================================================
# Plotting
# ============================================================

def save_barplot(values, labels, out_path, title):
    plt.figure(figsize=(10, 4))
    bars = plt.bar(np.arange(len(values)), values, width=0.6)
    plt.xticks(np.arange(len(values)), labels, rotation=45)
    plt.title(title)

    for bar in bars:
        height = bar.get_height()
        text = f"{height:.1f}" if isinstance(height, (float, np.floating)) else f"{int(height)}"
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            text,
            ha="center",
            va="bottom",
            fontsize=8
        )

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


# ============================================================
# Split selection
# ============================================================

def select_train_val_balanced(scene_pixel_counts_balance, train_percentage):
    """
    Select train/val splits using only WATER and ICE for balancing.

    Parameters
    ----------
    scene_pixel_counts_balance : np.ndarray
        Shape (n_scenes, 2), columns are [water_pixels, ice_pixels].
    train_percentage : int
        Percentage of scenes to assign to training.

    Returns
    -------
    train_ids : list[int]
    val_ids   : list[int]
    """
    n_scenes = scene_pixel_counts_balance.shape[0]

    if n_scenes == 0:
        return [], []

    n_train = max(1, n_scenes * train_percentage // 100)
    if n_scenes > 1:
        n_train = min(n_train, n_scenes - 1)
    else:
        n_train = 1

    all_ids = list(range(n_scenes))
    n_val = n_scenes - n_train

    if n_val <= 0:
        return sorted(all_ids), []

    has_water = scene_pixel_counts_balance[:, 0] > 0
    has_ice = scene_pixel_counts_balance[:, 1] > 0
    has_both = has_water & has_ice

    val_ids = []

    both_ids = np.where(has_both)[0].tolist()
    if both_ids:
        def balance_score(i):
            w, ice = scene_pixel_counts_balance[i]
            total = w + ice
            if total == 0:
                return 1.0
            return abs((w / total) - 0.5)

        best_both = min(both_ids, key=balance_score)
        val_ids.append(best_both)
    else:
        water_ids = np.where(has_water)[0].tolist()
        ice_ids = np.where(has_ice)[0].tolist()

        if water_ids:
            val_ids.append(water_ids[0])

        ice_candidates = [i for i in ice_ids if i not in val_ids]
        if ice_candidates and len(val_ids) < n_val:
            val_ids.append(ice_candidates[0])

    remaining = [i for i in all_ids if i not in val_ids]

    def imbalance_score(i):
        w, ice = scene_pixel_counts_balance[i]
        total = w + ice
        if total == 0:
            return 1.0
        return abs((w / total) - 0.5)

    remaining_sorted_for_val = sorted(remaining, key=imbalance_score, reverse=True)

    for sid in remaining_sorted_for_val:
        if len(val_ids) >= n_val:
            break
        val_ids.append(sid)

    val_ids = sorted(val_ids[:n_val])

    remaining_train = [i for i in all_ids if i not in val_ids]

    train_ids = []
    train_water = 0
    train_ice = 0

    while remaining_train and len(train_ids) < n_train:
        best_sid = None
        best_score = None
        best_total = None

        for sid in remaining_train:
            w, ice = scene_pixel_counts_balance[sid]
            new_w = train_water + w
            new_ice = train_ice + ice
            total = new_w + new_ice

            if total == 0:
                score = 1.0
            else:
                score = abs((new_w / total) - 0.5)

            current_total = w + ice

            if (
                best_score is None
                or score < best_score
                or (score == best_score and current_total > best_total)
            ):
                best_score = score
                best_sid = sid
                best_total = current_total

        train_ids.append(best_sid)
        train_water += scene_pixel_counts_balance[best_sid, 0]
        train_ice += scene_pixel_counts_balance[best_sid, 1]
        remaining_train.remove(best_sid)

    return sorted(train_ids), sorted(val_ids)


# ============================================================
# Args
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Create train/validation splits for annotated RCM products."
    )

    parser.add_argument(
        "--data_root",
        type=str,
        default="D:/Temp/CIS - Images - ArcticScope/03_dB_noise_sub_OK/",
        help="Root directory containing RCM product folders."
    )

    parser.add_argument(
        "--labels_root",
        type=str,
        default="D:/Temp/ArcticScope-executable_200m/Custom_Annotation/",
        help="Root directory containing annotation subfolders."
    )

    parser.add_argument(
        "--annotated_products_txt",
        type=str,
        default="annotated_products.txt",
        help="Text file listing annotated product folder names."
    )

    parser.add_argument(
        "--percentage",
        type=int,
        default=80,
        help="Percentage of scenes to assign to training."
    )

    parser.add_argument(
        "--validate_with_sar",
        action="store_true",
        help="Validate annotation size against the SAR raster size."
    )

    return parser.parse_args()


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    args = parse_args()

    data_root = args.data_root
    labels_root = args.labels_root
    annotated_products_txt = args.annotated_products_txt
    percentage = args.percentage
    validate_with_sar = args.validate_with_sar

    if percentage < 1 or percentage > 99:
        raise ValueError("--percentage must be between 1 and 99")

    print("========== CONFIG ==========")
    print(f"data_root: {data_root}")
    print(f"labels_root: {labels_root}")
    print(f"annotated_products_txt: {annotated_products_txt}")
    print(f"percentage (train): {percentage}")
    print(f"validate_with_sar: {validate_with_sar}")
    print("============================")

    product_dirs_rel = read_product_list(annotated_products_txt)
    product_dirs = find_existing_products(data_root, product_dirs_rel)

    n_scenes = len(product_dirs)
    if n_scenes == 0:
        raise RuntimeError("No valid annotated product folders were found.")

    n_train = max(1, n_scenes * percentage // 100)
    n_val = n_scenes - n_train

    print(f"Number of annotated scenes: {n_scenes}")
    print(f"Requested train scenes: {n_train}")
    print(f"Requested val scenes: {n_val}")

    start_time = time.time()

    iterator = [
        (product_dir, idx, labels_root, validate_with_sar, n_scenes)
        for idx, product_dir in enumerate(product_dirs)
    ]

    if len(product_dirs) > 1:
        hist = Parallel(Obtain_histograms, iterator)
    else:
        hist = [Obtain_histograms(item) for item in iterator]

    end_time = time.time()
    print(f"Processing time: {end_time - start_time:.2f} seconds")

    hist_pixel = np.zeros(N_CLASSES, dtype=np.int64)
    scene_count = np.zeros((N_CLASSES, n_scenes), dtype=np.int64)
    scene_pixel_counts_balance = np.zeros((n_scenes, 2), dtype=np.int64)

    for idx, item in enumerate(hist):
        h_p, s_c, s_p_balance = item
        hist_pixel += h_p
        scene_count += s_c
        scene_pixel_counts_balance[idx, :] = s_p_balance

    hist_scene = np.sum(scene_count, axis=1)

    class_names = [LABELS[i] for i in range(N_CLASSES)]

    if hist_pixel.sum() > 0:
        pixel_percent = hist_pixel / hist_pixel.sum() * 100.0
    else:
        pixel_percent = np.zeros_like(hist_pixel, dtype=np.float64)

    save_barplot(
        pixel_percent,
        class_names,
        "./data_info/classes_stats_pixel.png",
        "Pixel distribution (%)"
    )

    save_barplot(
        hist_scene,
        class_names,
        "./data_info/classes_stats_scene.png",
        "Scene coverage"
    )

    train_ids, val_ids = select_train_val_balanced(scene_pixel_counts_balance, percentage)

    # summarize train/val split using WATER and ICE only
    train_water = int(scene_pixel_counts_balance[train_ids, 0].sum()) if train_ids else 0
    train_ice = int(scene_pixel_counts_balance[train_ids, 1].sum()) if train_ids else 0
    val_water = int(scene_pixel_counts_balance[val_ids, 0].sum()) if val_ids else 0
    val_ice = int(scene_pixel_counts_balance[val_ids, 1].sum()) if val_ids else 0

    train_total_pixels_balance = train_water + train_ice
    val_total_pixels_balance = val_water + val_ice
    total_pixels_balance = int(hist_pixel[0] + hist_pixel[1])

    train_water_pct = 100.0 * train_water / train_total_pixels_balance if train_total_pixels_balance > 0 else 0.0
    train_ice_pct = 100.0 * train_ice / train_total_pixels_balance if train_total_pixels_balance > 0 else 0.0
    val_water_pct = 100.0 * val_water / val_total_pixels_balance if val_total_pixels_balance > 0 else 0.0
    val_ice_pct = 100.0 * val_ice / val_total_pixels_balance if val_total_pixels_balance > 0 else 0.0

    total_water_pct = 100.0 * int(hist_pixel[0]) / total_pixels_balance if total_pixels_balance > 0 else 0.0
    total_ice_pct = 100.0 * int(hist_pixel[1]) / total_pixels_balance if total_pixels_balance > 0 else 0.0

    train_file = f"./data_info/train_{percentage}.txt"
    val_file = f"./data_info/val_{100 - percentage}.txt"

    with open(train_file, "w", encoding="utf-8") as f:
        for sid in train_ids:
            f.write(os.path.basename(product_dirs[sid]) + "\n")

    with open(val_file, "w", encoding="utf-8") as f:
        for sid in val_ids:
            f.write(os.path.basename(product_dirs[sid]) + "\n")

    summary = {
        "n_scenes": n_scenes,
        "n_train": len(train_ids),
        "n_val": len(val_ids),
        "train_percentage": percentage,
        "val_percentage": 100 - percentage,
        "total_pixel_percentages_all_classes": {
            LABELS[i]: round(float(pixel_percent[i]), 4) for i in range(N_CLASSES)
        },
        "scene_counts_all_classes": {
            LABELS[i]: int(hist_scene[i]) for i in range(N_CLASSES)
        },
        "balance_classes": ["WATER", "ICE"],
        "total_pixel_percentages_balance_classes": {
            "WATER": round(total_water_pct, 2),
            "ICE": round(total_ice_pct, 2),
        },
        "train_pixel_percentages_balance_classes": {
            "WATER": round(train_water_pct, 2),
            "ICE": round(train_ice_pct, 2),
        },
        "val_pixel_percentages_balance_classes": {
            "WATER": round(val_water_pct, 2),
            "ICE": round(val_ice_pct, 2),
        },
        "train_scenes": [os.path.basename(product_dirs[sid]) for sid in train_ids],
        "val_scenes": [os.path.basename(product_dirs[sid]) for sid in val_ids],
    }

    with open("./data_info/train_val_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"Saved train split to: {train_file}")
    print(f"Saved val split to: {val_file}")
    print("Saved summary to: ./data_info/train_val_summary.json")
    print(f"train scenes: {len(train_ids)}")
    print(f"val scenes: {len(val_ids)}")
    print(f"train WATER: {train_water_pct:.2f}%")
    print(f"train ICE: {train_ice_pct:.2f}%")
    print(f"val WATER: {val_water_pct:.2f}%")
    print(f"val ICE: {val_ice_pct:.2f}%")