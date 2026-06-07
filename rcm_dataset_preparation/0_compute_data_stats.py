'''
No@
March 19, 2026

Compute global statistics (mean, standard deviation, minimum, and maximum) for
HH and HV channels across a dataset of SAR products at multiple pixel spacings.

For each product, the script:
1. Loads HH/HV data and metadata (including geometry type).
2. Resamples the data to specified pixel spacings (e.g., 50, 100, 200 m)
   using average resampling in either Earth or sensor geometry.
3. Accumulates partial statistics (count, sum, sum of squares, min, max)
   for each channel and spacing.

Processing is parallelized across products. Final global statistics are computed
using aggregated sums and saved to a JSON file with the hierarchy:
    channel (HH/HV) → pixel_spacing → {mean, std, min, max}.

'''

import os
import glob
import json
import math
import traceback

import numpy as np

import argparse
from parallel_stuff import Parallel
from rcm_functions import load_rcm_product, scale_hh_hv

from datetime import datetime
import matplotlib.pyplot as plt


# conda create -n rcm_dataset python numpy matplotlib rasterio affine lxml

# ============================================================
# Stats helpers
# ============================================================

def init_partial_stats():
    return {
        "count": 0,
        "sum": 0.0,
        "sum_sq": 0.0,
        "min": np.inf,
        "max": -np.inf,
    }


def update_partial_stats(stats, arr):
    if arr is None:
        return stats

    valid = np.isfinite(arr)
    if not np.any(valid):
        return stats

    vals = arr[valid].astype(np.float64, copy=False)

    stats["count"] += vals.size
    stats["sum"] += vals.sum()
    stats["sum_sq"] += np.square(vals).sum()

    vmin = vals.min()
    vmax = vals.max()

    if vmin < stats["min"]:
        stats["min"] = float(vmin)
    if vmax > stats["max"]:
        stats["max"] = float(vmax)

    return stats


def merge_partial_stats(a, b):
    if b["count"] == 0:
        return a
    if a["count"] == 0:
        return b.copy()

    return {
        "count": a["count"] + b["count"],
        "sum": a["sum"] + b["sum"],
        "sum_sq": a["sum_sq"] + b["sum_sq"],
        "min": min(a["min"], b["min"]),
        "max": max(a["max"], b["max"]),
    }


def finalize_stats(stats):
    if stats["count"] == 0:
        return {
            "mean": None,
            "std": None,
            "min": None,
            "max": None,
        }

    mean = stats["sum"] / stats["count"]
    var = (stats["sum_sq"] / stats["count"]) - (mean ** 2)
    var = max(0.0, var)  # numerical safety
    std = math.sqrt(var)

    return {
        "mean": float(mean),
        "std": float(std),
        "min": float(stats["min"]),
        "max": float(stats["max"]),
    }

def extract_month_from_acquisition_time(acquisition_time):
    """
    Convert an acquisition timestamp string into month number [1..12].
    Returns None if missing or invalid.
    """
    if acquisition_time is None:
        return None
    
    return acquisition_time.month


def save_month_histogram(worker_results, out_json_path):
    """
    Count successful products by acquisition month and save a bar plot.
    """
    month_counts = np.zeros(12, dtype=np.int64)

    for res in worker_results:
        if not res["ok"]:
            continue

        month = res.get("acquisition_month", None)
        if month is None:
            continue

        if 1 <= month <= 12:
            month_counts[month - 1] += 1

    month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    plt.figure(figsize=(10, 4))
    bars = plt.bar(np.arange(12), month_counts, width=0.7)
    plt.xticks(np.arange(12), month_labels)
    plt.ylabel("Number of scenes")
    plt.title("RCM scene count by acquisition month")

    for bar in bars:
        height = int(bar.get_height())
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            f"{height}",
            ha="center",
            va="bottom",
            fontsize=8
        )

    plt.tight_layout()

    plt.savefig("./data_info/month_hist_all_data.png", dpi=300, bbox_inches="tight")
    plt.close()

    return month_counts

# ============================================================
# Worker
# ============================================================


def process_one_product(spacings, data_dir):
    """
    Returns per-product partial stats for all requested spacings:
    {
        "HH": {"50": partial_stats, "100": partial_stats, ...},
        "HV": {"50": partial_stats, "100": partial_stats, ...},
        "product": "...",
        "ok": True/False,
        "error": ...
    }
    """
    result = {
        "HH": {str(s): init_partial_stats() for s in spacings},
        "HV": {str(s): init_partial_stats() for s in spacings},
        "product": str(data_dir),
        "ok": False,
        "error": None,
        "acquisition_month": None,
    }

    try:
        rcm_product = load_rcm_product(data_dir)
        if rcm_product is None:
            result["error"] = "load_rcm_product returned None"
            return result

        result["acquisition_month"] = extract_month_from_acquisition_time(
            rcm_product.get("acquisition_time", None)
            )
        
        for spacing in spacings:
            scaled = scale_hh_hv(rcm_product, target_spacing_m=spacing)

            # plt.imshow(scaled["hv"], cmap='gray'); plt.colorbar(); plt.show()
            result["HH"][str(spacing)] = update_partial_stats(
                result["HH"][str(spacing)],
                scaled["hh"]
            )
            result["HV"][str(spacing)] = update_partial_stats(
                result["HV"][str(spacing)],
                scaled["hv"]
            )

        result["ok"] = True
        return result

    except Exception as e:
        result["error"] = f"{type(e).__name__}: {e}\n{traceback.format_exc()}"
        return result


# ============================================================
# Reducer
# ============================================================

def reduce_results(worker_results, spacings):
    global_stats = {
        "HH": {str(s): init_partial_stats() for s in spacings},
        "HV": {str(s): init_partial_stats() for s in spacings},
    }

    failures = []

    for res in worker_results:
        if not res["ok"]:
            failures.append({
                "product": res["product"],
                "error": res["error"],
            })
            continue

        for spacing in spacings:
            s = str(spacing)
            global_stats["HH"][s] = merge_partial_stats(global_stats["HH"][s], res["HH"][s])
            global_stats["HV"][s] = merge_partial_stats(global_stats["HV"][s], res["HV"][s])

    final_json = {
        "HH": {},
        "HV": {}
    }

    for channel in ["HH", "HV"]:
        for spacing in spacings:
            s = str(spacing)
            final_json[channel][s] = finalize_stats(global_stats[channel][s])

    return final_json, failures


# ============================================================
# Main entry point
# ============================================================


def load_product_dirs_from_txt(data_root, txt_path):
    """
    Load product directories from a txt file.

    Each line should contain either:
        - folder name
        - or relative path from data_root

    Returns list of full paths that exist.
    """
    with open(txt_path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]

    data_dirs = []
    for line in lines:
        full_path = os.path.join(data_root, line)
        if os.path.isdir(full_path):
            data_dirs.append(full_path)
        else:
            print(f"[WARNING] Skipping missing folder: {full_path}")

    return sorted(data_dirs)

def compute_dataset_stats_json(
    data_root,
    out_json_path,
    spacings=(50, 100, 200),
    pattern="*",
):
    data_dirs = sorted(set(glob.glob(os.path.join(data_root, pattern))))
    # data_dirs = load_product_dirs_from_txt(data_root, "annotated_products.txt")
    # data_dirs = load_product_dirs_from_txt(data_root, "train_80.txt")

    data_dirs = [d for d in data_dirs if os.path.isdir(d)]

    print(f"Found {len(data_dirs)} product folders")

    if len(data_dirs) > 1:
        worker_results = Parallel(process_one_product, data_dirs, 
                                  list(spacings), n_cores=6)
    else:
        # test line
        worker_results = [process_one_product(list(spacings), path_) for path_ in data_dirs]

    final_json, failures = reduce_results(worker_results, spacings)

    # register acquisition month histogram in final JSON and save bar plot
    month_counts = save_month_histogram(worker_results, out_json_path)
    final_json["acquisition_month_histogram"] = {
        "Jan": int(month_counts[0]),
        "Feb": int(month_counts[1]),
        "Mar": int(month_counts[2]),
        "Apr": int(month_counts[3]),
        "May": int(month_counts[4]),
        "Jun": int(month_counts[5]),
        "Jul": int(month_counts[6]),
        "Aug": int(month_counts[7]),
        "Sep": int(month_counts[8]),
        "Oct": int(month_counts[9]),
        "Nov": int(month_counts[10]),
        "Dec": int(month_counts[11]),
    }
    
    with open(out_json_path, "w", encoding="utf-8") as f:
        json.dump(final_json, f, indent=2)

    print(f"Saved stats JSON to: {out_json_path}")
    print(f"Successful products: {sum(r['ok'] for r in worker_results)}")
    print(f"Failed products: {len(failures)}")

    return final_json, failures


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute SAR HH/HV statistics at multiple pixel spacings"
    )

    parser.add_argument(
        "--data_root",
        type=str,
        default="D:/Temp/CIS - Images - ArcticScope/03_dB_noise_sub_OK",
        help="Root directory containing SAR product folders"
    )

    parser.add_argument(
        "--out_json",
        type=str,
        default="./data_info/rcm_stats_by_pixel_spacing.json",
        help="Output JSON file path"
    )

    parser.add_argument(
        "--spacings",
        type=int,
        nargs="+",
        default=[50, 100, 200],    #   default=[200],  #   
        help="List of pixel spacings (meters), e.g. --spacings 50 100 200"
    )

    parser.add_argument(
        "--pattern",
        type=str,
        default="*",
        help="Glob pattern to filter subfolders"
    )

    return parser.parse_args()

# ============================================================
# Example use
# ============================================================

if __name__ == "__main__":
    args = parse_args()

    stats_json, failures = compute_dataset_stats_json(
        data_root=args.data_root,
        out_json_path=args.out_json,
        spacings=tuple(args.spacings),
        pattern=args.pattern,
    )

    if failures:
        fail_log = args.out_json.replace(".json", "_failures.json")
        with open(fail_log, "w", encoding="utf-8") as f:
            json.dump(failures, f, indent=2)
        print(f"Saved failure log to: {fail_log}")