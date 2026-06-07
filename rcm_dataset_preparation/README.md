# RCM-sea-ice-mmseg
This repository includes the different steps to follow towards updating an ice/water classification model on RCM images containing pixel level annotations

## Dataset preparation

Before running anything, make sure that:
- The SAR product folders (`--data_root`) and their corresponding annotation folders (`--labels_root`) are aligned, i.e. each annotated product has a matching SAR product folder with the same name.
- [`annotated_products.txt`](annotated_products.txt) is up to date and lists exactly the product folders that have annotations.

Once that is verified, run the following scripts **in sequence**:

### 1. `0_compute_data_stats.py`
Computes SAR HH/HV statistics at multiple pixel spacings over the product folders.

```bash
python 0_compute_data_stats.py ^
    --data_root "D:/Temp/CIS - Images - ArcticScope/03_dB_noise_sub_OK" ^
    --out_json ./data_info/rcm_stats_by_pixel_spacing.json ^
    --spacings 50 100 200
```

### 2. `1_class_stats_data_split.py`
Creates the train/validation split for the annotated products listed in `annotated_products.txt`.

```bash
python 1_class_stats_data_split.py ^
    --data_root "D:/Temp/CIS - Images - ArcticScope/03_dB_noise_sub_OK/" ^
    --labels_root "D:/Temp/ArcticScope-executable_200m/Custom_Annotation/" ^
    --annotated_products_txt annotated_products.txt ^
    --percentage 80 ^
    --validate_with_sar
```

### 3. `2_patch_dataset_creator.py`
Extracts the image/label patches that will be used to train the model.

```bash
python 2_patch_dataset_creator.py ^
    --data_root "D:/Temp/CIS - Images - ArcticScope/03_dB_noise_sub_OK/" ^
    --labels_root "D:/Temp/ArcticScope-executable_200m/Custom_Annotation/" ^
    --annotated_products_txt annotated_products.txt ^
    --output ./patches/ ^
    --shp_path ./StatCan_ocean/StatCan_ocean.shp ^
    --pixel_spacing 200 ^
    --labels_pixel_spacing 200 ^
    --patch_size 384 ^
    --validate_label_shape
```

> Adjust the paths and parameters above (e.g. `--pixel_spacing`, `--patch_size`, `--overlap`) to match your data and target configuration. Run `python <script>.py --help` for the full list of options.

## Loading a trained model

[`4_load_model.py`](4_load_model.py) loads a trained mmsegmentation model from its config and checkpoint files and runs inference on an input image.

```bash
python 4_load_model.py ^
    --config path/to/config.py ^
    --checkpoint path/to/checkpoint.pth ^
    --device cpu
```
