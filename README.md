# RCM-sea-ice-mmseg

This repository is used to train and update the segmentation model behind the active learning pipeline of the [ArcticScope](https://github.com/jnoat92/ArcticScope.git) project. The model currently classifies ice and water on SAR (RCM) images, but the pipeline is designed to be extended to other tasks, model architectures, and sensors.

## Structure

- [`rcm_dataset_preparation/`](rcm_dataset_preparation/) — Scripts to prepare the training dataset once images and annotations have been collected. See its [README](rcm_dataset_preparation/README.md) for usage.
- [`sea-ice-mmseg/`](sea-ice-mmseg/) — An adaptation of [mmsegmentation](https://github.com/open-mmlab/mmsegmentation) containing the complete pipeline to train the models. See its [README](sea-ice-mmseg/README.md) for usage.
- [`windows_create_env.sh`](windows_create_env.sh) — Script to set up the environment on Windows.
