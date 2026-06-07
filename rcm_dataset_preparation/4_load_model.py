"""
No@
Jun 7, 2026
Script to load an mmsegmentation model from a config and checkpoint, 
and run inference.
Requirements: mmsegmentation, mmengine, torch
arguments:
    --config: Path to the mmseg .py config file.
    --checkpoint: Path to the .pth checkpoint file.
    --device: Device to load model on (default: cpu)
"""
import argparse
import functools

import numpy as np
import torch

def Normalize_mean_std(
    img: np.ndarray,
    mean=(-16.849971303873264, -27.80859960890284),
    std=(6.601613867039019, 6.390188051517371),
    eps: float = 1e-12,
):
    """
    Standardize image using external mean/std.


    - NaNs are replaced by mean (=> become 0 after normalization)
    - Fully vectorized
    """
    assert img.ndim in (2, 3), "img must be (H,W) or (H,W,C)"
    img = img.astype(np.float32, copy=True)
    squeeze = (img.ndim == 2)
    img_ = img[..., None] if squeeze else img
    H, W, C = img_.shape
    mean = np.asarray(mean, dtype=np.float32)
    std = np.asarray(std, dtype=np.float32)
    assert mean.shape == (C,), f"mean must have shape ({C},)"
    assert std.shape == (C,), f"std must have shape ({C},)"
    assert np.all(std > 0), "std must be > 0"
    # Replace NaNs with mean (so normalized value becomes 0)
    img_ = np.where(
        np.isnan(img_),
        mean[None, None, :],
        img_
    )
    # Normalize
    img_norm = (img_ - mean[None, None, :]) / (std[None, None, :] + eps)
    return img_norm[..., 0] if squeeze else img_norm

def load_model(config_path, checkpoint_path, device="cpu"):
    """
    Load an mmsegmentation model from a config and MMEngine checkpoint.

    Args:
        config_path (str): Path to the mmseg .py config file.
        checkpoint_path (str): Path to the .pth checkpoint file.
        device (str): Device string, e.g. "cpu" or "cuda:0".

    Returns:
        model: Loaded mmseg model in eval mode.
    """
    from mmseg.apis import init_model

    # PyTorch 2.6 changed weights_only default to True, but MMEngine checkpoints
    # embed HistoryBuffer and numpy objects not in the safe-globals allowlist.
    _original_load = torch.load
    torch.load = functools.partial(_original_load, weights_only=False)
    try:
        model = init_model(config_path, checkpoint_path, device=device)
    finally:
        torch.load = _original_load

    return model

def get_input_divisor(model):
    """
    Derive the minimum spatial divisor required by the loaded backbone.

    Covers:
      - UNet-style (num_stages + strides + downsamples attributes)
      - ResNet/generic (strides attribute only)
      - Unknown: falls back to 32
    """
    backbone = model.backbone
    # UNet: mirrors _check_input_divisible logic
    if hasattr(backbone, 'num_stages') and hasattr(backbone, 'strides') and hasattr(backbone, 'downsamples'):
        rate = 1
        for i in range(1, backbone.num_stages):
            if backbone.strides[i] == 2 or backbone.downsamples[i - 1]:
                rate *= 2
        return rate
    # ResNet / generic backbones that expose per-stage strides
    if hasattr(backbone, 'strides'):
        return int(np.prod(backbone.strides))
    return 32

def inference(model, img, device):
    """
    Run inference on a normalized (H, W, C) image.

    Args:
        model: Loaded mmseg model (EncoderDecoder subclass).
        img (np.ndarray): Normalized float32 array of shape (H, W, C).
        device: Torch device to run inference on.

    Returns:
        numpy.ndarray: Predicted segmentation map of shape (H, W).
    """
    height, width = img.shape[:2]

    # (H, W, C) -> (1, C, H, W)
    tensor = torch.from_numpy(img).permute(2, 0, 1).unsqueeze(0).to(device)

    # Pad to meet backbone input divisor requirements (e.g. 32 for ResNet, 64 for UNet)
    divisor = get_input_divisor(model)
    pad_h = (divisor - height % divisor) % divisor
    pad_w = (divisor - width % divisor) % divisor
    if pad_h > 0 or pad_w > 0:
        tensor = torch.nn.functional.pad(tensor, (0, pad_w, 0, pad_h), mode="reflect")

    with torch.no_grad():
        # img_shape tells predict_by_feat to resize logits back to original size
        seg_logits = model.encode_decode(tensor, batch_img_metas=[{"img_shape": (height, width)}])

    seg_map = seg_logits.argmax(dim=1).squeeze().cpu().numpy()
    return seg_map


def main():
    parser = argparse.ArgumentParser(description="Load and verify an mmsegmentation checkpoint")
    parser.add_argument("--config", default="D:/Temp/Model_update_code/sea-ice-mmseg/configs/arcticscope/unet_1xb8-amp-coslr-30ki_rcm.py", type=str, help="Path to the mmseg config .py file")
    parser.add_argument("--checkpoint", default="D:/Temp/Model_update_code/sea-ice-mmseg/work_dirs/unet_1xb8-amp-coslr-30ki_rcm/best_mFscore_iter_3500.pth", type=str, help="Path to the .pth checkpoint file")
    parser.add_argument("--device", type=str, default="cpu", help="Device to load model on (default: cpu)")
    args = parser.parse_args()

    print(f"Config:     {args.config}")
    print(f"Checkpoint: {args.checkpoint}")
    print(f"Device:     {args.device}")

    img = np.random.rand(500, 500, 2).astype(np.float32)
    img = Normalize_mean_std(img)

    model = load_model(args.config, args.checkpoint, device=args.device)
    device = next(model.parameters()).device
    print(f"Model loaded: {type(model).__name__}  |  Parameters: {sum(p.numel() for p in model.parameters()):,}")

    print(f"Input divisor:  {get_input_divisor(model)}")
    print("Running inference on 500x500x2 toy image...")
    seg_map = inference(model, img, device)
    print(f"Output shape:   {seg_map.shape}")
    print(f"Unique classes: {np.unique(seg_map)}")


if __name__ == "__main__":
    main()
