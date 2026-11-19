"""Optional AI second opinion: GelGenie's band-segmentation network (U-Net), run with onnxruntime.

GelGenie (Aquilina et al., Nature Communications 2025; code Apache-2.0, models Apache-2.0 on
HuggingFace) segments band pixels. autogel uses it only on request, as a second opinion:
* bands autogel detected but the network sees no band in  -> flagged (``ai_support`` = False)
* band regions the network sees inside a lane where autogel has no band -> offered as new bands

Install the extra:  pip install autogel[ai]   (onnxruntime, ~17 MB). The model file (57 MB) is
downloaded once, on request, into ~/.autogel/models and checked against its SHA-256.
Runs on the CPU; large images are scaled to at most MAX_SIDE pixels on the long side.
"""
from __future__ import annotations

import hashlib
import os

import numpy as np
from scipy import ndimage

MODEL = {
    "name": "GelGenie-Universal-Dec-2023",
    "url": "https://huggingface.co/mattaq/GelGenie-Universal-Dec-2023/resolve/main/onnx_checkpoints/"
           "unet_dec_21_epoch_579.onnx",
    "file": "gelgenie_universal_dec2023.onnx",
    "sha256": "10ea972a40625f03b5a2987e1a98b532c6d4ded1daae32f1eb91aa24db8d9c48",
    "size_mb": 57.3,
    "licence": "Apache-2.0",
    "source": "https://huggingface.co/mattaq/GelGenie-Universal-Dec-2023",
    "trained_on": "GelGenie training split (Zenodo 10.5281/zenodo.14641949); not the test images or external_gels",
}
MAX_SIDE = 1536
_session = None


def model_path() -> str:
    base = os.environ.get("AUTOGEL_HOME") or os.path.join(os.path.expanduser("~"), ".autogel")
    return os.path.join(base, "models", MODEL["file"])


def available() -> tuple[bool, str]:
    """(usable, reason) - onnxruntime installed and model downloaded?"""
    try:
        import onnxruntime  # noqa: F401
    except ImportError:
        return False, "the AI extra is not installed (pip install autogel[ai])"
    if not os.path.exists(model_path()):
        return False, f"the model ({MODEL['size_mb']:.0f} MB) is not downloaded yet"
    return True, "ok"


def download_model(progress=None) -> str:
    """Download the model once (only when the user asks); verified by SHA-256."""
    import ssl
    import urllib.request
    try:
        import certifi
        ctx = ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        ctx = ssl.create_default_context()
    dst = model_path()
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    h = hashlib.sha256()
    req = urllib.request.Request(MODEL["url"], headers={"User-Agent": "autogel"})
    with urllib.request.urlopen(req, context=ctx, timeout=600) as r, open(dst + ".part", "wb") as f:
        total = int(r.headers.get("Content-Length") or 0)
        done = 0
        while chunk := r.read(1 << 20):
            f.write(chunk)
            h.update(chunk)
            done += len(chunk)
            if progress:
                progress(done, total)
    if h.hexdigest() != MODEL["sha256"]:
        os.remove(dst + ".part")
        raise ValueError("downloaded model failed the checksum test - not used")
    os.replace(dst + ".part", dst)
    return dst


def _get_session():
    global _session
    if _session is None:
        import onnxruntime as ort
        opts = ort.SessionOptions()
        opts.intra_op_num_threads = max(1, (os.cpu_count() or 2) - 1)
        _session = ort.InferenceSession(model_path(), opts, providers=["CPUExecutionProvider"])
    return _session


def segment(raw: np.ndarray) -> np.ndarray:
    """Band-pixel probability map (0..1) of an original (un-oriented) grey image, same size.

    Preprocessing as GelGenie: grey, scaled by the bit depth to 0..1, padded to a multiple of 32."""
    a = np.asarray(raw)
    if a.ndim == 3:
        a = a[..., :3].mean(-1)
    if np.issubdtype(a.dtype, np.integer):
        scale = 255.0 if a.max() <= 255 else 65535.0
    else:
        scale = float(max(a.max(), 1e-9))
    x = a.astype(np.float32) / scale
    H, W = x.shape
    f = max(1.0, max(H, W) / MAX_SIDE)
    if f > 1:
        x = ndimage.zoom(x, 1 / f, order=1)
    h, w = x.shape
    ph, pw = 32 * (h // 32 + 1) - h, 32 * (w // 32 + 1) - w
    x = np.pad(x, ((ph // 2, ph - ph // 2), (pw // 2, pw - pw // 2)))
    sess = _get_session()
    out = sess.run(None, {sess.get_inputs()[0].name: x[None, None]})[0][0]       # (2, h, w) logits
    e = np.exp(out - out.max(0, keepdims=True))
    prob = (e[1] / e.sum(0))[ph // 2: ph // 2 + h, pw // 2: pw // 2 + w]
    if f > 1:
        prob = ndimage.zoom(prob, (H / prob.shape[0], W / prob.shape[1]), order=1)
    return prob.astype(np.float32)


def oriented_mask(ga, threshold: float = 0.5) -> np.ndarray:
    """The network's band mask in autogel's oriented frame (cached on the analysis)."""
    cache = getattr(ga, "_ai_mask", None)
    if cache is not None and cache[0] == threshold:
        return cache[1]
    from .imageio import apply_orientation, load_raw, reduce
    raw, _, _ = load_raw(ga.path)
    raw = reduce(raw, ga.gel.reduced)
    prob = segment(raw)
    o = {**ga.gel.orientation, "invert": False}
    m = apply_orientation((prob >= threshold).astype(np.float32), o, order=0) > 0.5
    ga._ai_mask = (threshold, m)
    return m


def review(ga, threshold: float = 0.5, min_px: int = 15) -> dict:
    """Compare autogel's bands with the network's band regions.

    Returns {"support": [bool per band in ga.bands], "missing": [Band, ...]}: a band is supported
    when band pixels of the network fall inside it; "missing" are network regions inside a lane
    that overlap no autogel band (offered as new bands, source "ai")."""
    m = oriented_mask(ga, threshold)
    lab, n = ndimage.label(m)
    support = []
    for b in ga.bands:
        lane = ga.lanes[b.lane]
        d = lane.dx(b.y)
        y0, y1 = int(max(0, b.y0)), int(min(m.shape[0], b.y1 + 1))
        x0, x1 = int(max(0, lane.x0 + d)), int(min(m.shape[1], lane.x1 + d))
        support.append(bool(m[y0:y1, x0:x1].sum() >= min_px))
    missing = []
    for k, sl in enumerate(ndimage.find_objects(lab), start=1):
        if sl is None or (lab[sl] == k).sum() < min_px:
            continue
        cy, cx = (sl[0].start + sl[0].stop) / 2, (sl[1].start + sl[1].stop) / 2
        lanes = [i for i, l in enumerate(ga.lanes) if l.contains(cx, cy) and l.y0 + 10 < cy < l.y1]
        if not lanes:
            continue
        i = lanes[0]
        if any(bb.lane == i and bb.y0 - 2 <= cy <= bb.y1 + 2 for bb in ga.bands):
            continue
        missing.append(ga.band_for(i, float(sl[0].start), float(sl[0].stop - 1), source="ai"))
    return {"support": support, "missing": missing}
