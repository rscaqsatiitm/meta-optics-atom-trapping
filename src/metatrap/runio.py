"""Run bookkeeping: config loading and provenance records for every result."""

from __future__ import annotations

import datetime as _dt
import json
import platform
import subprocess
from importlib import metadata
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
_PACKAGES = ("numpy", "scipy", "jax", "fmmax", "grcwa", "h5py", "matplotlib", "metatrap")


def load_config(path: str | Path) -> dict[str, Any]:
    with open(path) as f:
        cfg = yaml.safe_load(f)
    cfg.setdefault("run_id", Path(path).stem)
    cfg["_config_path"] = str(path)
    return cfg


def git_state() -> dict[str, str]:
    def _git(*args: str) -> str:
        try:
            return subprocess.check_output(["git", *args], cwd=REPO_ROOT, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            return "unknown"

    return {"commit": _git("rev-parse", "HEAD"), "dirty": str(bool(_git("status", "--porcelain")))}


def versions() -> dict[str, str]:
    out = {"python": platform.python_version()}
    for p in _PACKAGES:
        try:
            out[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            out[p] = "not installed"
    return out


def run_dir(run_id: str, root: str | Path | None = None) -> Path:
    d = Path(root) if root else REPO_ROOT / "results" / run_id
    d.mkdir(parents=True, exist_ok=True)
    return d


def _jsonable(x: Any) -> Any:
    import numpy as np

    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, complex):
        return {"re": x.real, "im": x.imag}
    return x


def write_json(path: str | Path, obj: Any) -> None:
    with open(path, "w") as f:
        json.dump(_jsonable(obj), f, indent=2, sort_keys=False)


def write_params(directory: Path, cfg: dict[str, Any], extra: dict[str, Any] | None = None) -> None:
    rec = {
        "config": cfg,
        "git": git_state(),
        "versions": versions(),
        "host": {"machine": platform.machine(), "system": platform.system()},
        "utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
    }
    if extra:
        rec.update(extra)
    write_json(directory / "params.json", rec)
