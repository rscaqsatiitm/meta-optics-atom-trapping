"""Build the meta-atom library. Usage: python scripts/run_library.py configs/library_1075.yaml"""

import sys
import time

from metatrap.rcwa.library import LibrarySpec, run_sweep
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_params


def main(cfg_path):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    spec = LibrarySpec.from_config(cfg)
    h5 = REPO_ROOT / cfg.get("output", f"data/{cfg['run_id']}.h5")
    h5.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    run_sweep(spec, h5, log=lambda s: print(s, flush=True))
    print(f"done in {(time.time()-t0)/60:.1f} min -> {h5}")


if __name__ == "__main__":
    main(sys.argv[1])
