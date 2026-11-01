"""Build the meta-atom library.

Usage: python scripts/run_library.py configs/library_main.yaml [--max-batches N]

With --max-batches the process exits after N new batches (exit code 3 if work remains),
so a shell loop can restart it with a fresh JAX compilation cache:
    until python scripts/run_library.py configs/library_main.yaml --max-batches 40; do :; done
"""

import sys
import time

from metatrap.rcwa.library import LibrarySpec, is_complete, run_sweep
from metatrap.runio import REPO_ROOT, load_config, run_dir, write_params


def main(cfg_path, max_batches=None):
    cfg = load_config(cfg_path)
    out = run_dir(cfg["run_id"])
    write_params(out, cfg)
    spec = LibrarySpec.from_config(cfg)
    h5 = REPO_ROOT / cfg.get("output", f"data/{cfg['run_id']}.h5")
    h5.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    run_sweep(spec, h5, log=lambda s: print(s, flush=True), max_batches=max_batches)
    if not is_complete(h5):
        print(f"partial ({(time.time()-t0)/60:.1f} min), work remains", flush=True)
        sys.exit(3)
    print(f"done in {(time.time()-t0)/60:.1f} min -> {h5}")


if __name__ == "__main__":
    mb = int(sys.argv[sys.argv.index("--max-batches") + 1]) if "--max-batches" in sys.argv else None
    main(sys.argv[1], mb)
