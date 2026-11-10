# Running the 10^6-site design on a desktop GPU

This note covers the full-aperture design and simulation of the 1000 x 1000-site direct-projection hologram (`configs/fullaperture_1M.yaml`) on a Windows desktop with an NVIDIA GPU, for example an i5-13500 with 16 GB RAM and an RTX 3050 (6 GB). The code was written and tested on a CPU-only Linux container. The GPU kernels were run there under the numba CUDA simulator and compiled to PTX, but they have **not yet been run on a GPU**. Steps 3 and 4 below are the first real-GPU checks.

## 1. What runs where

| Stage | Cost for 10^6 sites (DERIVED) | Where |
|---|---|---|
| WGS design, 60 + 2 x 30 iterations, two kernel products each | 5.2e12 kernel evaluations per product | GPU (`--gpu`) |
| Render of the 18.4 mm metasurface, 1.07e9 pixels inside the aperture | windowed-sinc lookup, 256 taps per pixel; 1.1 M pixel/s on a shared 4-vCPU container (about 16 min); expected faster on 14 cores | CPU (numba) |
| Folded propagation, 18432^2 (9.2 mm window) | complex64, folded and propagated in place, about 2 copies of 2.7 GB | CPU |
| Tiled focal metrics, 10^6 sites | 1024^2 tiles upsampled x4, per-site loop | CPU |

At the 100k design's measured CPU rate (6.9e8 evaluations/s on 4 vCPU), one 10^6 product would take about 2 h. With about 240 products, the CPU-only design would need several weeks on 4 vCPU and roughly a week on 14 cores. The GPU is therefore needed for the design step only.

Expected GPU rate (DERIVED from the instruction mix, not measured): about 3 special-function operations (sqrt, sin, cos) and about 25 single-precision operations per evaluation. The upper bound is about 1e11 evaluations/s, and 3e10 to 6e10 is likely in practice. That gives 1.5 to 3 min per product and 6 to 12 h for the WGS. `scripts/bench_gpu.py` measures the actual rate (step 4).

## 2. Setup (Windows, WSL2)

WSL2 with Ubuntu is the simplest route. The Windows NVIDIA driver provides CUDA inside WSL2, so no Linux driver is installed.

1. In PowerShell (administrator): `wsl --install -d Ubuntu`, then reboot.
2. WSL2 gets half of the host RAM by default (8 GB here). The 10^6 run peaks at a few GB, so this is enough, but 12 GB leaves margin. Create `C:\Users\<you>\.wslconfig` containing:
   ```
   [wsl2]
   memory=12GB
   ```
   Then run `wsl --shutdown` and reopen Ubuntu.
3. In Ubuntu:
   ```bash
   sudo apt update && sudo apt install -y python3-venv git
   git clone https://github.com/rscaqsatiitm/meta-optics-atom-trapping.git
   cd meta-optics-atom-trapping && git checkout claude/pensive-euler-h93h6y
   python3 -m venv .venv && . .venv/bin/activate
   pip install -e ".[test,gpu]"
   ```
   The `gpu` extra installs `numba-cuda[cu13]` with the CUDA 13 runtime libraries from PyPI. It also pins `numpy<2.4`, because numba-cuda 0.30.4 still imports `np.row_stack`. Use `cu12` instead if the driver is older than CUDA 13.
4. Check that the GPU is visible: `nvidia-smi`, then `python -c "from numba import cuda; cuda.detect()"`.

Keep the clone inside the Linux file system (`~/...`), not under `/mnt/c`. File access across the boundary is slow.

## 3. Checks before the long run

```bash
python -m pytest -q
python scripts/check_gpu_ptx.py
```

The first command runs the test suite. On a GPU machine `tests/test_lazy_gpu.py` now runs on the real device instead of being skipped. The second compiles the kernels and confirms that the per-pair arithmetic is single precision. Float64 work would run at 1/64 rate on this GPU. An implicit int64-to-float64 promotion of exactly this kind was found and removed while preparing this note.

## 4. Benchmark and validation

```bash
python scripts/bench_gpu.py          # 100k geometry: accuracy vs CPU, evaluations/s, 10^6 extrapolation
python scripts/fullaperture_lazy.py configs/fullaperture_100k_gpu_check.yaml --gpu --tiled
```

The second command re-evaluates the finished 100k design (`data/designs/fullaperture_100k_round2_state.npz`). It renders, propagates and measures the design with ideal and library pillars. Its `results/fullaperture_100k_gpu_check/metrics.json` should reproduce `results/fullaperture_100k_tiled/metrics.json` from the container run (CPU kernels) to about 1e-3, the single-precision tolerance. This is the end-to-end GPU check. It takes roughly 1 h, dominated by the two CPU renders.

## 5. The 10^6-site run

```bash
nohup python -u scripts/fullaperture_lazy.py configs/fullaperture_1M.yaml --gpu --tiled > fa1M.log 2>&1 &
tail -f fa1M.log
```

The run checkpoints the WGS every 10 iterations and saves the design state before and after each feedback round (`results/fullaperture_1M/`). After an interruption (reboot, Windows update, closed terminal), the same command resumes from the last checkpoint.

To keep Windows from sleeping, set the power plan to never sleep while the run is active.

Each product is split into GPU launches of about 2e9 evaluations (`pairs_per_launch`). This keeps every launch well under the 2 s Windows display-driver timeout (TDR), which also applies under WSL2 when the RTX 3050 drives a monitor. If the log shows `CUDA_ERROR_LAUNCH_TIMEOUT`, halve `pairs_per_launch` in `LazyKGGPU`. Do not raise the TDR delay in the registry.

## 6. Outputs to send back

`results/fullaperture_1M/metrics.json` and `params.json` are small and can be committed. The `per_site_*.npz` files are about 45 MB each for 10^6 sites. The `.npz` files are git-ignored, so send them separately if the per-site maps are wanted.
