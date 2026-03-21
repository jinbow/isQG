# Python → Julia port (IsQG)

This repository now contains a Julia package port of the Python implementation in `src/`.

## Layout

- Python (reference): `src/isqg.py`, `src/invert.py`, `src/test_matrix.py`
- Julia package: `julia/IsQG/`
  - `julia/IsQG/src/types.jl`: `ISQGData` type (Python `isqg_data`)
  - `julia/IsQG/src/core.jl`: `precondition!`, `solve_sqg!`, `solve_psii!`, `solve!`
  - `julia/IsQG/src/io.jl`: `invertit(filein, fileout)` (`.mat` IO via `MAT.jl`)
  - `julia/IsQG/src/matrices.jl`: `construct_matrices` (port of `src/test_matrix.py`)
  - `julia/IsQG/src/utils.jl`: `lonlat2xy`, `lat2f`, `gradxy`, `psi2uv`, `vorticity`, `twopave`, `N2rho`, `anomaly`, fitting helpers

## Julia dependencies

The Julia package uses:

- `LinearAlgebra`, `SparseArrays` (stdlib)
- `FFTW` (FFT acceleration; optional via `IsQGFFTWExt`, otherwise a slow pure-Julia fallback is used)
- `MAT` (MATLAB `.mat` IO; optional via `IsQGMATExt`)
- `Plots` (diagnostic plots; optional via `IsQGPlotsExt`)

## Run tests

From the repo root:

```sh
cd julia/IsQG
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

## Generate/update Python reference data

The Julia tests can compare against a small embedded reference generated from the Python implementation:

```sh
python julia/IsQG/tools/gen_reference_jl.py
```

This writes `julia/IsQG/test/reference/reference_data.jl`.

## Minimal example

```sh
cd julia/IsQG
julia --project=. examples/minimal.jl
```

If `Plots.jl` is available, the example writes horizontal + vertical diagnostic plots to `julia/IsQG/examples/output/`.

## MATLAB IO / invertit

```sh
cd julia/IsQG
julia --project=. examples/invertit_cli.jl input.mat output.mat
```
