# HWG Pipeline

`hwg-pipeline` is a SageMath-only scientific codebase. Essential code is run
with `./scripts/sage-python`, ensuring that Sage's exact arithmetic and
representation-theory facilities are available even when `sage` is not on the
noninteractive shell's `PATH`.

The implemented stages expand structured HWGs, restore irreducible Sage
characters, and compute exact dimension-refined and unrefined Hilbert series.
Later stages such as Adams operations, plethystic logarithms, branching, charge
maps, and monopole-formula calculations remain unimplemented.

## Environment

Create and activate the conda-forge environment:

```console
conda env create -f environment.yml
conda activate hwg-pipeline
```

Install the package in editable mode and run its checks:

```console
./scripts/sage-python -m pip install -e .
make test
make environment-check
```

The launcher checks `HWG_SAGE_EXECUTABLE` first, then `sage` on `PATH`, and
finally `$HOME/miniforge3/envs/hwg-sage/bin/sage`. For example:

```console
HWG_SAGE_EXECUTABLE=/opt/sage/bin/sage ./scripts/sage --version
./scripts/sage-python -m hwg_pipeline characters su3_nf5_k3o2_infinite --order 10
```

## Repository layout

- `src/hwg_pipeline/`: generic exact-arithmetic pipeline implementation
- `scripts/`: noninteractive Sage and Sage-Python launchers
- `tests/`: environment and mathematical tests
- `theories/`: theory-specific data
- `references/overleaf/`: source LaTeX
- `generated/`: generated outputs
- `docs/`: documentation
- `notebooks/`: non-essential exploratory notebooks

See [`SPEC.md`](SPEC.md) for the planned pipeline and [`AGENTS.md`](AGENTS.md)
for project-specific contribution instructions.
