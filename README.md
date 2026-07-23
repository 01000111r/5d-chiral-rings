# HWG Pipeline

`hwg-pipeline` is the repository scaffold for a SageMath-only scientific
codebase. Essential code is run with `sage -python`, ensuring that Sage's exact
arithmetic and representation-theory facilities are available.

This initial version contains infrastructure and environment checks only. It
does **not** implement HWG expansion, representation conversion, characters,
plethystic functions, branching, charge maps, monopole formula calculations,
or any other stage of the planned mathematical pipeline.

## Environment

Create and activate the conda-forge environment:

```console
conda env create -f environment.yml
conda activate hwg-pipeline
```

Install the package in editable mode and run its checks:

```console
sage -python -m pip install -e .
make test
make environment-check
```

## Repository layout

- `src/hwg_pipeline/`: generic implementation code (currently only the package
  scaffold)
- `tests/`: environment and, in future milestones, mathematical tests
- `theories/`: theory-specific data
- `references/overleaf/`: source LaTeX
- `generated/`: generated outputs
- `docs/`: documentation
- `notebooks/`: non-essential exploratory notebooks

See [`SPEC.md`](SPEC.md) for the planned pipeline and [`AGENTS.md`](AGENTS.md)
for project-specific contribution instructions.
