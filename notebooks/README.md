# HWG pipeline notebooks

## Project walkthrough

[`hwg_pipeline_project_walkthrough.ipynb`](hwg_pipeline_project_walkthrough.ipynb) is the deterministic, stored-results audit notebook for `su3_nf5_k3o2_infinite` through order 10 and manifest SU(5) branching.

Generate it from the repository root with:

```bash
./scripts/sage-python -m hwg_pipeline project-notebook \
  su3_nf5_k3o2_infinite --order 10 \
  --branching su3_nf5_k3o2_to_manifest \
  --through branching --strict
```

The notebook embeds all completed mathematical results in Markdown. Its code cells only load stored files and run small exact-arithmetic or Sage convention demonstrations; they never invoke the physical pipeline.
