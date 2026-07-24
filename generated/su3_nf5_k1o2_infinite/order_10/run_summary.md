# Run summary

All seven calculation stages passed. Refined, q-refined dimension, and independent scalar reconstructions passed; the refined difference is empty. The six-stage deterministic rerun was byte-for-byte identical for mathematical JSON and TeX outputs. Sage 10.7; Python 3.11.15.

## Stage elapsed time and peak RSS

- `./scripts/sage-python -m hwg_pipeline expand su3_nf5_k1o2_infinite --order 10`: 2.813 s; peak RSS unavailable.
- `./scripts/sage-python -m hwg_pipeline characters su3_nf5_k1o2_infinite --order 10`: 3.787 s; peak RSS unavailable.
- `./scripts/sage-python -m hwg_pipeline plethystic-log su3_nf5_k1o2_infinite --order 10 --formal-log direct`: 18.191 s; peak RSS unavailable.
- `./scripts/sage-python -m hwg_pipeline reconstruct su3_nf5_k1o2_infinite --order 10`: 46.173 s; peak RSS unavailable.
- `./scripts/sage-python -m hwg_pipeline analyze-pl su3_nf5_k1o2_infinite --order 10`: 3.852 s; peak RSS unavailable.
- `./scripts/sage-python -m hwg_pipeline compact-latex su3_nf5_k1o2_infinite --order 10`: 2.739 s; peak RSS unavailable.

PDF compilation and peak RSS are unavailable because no LaTeX compiler and no `/usr/bin/time` executable are installed. Excluded calculations are not applicable. Completed campaign outputs and the n=3 benchmark are unchanged.
