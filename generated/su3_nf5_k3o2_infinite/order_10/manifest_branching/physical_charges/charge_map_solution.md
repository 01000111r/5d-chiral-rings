# Exact charge map: `su3_nf5_k3o2_manifest_physical`

## Defining equations

- `(x,q)=(6,0) -> (B,I)=(0,1)`.
- `(x,q)=(2,1) -> (B,I)=(-3,0)`.

## Solution

`A = [[0,-3],[1/6,-1/3]]`, determinant `1/2`, rank 2.

`B = -3q`; `I = (x-2q)/6`.

Inverse: `q=-B/3`; `x=6I-2B/3`.

All defining residuals are `(0,0)` and both redundant conjugate validation anchors pass.

All 769 transformed character/PL sectors have integral physical charges.

| degree | SU(5) labels | x | q | B | I | multiplicity | classification |
|---:|---|---:|---:|---:|---:|---:|---|
| 2 | [0, 0, 0, 0] | {'numerator': 0, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 1, 'denominator': 1} | neutral_singlet_candidate |
| 2 | [0, 0, 0, 0] | {'numerator': 0, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 1, 'denominator': 1} | neutral_singlet_candidate |
| 2 | [0, 0, 0, 1] | {'numerator': -6, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': -1, 'denominator': 1} | {'numerator': 1, 'denominator': 1} | anti_instanton_generator_candidate |
| 2 | [1, 0, 0, 0] | {'numerator': 6, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 1, 'denominator': 1} | {'numerator': 1, 'denominator': 1} | instanton_generator_candidate |
| 2 | [1, 0, 0, 1] | {'numerator': 0, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 1, 'denominator': 1} | meson_adjoint_candidate |
| 3 | [0, 0, 0, 1] | {'numerator': 4, 'denominator': 1} | {'numerator': -1, 'denominator': 1} | {'numerator': 3, 'denominator': 1} | {'numerator': 1, 'denominator': 1} | {'numerator': 1, 'denominator': 1} | mixed_baryon_instanton_generator_candidate |
| 3 | [0, 0, 1, 0] | {'numerator': -2, 'denominator': 1} | {'numerator': -1, 'denominator': 1} | {'numerator': 3, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 1, 'denominator': 1} | baryon_generator_candidate |
| 3 | [0, 1, 0, 0] | {'numerator': 2, 'denominator': 1} | {'numerator': 1, 'denominator': 1} | {'numerator': -3, 'denominator': 1} | {'numerator': 0, 'denominator': 1} | {'numerator': 1, 'denominator': 1} | antibaryon_generator_candidate |
| 3 | [1, 0, 0, 0] | {'numerator': -4, 'denominator': 1} | {'numerator': 1, 'denominator': 1} | {'numerator': -3, 'denominator': 1} | {'numerator': -1, 'denominator': 1} | {'numerator': 1, 'denominator': 1} | mixed_baryon_instanton_generator_candidate |

## Audit status

- **Verified:** the input branching data and exact preservation checks.
- **Manually supplied:** the physical charge anchors and their convention.
- **Computationally derived:** the exact rational matrix and inverse.
- **Conservative:** candidate operator names use only representation and solved-charge rules.
- The charge map is convention-dependent; reversing both instanton and baryon orientations gives an equivalent alternative convention.
- The program did not infer physical charge meanings without anchors.
- The two neutral singlets have not yet been microscopically distinguished.
- Mixed-charge generators have not yet been assigned explicit composite formulas.
- No explicit polynomial relations have been constructed.
