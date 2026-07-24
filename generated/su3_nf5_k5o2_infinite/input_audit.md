# Input audit: `su3_nf5_k5o2_infinite`

## Source and metadata

- Source: `references/overleaf/su3_5f_6f_hwg_results.tex`; paper equation `10.3`.
- Theory: SU(3)+5F, exact $|k|=5/2$, infinite coupling.
- Enhanced symmetry: $SO(10)\times U(1)_q$; D5 coroot Dynkin convention.
- Source SHA-256: `2e79f9e66dfb0a2a7d2086f642aa94b45b5284287ce41ade80a128dc0439b359`.
- Fixture SHA-256: `6d82c13ca6ebe4b89b1d44f19363a0e6d7d2c27292f3caec8505b3769ccb66ef`.

## Original PE expression
```tex
\PE\!\left[(\mu_2+1)t^2+(q\mu_4+q^{-1}\mu_5)t^3
+\mu_4\mu_5t^4-\mu_4\mu_5t^6\right]
```

Normalized: `\operatorname{PE}\!\left[mu_2 t^{2} + t^{2} + mu_4 q t^{3} + mu_5 q^{-1} t^{3} + mu_4 mu_5 t^{4} - mu_4 mu_5 t^{6}\right]`

| # | coefficient | degree | D5 labels | q charge |
|---:|---:|---:|---|---:|
| 1 | +1 | 2 | `[0, 0, 0, 0, 0]` | +0 |
| 2 | +1 | 2 | `[0, 1, 0, 0, 0]` | +0 |
| 3 | +1 | 3 | `[0, 0, 0, 1, 0]` | +1 |
| 4 | +1 | 3 | `[0, 0, 0, 0, 1]` | -1 |
| 5 | +1 | 4 | `[0, 0, 0, 1, 1]` | +0 |
| 6 | -1 | 6 | `[0, 0, 0, 1, 1]` | +0 |

## Original rational product
```tex
\frac{1-\mu_4\mu_5t^6}
{(1-t^2)(1-\mu_2t^2)(1-q\mu_4t^3)(1-q^{-1}\mu_5t^3)(1-\mu_4\mu_5t^4)}
```

Normalized: `\left(1 - mu_2 t^{2}\right)^{-1} \left(1 - mu_4 mu_5 t^{4}\right)^{-1} \left(1 - mu_4 q t^{3}\right)^{-1} \left(1 - mu_5 q^{-1} t^{3}\right)^{-1} \left(1 - t^{2}\right)^{-1} \left(1 - mu_4 mu_5 t^{6}\right)`

| # | power | location | degree | D5 labels | q charge |
|---:|---:|---|---:|---|---:|
| 1 | +1 | numerator | 6 | `[0, 0, 0, 1, 1]` | +0 |
| 2 | -1 | denominator | 2 | `[0, 0, 0, 0, 0]` | +0 |
| 3 | -1 | denominator | 2 | `[0, 1, 0, 0, 0]` | +0 |
| 4 | -1 | denominator | 3 | `[0, 0, 0, 1, 0]` | +1 |
| 5 | -1 | denominator | 3 | `[0, 0, 0, 0, 1]` | -1 |
| 6 | -1 | denominator | 4 | `[0, 0, 0, 1, 1]` | +0 |

## Validation

- **PASS:** fixture matches source lines 127–133 term by term, including signs and charges.
- **PASS:** exact rational arithmetic; no YAML floats.
- **PASS:** PE coefficients are negatives of rational-product powers for each monomial.
- **PASS:** D5 convention is `WeylCharacterRing("D5", style="coroots")`.
- Unresolved ambiguities: none.
