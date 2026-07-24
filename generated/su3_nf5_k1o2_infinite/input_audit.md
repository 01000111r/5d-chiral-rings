# Input audit: `su3_nf5_k1o2_infinite`

**PASS.** The calculation input was transcribed term-by-term from equation 12.3 and its explicit $N=3$ specialization in `references/overleaf/su3_5f_6f_hwg_results.tex` before expansion. No ambiguity remains.

- Ordered factors: `A4 (SU(5))`, `A1 (SU(2))`; Sage `coroots` Dynkin convention.
- Abelian convention: integer exponent of `q`; conjugation reverses charge.
- Source SHA-256: `2e79f9e66dfb0a2a7d2086f642aa94b45b5284287ce41ade80a128dc0439b359`
- Fixture SHA-256: `486338e780eb512a60816963b6355494ab736f9f8bcff2e91542c9da72d31f7d`

## Original PE
```tex
\PE\!\left[(\mu_1\mu_4+\nu^2+1)t^2
+\nu(q\mu_2+q^{-1}\mu_3)t^3
+\mu_2\mu_3t^4-\nu^2\mu_2\mu_3t^6\right]
```

## Normalized PE
`PE[t^2 + nu^2 t^2 + mu_1 mu_4 t^2 + nu q mu_2 t^3 + nu q^(-1) mu_3 t^3 + mu_2 mu_3 t^4 - nu^2 mu_2 mu_3 t^6]`

| # | coefficient | degree | A4 labels | A1 labels | q |
|---:|---:|---:|---|---|---:|
| 1 | +1 | 2 | `[0, 0, 0, 0]` | `[0]` | 0 |
| 2 | +1 | 2 | `[0, 0, 0, 0]` | `[2]` | 0 |
| 3 | +1 | 2 | `[1, 0, 0, 1]` | `[0]` | 0 |
| 4 | +1 | 3 | `[0, 1, 0, 0]` | `[1]` | 1 |
| 5 | +1 | 3 | `[0, 0, 1, 0]` | `[1]` | -1 |
| 6 | +1 | 4 | `[0, 1, 1, 0]` | `[0]` | 0 |
| 7 | -1 | 6 | `[0, 1, 1, 0]` | `[2]` | 0 |

## Original rational product
```tex
\frac{1-\nu^2\mu_2\mu_3t^6}
{(1-t^2)(1-\nu^2t^2)(1-\mu_1\mu_4t^2)
 (1-\nu q\mu_2t^3)(1-\nu q^{-1}\mu_3t^3)(1-\mu_2\mu_3t^4)}
```

## Normalized rational product
`(1-nu^2 mu_2 mu_3 t^6)/((1-t^2)(1-nu^2 t^2)(1-mu_1 mu_4 t^2)(1-nu q mu_2 t^3)(1-nu q^(-1) mu_3 t^3)(1-mu_2 mu_3 t^4))`

The complete numerator/denominator factor list is stored in `input_audit.json`. Source-versus-fixture comparison: **PASS**. Unresolved ambiguities: **none**.
