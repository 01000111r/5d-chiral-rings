# Input audit: `su3_nf5_k3o2_infinite`

## Source metadata

| Field | Transcribed value |
|---|---|
| Theory | `SU(3)+5F` at infinite coupling |
| Theory ID | `su3_nf5_k3o2_infinite` |
| Gauge algebra / display name | `A2` / `SU(3)` |
| Number of flavours | `5` |
| Chern--Simons level | `\lvert k\rvert=3/2` (absolute-value convention) |
| Enhanced symmetry | `A5` / `SU(6)`, with abelian fugacity `q` |
| Grading variable | `t` |
| Paper equation | `11.3` |
| Committed source | `references/overleaf/su3_5f_6f_hwg_results.tex` |

## Plethystic-exponential transcription

### Original explicit formula (verbatim source metadata)

```tex
\PE\!\left[(\mu_1\mu_5+1)t^2+(q\mu_2+q^{-1}\mu_4)t^3
+\mu_2\mu_4t^4-\mu_2\mu_4t^6\right]
```

### Normalized generated formula

```tex
\operatorname{PE}\!\left[mu_1 mu_5 t^{2} + t^{2} + mu_2 q t^{3} + mu_4 q^{-1} t^{3} + mu_2 mu_4 t^{4} - mu_2 mu_4 t^{6}\right]
```

### Term-by-term table

| # | coefficient | degree | A5 Dynkin labels | q charge |
|---:|---:|---:|---|---:|
| 1 | +1 | 2 | `[1,0,0,0,1]` | 0 |
| 2 | +1 | 2 | `[0,0,0,0,0]` | 0 |
| 3 | +1 | 3 | `[0,1,0,0,0]` | +1 |
| 4 | +1 | 3 | `[0,0,0,1,0]` | -1 |
| 5 | +1 | 4 | `[0,1,0,1,0]` | 0 |
| 6 | -1 | 6 | `[0,1,0,1,0]` | 0 |

## Rational-product transcription

### Original explicit formula (verbatim source metadata)

```tex
\frac{1-\mu_2\mu_4t^6}
{(1-t^2)(1-\mu_1\mu_5t^2)(1-q\mu_2t^3)(1-q^{-1}\mu_4t^3)(1-\mu_2\mu_4t^4)}
```

### Normalized generated rational product

```tex
\left(1 - mu_1 mu_5 t^{2}\right)^{-1} \left(1 - mu_2 mu_4 t^{4}\right)^{-1} \left(1 - mu_2 q t^{3}\right)^{-1} \left(1 - mu_4 q^{-1} t^{3}\right)^{-1} \left(1 - t^{2}\right)^{-1} \left(1 - mu_2 mu_4 t^{6}\right)
```

### Product-factor table

| # | factor | power | location |
|---:|---|---:|---|
| 1 | `(1 - mu_2 mu_4 t^6)` | +1 | numerator |
| 2 | `(1 - t^2)` | -1 | denominator |
| 3 | `(1 - mu_1 mu_5 t^2)` | -1 | denominator |
| 4 | `(1 - q mu_2 t^3)` | -1 | denominator |
| 5 | `(1 - q^(-1) mu_4 t^3)` | -1 | denominator |
| 6 | `(1 - mu_2 mu_4 t^4)` | -1 | denominator |

## Validation results

- **PASS — source transcription:** both original strings are retained literally in the fixture.
- **PASS — exactness:** the level and every charge are exact rational values; no floating-point values occur.
- **PASS — PE structure:** exactly six structured terms match the source monomials and signs.
- **PASS — product structure:** exactly six factors are present; the numerator has power `+1` and all denominator factors have power `-1`.
- **PASS — mathematical equivalence:** for every monomial, the product power is the negative of its PE coefficient, as required by `PE[sum(c m)] = product(1-m)^(-c)`.
- Normalization changes notation, spacing, line wrapping, and factor order only; these are not mathematical differences.
