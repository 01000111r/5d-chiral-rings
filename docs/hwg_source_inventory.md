# HWG LaTeX source inventory

## Scope and transcription policy

This inventory audits `references/overleaf/su3_5f_6f_hwg_results.tex`. It is a
transcription of that source, not a derivation or a correction. In particular,
the PE exponents below retain the source's ordering, signs, indices, and
parenthesisation. A paper equation number is recorded only where the source
actually cites one.

The source says that the paper's gauge-rank variable `n` has been rewritten as
`N`, that `\mu_i` are highest-weight fugacities, that `q` is a `U(1)`
fugacity, and that `\nu`, `\nu_1`, and `\nu_2` are `SU(2)` highest-weight
fugacities. It also says that the results depend on `|k|`. The grading variable
`t` is not described there as a highest-weight fugacity.

For stable proposed IDs, this inventory uses
`su<rank>_nf<flavours>_k<nonnegative-level>`, with `p` standing for the slash
in a half-integer (for example, `3p2` means `3/2`). General families use `sun`
and retain `N` in the ID.

## General \(N_f=2N\) formulas

### `sun_nf2n_k2`

- **Five-dimensional gauge theory:** `SU(N)+2N F` (the general family in the
  source's `N_f=2N` section, with `N` identified there as the gauge rank).
- **\(N_f\):** `2N`.
- **Level:** `|k|=2`.
- **Enhanced symmetry:** `SO(4N)\times U(1)`.
- **Simple factors and ranks:** `SO(4N)` has rank `2N`.
- **Abelian factors:** one `U(1)`.
- **Highest-weight fugacities:** `\mu_{2},\mu_{4},\ldots,\mu_{2N-2},\mu_{2N}`;
  `q` is the `U(1)` fugacity, not a highest-weight fugacity.
- **Paper equation cited in the source:** `(7.3)`.
- **Exact PE exponent copied from the source:**

  ```tex
  \sum_{i=1}^{N-1}\mu_{2i}t^{2i}+t^2+
  \left(q+q^{-1}\right)\mu_{2N}t^N
  ```

- **Explicit rational-product form present:** no.
- **Positive exponent terms:** `\sum_{i=1}^{N-1}\mu_{2i}t^{2i}`; `t^2`;
  `\left(q+q^{-1}\right)\mu_{2N}t^N`.
- **Negative exponent terms:** none.
- **Notation/convention care:** the source labels the level by `|k|`; the
  indexed sum is left unexpanded; `q+q^{-1}` contains the two `U(1)` charges.

### `sun_nf2n_k1`

- **Five-dimensional gauge theory:** `SU(N)+2N F`.
- **\(N_f\):** `2N`.
- **Level:** `|k|=1`.
- **Enhanced symmetry:** `SU(2N+1)\times U(1)`.
- **Simple factors and ranks:** `SU(2N+1)` has rank `2N`.
- **Abelian factors:** one `U(1)`.
- **Highest-weight fugacities:** `\mu_i` with indices appearing as `i`,
  `2N-i+1`, `N`, and `N+1`; `q` is the `U(1)` fugacity.
- **Paper equation cited in the source:** `(8.3)`.
- **Exact PE exponent copied from the source:**

  ```tex
  \sum_{i=1}^{N-1}\mu_i\mu_{2N-i+1}t^{2i}+t^2+
  \left(q\mu_N+q^{-1}\mu_{N+1}\right)t^N
  ```

- **Explicit rational-product form present:** no.
- **Positive exponent terms:**
  `\sum_{i=1}^{N-1}\mu_i\mu_{2N-i+1}t^{2i}`; `t^2`;
  `\left(q\mu_N+q^{-1}\mu_{N+1}\right)t^N`.
- **Negative exponent terms:** none.
- **Notation/convention care:** the source labels the level by `|k|`; the
  indexed sum is not expanded; `q` and `q^{-1}` distinguish the two terms at
  order `t^N`.

### `sun_nf2n_k0`

- **Five-dimensional gauge theory:** `SU(N)+2N F`.
- **\(N_f\):** `2N`.
- **Level:** `k=0` (the source does not write `|k|=0` in this heading).
- **Enhanced symmetry:** `SU(2N)\times SU(2)\times SU(2)`.
- **Simple factors and ranks:** `SU(2N)` has rank `2N-1`; each `SU(2)` has
  rank `1`.
- **Abelian factors:** none.
- **Highest-weight fugacities:** `\mu_i` with indices appearing as `i`,
  `2N-i`, and `N`; `\nu_1` and `\nu_2` are the two `SU(2)` highest-weight
  fugacities.
- **Paper equation cited in the source:** `(9.5)`.
- **Exact PE exponent copied from the source:**

  ```tex
  \sum_{i=1}^{N}\mu_i\mu_{2N-i}t^{2i}
  +(\nu_1^2+\nu_2^2)t^2+t^4
  +\nu_1\nu_2\mu_N\left(t^N+t^{N+2}\right)
  -\nu_1^2\nu_2^2\mu_N^2t^{2N+4}
  ```

- **Explicit rational-product form present:** no.
- **Positive exponent terms:**
  `\sum_{i=1}^{N}\mu_i\mu_{2N-i}t^{2i}`;
  `(\nu_1^2+\nu_2^2)t^2`; `t^4`;
  `\nu_1\nu_2\mu_N\left(t^N+t^{N+2}\right)`.
- **Negative exponent terms:**
  `-\nu_1^2\nu_2^2\mu_N^2t^{2N+4}`.
- **Notation/convention care:** this is the only `N_f=2N` heading written
  with `k=0`; the sum includes its `i=N` endpoint; the two `\nu` variables
  refer to distinct `SU(2)` factors.

## Explicit \(SU(3)+6F\) specialisations

The explicit section states the convention
`\PE[\sum_a x_a-\sum_b y_b]=\prod_b(1-y_b)/\prod_a(1-x_a)`.

### `su3_nf6_k2`

- **Five-dimensional gauge theory:** `SU(3)+6F`.
- **\(N_f\):** `6`.
- **Level:** `|k|=2`.
- **Enhanced symmetry:** `SO(12)\times U(1)`.
- **Simple factors and ranks:** `SO(12)` has rank `6`.
- **Abelian factors:** one `U(1)`.
- **Highest-weight fugacities:** `\mu_2`, `\mu_4`, `\mu_6`; `q` is the
  `U(1)` fugacity.
- **Paper equation cited in the source:** none in the explicit entry; its
  corresponding general formula cites `(7.3)`.
- **Exact PE exponent copied from the source:**

  ```tex
  (\mu_2+1)t^2+(q+q^{-1})\mu_6t^3+\mu_4t^4
  ```

- **Explicit rational-product form present:** yes:

  ```tex
  \frac{1}{(1-t^2)(1-\mu_2t^2)(1-q\mu_6t^3)(1-q^{-1}\mu_6t^3)(1-\mu_4t^4)}
  ```

- **Positive exponent terms:** `(\mu_2+1)t^2`;
  `(q+q^{-1})\mu_6t^3`; `\mu_4t^4`.
- **Negative exponent terms:** none.
- **Notation/convention care:** the exponent groups multiple monomials inside
  parentheses; the source gives both PE and rational-product forms.

### `su3_nf6_k1`

- **Five-dimensional gauge theory:** `SU(3)+6F`.
- **\(N_f\):** `6`.
- **Level:** `|k|=1`.
- **Enhanced symmetry:** `SU(7)\times U(1)`.
- **Simple factors and ranks:** `SU(7)` has rank `6`.
- **Abelian factors:** one `U(1)`.
- **Highest-weight fugacities:** `\mu_1`, `\mu_2`, `\mu_3`, `\mu_4`,
  `\mu_5`, `\mu_6`; `q` is the `U(1)` fugacity.
- **Paper equation cited in the source:** none in the explicit entry; its
  corresponding general formula cites `(8.3)`.
- **Exact PE exponent copied from the source:**

  ```tex
  (\mu_1\mu_6+1)t^2+(q\mu_3+q^{-1}\mu_4)t^3+\mu_2\mu_5t^4
  ```

- **Explicit rational-product form present:** yes:

  ```tex
  \frac{1}{(1-t^2)(1-\mu_1\mu_6t^2)(1-q\mu_3t^3)(1-q^{-1}\mu_4t^3)(1-\mu_2\mu_5t^4)}
  ```

- **Positive exponent terms:** `(\mu_1\mu_6+1)t^2`;
  `(q\mu_3+q^{-1}\mu_4)t^3`; `\mu_2\mu_5t^4`.
- **Negative exponent terms:** none.
- **Notation/convention care:** `q` and `q^{-1}` multiply different
  highest-weight fugacities.

### `su3_nf6_k0`

- **Five-dimensional gauge theory:** `SU(3)+6F`.
- **\(N_f\):** `6`.
- **Level:** `k=0`.
- **Enhanced symmetry:** `SU(6)\times SU(2)\times SU(2)`.
- **Simple factors and ranks:** `SU(6)` has rank `5`; each `SU(2)` has rank
  `1`.
- **Abelian factors:** none.
- **Highest-weight fugacities:** `\mu_1`, `\mu_2`, `\mu_3`, `\mu_4`,
  `\mu_5`; `\nu_1` and `\nu_2` are the two `SU(2)` highest-weight
  fugacities.
- **Paper equation cited in the source:** none in the explicit entry; its
  corresponding general formula cites `(9.5)`.
- **Exact PE exponent copied from the source:**

  ```tex
  (\mu_1\mu_5+\nu_1^2+\nu_2^2)t^2
  +\nu_1\nu_2\mu_3t^3
  +(\mu_2\mu_4+1)t^4
  +\nu_1\nu_2\mu_3t^5
  +\mu_3^2t^6
  -\nu_1^2\nu_2^2\mu_3^2t^{10}
  ```

- **Explicit rational-product form present:** yes:

  ```tex
  \frac{1-\nu_1^2\nu_2^2\mu_3^2t^{10}}
  {(1-t^4)(1-\mu_1\mu_5t^2)(1-\nu_1^2t^2)(1-\nu_2^2t^2)
   (1-\nu_1\nu_2\mu_3t^3)(1-\mu_2\mu_4t^4)
   (1-\nu_1\nu_2\mu_3t^5)(1-\mu_3^2t^6)}
  ```

- **Positive exponent terms:**
  `(\mu_1\mu_5+\nu_1^2+\nu_2^2)t^2`;
  `\nu_1\nu_2\mu_3t^3`; `(\mu_2\mu_4+1)t^4`;
  `\nu_1\nu_2\mu_3t^5`; `\mu_3^2t^6`.
- **Negative exponent terms:**
  `-\nu_1^2\nu_2^2\mu_3^2t^{10}`.
- **Notation/convention care:** `k=0` is not written with absolute-value
  bars; the two `\nu` variables represent distinct `SU(2)` factors; the
  negative PE monomial appears as the rational numerator.

## General \(N_f=2N-1\) formulas

### `sun_nf2nminus1_k5p2`

- **Five-dimensional gauge theory:** `SU(N)+(2N-1)F`.
- **\(N_f\):** `2N-1`.
- **Level:** `|k|=\frac52`.
- **Enhanced symmetry:** `SO(4N-2)\times U(1)`.
- **Simple factors and ranks:** `SO(4N-2)` has rank `2N-1`.
- **Abelian factors:** one `U(1)`.
- **Highest-weight fugacities:** `\mu_{2},\mu_{4},\ldots,\mu_{2N-4}`,
  `\mu_{2N-2}`, and `\mu_{2N-1}`; `q` is the `U(1)` fugacity.
- **Paper equation cited in the source:** `(10.3)`.
- **Exact PE exponent copied from the source:**

  ```tex
  \sum_{i=1}^{N-2}\mu_{2i}t^{2i}+t^2
  +t^N\left(q\mu_{2N-2}+q^{-1}\mu_{2N-1}\right)
  +\mu_{2N-2}\mu_{2N-1}\left(t^{2N-2}-t^{2N}\right)
  ```

- **Explicit rational-product form present:** no.
- **Positive exponent terms:** `\sum_{i=1}^{N-2}\mu_{2i}t^{2i}`; `t^2`;
  `t^N\left(q\mu_{2N-2}+q^{-1}\mu_{2N-1}\right)`;
  `+\mu_{2N-2}\mu_{2N-1}t^{2N-2}` (the positive part of the source's
  final parenthesis).
- **Negative exponent terms:**
  `-\mu_{2N-2}\mu_{2N-1}t^{2N}` (the negative part of that parenthesis).
- **Notation/convention care:** the last source term contains both a positive
  and a negative monomial in one parenthesis; this inventory does not simplify
  it in the exact transcription.

### `sun_nf2nminus1_k3p2`

- **Five-dimensional gauge theory:** `SU(N)+(2N-1)F`.
- **\(N_f\):** `2N-1`.
- **Level:** `|k|=\frac32`.
- **Enhanced symmetry:** `SU(2N)\times U(1)`.
- **Simple factors and ranks:** `SU(2N)` has rank `2N-1`.
- **Abelian factors:** one `U(1)`.
- **Highest-weight fugacities:** `\mu_i` with indices appearing as `i`,
  `2N-i`, `N-1`, and `N+1`; `q` is the `U(1)` fugacity.
- **Paper equation cited in the source:** `(11.3)`.
- **Exact PE exponent copied from the source:**

  ```tex
  \sum_{i=1}^{N-1}\mu_i\mu_{2N-i}t^{2i}+t^2
  +\left(q\mu_{N-1}+q^{-1}\mu_{N+1}\right)t^N
  -\mu_{N-1}\mu_{N+1}t^{2N}
  ```

- **Explicit rational-product form present:** no.
- **Positive exponent terms:**
  `\sum_{i=1}^{N-1}\mu_i\mu_{2N-i}t^{2i}`; `t^2`;
  `\left(q\mu_{N-1}+q^{-1}\mu_{N+1}\right)t^N`.
- **Negative exponent terms:** `-\mu_{N-1}\mu_{N+1}t^{2N}`.
- **Notation/convention care:** `q` and `q^{-1}` multiply different
  highest-weight fugacities; the sum remains indexed.

### `sun_nf2nminus1_k1p2`

- **Five-dimensional gauge theory:** `SU(N)+(2N-1)F`.
- **\(N_f\):** `2N-1`.
- **Level:** `|k|=\frac12`.
- **Enhanced symmetry:** `SU(2N-1)\times SU(2)\times U(1)`.
- **Simple factors and ranks:** `SU(2N-1)` has rank `2N-2`; `SU(2)` has
  rank `1`.
- **Abelian factors:** one `U(1)`.
- **Highest-weight fugacities:** `\mu_i` with indices appearing as `i`,
  `2N-i-1`, `N-1`, and `N`; `\nu` is the `SU(2)` highest-weight fugacity;
  `q` is the `U(1)` fugacity.
- **Paper equation cited in the source:** `(12.3)`.
- **Exact PE exponent copied from the source:**

  ```tex
  \sum_{i=1}^{N-1}\mu_i\mu_{2N-i-1}t^{2i}
  +(\nu^2+1)t^2
  +\nu\left(q\mu_{N-1}+q^{-1}\mu_N\right)t^N
  -\nu^2\mu_{N-1}\mu_Nt^{2N}
  ```

- **Explicit rational-product form present:** no.
- **Positive exponent terms:**
  `\sum_{i=1}^{N-1}\mu_i\mu_{2N-i-1}t^{2i}`;
  `(\nu^2+1)t^2`;
  `\nu\left(q\mu_{N-1}+q^{-1}\mu_N\right)t^N`.
- **Negative exponent terms:** `-\nu^2\mu_{N-1}\mu_Nt^{2N}`.
- **Notation/convention care:** `\nu` is a highest-weight fugacity for the
  simple `SU(2)` factor, while `q` is the abelian fugacity.

## Explicit \(SU(3)+5F\) specialisations

### `su3_nf5_k5p2`

- **Five-dimensional gauge theory:** `SU(3)+5F`.
- **\(N_f\):** `5`.
- **Level:** `|k|=\frac52`.
- **Enhanced symmetry:** `SO(10)\times U(1)`.
- **Simple factors and ranks:** `SO(10)` has rank `5`.
- **Abelian factors:** one `U(1)`.
- **Highest-weight fugacities:** `\mu_2`, `\mu_4`, `\mu_5`; `q` is the
  `U(1)` fugacity.
- **Paper equation cited in the source:** none in the explicit entry; its
  corresponding general formula cites `(10.3)`.
- **Exact PE exponent copied from the source:**

  ```tex
  (\mu_2+1)t^2+(q\mu_4+q^{-1}\mu_5)t^3
  +\mu_4\mu_5t^4-\mu_4\mu_5t^6
  ```

- **Explicit rational-product form present:** yes:

  ```tex
  \frac{1-\mu_4\mu_5t^6}
  {(1-t^2)(1-\mu_2t^2)(1-q\mu_4t^3)(1-q^{-1}\mu_5t^3)(1-\mu_4\mu_5t^4)}
  ```

- **Positive exponent terms:** `(\mu_2+1)t^2`;
  `(q\mu_4+q^{-1}\mu_5)t^3`; `+\mu_4\mu_5t^4`.
- **Negative exponent terms:** `-\mu_4\mu_5t^6`.
- **Notation/convention care:** `q` and `q^{-1}` multiply different spinor
  fugacities as written; the negative monomial appears in the rational
  numerator.

### `su3_nf5_k3p2`

- **Five-dimensional gauge theory:** `SU(3)+5F`.
- **\(N_f\):** `5`.
- **Level:** `|k|=\frac32`.
- **Enhanced symmetry:** `SU(6)\times U(1)`.
- **Simple factors and ranks:** `SU(6)` has rank `5`.
- **Abelian factors:** one `U(1)`.
- **Highest-weight fugacities:** `\mu_1`, `\mu_2`, `\mu_4`, `\mu_5`; `q`
  is the `U(1)` fugacity.
- **Paper equation cited in the source:** none in the explicit entry; its
  corresponding general formula cites `(11.3)`.
- **Exact PE exponent copied from the source:**

  ```tex
  (\mu_1\mu_5+1)t^2+(q\mu_2+q^{-1}\mu_4)t^3
  +\mu_2\mu_4t^4-\mu_2\mu_4t^6
  ```

- **Explicit rational-product form present:** yes:

  ```tex
  \frac{1-\mu_2\mu_4t^6}
  {(1-t^2)(1-\mu_1\mu_5t^2)(1-q\mu_2t^3)(1-q^{-1}\mu_4t^3)(1-\mu_2\mu_4t^4)}
  ```

- **Positive exponent terms:** `(\mu_1\mu_5+1)t^2`;
  `(q\mu_2+q^{-1}\mu_4)t^3`; `+\mu_2\mu_4t^4`.
- **Negative exponent terms:** `-\mu_2\mu_4t^6`.
- **Notation/convention care:** `q` and `q^{-1}` multiply different
  highest-weight fugacities; the negative monomial appears in the rational
  numerator.

### `su3_nf5_k1p2`

- **Five-dimensional gauge theory:** `SU(3)+5F`.
- **\(N_f\):** `5`.
- **Level:** `|k|=\frac12`.
- **Enhanced symmetry:** `SU(5)\times SU(2)\times U(1)`.
- **Simple factors and ranks:** `SU(5)` has rank `4`; `SU(2)` has rank `1`.
- **Abelian factors:** one `U(1)`.
- **Highest-weight fugacities:** `\mu_1`, `\mu_2`, `\mu_3`, `\mu_4`; `\nu`
  is the `SU(2)` highest-weight fugacity; `q` is the `U(1)` fugacity.
- **Paper equation cited in the source:** none in the explicit entry; its
  corresponding general formula cites `(12.3)`.
- **Exact PE exponent copied from the source:**

  ```tex
  (\mu_1\mu_4+\nu^2+1)t^2
  +\nu(q\mu_2+q^{-1}\mu_3)t^3
  +\mu_2\mu_3t^4-\nu^2\mu_2\mu_3t^6
  ```

- **Explicit rational-product form present:** yes:

  ```tex
  \frac{1-\nu^2\mu_2\mu_3t^6}
  {(1-t^2)(1-\nu^2t^2)(1-\mu_1\mu_4t^2)
   (1-\nu q\mu_2t^3)(1-\nu q^{-1}\mu_3t^3)(1-\mu_2\mu_3t^4)}
  ```

- **Positive exponent terms:** `(\mu_1\mu_4+\nu^2+1)t^2`;
  `+\nu(q\mu_2+q^{-1}\mu_3)t^3`; `+\mu_2\mu_3t^4`.
- **Negative exponent terms:** `-\nu^2\mu_2\mu_3t^6`.
- **Notation/convention care:** `\nu` and `q` belong to different symmetry
  factors; the negative monomial appears in the rational numerator.

## Explicit-case presence check

All six requested explicit cases are present in the source:

- `SU(3)+6F`, `|k|=2` — `su3_nf6_k2`.
- `SU(3)+6F`, `|k|=1` — `su3_nf6_k1`.
- `SU(3)+6F`, `k=0` — `su3_nf6_k0`.
- `SU(3)+5F`, `|k|=5/2` — `su3_nf5_k5p2`.
- `SU(3)+5F`, `|k|=3/2` — `su3_nf5_k3p2`.
- `SU(3)+5F`, `|k|=1/2` — `su3_nf5_k1p2`.

## Audit observations (not corrections)

- The source cites paper equation numbers on the six general formulas, but not
  on the six explicit `N=3` entries. This inventory therefore distinguishes
  “none in the explicit entry” from the citation on the corresponding general
  formula rather than silently assigning that number to the explicit display.
- The source says generally that results depend on `|k|`, while both the
  general and explicit zero-level headings use `k=0` rather than `|k|=0`.
- The `SU(3)+6F`, `k=0` explicit exponent contains `+\mu_3^2t^6`; it has been
  retained exactly as printed.
- No mismatch is apparent between any printed explicit PE exponent and the
  result of substituting `N=3` into its corresponding printed general exponent.
  This is only a source-internal transcription check, not a validation of the
  physics formulas against the cited paper.

## Recommended first explicit entry

Enter **`su3_nf5_k3p2` (`SU(3)+5F`, `|k|=3/2`) first**. The source reveals no
reason to depart from the requested preference: this case has one simple
non-abelian factor and one abelian factor, an explicit PE exponent, an explicit
rational-product form, and both positive and negative exponent terms. It is
therefore a compact first transcription target while still exercising the
sign and numerator convention visible in the source.
