# Input audit: `su3_nf5_finite`

Source: `references/overleaf/su3_nf5_nf6_finite_hwg_results.tex` (`44759e2f6ffbc7004a7b4a1d370bb8a5395a2da1fe0ce2675348410e13bd30c1`)

The committed local LaTeX extract is accepted as sufficient source authority; no full paper PDF is required.

General formula: `PE[t^2 + (mu_Nc beta + mu_(Nf-Nc) beta^(-1)) t^Nc + sum_(j=1)^Nc mu_j mu_(Nf-j) t^(2j) - mu_Nc mu_(Nf-Nc) t^(2Nc)]`

Unsimplified specialization: `t^2 + (mu_3 beta + mu_2 beta^(-1)) t^3 + mu_1 mu_4 t^2 + mu_2 mu_3 t^4 + mu_3 mu_2 t^6 - mu_3 mu_2 t^6`

Exact cancellation: `+mu_3 mu_2 t^6 - mu_3 mu_2 t^6 = 0`

Simplified HWG: `PE[(1 + mu_1 mu_4)t^2 + (beta mu_3 + beta^(-1) mu_2)t^3 + mu_2 mu_3 t^4]`

Product: `1/[(1-t^2)(1-mu_1 mu_4 t^2)(1-beta mu_3 t^3)(1-beta^(-1) mu_2 t^3)(1-mu_2 mu_3 t^4)]`

## Checks

- **PASS — local_source_exists**
- **PASS — source_says_finite_coupling**
- **PASS — source_identifies_equation_5_52**
- **PASS — domain_5_ge_5**
- **PASS — degree_six_terms_identical_apart_from_sign**
- **PASS — degree_six_sum_zero**
- **PASS — five_positive_simplified_terms**
- **PASS — five_denominator_no_numerator**
- **PASS — all_A4_labels_length_four**
- **PASS — beta_only_abelian_fugacity**
- **PASS — no_instanton_fugacity**
- **PASS — all_passed**
