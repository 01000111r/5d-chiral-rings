# Source/input audit: `su4_nf7_k5o2_infinite`

- **PASS — authoritative source:** Eq. (10.3) in `references/overleaf/su3_5f_6f_hwg_results.tex`.
- **PASS — specialization:** set `N=4`, hence `N_f=2N-1=7` and enhanced symmetry `SO(14) x U(1)_q` (`D7`).
- **PASS — finite sum:** `i=1,2` gives `mu_2 t^2 + mu_4 t^4`.
- **PASS — singlet:** `t^2` is retained.
- **PASS — charged term:** `t^N(q mu_(2N-2)+q^-1 mu_(2N-1))` gives `(q mu_6+q^-1 mu_7)t^4`.
- **PASS — spinor-product term:** `mu_(2N-2)mu_(2N-1)(t^(2N-2)-t^(2N))` gives `mu_6 mu_7(t^6-t^8)`.
- **PASS — fixture equality:** the seven signed structured terms equal this specialization exactly.

Thus

`HWG = PE[(1+mu_2)t^2+(mu_4+q mu_6+q^-1 mu_7)t^4+mu_6 mu_7 t^6-mu_6 mu_7 t^8]`.

The rational product stored in the fixture is an exact algebraic conversion of this specialized PE. The source prints the general equation, not the `N=4` expression literally.
