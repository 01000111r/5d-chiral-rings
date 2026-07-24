# Input audit: `su3_nf6_finite`

- **Source authority:** `references/overleaf/su3_nf5_nf6_finite_hwg_results.tex` (`SHA-256 44759e2f6ffbc7004a7b4a1d370bb8a5395a2da1fe0ce2675348410e13bd30c1`)
- **Fixture:** `theories/su3_nf6_finite.yaml` (`SHA-256 a7aac0d923723d12c9b10e6db988e15a2890731e915597ea0ba83cb9145ead68`)
- **Paper:** A. Bourget et al., *Brane webs and magnetic quivers for SQCD*, equation (5.52).
- **Substitution:** $N_c=3$, $N_f=6$; $6\geq5$.
- **Unsimplified exponent:** $t^2+(\mu_3\beta+\mu_3\beta^{-1})t^3+\mu_1\mu_5t^2+\mu_2\mu_4t^4+\mu_3^2t^6-\mu_3^2t^6$.
- **Exact cancellation:** $+\mu_3^2t^6-\mu_3^2t^6=0$.
- **Simplified PE:** `PE[(1 + mu_1 mu_5)t^2 + mu_3(beta + beta^(-1))t^3 + mu_2 mu_4 t^4]`
- **Product:** `1/[(1-t^2)(1-mu_1 mu_5 t^2)(1-beta mu_3 t^3)(1-beta^(-1) mu_3 t^3)(1-mu_2 mu_4 t^4)]`
- **Convention:** A5 coroot Dynkin labels; $B_\beta$ is the exponent of $\beta$, $B=3B_\beta$, and $I=0$.
- **Self-conjugacy:** $[0,0,1,0,0]$ and $[0,1,0,1,0]$ are self-conjugate. Baryon and antibaryon remain distinct at opposite beta charge.
- **Checks:** all twelve required input checks passed.
- **Unresolved ambiguities:** none.
