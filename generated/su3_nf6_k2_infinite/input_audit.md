# Input audit: `su3_nf6_k2_infinite`

- **PASS** — stored source `references/overleaf/su3_5f_6f_hwg_results.tex`, equation (7.3), was read before expansion.
- **PASS** — fixture and explicit source agree term by term: five positive PE terms and five denominator factors.
- Source SHA-256: `2e79f9e66dfb0a2a7d2086f642aa94b45b5284287ce41ade80a128dc0439b359`
- Fixture SHA-256: `1f99f98c4c48ee59dcd88153c8bf7acc37c7a986f0a24ba535e09b409ccbdb58`
- Convention: `WeylCharacterRing("D6", style="coroots")`; source uses `mu_6 = [0,0,0,0,0,1]`.
- Charges are exact integer exponents of `q`; no unresolved ambiguities.

## Source PE

```latex
\PE\!\left[(\mu_2+1)t^2+(q+q^{-1})\mu_6t^3+\mu_4t^4\right]
```

## Source rational product

```latex
\frac{1}{(1-t^2)(1-\mu_2t^2)(1-q\mu_6t^3)(1-q^{-1}\mu_6t^3)(1-\mu_4t^4)}
```
