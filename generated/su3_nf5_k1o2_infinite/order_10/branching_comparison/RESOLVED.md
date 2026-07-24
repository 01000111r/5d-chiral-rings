# Product-group branching blocker resolved

## Original blocker

The shared `_terms` parser retained only the first stored irrep and the generator imposed one length-five label check.  This discarded the UV `A1` representation belonging to the ordered `A4 x A1` parent.

## Shared fix

The resolving commit introduces ordered `ProductIrrep`/`FactorIrrep` values, validates each factor independently, preserves `A4`, and branches `A1` to its exact weights while retaining external `q`.  JSON and LaTeX preserve parentage and render product labels with a semicolon.

## Verification

- Focused comparison tests: `35 passed`.
- Full suite: `187 passed`.
- The strict order-10 branching-comparison command completed twice with deterministic mathematical outputs.
- Existing `su3_nf5_k3o2_infinite` and `su3_nf5_k5o2_infinite` reports remained byte-for-byte unchanged.
- The resolving commit is the commit containing this file (see repository history).
