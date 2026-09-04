"""Complete exact UV checks for SU(4)+7F at |k|=5/2 through t^10."""

import json
from pathlib import Path

import pytest
from sage.all import QQ, ZZ

from hwg_pipeline import load_theory
from hwg_pipeline.characters import restore_characters, unrefine
from hwg_pipeline.expansion import expand_pe, expand_rational_product
from hwg_pipeline.plethystic import (
    dimension_refine_virtual,
    plethystic_logarithm,
    scalar_plethystic_logarithm,
    unrefine_virtual,
)
from hwg_pipeline.sage_backend import irrep, irrep_dimension


THEORY_PATH = Path("theories/su4_nf7_k5o2_infinite.yaml")
SOURCE_PATH = Path("references/overleaf/su3_5f_6f_hwg_results.tex")
OUTPUT = Path("generated/su4_nf7_k5o2_infinite/order_10")


def signature(item):
    monomial = item.monomial
    value = getattr(item, "coefficient", getattr(item, "power", None))
    return (int(value), int(monomial.t_degree),
            tuple(map(int, monomial.representations[0].dynkin_labels)),
            int(monomial.abelian_charges[0][1]))


@pytest.fixture(scope="module")
def calculation():
    theory = load_theory(THEORY_PATH)
    hwg = expand_pe(theory, 10)
    characters = restore_characters(theory, hwg)
    pl = plethystic_logarithm(characters, 10)
    return theory, hwg, characters, pl


def test_eq_10_3_and_exact_n4_specialization():
    source = SOURCE_PATH.read_text(encoding="utf-8")
    assert r"(paper eq.~(10.3))" in source
    assert r"\sum_{i=1}^{N-2}\mu_{2i}t^{2i}+t^2" in source
    assert r"t^N\left(q\mu_{2N-2}+q^{-1}\mu_{2N-1}\right)" in source
    assert r"\mu_{2N-2}\mu_{2N-1}\left(t^{2N-2}-t^{2N}\right)" in source
    audit = json.loads(Path("generated/su4_nf7_k5o2_infinite/input_audit.json").read_text())
    assert audit["exact_match"] is True
    assert audit["specialization"]["terms"] == [
        {"origin": "sum, i=1", "result": "mu_2 t^2"},
        {"origin": "sum, i=2", "result": "mu_4 t^4"},
        {"origin": "singlet", "result": "t^2"},
        {"origin": "charged spinors", "result": "(q mu_6 + q^-1 mu_7) t^4"},
        {"origin": "spinor product", "result": "mu_6 mu_7 (t^6 - t^8)"},
    ]


def test_fixture_has_exact_metadata_and_seven_signed_terms(calculation):
    theory = calculation[0]
    assert (theory.id, theory.gauge_algebra, theory.gauge_display_name) == (
        "su4_nf7_k5o2_infinite", "A3", "SU(4)")
    assert theory.number_of_flavours == 7
    assert theory.chern_simons_level == QQ(5) / 2
    assert theory.chern_simons_convention == "absolute value"
    assert theory.coupling == "infinite"
    assert theory.simple_factors[0].cartan_name == "D7"
    assert theory.simple_factors[0].display_name == "SO(14)"
    assert theory.source_references[0].equation == "10.3"
    assert [signature(term) for term in theory.pe.terms] == [
        (1, 2, (0, 0, 0, 0, 0, 0, 0), 0),
        (1, 2, (0, 1, 0, 0, 0, 0, 0), 0),
        (1, 4, (0, 0, 0, 1, 0, 0, 0), 0),
        (1, 4, (0, 0, 0, 0, 0, 1, 0), 1),
        (1, 4, (0, 0, 0, 0, 0, 0, 1), -1),
        (1, 6, (0, 0, 0, 0, 0, 1, 1), 0),
        (-1, 8, (0, 0, 0, 0, 0, 1, 1), 0),
    ]


def test_exact_rational_product_and_independent_expansion(calculation):
    theory, hwg, _, _ = calculation
    assert [signature(term) for term in theory.rational_product.factors] == [
        (1, 8, (0, 0, 0, 0, 0, 1, 1), 0),
        (-1, 2, (0, 0, 0, 0, 0, 0, 0), 0),
        (-1, 2, (0, 1, 0, 0, 0, 0, 0), 0),
        (-1, 4, (0, 0, 0, 1, 0, 0, 0), 0),
        (-1, 4, (0, 0, 0, 0, 0, 1, 0), 1),
        (-1, 4, (0, 0, 0, 0, 0, 0, 1), -1),
        (-1, 6, (0, 0, 0, 0, 0, 1, 1), 0),
    ]
    assert hwg == expand_rational_product(theory, 10)
    assert all(term.monomial.t_degree <= 10 for term in theory.pe.terms)
    assert not any(monomial.t_degree % 2 for monomial, _ in hwg)


def test_exact_d7_dimensions_and_spinor_conjugation(calculation):
    factor = calculation[0].simple_factors[0]
    expected = {
        (0, 0, 0, 0, 0, 0, 0): 1,
        (0, 1, 0, 0, 0, 0, 0): 91,
        (0, 0, 0, 1, 0, 0, 0): 1001,
        (0, 0, 0, 0, 0, 1, 0): 64,
        (0, 0, 0, 0, 0, 0, 1): 64,
        (0, 0, 0, 0, 0, 1, 1): 3003,
        (2, 0, 0, 0, 0, 0, 0): 104,
    }
    assert {labels: int(irrep_dimension(factor, labels)) for labels in expected} == expected
    assert irrep(factor, (0, 0, 0, 0, 0, 1, 0)).dual() == irrep(
        factor, (0, 0, 0, 0, 0, 0, 1))


def test_complete_low_degree_native_pl_and_scalar_benchmarks(calculation):
    _, _, characters, pl = calculation
    sectors = {}
    for (degree, charges), content in pl:
        for labels, coefficient in content:
            sectors.setdefault(int(degree), {})[(int(dict(charges)["q"]), labels[0])] = int(coefficient)
    assert sectors[2] == {
        (0, (0, 0, 0, 0, 0, 0, 0)): 1,
        (0, (0, 1, 0, 0, 0, 0, 0)): 1,
    }
    assert sectors[4] == {
        (1, (0, 0, 0, 0, 0, 1, 0)): 1,
        (-1, (0, 0, 0, 0, 0, 0, 1)): 1,
        (0, (0, 0, 0, 0, 0, 0, 0)): -1,
        (0, (2, 0, 0, 0, 0, 0, 0)): -1,
    }
    dimensions = {(int(degree), int(dict(charges)["q"])): int(coefficient)
                  for (degree, charges), coefficient in dimension_refine_virtual(pl)}
    assert dimensions[(2, 0)] == 92
    assert {charge: dimensions[(4, charge)] for charge in (-1, 0, 1)} == {-1: 64, 0: -105, 1: 64}
    assert dict(unrefine_virtual(pl))[ZZ(4)] == 23
    assert scalar_plethystic_logarithm(unrefine(characters), 10) == unrefine_virtual(pl)
    assert not any(degree % 2 for (degree, _), _ in pl)


def test_complete_persisted_pl_and_exact_reconstruction():
    expected_dim = {
        "2": {0: 92},
        "4": {-1: 64, 0: -105, 1: 64},
        "6": {-1: -896, 0: 104, 1: -896},
        "8": {-2: -364, -1: 7552, 0: -7190, 1: 7552, 2: -364},
        "10": {-2: 17018, -1: -48320, 0: 186851, 1: -48320, 2: 17018},
    }
    payload = json.loads((OUTPUT / "q_refined_dimension_pl.json").read_text())
    actual = {degree: {int(entry["abelian_charges"]["q"]): entry["coefficient"] for entry in entries}
              for degree, entries in payload["coefficients_by_t_degree"].items()}
    assert actual == expected_dim
    unrefined = json.loads((OUTPUT / "unrefined_plethystic_logarithm.json").read_text())
    assert unrefined["coefficients_by_t_degree"] == {
        "2": 92, "4": 23, "6": -1688, "8": 7186, "10": 124247}
    difference = json.loads((OUTPUT / "reconstruction_difference.json").read_text())
    checks = json.loads((OUTPUT / "reconstruction_checks.json").read_text())
    assert difference == {"maximum_t_degree": 10, "mismatch_count": 0, "mismatches": [],
                          "theory_id": "su4_nf7_k5o2_infinite"}
    assert checks["validation_results"]["all_passed"] is True
