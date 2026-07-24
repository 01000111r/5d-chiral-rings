import hashlib
import json
from pathlib import Path
import re

import pytest

from sage.all import WeylCharacterRing

from hwg_pipeline.characters import dimension_refine, restore_characters, unrefine
from hwg_pipeline.expansion import expand_pe, expand_rational_product
from hwg_pipeline.io import load_theory
from hwg_pipeline.plethystic import plethystic_exponential, plethystic_logarithm, VirtualCharacterSeries


@pytest.fixture
def root_dir():
    return Path(__file__).resolve().parents[1]


def fixture(root_dir):
    return load_theory(root_dir / "theories/su3_nf5_finite.yaml")


def test_authoritative_source_and_exact_specialization(root_dir):
    source = (root_dir / "references/overleaf/su3_nf5_nf6_finite_hwg_results.tex").read_text()
    assert "Equation~(5.52)" in source and "finite-coupling" in source
    assert "+\\mu_3\\mu_2t^6\n-\\mu_3\\mu_2t^6" in source
    theory = fixture(root_dir)
    assert len(theory.pe.terms) == 5
    assert sum((x.coefficient for x in (type(theory.pe.terms[0])(1, theory.pe.terms[-1].monomial),
                                        type(theory.pe.terms[0])(-1, theory.pe.terms[-1].monomial)))) == 0


def test_three_exact_routes_and_beta_preservation(root_dir):
    theory = fixture(root_dir)
    pe = expand_pe(theory, 10)
    assert pe == expand_rational_product(theory, 10)
    assert all(set(dict(m.abelian_charges)) == {"beta"} for m, _ in pe)
    assert all(len(m.representations[0].dynkin_labels) == 4 for m, _ in pe)


def test_a4_conventions_dimensions_and_conjugation():
    ring = WeylCharacterRing("A4", style="coroots")
    expected = {(0,0,0,0):1, (1,0,0,0):5, (0,0,0,1):5,
                (1,0,0,1):24, (0,1,0,0):10, (0,0,1,0):10,
                (0,1,1,0):75}
    assert {labels: ring(labels).degree() for labels in expected} == expected
    assert tuple(reversed((0,0,1,0))) == (0,1,0,0)
    assert tuple(reversed((0,1,1,0))) == (0,1,1,0)


def test_dimension_unrefinement_pl_and_refined_reconstruction(root_dir):
    theory = fixture(root_dir)
    characters = restore_characters(theory, expand_pe(theory, 10))
    dimensions = dimension_refine(characters)
    assert dict(unrefine(characters)) == {
        d: sum(c for (sd, _), c in dimensions if sd == d) for d, _ in unrefine(characters)}
    pl = plethystic_logarithm(characters, 10)
    expected = VirtualCharacterSeries.from_character_series(characters, 10)
    assert plethystic_exponential(pl, 10) == expected


def test_audit_normalization_no_instanton_and_deterministic_report(root_dir):
    audit = json.loads((root_dir / "generated/su3_nf5_finite/input_audit.json").read_text())
    assert audit["physical_baryon_conversion"] == "B = 3 B_beta"
    assert audit["finite_coupling_instanton_charge"] == "I = 0"
    assert audit["input_checks"]["all_passed"]
    report = root_dir / "generated/su3_nf5_finite/order_10/compact_report/compact_results.tex"
    before = hashlib.sha256(report.read_bytes()).hexdigest()
    text = report.read_text()
    assert "beta" in text and "U(1)_q" not in text and "t,q" not in text
    assert r"\beta" in text
    assert "{ABELIAN_ID}" not in text
    assert not re.search(r"(?<![A-Za-z\\])beta(?:\^|\b)", text)
    assert hashlib.sha256(report.read_bytes()).hexdigest() == before
