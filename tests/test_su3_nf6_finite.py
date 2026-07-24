import hashlib
import json
from pathlib import Path

from sage.all import WeylCharacterRing

from hwg_pipeline.characters import dimension_refine, restore_characters, unrefine
from hwg_pipeline.expansion import expand_pe, expand_rational_product
from hwg_pipeline.io import load_theory
from hwg_pipeline.plethystic import plethystic_exponential, plethystic_logarithm, VirtualCharacterSeries

ROOT = Path(__file__).resolve().parents[1]


def theory():
    return load_theory(ROOT / "theories/su3_nf6_finite.yaml")


def test_source_specialization_and_exact_cancellation():
    source = (ROOT / "references/overleaf/su3_nf5_nf6_finite_hwg_results.tex").read_text()
    assert "Equation~(5.52)" in source and "finite-coupling" in source
    assert "+\\mu_3^2t^6\n-\\mu_3^2t^6" in source
    audit = json.loads((ROOT / "generated/su3_nf6_finite/input_audit.json").read_text())
    assert audit["input_checks"]["all_passed"]
    assert sum(x["coefficient"] for x in audit["cancelled_source_terms"]) == 0


def test_three_routes_beta_distinction_and_unrefinement():
    data = theory()
    expanded = expand_pe(data, 10)
    assert expanded == expand_rational_product(data, 10)
    chars = restore_characters(data, expanded)
    degree_three = {(dict(q)["beta"], labels) for (d, q), content in chars if d == 3 for labels, _ in content}
    assert (1, ((0, 0, 1, 0, 0),)) in degree_three
    assert (-1, ((0, 0, 1, 0, 0),)) in degree_three
    dimensions = dimension_refine(chars)
    assert dict(unrefine(chars)) == {d: sum(c for (sd, _), c in dimensions if sd == d) for d, _ in unrefine(chars)}


def test_a5_exact_conventions_and_conjugation():
    ring = WeylCharacterRing("A5", style="coroots")
    expected = {(0,0,0,0,0):1, (1,0,0,0,0):6, (0,0,0,0,1):6,
                (1,0,0,0,1):35, (0,0,1,0,0):20, (0,1,0,1,0):189}
    assert {labels: ring(labels).degree() for labels in expected} == expected
    assert tuple(reversed((0,0,1,0,0))) == (0,0,1,0,0)
    assert tuple(reversed((0,1,0,1,0))) == (0,1,0,1,0)


def test_refined_pl_reconstructs_and_normalizations():
    data = theory()
    chars = restore_characters(data, expand_pe(data, 10))
    pl = plethystic_logarithm(chars, 10)
    assert plethystic_exponential(pl, 10) == VirtualCharacterSeries.from_character_series(chars, 10)
    audit = json.loads((ROOT / "generated/su3_nf6_finite/input_audit.json").read_text())
    assert audit["physical_baryon_conversion"] == "B = 3 B_beta"
    assert audit["finite_coupling_instanton_charge"] == "I = 0"


def test_deterministic_beta_report_without_instanton_fugacity():
    report = ROOT / "generated/su3_nf6_finite/order_10/compact_report/compact_results.tex"
    before = hashlib.sha256(report.read_bytes()).hexdigest()
    text = report.read_text()
    assert r"\beta" in text and "U(1)_q" not in text and "instanton fugacity" not in text
    assert hashlib.sha256(report.read_bytes()).hexdigest() == before
