"""Authoritative exact benchmarks for the finite SU(4)+7F Higgs branch."""
import hashlib
import json
from pathlib import Path

from sage.all import WeylCharacterRing

from hwg_pipeline.characters import dimension_refine, restore_characters, unrefine
from hwg_pipeline.expansion import expand_pe, expand_rational_product
from hwg_pipeline.io import load_theory
from hwg_pipeline.plethystic import plethystic_exponential, plethystic_logarithm, VirtualCharacterSeries

ROOT = Path(__file__).resolve().parents[1]
THEORY = ROOT / "theories/su4_nf7_finite.yaml"
OUT = ROOT / "generated/su4_nf7_finite/order_10"


def fixture():
    return load_theory(THEORY)


def charged(payload):
    return {int(d): {int(e["abelian_charges"]["beta"]): e["coefficient"] for e in entries}
            for d, entries in payload["coefficients_by_t_degree"].items()}


def test_authoritative_literal_specialization_and_audit():
    source = (ROOT / "references/overleaf/su4_nf7_finite_hwg_results.tex").read_text()
    assert "equation~(5.52)" in source and "N_c=4,N_f=7" in source
    assert "+\\mu_4\\mu_3t^8-\\mu_4\\mu_3t^8=0" in source
    audit = json.loads((ROOT / "generated/su4_nf7_finite/input_audit.json").read_text())
    assert audit["substitution"] == {"N_c": 4, "N_f": 7}
    assert audit["domain_check"] == "7 >= 2*4-1"
    assert audit["exact_cancellation"] == "+mu_4 mu_3 t^8 - mu_4 mu_3 t^8 = 0"
    assert audit["source_sha256"] == hashlib.sha256((ROOT / audit["source_file"]).read_bytes()).hexdigest()
    assert audit["fixture_sha256"] == hashlib.sha256(THEORY.read_bytes()).hexdigest()
    assert audit["unresolved_ambiguities"] == [] and audit["input_checks"]["all_passed"]


def test_exact_six_term_pe_product_and_charge_conventions():
    theory = fixture()
    assert len(theory.pe.terms) == len(theory.rational_product.factors) == 6
    assert all(term.coefficient == 1 for term in theory.pe.terms)
    assert all(factor.power == -1 for factor in theory.rational_product.factors)
    assert expand_pe(theory, 10) == expand_rational_product(theory, 10)
    expected = [(2,(0,0,0,0,0,0),0),(2,(1,0,0,0,0,1),0),
                (4,(0,1,0,0,1,0),0),(4,(0,0,0,1,0,0),1),
                (4,(0,0,1,0,0,0),-1),(6,(0,0,1,1,0,0),0)]
    actual = [(int(x.monomial.t_degree), tuple(x.monomial.representations[0].dynkin_labels),
               int(dict(x.monomial.abelian_charges)["beta"])) for x in theory.pe.terms]
    assert actual == expected
    assert [x.id for x in theory.abelian_factors] == ["beta"]
    audit = json.loads((ROOT / "generated/su4_nf7_finite/input_audit.json").read_text())
    assert audit["physical_baryon_conversion"] == "B = 4 B_beta"
    assert audit["finite_coupling_instanton_charge"] == "I = 0"


def test_a6_dimensions_and_conjugation():
    ring = WeylCharacterRing("A6", style="coroots")
    expected = {(0,0,0,0,0,0):1, (1,0,0,0,0,1):48,
                (0,0,0,1,0,0):35, (0,0,1,0,0,0):35,
                (0,1,0,0,1,0):392, (0,0,1,1,0,0):784}
    assert {labels: ring(labels).degree() for labels in expected} == expected
    assert tuple(reversed((0,0,0,1,0,0))) == (0,0,1,0,0,0)


def test_independent_dimension_series_benchmarks():
    refined = charged(json.loads((OUT / "beta_refined_dimension_series.json").read_text()))
    assert refined == {0:{0:1}, 2:{0:49}, 4:{-1:35,0:1176,1:35},
        6:{-1:1358,0:18472,1:1358},
        8:{-2:490,-1:26166,0:214963,1:26166,2:490},
        10:{-2:16170,-1:335223,0:1987047,1:335223,2:16170}}
    assert all(v.get(q) == v.get(-q) for v in refined.values() for q in v)
    plain = json.loads((OUT / "unrefined_hilbert_series.json").read_text())["coefficients_by_t_degree"]
    assert plain == {"0":1,"1":0,"2":49,"3":0,"4":1246,"5":0,
                     "6":21188,"7":0,"8":268275,"9":0,"10":2689833}


def test_complete_low_degree_and_dimension_pl_benchmarks():
    raw = json.loads((OUT / "refined_plethystic_logarithm.json").read_text())["coefficients_by_t_degree"]
    def sector(d):
        return {(int(e["abelian_charges"]["beta"]), tuple(e["irreducible_representations"][0]["dynkin_labels"])): e["coefficient"] for e in raw[str(d)]}
    assert sector(2) == {(0,(0,0,0,0,0,0)):1,(0,(1,0,0,0,0,1)):1}
    assert sector(4) == {(1,(0,0,0,1,0,0)):1,(-1,(0,0,1,0,0,0)):1,
                         (0,(0,0,0,0,0,0)):-1,(0,(1,0,0,0,0,1)):-1}
    assert not any(int(d) <= 9 and int(d) % 2 for d in raw)
    refined = charged(json.loads((OUT / "beta_refined_dimension_pl.json").read_text()))
    assert refined == {2:{0:49},4:{-1:35,0:-49,1:35},6:{-1:-357,0:48,1:-357},
        8:{-2:-140,-1:2499,0:-490,1:2499,2:-140},
        10:{-2:4655,-1:-13916,0:12690,1:-13916,2:4655}}
    assert all(v.get(q) == v.get(-q) for v in refined.values() for q in v)
    plain = json.loads((OUT / "unrefined_plethystic_logarithm.json").read_text())["coefficients_by_t_degree"]
    assert plain == {"2":49,"4":21,"6":-666,"8":4228,"10":-5832}


def test_exact_refined_reconstruction_and_persisted_checks():
    theory = fixture()
    characters = restore_characters(theory, expand_pe(theory, 10))
    pl = plethystic_logarithm(characters, 10)
    assert plethystic_exponential(pl, 10) == VirtualCharacterSeries.from_character_series(characters, 10)
    dimensions = dimension_refine(characters)
    assert dict(unrefine(characters)) == {d: sum(c for (sd,_),c in dimensions if sd == d) for d,_ in unrefine(characters)}
    difference = json.loads((OUT / "reconstruction_difference.json").read_text())
    checks = json.loads((OUT / "reconstruction_checks.json").read_text())["validation_results"]
    assert difference["mismatch_count"] == 0 and difference["mismatches"] == []
    assert all(checks.values())
