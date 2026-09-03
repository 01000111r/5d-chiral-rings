"""Verification gates for SU(4)+7F, |k|=1/2 through physical branching."""

import json
from pathlib import Path

import yaml
from sage.all import QQ, ZZ, matrix

from hwg_pipeline import load_theory
from hwg_pipeline.branching import branch_a1_to_u1
from hwg_pipeline.branching_conventions import solve_two_anchor_map
from hwg_pipeline.charge_maps import load_charge_map_spec, solve_charge_map
from hwg_pipeline.sage_backend import irrep_dimension


ROOT = Path(__file__).resolve().parents[1]
THEORY = ROOT / "theories/su4_nf7_k1o2_infinite.yaml"
OUT = ROOT / "generated/su4_nf7_k1o2_infinite/order_10"
RAW = OUT / "manifest_branching"
PHYSICAL = RAW / "physical_charges"


def _json(path):
    return json.loads(path.read_text())


def _q(value):
    return QQ(value["numerator"]) / QQ(value["denominator"]) if isinstance(value, dict) else QQ(value)


def _index(path, physical=False):
    data = _json(path)["coefficients_by_t_degree"]
    result = {}
    for degree, entries in data.items():
        for entry in entries:
            charges = entry["physical_charges"] if physical else entry["raw_charges"]
            multiplicity = _q(entry.get("coefficient", entry.get("multiplicity")))
            names = ("B", "I") if physical else ("x", "q")
            result[(int(degree), tuple(entry["child_dynkin_labels"]), *(_q(charges[n]) for n in names))] = multiplicity
    return result


def test_authoritative_source_and_exact_fixture():
    source = (ROOT / "references/overleaf/su3_5f_6f_hwg_results.tex").read_text()
    assert "paper eq.~(12.3)" in source
    assert r"\sum_{i=1}^{N-1}\mu_i\mu_{2N-i-1}t^{2i}" in source
    assert r"\nu\left(q\mu_{N-1}+q^{-1}\mu_N\right)t^N" in source
    assert r"-\nu^2\mu_{N-1}\mu_Nt^{2N}" in source
    theory = load_theory(THEORY)
    assert theory.source_references[0].equation == "12.3"
    assert [factor.cartan_name for factor in theory.simple_factors] == ["A6", "A1"]
    assert len(theory.pe.terms) == len(theory.rational_product.factors) == 8
    assert [int(x.coefficient) for x in theory.pe.terms] == [1] * 7 + [-1]


def test_exact_dimensions_and_conjugations():
    theory = load_theory(THEORY)
    a6, a1 = theory.simple_factors
    dimensions = {
        (0, 0, 0, 0, 0, 0): 1, (1, 0, 0, 0, 0, 1): 48,
        (0, 0, 1, 0, 0, 0): 35, (0, 0, 0, 1, 0, 0): 35,
        (0, 1, 0, 0, 1, 0): 392, (0, 0, 1, 1, 0, 0): 784,
    }
    assert {labels: int(irrep_dimension(a6, labels)) for labels in dimensions} == dimensions
    assert [int(irrep_dimension(a1, (m,))) for m in range(3)] == [1, 2, 3]
    assert tuple(reversed((0, 0, 1, 0, 0, 0))) == (0, 0, 0, 1, 0, 0)
    assert all(tuple(reversed((m,))) == (m,) for m in range(3))


def test_expansion_pl_and_reconstruction_gates():
    assert _json(OUT / "checks.json")["validation_results"]["pe_equals_rational_product"]
    hwg = _json(OUT / "hwg_expansion.json")["coefficients_by_t_degree"]
    assert set(hwg) == {"0", "2", "4", "6", "8", "10"}
    pl = _json(OUT / "refined_plethystic_logarithm.json")["coefficients_by_t_degree"]
    assert set(pl) == {"2", "4", "6", "8", "10"}
    def native(entries):
        return {(int(e["abelian_charges"]["q"]),
                 tuple(tuple(r["dynkin_labels"]) for r in e["irreducible_representations"])): QQ(e["coefficient"])
                for e in entries}
    z6, z1 = (0,) * 6, (0,)
    assert native(pl["2"]) == {(0, (z6, z1)): 1, (0, ((1, 0, 0, 0, 0, 1), z1)): 1,
                                (0, (z6, (2,))): 1}
    assert native(pl["4"]) == {
        (1, ((0, 0, 1, 0, 0, 0), (1,))): 1,
        (-1, ((0, 0, 0, 1, 0, 0), (1,))): 1,
        (0, (z6, z1)): -2, (0, ((1, 0, 0, 0, 0, 1), z1)): -1,
    }
    assert _json(OUT / "plethystic_logarithm_checks.json")["validation_results"]["direct_scalar_matches_refined_unrefinement"]
    qpl = _json(OUT / "q_refined_dimension_pl.json")["coefficients_by_t_degree"]
    assert qpl["2"] == [{"abelian_charges": {"q": "0"}, "coefficient": 52}]
    assert {(int(e["abelian_charges"]["q"]), e["coefficient"]) for e in qpl["4"]} == {(-1, 70), (0, -50), (1, 70)}
    assert _json(OUT / "unrefined_plethystic_logarithm.json")["coefficients_by_t_degree"] == {"2": 52, "4": 90, "6": -1520, "8": 4585, "10": 102164}
    assert _json(OUT / "reconstruction_difference.json")["mismatch_count"] == 0
    assert _json(OUT / "reconstruction_checks.json")["validation_results"]["all_passed"]


def test_raw_branching_exact_low_degrees_and_all_checks():
    assert branch_a1_to_u1((0,)) == (0,)
    assert branch_a1_to_u1((1,)) == (1, -1)
    assert branch_a1_to_u1((2,)) == (2, 0, -2)
    index = _index(RAW / "branched_refined_plethystic_logarithm.json")
    z = (0,) * 6
    assert {k: v for k, v in index.items() if k[0] == 2} == {
        (2, z, -2, 0): 1, (2, z, 0, 0): 2, (2, z, 2, 0): 1,
        (2, (1, 0, 0, 0, 0, 1), 0, 0): 1,
    }
    assert {k: v for k, v in index.items() if k[0] == 4} == {
        (4, (0, 0, 1, 0, 0, 0), 1, 1): 1, (4, (0, 0, 1, 0, 0, 0), -1, 1): 1,
        (4, (0, 0, 0, 1, 0, 0), 1, -1): 1, (4, (0, 0, 0, 1, 0, 0), -1, -1): 1,
        (4, z, 0, 0): -2, (4, (1, 0, 0, 0, 0, 1), 0, 0): -1,
    }
    checks = _json(RAW / "branching_checks.json")["validation_results"]
    assert checks["all_passed"] and checks["all_parent_dimensions_preserved"]


def test_exact_charge_solution_anchors_lattice_and_physical_branching():
    spec = load_charge_map_spec(ROOT / "theories/charge_maps/su4_nf7_k1o2_manifest_canonical.yaml")
    solution = solve_charge_map(spec)
    expected = matrix(QQ, [[-QQ(7)/4, -QQ(9)/4], [QQ(1)/2, -QQ(1)/2]])
    assert solution.matrix == expected and solution.matrix.rank() == 2 and solution.matrix.det() == 2
    assert solution.inverse_matrix == matrix(QQ, [[-QQ(1)/4, QQ(9)/8], [-QQ(1)/4, -QQ(7)/8]])
    standalone = yaml.safe_load((ROOT / "theories/branching/su4_nf7_k1o2_to_finite.yaml").read_text())
    assert solve_two_anchor_map(standalone["raw_charge_order"], standalone["classical_anchor"], standalone["instanton_anchor"])[0] == expected
    checks = _json(PHYSICAL / "charge_map_checks.json")["validation_results"]
    assert checks["all_passed"] and checks["all_defining_anchors_exist_in_raw_data"]
    index = _index(PHYSICAL / "physical_branched_refined_plethystic_logarithm.json", True)
    z = (0,) * 6
    assert {k: v for k, v in index.items() if k[0] == 2} == {
        (2, z, -QQ(7)/2, 1): 1, (2, z, 0, 0): 2, (2, z, QQ(7)/2, -1): 1,
        (2, (1, 0, 0, 0, 0, 1), 0, 0): 1,
    }
    assert all((2 * B) in ZZ and I in ZZ for (_, _, B, I) in index)
    for (degree, labels, B, I), multiplicity in index.items():
        assert index[(degree, tuple(reversed(labels)), -B, -I)] == multiplicity


def test_frozen_references_are_still_regression_cases():
    assert _json(ROOT / "generated/su4_nf7_finite/order_10/reconstruction_checks.json")["validation_results"]["all_passed"]
    assert _json(ROOT / "generated/su4_nf7_k3o2_infinite/order_10/reconstruction_checks.json")["validation_results"]["all_passed"]
    preflight = _json(ROOT / "generated/branching_convention_preflight/su4_nf7_signed_branching/preflight_results.json")
    assert preflight["all_passed"] and {x["theory_id"] for x in preflight["theories"]} == {"su4_nf7_k1o2_infinite", "su4_nf7_k3o2_infinite"}
