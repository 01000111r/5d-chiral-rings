"""Exact benchmarks for the raw SU(8) -> SU(7) x U(1)_x branching."""

import json
from pathlib import Path

from sage.all import QQ, ZZ

from hwg_pipeline.branching import branch_irrep, load_branching_spec
from hwg_pipeline.model import SimpleGroupSpec
from hwg_pipeline.sage_backend import irrep_dimension


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "generated/su4_nf7_k3o2_infinite/order_10"
SPEC = load_branching_spec(ROOT / "theories/branchings/su4_nf7_k3o2_to_manifest.yaml")
A7 = SimpleGroupSpec("uv", "A", 7, "SU(8)", tuple(f"u{i}" for i in range(7)))
A6 = SPEC.child_group


def pieces(labels):
    return {(tuple(map(int, p.child_dynkin_labels)), int(p.x_charge)): int(p.multiplicity)
            for p in branch_irrep(A7, A6, labels)}


def test_specification_and_independent_low_irrep_anchors():
    assert (SPEC.parent_simple_factor, SPEC.child_simple_factor, SPEC.child_rank) == ("A7", "A6", 6)
    assert SPEC.raw_branching_u1_name == "x" and SPEC.preserved_abelian_factors == ("q",)
    expected = {
        (1, 0, 0, 0, 0, 0, 0): {((1, 0, 0, 0, 0, 0), 1): 1, ((0,)*6, -7): 1},
        (1, 0, 0, 0, 0, 0, 1): {((1, 0, 0, 0, 0, 1), 0): 1, ((1, 0, 0, 0, 0, 0), 8): 1,
                                      ((0, 0, 0, 0, 0, 1), -8): 1, ((0,)*6, 0): 1},
        (0, 0, 1, 0, 0, 0, 0): {((0, 0, 1, 0, 0, 0), 3): 1, ((0, 1, 0, 0, 0, 0), -5): 1},
        (0, 0, 0, 0, 1, 0, 0): {((0, 0, 0, 0, 1, 0), 5): 1, ((0, 0, 0, 1, 0, 0), -3): 1},
    }
    expected_dimensions = {tuple(expected): total for expected, total in zip(expected, (8, 63, 56, 56))}
    for parent_labels, children in expected.items():
        assert pieces(parent_labels) == children
        assert irrep_dimension(A7, parent_labels) == expected_dimensions[parent_labels]
        assert sum(m * irrep_dimension(A6, labels) for (labels, _), m in children.items()) == expected_dimensions[parent_labels]


def _branched_pl():
    payload = json.loads((SOURCE / "manifest_branching/branched_refined_plethystic_logarithm.json").read_text())
    return payload, payload["coefficients_by_t_degree"]


def _key(entry, coefficient="coefficient"):
    return ((tuple(entry["child_dynkin_labels"]), int(entry["raw_charges"]["x"]),
             int(entry["raw_charges"]["q"])), QQ(entry[coefficient]))


def test_complete_t2_and_t4_pl_benchmarks():
    _, pl = _branched_pl()
    assert dict(_key(e) for e in pl["2"]) == {
        ((0,)*6, 0, 0): 2, ((1, 0, 0, 0, 0, 1), 0, 0): 1,
        ((1, 0, 0, 0, 0, 0), 8, 0): 1, ((0, 0, 0, 0, 0, 1), -8, 0): 1}
    assert dict(_key(e) for e in pl["4"]) == {
        ((0, 0, 1, 0, 0, 0), 3, 1): 1, ((0, 1, 0, 0, 0, 0), -5, 1): 1,
        ((0, 0, 0, 0, 1, 0), 5, -1): 1, ((0, 0, 0, 1, 0, 0), -3, -1): 1,
        ((0,)*6, 0, 0): -2, ((1, 0, 0, 0, 0, 1), 0, 0): -1,
        ((1, 0, 0, 0, 0, 0), 8, 0): -1, ((0, 0, 0, 0, 0, 1), -8, 0): -1}


def test_complete_outputs_and_independent_consistency_checks():
    payload, by_degree = _branched_pl()
    checks = json.loads((SOURCE / "manifest_branching/branching_checks.json").read_text())["validation_results"]
    assert payload["raw_charge_basis"] == "(x,q)" and payload["physical_charge_map_assumed"] is False
    assert set(by_degree) == {"2", "4", "6", "8", "10"}
    assert checks["all_passed"] and checks["character_unrefinement_preserved"]
    assert checks["plethystic_log_unrefinement_preserved"] and checks["plethystic_log_conjugation_consistent"]
    assert checks["all_parent_dimensions_preserved"] and checks["x_charges_are_exact_integers"]
    for entries in by_degree.values():
        index = dict(_key(e) for e in entries)
        assert all(ZZ(e["raw_charges"]["x"]) in ZZ for e in entries)
        for (labels, x, q), coefficient in index.items():
            assert index[(tuple(reversed(labels)), -x, -q)] == coefficient
