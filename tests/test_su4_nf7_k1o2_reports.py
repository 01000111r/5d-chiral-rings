"""Semantic checks for the frozen SU(4)+7F comparison reports."""

import json
from pathlib import Path


ROOT = Path(__file__).parents[1]
BASE = ROOT / "generated/su4_nf7_k1o2_infinite/order_10"


def load(name):
    return json.loads((BASE / name).read_text())


def test_compact_report_is_native_and_matches_persisted_data():
    tex = (BASE / "compact_report/compact_results.tex").read_text()
    checks = load("compact_report/compact_results_checks.json")["validation_results"]
    stored = load("refined_plethystic_logarithm.json")["coefficients_by_t_degree"]
    assert "Rational-product form" in tex
    assert "[0,0,0,0,0,0;2]" in tex
    assert "_{A6}" in tex and "_{A1" in tex and "U}(1)_{q}" in tex
    assert checks["all_passed"] and checks["negative_pl_coefficients_retained"]
    assert sum(len(v) for v in stored.values()) == load("compact_report/compact_results_manifest.json")["term_counts"]["refined_character_plethystic_logarithm"]


def test_comparison_partitions_every_channel_and_preserves_signs():
    comparison = load("branching_comparison/finite_uv_comparison.json")
    physical = load("branching_comparison/physical_branching.json")["combined_by_degree"]
    finite = comparison["finite_terms"]
    assert len(finite) == sum(row["native_finite_terms"] for row in comparison["degree_summary"])
    assert len(physical) == sum(row["combined_physical_uv_sectors"] for row in comparison["degree_summary"])
    assert len(comparison["uv_only_channels"]) == sum(row["uv_only_channels"] for row in comparison["degree_summary"])
    assert all(term["status"] in {"exact-match", "representation-match-different-multiplicity", "representation-match-different-sign", "absent"} for term in finite)
    assert any(term["signed_multiplicity"] < 0 for term in finite)
    categories = comparison["physical_sector_counts"]
    assert sum(categories[name] for name in ("neutral", "classical_baryonic", "pure_instanton", "mixed_baryon_instanton")) == len(physical)


def test_report_charge_map_and_frozen_raw_branch_agree_exactly():
    comparison = load("branching_comparison/finite_uv_comparison.json")
    charge_map = load("branching_comparison/charge_map.json")
    assert charge_map["solution_matrix"] == [["-7/4", "-9/4"], ["1/2", "-1/2"]]
    assert charge_map["inverse_matrix"] == [["-1/4", "9/8"], ["-1/4", "-7/8"]]
    checks = load("branching_comparison/branching_checks.json")["checks"]
    statuses = [value for group in checks.values() for value in group.values()
                if value in {"pass", "fail", "pending", "unavailable"}]
    assert "fail" not in statuses and "pending" not in statuses
    assert checks["presentation"]["latex_compile"] in {"pass", "unavailable"}
    statement = comparison["statement"]
    assert statement == "This compares representation channels in two different coordinate rings; it does not assert that the finite plethystic logarithm equals the I=0 sector of the UV plethystic logarithm."
