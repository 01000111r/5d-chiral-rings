"""Read-only verification of the frozen-input k=5/2 report layer."""
import json
from pathlib import Path

ROOT = Path(__file__).parents[1]
BASE = ROOT / "generated/su4_nf7_k5o2_infinite/order_10"


def load(relative):
    return json.loads((BASE / relative).read_text())


def test_compact_report_is_native_and_frozen_exact():
    tex = (BASE / "compact_report/compact_results.tex").read_text()
    checks = load("compact_report/compact_results_checks.json")["validation_results"]
    assert all(checks.values())
    assert r"\mathrm{SO}(14)\times\mathrm{U}(1)_{q}" in tex
    assert "Rational-product form" in tex and "q" in tex
    assert "SU(7)" not in tex and "B,I" not in tex
    assert "124247t^{10}" in tex


def test_comparison_is_complete_and_accounted_degree_by_degree():
    comparison = load("branching_comparison/finite_uv_comparison.json")
    assert comparison["statement"] == (
        "This compares representation channels in two different coordinate rings; "
        "it does not assert that the finite plethystic logarithm equals the I=0 "
        "sector of the UV plethystic logarithm.")
    assert [row["degree"] for row in comparison["degree_summary"]] == [2, 4, 6, 8, 10]
    for row in comparison["degree_summary"]:
        assert row["native_finite_terms"] == sum(row[k] for k in (
            "finite_exact_matches", "finite_different_multiplicity",
            "finite_different_sign", "finite_absent"))
        assert row["combined_physical_uv_sectors"] == sum(row["physical_category_counts"].values())
    assert sum(comparison["summary"].values()) == len(comparison["finite_terms"])
    finite_keys = {(x["degree"], tuple(x["labels"]), str(x["B"]), x["I"])
                   for x in comparison["finite_terms"]}
    assert all((x["degree"], tuple(x["labels"]), str(x["B"]), x["I"])
               not in finite_keys for x in comparison["uv_only_channels"])


def test_low_degrees_and_required_conventions_are_rendered():
    comparison = load("branching_comparison/finite_uv_comparison.json")
    by = {(x["degree"], tuple(x["labels"]), str(x["B"])): x
          for x in comparison["finite_terms"]}
    assert by[(2, (1, 0, 0, 0, 0, 1), "0")]["status"] == "exact-match"
    assert by[(2, (0, 0, 0, 0, 0, 0), "0")]["status"] == "representation-match-different-multiplicity"
    for labels, baryon in (((0, 0, 1, 0, 0, 0), "-4"),
                            ((0, 0, 0, 1, 0, 0), "4"),
                            ((1, 0, 0, 0, 0, 1), "0"),
                            ((0, 0, 0, 0, 0, 0), "0")):
        assert by[(4, labels, baryon)]["status"] == "exact-match"
    tex = (BASE / "branching_comparison/branching_comparison.tex").read_text()
    assert r"D_7\to A_6\times U(1)_x" in tex
    assert r"x=-2\sum_iw_i" in tex
    assert "Four exact D7 branching benchmarks" in tex
    assert r"91=21+48+1+21" in tex and r"104=28+48+28" in tex
    assert r"B(E_-)=-N-k=-3/2" in tex and r"r=2" in tex
    assert r"\frac38&\frac{35}{8}" in tex
    assert "rank $2$ and determinant $1$" in tex
