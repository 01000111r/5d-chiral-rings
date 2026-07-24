"""Focused formatting-only tests for compact stored-result reports."""

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path

from hwg_pipeline.compact_latex import (
    _render_dimension_series, _render_grouped, _render_scalar, _term_chunks,
    generate_compact_latex, render_dynkin, render_exact,
    render_highest_weight, render_laurent_coefficient, render_q_power,
)


def _rep(labels, q, coefficient=1, key="multiplicity"):
    return {"abelian_charges": {"q": str(q)},
            "irreducible_representations": [
                {"cartan_factor_id": "enhanced", "dynkin_labels": labels}],
            key: coefficient}


def test_exact_integer_and_rational_rendering():
    assert render_exact(-3) == "-3"
    assert render_exact({"numerator": 3, "denominator": 2}) == r"\frac{3}{2}"


def test_highest_weight_monomial_rendering():
    assert render_highest_weight([2, 0, 0, 1, 0], -1) == r"q^{-1}\mu_{1}^{2}\mu_{4}"


def test_dynkin_label_and_q_power_rendering():
    assert render_dynkin([1, 0, 2, 0, 1]) == "[1,0,2,0,1]"
    assert (render_q_power(-2), render_q_power(0), render_q_power(1)) == (r"q^{-2}", "", "q")


def _dimension(q, coefficient):
    return {"abelian_charges": {"q": str(q)}, "coefficient": coefficient}


def test_equal_dimension_conjugates_retain_q_and_factor_compactly():
    # Synthetic dimensions of [1,0]q and [0,1]q^-1 are both d=3.
    assert render_laurent_coefficient([_dimension(1, 3), _dimension(-1, 3)]) == (
        r"3\bigl(q+q^{-1}\bigr)")


def test_asymmetric_dimension_sectors_are_not_incorrectly_combined():
    rendered = render_laurent_coefficient([_dimension(-1, 2), _dimension(1, 3)])
    assert rendered == r"2q^{-1}\,+3q"


def test_positive_and_negative_higher_q_powers_survive():
    assert render_laurent_coefficient([_dimension(-2, 5), _dimension(2, 5)]) == (
        r"5\bigl(q^{2}+q^{-2}\bigr)")


def test_full_unrefinement_signed_polynomial():
    rendered = _render_scalar("P(t)", {"0": 1, "1": 0, "2": -3, "3": 2}, 3)
    assert r"1\,-3t^{2}\,+2t^{3}" in rendered
    assert r"O(t^{4})" in rendered


def test_dimension_series_contains_no_dynkin_labels_and_does_not_set_q_to_one():
    rendered = _render_dimension_series("D", {"3": [_dimension(-1, 15), _dimension(1, 15)]}, 3)
    assert r"15\bigl(q+q^{-1}\bigr)\Bigr)t^{3}" in rendered
    assert "[" not in rendered and "q=1" not in rendered


def test_signed_dimension_pl_survives():
    rendered = _render_dimension_series("P", {"4": [_dimension(0, -36)]}, 4)
    assert r"-36t^{4}" in rendered


def test_signed_pl_terms_are_retained():
    groups = {"2": [_rep([0, 0, 0, 0, 0], 0, -2, "coefficient"),
                     _rep([1, 0, 0, 0, 1], 1, 1, "coefficient")]}
    rendered = _render_grouped("P", groups, 2, "character")
    assert "-2[0,0,0,0,0]" in rendered and "+q[1,0,0,0,1]" in rendered


def test_degree_grouping_and_complete_cutoff():
    groups = {"0": [_rep([0, 0, 0, 0, 0], 0)], "3": [_rep([0, 1, 0, 0, 0], 1)]}
    rendered = _render_grouped("H", groups, 3, "character")
    assert r"t^{3}" in rendered and r"O(t^{4})" in rendered


def test_long_coefficients_wrap_only_between_terms():
    groups = {"2": [_rep([i, 0, 0, 0, 0], 0) for i in range(6)]}
    rendered = _render_grouped("H", groups, 2, "character", chunk=2)
    assert rendered.count(r"&{}") == 2
    assert rendered.count(r"t^{2}") == 1


def test_width_aware_wrapping_never_splits_a_mathematical_term():
    terms = [r"q^{-1}[1,0,0,0,1]", r"+12q[0,1,0,1,0]", r"-7[2,0,0,0,2]"]
    chunks = _term_chunks(terms, maximum_terms=9, maximum_width=32)
    assert all(len(chunk) <= 32 for chunk in chunks)
    assert r"\,".join(chunks) == r"\,".join(terms)


def test_scalar_polynomial_wraps_within_configured_width():
    rendered = _render_scalar("H(t)", {str(d): d * 12345 for d in range(1, 11)}, 10,
                              maximum_width=38)
    assert rendered.count(r"\\") >= 3
    assert r"O(t^{11})" in rendered


def test_synthetic_rendering_is_deterministic_and_does_not_mutate():
    groups = {"2": [_rep([1, 0, 0, 0, 1], 0)]}
    original = deepcopy(groups)
    assert _render_grouped("H", groups, 2, "character") == _render_grouped("H", groups, 2, "character")
    assert groups == original


def test_generated_document_has_no_tables_and_all_manifest_terms(repo_root):
    generate_compact_latex(repo_root, "su3_nf5_k3o2_infinite", 10)
    report = repo_root / "generated/su3_nf5_k3o2_infinite/order_10/compact_report"
    tex = (report / "compact_results.tex").read_text()
    manifest = json.loads((report / "compact_results_manifest.json").read_text())
    assert "longtable" not in tex and "tabular" not in tex and "booktabs" not in tex
    source = json.loads((repo_root / "generated/su3_nf5_k3o2_infinite/order_10/character_series.json").read_text())
    assert manifest["term_counts"]["refined_character_hilbert_series"] == sum(
        len(x) for x in source["coefficients_by_t_degree"].values())
    assert r"O(t^{11})" in tex
    assert r"H_{\mathrm{dim}}(t,q)" in tex and r"\operatorname{PL}_{\mathrm{dim}}(t,q)" in tex
    assert r"\resizebox{0.98\linewidth}{!}" in tex
    assert "q=1" not in tex
    assert "$q$-unrefined" not in tex


def test_generation_does_not_mutate_stored_results(repo_root):
    directory = repo_root / "generated/su3_nf5_k3o2_infinite/order_10"
    paths = sorted(directory.glob("*.json"))
    before = {p: sha256(p.read_bytes()).hexdigest() for p in paths}
    generate_compact_latex(repo_root, "su3_nf5_k3o2_infinite", 10)
    assert before == {p: sha256(p.read_bytes()).hexdigest() for p in paths}


def pytest_generate_tests(metafunc):
    if "repo_root" in metafunc.fixturenames:
        metafunc.parametrize("repo_root", [Path(__file__).resolve().parents[1]])
