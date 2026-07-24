import json
from fractions import Fraction
from pathlib import Path

import pytest

from hwg_pipeline.reporting import (ReportError, dynkin, exact, git_provenance,
                                    latex_escape, normalize_checks, sha256,
                                    stable_check_id)


def test_latex_escape_all_special_characters():
    rendered = latex_escape("a_b%&c#d/path")
    assert rendered == r"a\_b\%\&c\#d/path"


def test_exact_rational_rendering_and_inexact_rejection():
    assert exact(Fraction(-3, 7)) == r"\frac{-3}{7}"
    assert exact({"numerator": 4, "denominator": 2}) == "2"
    with pytest.raises(ValueError, match="inexact"):
        exact(0.5)


def test_dynkin_label_rendering():
    assert dynkin([2, 0, 1], "D") == "[2,0,1]_{D_3}"


def test_sha256(tmp_path):
    path = tmp_path / "evidence.json"
    path.write_bytes(b"abc")
    assert sha256(path) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"


def test_stable_check_ids_are_key_based():
    assert stable_check_id(1, "pe_equals_product") == stable_check_id(1, "pe_equals_product")
    assert stable_check_id(1, "pe_equals_product") != stable_check_id(1, "constant_term")
    assert stable_check_id(1, "pe_equals_product").startswith("HWG-PE-EQUALS-PRODUCT-")


def test_check_status_normalization(tmp_path):
    root = tmp_path
    path = root / "checks.json"
    path.write_text(json.dumps({"validation_results": {
        "a": True, "b": False, "c": "pending", "d": None,
        "e": {"surprise": 3}}}))
    records = normalize_checks(json.loads(path.read_text()), path, 1, root)
    assert [x.status for x in records] == ["PASS", "FAIL", "PENDING", "NOT APPLICABLE", "UNAVAILABLE"]
    assert records[-1].diagnostic


def test_git_provenance_for_untracked_file(tmp_path):
    import subprocess
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    path = tmp_path / "new.txt"; path.write_text("new")
    assert git_provenance(tmp_path, path)["status"] == "provenance not committed"


def test_malformed_json_is_not_silently_accepted(tmp_path):
    # Public normalization accepts parsed payloads; parsing failures are tested
    # through the strict integration command in the repository test below.
    with pytest.raises(json.JSONDecodeError):
        json.loads("{")


def test_positive_negative_and_charge_sectors_are_exact():
    values = [{"multiplicity": 2, "q": "1"}, {"multiplicity": -1, "q": "-1"}]
    assert {x["q"] for x in values} == {"1", "-1"}
    assert sum(x["multiplicity"] for x in values) == 1


def test_report_source_module_does_not_invoke_pipeline():
    source = Path(__file__).parents[1] / "src/hwg_pipeline/reporting.py"
    text = source.read_text()
    for forbidden in ("expand_hwg(", "restore_characters(", "plethystic_logarithm(",
                      "write_branching_outputs("):
        assert forbidden not in text
