import json
from fractions import Fraction
from pathlib import Path

import pytest

from hwg_pipeline.notebook_report import (NotebookError, _series, charge,
    dynkin_labels, exact_rational, git_provenance, markdown_table,
    normalize_checks, sha256, stable_check_id)


def test_exact_rational_rendering():
    assert exact_rational(Fraction(3, 2)) == "3/2"
    assert exact_rational({"numerator": -1, "denominator": 3}) == "-1/3"
    with pytest.raises(ValueError): exact_rational(0.5)


def test_dynkin_and_charge_rendering():
    assert dynkin_labels([1, 0, 1]) == "[1,0,1]_{A3}"
    assert charge(-1) == "q^{-1}"
    assert charge(0) == "1"


def test_markdown_escaping():
    rendered = markdown_table(("a",), (("x|y\nz",),))
    assert "x\\|y<br>z" in rendered


def test_complete_degree_grouping_and_signed_multiplicities():
    payload = {"coefficients_by_t_degree": {"4": [
        {"dynkin_labels": [1, 0], "coefficient": -2}], "2": [
        {"dynkin_labels": [0, 0], "multiplicity": 1},
        {"dynkin_labels": [1, 1], "multiplicity": 3}]}}
    table, count = _series(payload)
    assert count == 3
    assert table.index("| 2 |") < table.index("| 4 |")
    assert "-2" in table and "3" in table


def test_multiple_simple_factors():
    payload={"coefficients_by_t_degree":{"0":[{"irreducible_representations":[
        {"dynkin_labels":[1]},{"dynkin_labels":[0,1]}],"multiplicity":1}]}}
    table,count=_series(payload)
    assert count == 1 and "[1]_{A1} × [0,1]_{A2}" in table


def test_check_normalization_and_stable_ids(tmp_path):
    root=tmp_path; (root/"x.json").write_text(json.dumps({"validation_results":{
        "good":True,"bad":False,"unknown":[1]}}))
    records=normalize_checks(json.loads((root/"x.json").read_text()),root/"x.json",1,root)
    assert [x.status for x in records] == ["FAIL","PASS","UNAVAILABLE"]
    assert records[1].check_id == stable_check_id(1,"good")
    assert records[1].check_id.startswith("NB-HWG-")


def test_hashes(tmp_path):
    path=tmp_path/"x"; path.write_bytes(b"abc")
    assert sha256(path)=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"


def test_untracked_provenance(tmp_path):
    # The repository root used here is real; this test file may itself be
    # uncommitted while under development, which is the exact required case.
    root=Path(__file__).resolve().parents[1]
    path=root/"untracked-provenance-test.tmp"; path.write_text("x")
    try: assert git_provenance(root,path)["status"] == "not committed"
    finally: path.unlink()


def test_malformed_json_contract(tmp_path):
    path=tmp_path/"bad.json"; path.write_text("{")
    with pytest.raises(json.JSONDecodeError): json.loads(path.read_text())


def test_notebook_artifact_is_valid_deterministic_and_complete():
    root=Path(__file__).resolve().parents[1]
    path=root/"notebooks/hwg_pipeline_project_walkthrough.ipynb"
    first=path.read_bytes(); nb=json.loads(first)
    assert nb["nbformat"]==4 and nb["nbformat_minor"]>=0
    assert all("id" not in cell for cell in nb["cells"])
    assert [c["cell_type"] for c in nb["cells"]][:3]==["markdown"]*3
    source="".join("".join(c["source"]) for c in nb["cells"])
    assert "Appendix J" in source and "Complete branched refined PL" in source
    assert not any(ch in source for ch in ("\b","\f","\v"))
    assert path.read_bytes()==first


@pytest.mark.parametrize("stage", ("input","hwg","characters","dimensions",
    "plethystic-log","reconstruction","operator-analysis","branching","charge-map"))
def test_supported_stage_names(stage):
    from hwg_pipeline.notebook_report import STAGES
    assert stage in STAGES

