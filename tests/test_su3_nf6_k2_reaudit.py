"""Exact regression tests for the independently reaudited D6 report."""
import json
from pathlib import Path

import pytest
from sage.all import QQ, WeylCharacterRing

from hwg_pipeline.branching import D5_EMBEDDING, D6_EMBEDDING, branch_irrep, validate_d_type_branching
from hwg_pipeline.model import SimpleGroupSpec

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "generated/su3_nf6_k2_infinite/order_10/branching_comparison"

def group(cartan, rank):
    return SimpleGroupSpec(cartan.lower(), cartan[0], rank, cartan, tuple(f"w{i}" for i in range(rank)))

def pieces(labels):
    return {(tuple(p.child_dynkin_labels), int(p.x_charge), int(p.multiplicity))
            for p in branch_irrep(group("D6", 6), group("A5", 5), labels, D6_EMBEDDING)}

def test_adjoint_and_both_terminal_spinors_exactly_distinguished():
    assert pieces((0,1,0,0,0,0)) == {
        ((0,1,0,0,0),2,1), ((1,0,0,0,1),0,1),
        ((0,0,0,0,0),0,1), ((0,0,0,1,0),-2,1)}
    assert pieces((0,0,0,0,0,1)) == {
        ((1,0,0,0,0),-2,1), ((0,0,1,0,0),0,1), ((0,0,0,0,1),2,1)}
    assert pieces((0,0,0,0,1,0)) == {
        ((0,0,0,0,0),-3,1), ((0,1,0,0,0),-1,1),
        ((0,0,0,1,0),1,1), ((0,0,0,0,0),3,1)}

def test_dimension_compatible_old_fixed_x_conjugation_is_rejected():
    wrong=[((0,0,0,0,1),-2,1),((0,0,1,0,0),0,1),((1,0,0,0,0),2,1)]
    with pytest.raises(ValueError, match="fixed-charge restricted-weight evidence differs"):
        validate_d_type_branching("D6", (0,0,0,0,0,1), wrong)
    A5=WeylCharacterRing("A5",style="coroots")
    assert A5(1,0,0,0,0).degree() == A5(0,0,0,0,1).degree() == 6

def test_all_rules_anchors_map_lattice_and_rendering():
    audit=json.loads((OUT/"branching_accuracy_reaudit.json").read_text())
    assert (audit["rule_count"],audit["failure_count"]) == (96,28)
    assert audit["after_correction"]["pass_count"] == 96
    anchors=json.loads((OUT/"anchor_status_justification.json").read_text())
    assert all(a["verdict"] == "retain" and not a["change_justified"] for a in anchors["anchors"].values())
    cmap=json.loads((OUT/"charge_map.json").read_text())
    assert cmap["solution_matrix"] == [["-1/2","3"],["1/2","0"]]
    assert cmap["inverse"] == {"x":"2*I","q":"B/3+I/3"}
    physical=json.loads((OUT/"physical_branching.json").read_text())
    for parent in physical["parents"]:
        for child in parent["children"]:
            assert QQ(child["physical_charges"]["B"]) in __import__('sage.all',fromlist=['ZZ']).ZZ
            assert QQ(child["physical_charges"]["I"]) in __import__('sage.all',fromlist=['ZZ']).ZZ
    tex=(OUT/"branching_comparison.tex").read_text()
    assert "]_{5}" not in tex and "]_{5;" not in tex

def test_d5_low_branch_is_unchanged_by_d6_fix():
    assert {(tuple(p.child_dynkin_labels),int(p.x_charge)) for p in
      branch_irrep(group("D5",5),group("A4",4),(0,0,0,0,1),D5_EMBEDDING)} == {
        ((0,0,0,0),-5),((0,0,1,0),-1),((1,0,0,0),3)}
