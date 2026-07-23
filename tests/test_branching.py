"""Understandable exact tests for generic and fixture branching."""

import pytest
from pathlib import Path
from sage.all import ZZ

from hwg_pipeline.branching import (EMBEDDING, branch_irrep,
    branch_representation_content, load_branching_spec)
from hwg_pipeline.characters import RepresentationContent
from hwg_pipeline.model import SimpleGroupSpec
from hwg_pipeline.sage_backend import irrep_dimension, tensor_product


def group(rank, name="g"):
    return SimpleGroupSpec(name, "A", rank, f"SU({rank+1})",
                           tuple(f"u{i}" for i in range(rank)))


def terms(rank, labels):
    return {(tuple(map(int, p.child_dynkin_labels)), int(p.x_charge)): int(p.multiplicity)
            for p in branch_irrep(group(rank), group(rank-1,"h"), labels, EMBEDDING)}


def test_a2_fundamental_and_conjugate():
    assert terms(2, (1,0)) == {((1,),1):1, ((0,),-2):1}
    assert terms(2, (0,1)) == {((1,),-1):1, ((0,),2):1}


@pytest.mark.parametrize("labels", [(1,0),(0,1),(1,1),(2,1),(0,3)])
def test_a2_dimensions_are_preserved(labels):
    pieces = branch_irrep(group(2), group(1,"h"), labels)
    assert irrep_dimension(group(2), labels) == sum(
        p.multiplicity*irrep_dimension(group(1,"h"),p.child_dynkin_labels) for p in pieces)


def test_direct_sums_and_integer_scalars_preserve_external_charge():
    parent = group(2); child = group(1,"h")
    content = RepresentationContent((parent,), [(((1,0),),2), (((0,1),),3)])
    spec = type("Spec",(),{"child_group":child,"embedding_type":EMBEDDING})()
    branched = branch_representation_content(content,spec,7,(("q",ZZ(4)),))
    assert all(x.t_degree == 7 and x.abelian_charges == (("q",ZZ(4)),) for x in branched)
    assert sum(x.multiplicity*irrep_dimension(child,x.child_dynkin_labels) for x in branched) == 2*3+3*3


def _child_tensor(left, right, child):
    result={}
    for a in left:
        for b in right:
            for labels,m in tensor_product(child,a.child_dynkin_labels,b.child_dynkin_labels):
                key=tuple(labels),a.x_charge+b.x_charge
                result[key]=result.get(key,0)+a.multiplicity*b.multiplicity*m
    return result


def test_branching_commutes_with_tensor_products():
    parent,child=group(2),group(1,"h")
    left=branch_irrep(parent,child,(1,0)); right=branch_irrep(parent,child,(0,1))
    expected=_child_tensor(left,right,child)
    actual={}
    for labels,m in tensor_product(parent,(1,0),(0,1)):
        for p in branch_irrep(parent,child,labels):
            key=p.child_dynkin_labels,p.x_charge
            actual[key]=actual.get(key,0)+m*p.multiplicity
    assert actual == expected


def test_output_is_deterministic_and_invalid_embeddings_are_rejected():
    first=branch_irrep(group(2),group(1,"h"),(2,1)); assert first==branch_irrep(group(2),group(1,"h"),(2,1))
    assert first == tuple(sorted(first,key=lambda x:(x.x_charge,x.child_dynkin_labels)))
    with pytest.raises(ValueError): branch_irrep(group(2),group(1,"h"),(1,0),"invalid")
    with pytest.raises(ValueError): branch_irrep(group(1),SimpleGroupSpec("z","A",1,"",("z",)),(1,))


@pytest.mark.parametrize("labels, expected", [
 ((1,0,0,0,0),{((1,0,0,0),1):1,((0,0,0,0),-5):1}),
 ((0,0,0,0,1),{((0,0,0,1),-1):1,((0,0,0,0),5):1}),
 ((1,0,0,0,1),{((1,0,0,1),0):1,((0,0,0,0),0):1,((1,0,0,0),6):1,((0,0,0,1),-6):1}),
 ((0,1,0,0,0),{((0,1,0,0),2):1,((1,0,0,0),-4):1}),
 ((0,0,0,1,0),{((0,0,1,0),-2):1,((0,0,0,1),4):1}),])
def test_a5_conventions(labels,expected): assert terms(5,labels)==expected


def test_fixture_specification():
    root=Path(__file__).resolve().parents[1]
    spec=load_branching_spec(root/"theories/branchings/su3_nf5_k3o2_to_manifest.yaml")
    assert spec.raw_branching_u1_name=="x" and spec.preserved_abelian_factors==("q",)
