import pytest
from sage.all import QQ, matrix

from hwg_pipeline.charge_maps import (ChargeAnchor, ChargeMapSpec, ChargeVector,
    InconsistentChargeMapError, apply_charge_map, apply_charge_map_to_records,
    rational_json, solve_charge_map)


def cv(names, values): return ChargeVector(tuple(names), tuple(values))
def anchor(name, raw, physical, rn=("x","q"), pn=("B","I")):
    return ChargeAnchor(name, cv(rn,raw), cv(pn,physical))
def spec(*anchors, rn=("x","q"), pn=("B","I")):
    return ChargeMapSpec("test",tuple(rn),tuple(pn),tuple(anchors))


def test_unique_square_exact_solution_and_residuals():
    solved=solve_charge_map(spec(anchor("a",(1,0),(2,3)),anchor("b",(0,1),(4,5))))
    assert solved.matrix == matrix(QQ,[[2,4],[3,5]])
    assert solved.diagnostics.unique and all(not any(x) for x in solved.diagnostics.defining_residuals)


def test_underdetermined_has_nullspace_and_no_map():
    solved=solve_charge_map(spec(anchor("a",(1,0),(2,3))))
    assert solved.diagnostics.consistent and not solved.diagnostics.unique
    assert solved.diagnostics.nullspace and solved.matrix is None


def test_inconsistent_reports_ranks():
    with pytest.raises(InconsistentChargeMapError) as error:
        solve_charge_map(spec(anchor("a",(1,0),(2,3)),anchor("b",(1,0),(9,3))))
    assert error.value.diagnostics.coefficient_rank < error.value.diagnostics.augmented_rank


def test_consistent_overdetermined_retains_unique_solution():
    solved=solve_charge_map(spec(anchor("a",(1,0),(2,3)),anchor("b",(0,1),(4,5)),anchor("sum",(1,1),(6,8))))
    assert solved.diagnostics.unique and solved.matrix == matrix(QQ,[[2,4],[3,5]])


def test_rational_serialization_inverse_and_roundtrip():
    solved=solve_charge_map(spec(anchor("a",(2,0),(1,0)),anchor("b",(0,3),(0,1))))
    assert solved.matrix == matrix(QQ,[[QQ(1)/2,0],[0,QQ(1)/3]])
    assert rational_json(QQ(1)/3) == {"numerator":1,"denominator":3}
    raw=cv(("x","q"),(QQ(-3)/2,QQ(7)/4)); physical=apply_charge_map(solved,raw)
    from sage.all import vector
    assert tuple(solved.inverse_matrix*vector(QQ,physical.values)) == raw.values


def test_series_application_combines_and_retains_provenance_deterministically():
    solved=solve_charge_map(spec(anchor("a",(1,0),(1,0)),anchor("b",(0,1),(0,1))))
    base={"t_degree":2,"child_dynkin_labels":[1,0],"child_representation_dimension":3,"signed_multiplicity":1}
    records=[{**base,"raw_charges":{"x":"-1/2","q":"1/3"}}, {**base,"raw_charges":{"x":{"numerator":-1,"denominator":2},"q":"1/3"}}]
    result=apply_charge_map_to_records(solved,records)
    assert len(result)==1 and result[0]["t_degree"]==2 and result[0]["child_dynkin_labels"]==[1,0]
    assert result[0]["signed_multiplicity"]==rational_json(2) and len(result[0]["provenance"])==2
    assert result == apply_charge_map_to_records(solved,list(reversed(records)))


def test_more_than_two_charge_factors():
    rn=("a","b","c"); pn=("u","v","w")
    anchors=[anchor(str(i),tuple(1 if i==j else 0 for j in range(3)),tuple((i+1) if i==j else 0 for j in range(3)),rn,pn) for i in range(3)]
    solved=solve_charge_map(spec(*anchors,rn=rn,pn=pn))
    assert solved.matrix == matrix(QQ,[[1,0,0],[0,2,0],[0,0,3]])
