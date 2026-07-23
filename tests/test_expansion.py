"""Hand-checkable tests for exact sparse highest-weight expansions."""
from dataclasses import FrozenInstanceError, replace
import pytest
from sage.all import QQ, ZZ
from hwg_pipeline import (HWGTerm, HighestWeightMonomial, PlethysticExponentialSpec,
    RationalProductFactor, RationalProductSpec, RepresentationSpec, SparseSeries,
    expand_pe, expand_rational_product, load_theory, unit_monomial)

@pytest.fixture
def base(): return load_theory("theories/template.yaml")

def mono(base, degree, labels=(0, 0), charge=0):
    return HighestWeightMonomial(degree, (RepresentationSpec("flavor", labels),),
                                 (("topological", QQ(charge)),))

def with_pe(base, terms):
    return replace(base, pe=PlethysticExponentialSpec(tuple(terms), "test"))

def sig(series):
    return [(int(m.t_degree), tuple(int(x) for x in m.representations[0].dynkin_labels),
             m.abelian_charges[0][1], int(c)) for m, c in series]

def test_pe_t2(base):
    assert sig(expand_pe(with_pe(base, [HWGTerm(1, mono(base, 2))]), 6)) == [
        (0,(0,0),QQ(0),1), (2,(0,0),QQ(0),1), (4,(0,0),QQ(0),1), (6,(0,0),QQ(0),1)]

def test_pe_two_t2(base):
    assert [int(c) for _, c in expand_pe(with_pe(base, [HWGTerm(2, mono(base, 2))]), 6)] == [1,2,3,4]

def test_pe_t2_minus_t4(base):
    theory = with_pe(base, [HWGTerm(1, mono(base,2)), HWGTerm(-1, mono(base,4))])
    assert sig(expand_pe(theory,8)) == [(0,(0,0),QQ(0),1),(2,(0,0),QQ(0),1)]

def test_two_distinct_generators(base):
    theory=with_pe(base,[HWGTerm(1,mono(base,2,(1,0))), HWGTerm(1,mono(base,3,(0,1)))])
    assert sig(expand_pe(theory,5)) == [(0,(0,0),QQ(0),1),(2,(1,0),QQ(0),1),
        (3,(0,1),QQ(0),1),(4,(2,0),QQ(0),1),(5,(1,1),QQ(0),1)]

def test_generator_relation_cancels(base):
    x=mono(base,2,(1,0)); theory=with_pe(base,[HWGTerm(1,x),HWGTerm(-1,x**2)])
    assert sig(expand_pe(theory,8)) == [(0,(0,0),QQ(0),1),(2,(1,0),QQ(0),1)]

def test_multiplication_order_and_immediate_truncation(base):
    unit=unit_monomial(base)
    a=SparseSeries(((unit,1),(mono(base,6),2)),5)
    b=SparseSeries(((unit,1),(mono(base,2,(1,0)),1)),5)
    assert a*b == b*a == b
    assert all(m.t_degree <= 5 for m,_ in a*b)

def test_immutability(base):
    series=SparseSeries.unit(unit_monomial(base),4)
    with pytest.raises(FrozenInstanceError): series.max_degree=8
    assert series + series != series

def test_degree_stability(base):
    theory=with_pe(base,[HWGTerm(1,mono(base,2,(1,0))),HWGTerm(1,mono(base,3))])
    assert expand_pe(theory,7) == expand_pe(theory,8).truncate(7)

def test_reject_degree_zero_and_cancel_zero(base):
    with pytest.raises(ValueError,match="degree-zero"):
        expand_pe(with_pe(base,[HWGTerm(1,mono(base,0,(1,0)))]),4)
    m=mono(base,2)
    assert len(SparseSeries(((m,3),(m,-3)),4)) == 0

def test_pe_equals_product_physical_fixture():
    theory=load_theory("theories/su3_nf5_k3o2_infinite.yaml")
    pe=expand_pe(theory,10); product=expand_rational_product(theory,10)
    assert pe == product
    assert next(c for m,c in pe if m.t_degree == 0) == ZZ(1)
    assert all(m.t_degree <= 10 and c in ZZ and c >= 0 for m,c in pe)

def test_direct_signed_product_powers(base):
    m=mono(base,2)
    theory=replace(base,rational_product=RationalProductSpec((RationalProductFactor(m,-2),
        RationalProductFactor(m**2,1)),"test"))
    assert [int(c) for _,c in expand_rational_product(theory,6)] == [1,2,2,2]
