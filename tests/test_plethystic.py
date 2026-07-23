import pytest
from sage.all import QQ

from hwg_pipeline.io import load_theory
from hwg_pipeline.expansion import expand_pe
from hwg_pipeline.characters import restore_characters
from hwg_pipeline.model import SimpleGroupSpec
from hwg_pipeline.plethystic import (VirtualCharacterSeries, VirtualRepresentationContent,
    adams_series, formal_logarithm, plethystic_logarithm, unrefine_virtual)


def spec(letter, rank):
    return SimpleGroupSpec("g", letter, rank, "G", tuple(f"m{i}" for i in range(rank)))


def test_sage_adams_conventions_and_virtual_arithmetic():
    from hwg_pipeline.sage_backend import adams_decomposition, irrep_dimension
    a1, a5 = spec("A",1), spec("A",5)
    assert dict(adams_decomposition(a1,(1,),2)) == {(0,):-1,(2,):1}
    assert sum(m*irrep_dimension(a1,l) for l,m in adams_decomposition(a1,(1,),2)) == 2
    assert dict(adams_decomposition(a5,(1,0,0,0,0),2)) == {(0,1,0,0,0):-1,(2,0,0,0,0):1}
    assert sum(m*irrep_dimension(a5,l) for l,m in adams_decomposition(a5,(1,0,0,0,0),2)) == 6
    v=VirtualRepresentationContent.single_irrep((a1,),((1,),))
    one=VirtualRepresentationContent.trivial((a1,))
    assert v*one == v and (v-v).terms == () and (QQ(1)/2*v).terms[0][1] == QQ(1)/2


def fixture_series(order=10):
    theory=load_theory("theories/su3_nf5_k3o2_infinite.yaml")
    return theory, restore_characters(theory,expand_pe(theory,order))


def test_adams_series_scales_degree_charge_and_rejects_bad_k():
    theory, series=fixture_series(3)
    virtual=VirtualCharacterSeries.from_character_series(series,6)
    assert adams_series(virtual,1,6)==virtual
    doubled=adams_series(virtual,2,6)
    assert all(d%2==0 for (d,_),_ in doubled)
    assert {q for (d,c),_ in doubled if d==6 for _,q in c} == {-2,2}
    with pytest.raises(ValueError): adams_series(virtual,0,6)
    with pytest.raises(ValueError): adams_series(virtual,QQ(1)/2,6)


def test_formal_log_has_rationals_and_scalar_plethystic_examples():
    theory,_=fixture_series(2); z=tuple((a.id,0) for a in theory.abelian_factors)
    one=VirtualRepresentationContent.trivial(theory.simple_factors)
    geometric=VirtualCharacterSeries(theory,[((d,z),one) for d in (0,2,4,6)],6)
    log=formal_logarithm(geometric,6)
    assert any(c.denominator()>1 for _,x in log for _,c in x)
    assert unrefine_virtual(plethystic_logarithm(geometric,6)) == ((2,QQ(1)),)
    finite=VirtualCharacterSeries(theory,[((0,z),one),((2,z),one)],6)
    assert unrefine_virtual(plethystic_logarithm(finite,6)) == ((2,QQ(1)),(4,QQ(-1)))
    bad=VirtualCharacterSeries(theory,[((2,z),one)],6)
    with pytest.raises(ValueError): formal_logarithm(bad)


def test_physical_fixture_leading_terms_and_stability():
    _,series=fixture_series(10)
    pl=plethystic_logarithm(series,10)
    assert unrefine_virtual(pl)[:3] == ((2,QQ(36)),(3,QQ(30)),(4,QQ(-36)))
    _,series11=fixture_series(11)
    pl11=plethystic_logarithm(series11,11)
    assert tuple(pl11)[:len(tuple(pl))] == tuple(pl)
    assert all(c.denominator()==1 for _,x in pl for _,c in x)
    assert any(c<0 for _,x in pl for _,c in x)
