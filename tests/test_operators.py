import json
from pathlib import Path

import pytest
from sage.all import ZZ

from hwg_pipeline.io import load_theory
from hwg_pipeline.operators import (candidate_generators, enumerate_quadratic_channels,
    extract_operator_content, first_negative_degree, first_relation_candidates,
    symmetric_square_decomposition, tensor_decomposition)
from hwg_pipeline.plethystic import VirtualCharacterSeries, VirtualRepresentationContent
from hwg_pipeline.sage_backend import irrep_dimension, symmetric_power, tensor_product


ROOT = Path(__file__).parents[1]


def series(terms):
    theory = load_theory(ROOT / "theories/su3_nf5_k3o2_infinite.yaml")
    built = []
    for degree, charge, labels, multiplicity in terms:
        content = VirtualRepresentationContent.single_irrep(
            theory.simple_factors, (tuple(labels),), multiplicity)
        built.append(((ZZ(degree), (("q", ZZ(charge)),)), content))
    return VirtualCharacterSeries(theory, built, max(x[0] for x in terms))


def test_conservative_classification_mixed_multiplicity_and_order():
    pl = series([(5, 0, [2,0,0,0,0], 1), (2, 0, [1,0,0,0,1], 2),
                 (4, 0, [0,0,0,0,0], -1), (6, 0, [0,1,0,0,0], -1),
                 (4, 1, [0,0,0,1,0], 1)])
    assert first_negative_degree(pl) == 4
    assert [x.signed_multiplicity for x in candidate_generators(pl)] == [2]
    assert len(first_relation_candidates(pl)) == 1
    records = extract_operator_content(pl)
    assert [x.t_degree for x in records] == sorted(x.t_degree for x in records)
    assert {x.classification for x in records if x.t_degree == 5} == {"higher_positive_correction"}
    assert {x.classification for x in records if x.t_degree == 6} == {"higher_negative_correction"}
    assert all(x.mixed_degree for x in records if x.t_degree == 4)


def test_sage_product_decompositions_and_dimensions():
    theory = load_theory(ROOT / "theories/su3_nf5_k3o2_infinite.yaml")
    a5 = theory.simple_factors[0]
    assert tensor_product(a5, [1,0,0,0,0], [0,0,0,0,1]) == (
        ((0,0,0,0,0), 1), ((1,0,0,0,1), 1))
    assert sum(irrep_dimension(a5, x)*m for x,m in tensor_product(a5,
        [1,0,0,0,0], [0,0,0,0,1])) == 6*6
    from hwg_pipeline.model import SimpleGroupSpec
    a1 = SimpleGroupSpec("a", "A", 1, "SU(2)", ("u",))
    assert symmetric_power(a1, [1], 2) == (((2,), 1),)
    assert tensor_product(a1, [1], [1]) == (((0,), 1), ((2,), 1))


def test_physical_fixture_content_and_channels():
    from hwg_pipeline.__main__ import _load_virtual_pl
    theory = load_theory(ROOT / "theories/su3_nf5_k3o2_infinite.yaml")
    path = ROOT / "generated/su3_nf5_k3o2_infinite/order_10/refined_plethystic_logarithm.json"
    pl = _load_virtual_pl(theory, path, 10)
    generators, relations = candidate_generators(pl), first_relation_candidates(pl)
    assert first_negative_degree(pl) == 4
    assert [(int(x.t_degree), int(dict(x.abelian_charges)["q"]), x.dynkin_labels[0], int(x.representation_dimension)) for x in generators] == [
        (2,0,(0,0,0,0,0),1), (2,0,(1,0,0,0,1),35),
        (3,-1,(0,0,0,1,0),15), (3,1,(0,1,0,0,0),15)]
    assert [(x.dynkin_labels[0], int(x.signed_multiplicity)) for x in relations] == [
        ((0,0,0,0,0),-1), ((1,0,0,0,1),-1)]
    channels = enumerate_quadratic_channels(theory, generators, relations)
    assert [x["product_type"] for x in channels] == ["symmetric_square", "tensor", "symmetric_square"]
    free = {}
    for channel in channels:
        assert channel["t_degree"] == 4 and dict(channel["abelian_charges"])["q"] == 0
        for labels, mult in channel["decomposition"]:
            free[labels] = free.get(labels, 0) + mult
    assert free == {((0,0,0,0,0),):2, ((1,0,0,0,1),):2,
                    ((0,1,0,1,0),):1, ((2,0,0,0,2),):1}
    assert sum(irrep_dimension(theory.simple_factors[0], k[0])*v for k,v in free.items()) == 666


def test_invalid_product_data_rejected():
    theory = load_theory(ROOT / "theories/su3_nf5_k3o2_infinite.yaml")
    with pytest.raises(ValueError):
        tensor_decomposition(theory, (), ())
