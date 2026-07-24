"""Independent fixture and product-group checks for SU(3)+5F, |k|=1/2."""

from pathlib import Path

from sage.all import QQ, ZZ

from hwg_pipeline import (RepresentationContent, SimpleGroupSpec, expand_hwg,
                          irrep_dimension, load_theory, restore_characters)
from hwg_pipeline.plethystic import (VirtualCharacterSeries,
    VirtualRepresentationContent, adams_series, plethystic_exponential,
    plethystic_logarithm)


def _group(identifier, rank):
    return SimpleGroupSpec(identifier, "A", rank, f"A{rank}",
                           tuple(f"m{i}" for i in range(rank)))


def test_source_formula_and_fixture_are_termwise_identical():
    source = Path("references/overleaf/su3_5f_6f_hwg_results.tex").read_text()
    assert r"\HWG_{5F,\frac12}" in source
    assert r"(\mu_1\mu_4+\nu^2+1)t^2" in source
    theory = load_theory("theories/su3_nf5_k1o2_infinite.yaml")
    assert theory.source_references[0].equation == "12.3"
    assert len(theory.pe.terms) == 7 and len(theory.rational_product.factors) == 7
    assert tuple(x.cartan_name for x in theory.simple_factors) == ("A4", "A1")
    pe_terms = {term.monomial: term.coefficient for term in theory.pe.terms}
    product_terms = {factor.monomial: -factor.power
                     for factor in theory.rational_product.factors}
    assert pe_terms == product_terms
    assert expand_hwg(theory, 10) == expand_hwg(theory, 11).truncate(10)


def test_a4_a1_coroot_dimensions_and_product_dimensions():
    a4, a1 = _group("a4", 4), _group("a1", 1)
    expected_a4 = { (0,0,0,0): 1, (1,0,0,0): 5, (0,0,0,1): 5,
        (1,0,0,1): 24, (0,1,0,0): 10, (0,0,1,0): 10,
        (0,1,1,0): 75 }
    assert {labels: irrep_dimension(a4, labels) for labels in expected_a4} == expected_a4
    assert [irrep_dimension(a1, (n,)) for n in range(3)] == [1, 2, 3]
    products = [(((1,0,0,1),(0,)),24), (((0,0,0,0),(2,)),3),
                (((0,1,0,0),(1,)),20), (((0,0,1,0),(1,)),20),
                (((0,1,1,0),(0,)),75), (((0,1,1,0),(2,)),225)]
    for labels, dimension in products:
        assert RepresentationContent.single_irrep((a4, a1), labels).total_dimension() == dimension


def test_product_operations_are_factorwise_and_keep_trivial_factors():
    a2, a1 = _group("a2", 2), _group("a1", 1)
    specs = (a2, a1)
    left = RepresentationContent.single_irrep(specs, ((1,0),(1,)))
    right = RepresentationContent.single_irrep(specs, ((0,1),(1,)))
    product = left * right
    assert product.total_dimension() == left.total_dimension() * right.total_dimension()
    assert all(len(key) == 2 and len(key[0]) == 2 and len(key[1]) == 1 for key, _ in product)
    assert RepresentationContent.single_irrep(specs, ((1,0),(0,))) != left


def test_product_adams_and_small_refined_pe_pl_round_trip():
    a1x, a1y = _group("x", 1), _group("y", 1)
    theory = type("Theory", (), {"simple_factors": (a1x, a1y),
                                  "abelian_factors": ()})()
    content = VirtualRepresentationContent.single_irrep((a1x, a1y), ((1,), (1,)))
    source = VirtualCharacterSeries(theory, [((2, ()), content)], 6)
    doubled = adams_series(source, 2, 6)
    assert all(len(labels) == 2 for _, terms in doubled for labels, _ in terms)
    hilbert = plethystic_exponential(source, 6)
    assert plethystic_logarithm(hilbert, 6) == source


def test_physical_restoration_keeps_factor_order_and_q_sectors():
    theory = load_theory("theories/su3_nf5_k1o2_infinite.yaml")
    series = restore_characters(theory, expand_hwg(theory, 3))
    degree_three = [(charges, labels) for (degree, charges), content in series
                    if degree == 3 for labels, coefficient in content if coefficient]
    assert {labels for _, labels in degree_three} == {
        ((0,1,0,0),(1,)), ((0,0,1,0),(1,))}
    assert {dict(charges)["q"] for charges, _ in degree_three} == {QQ(-1), QQ(1)}
