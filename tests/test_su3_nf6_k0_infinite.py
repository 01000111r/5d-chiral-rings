"""Source, cutoff, and ordered triple-product checks for SU(3)+6F, k=0."""
from pathlib import Path

from sage.all import QQ, WeylCharacterRing

from hwg_pipeline import (RepresentationContent, SimpleGroupSpec, expand_hwg,
                          irrep_dimension, load_theory)
from hwg_pipeline.plethystic import (VirtualCharacterSeries,
    VirtualRepresentationContent, adams_series, plethystic_exponential,
    plethystic_logarithm)


PATH = "theories/su3_nf6_k0_infinite.yaml"


def _group(identifier, rank):
    return SimpleGroupSpec(identifier, "A", rank, f"A{rank}",
                           tuple(f"m{i}" for i in range(rank)))


def test_source_fixture_termwise_and_cutoff_boundary():
    source = Path("references/overleaf/su3_5f_6f_hwg_results.tex").read_text()
    assert r"\HWG_{6F,0}" in source
    assert r"-\nu_1^2\nu_2^2\mu_3^2t^{10}" in source
    theory = load_theory(PATH)
    assert theory.chern_simons_level == QQ(0)
    assert theory.chern_simons_convention == "signed zero"
    assert theory.source_references[0].equation == "9.5"
    assert tuple(x.id for x in theory.simple_factors) == ("a5", "a1_1", "a1_2")
    assert tuple(x.cartan_name for x in theory.simple_factors) == ("A5", "A1", "A1")
    assert not theory.abelian_factors and len(theory.pe.terms) == 9
    assert len(theory.rational_product.factors) == 9
    assert {x.monomial: x.coefficient for x in theory.pe.terms} == {
        x.monomial: -x.power for x in theory.rational_product.factors}
    boundary = [x for x in theory.pe.terms if x.monomial.t_degree == 10]
    assert len(boundary) == 1 and boundary[0].coefficient == -1
    assert expand_hwg(theory, 10) == expand_hwg(theory, 11).truncate(10)


def test_a5_and_separate_a1_coroot_conventions():
    a5 = WeylCharacterRing("A5", style="coroots")
    labels = [(0,0,0,0,0), (1,0,0,0,0), (0,0,0,0,1),
              (1,0,0,0,1), (0,0,1,0,0), (0,1,0,1,0), (0,0,2,0,0)]
    assert [a5(x).degree() for x in labels] == [1, 6, 6, 35, 20, 189, 175]
    assert a5(labels[4]).dual() == a5(labels[4])
    assert a5(labels[6]).dual() == a5(labels[6])
    a1_1 = WeylCharacterRing("A1", style="coroots")
    a1_2 = WeylCharacterRing("A1", style="coroots")
    for ring in (a1_1, a1_2):
        assert [ring((n,)).degree() for n in range(3)] == [1, 2, 3]
        assert all(ring((n,)).dual() == ring((n,)) for n in range(3))


def test_triple_dimensions_and_identical_factor_positions_are_distinct():
    specs = (_group("a5", 5), _group("a1_1", 1), _group("a1_2", 1))
    cases = [(((1,0,0,0,1),(0,),(0,)),35),
             (((0,0,0,0,0),(2,),(0,)),3),
             (((0,0,0,0,0),(0,),(2,)),3),
             (((0,0,1,0,0),(1,),(1,)),80),
             (((0,1,0,1,0),(0,),(0,)),189),
             (((0,0,2,0,0),(0,),(0,)),175),
             (((0,0,2,0,0),(2,),(2,)),1575)]
    contents = []
    for labels, dimension in cases:
        content = RepresentationContent.single_irrep(specs, labels)
        assert content.total_dimension() == dimension
        contents.append(content)
    assert contents[1] != contents[2]


def test_triple_product_operations_serialization_shape_and_round_trip():
    specs = (_group("a2", 2), _group("left", 1), _group("right", 1))
    left = RepresentationContent.single_irrep(specs, ((1,0),(1,),(0,)))
    right = RepresentationContent.single_irrep(specs, ((0,1),(0,),(1,)))
    product = left * right
    assert product.total_dimension() == left.total_dimension() * right.total_dimension()
    assert all(tuple(map(len, labels)) == (2,1,1) for labels, _ in product)
    theory = type("Theory", (), {"simple_factors": specs, "abelian_factors": ()})()
    source = VirtualCharacterSeries(theory, [((2, ()),
        VirtualRepresentationContent.single_irrep(specs, ((1,0),(1,),(0,))))], 6)
    doubled = adams_series(source, 2, 6)
    assert all(tuple(map(len, labels)) == (2,1,1)
               for _, terms in doubled for labels, _ in terms)
    assert plethystic_logarithm(plethystic_exponential(source, 6), 6) == source


def test_no_abelian_compact_report_has_only_two_level_hierarchy():
    report = Path("generated/su3_nf6_k0_infinite/order_10/compact_report/compact_results.tex")
    if report.exists():
        text = report.read_text()
        assert "q-refined" not in text and "U}(1)_q" not in text
        assert r"[a_1,\ldots,a_{5};b;c]:=" in text
        assert r"A1^{(1)}" in text and r"A1^{(2)}" in text
