"""Source-transcription and independent D6 convention checks."""
from pathlib import Path

from sage.all import QQ, WeylCharacterRing

from hwg_pipeline import load_theory


PATH = "theories/su3_nf6_k2_infinite.yaml"
PE = r"\PE\!\left[(\mu_2+1)t^2+(q+q^{-1})\mu_6t^3+\mu_4t^4\right]"
PRODUCT = r"\frac{1}{(1-t^2)(1-\mu_2t^2)(1-q\mu_6t^3)(1-q^{-1}\mu_6t^3)(1-\mu_4t^4)}"


def signature(item):
    monomial = item.monomial
    return (int(getattr(item, "coefficient", getattr(item, "power", 0))),
            int(monomial.t_degree),
            tuple(map(int, monomial.representations[0].dynkin_labels)),
            monomial.abelian_charges[0][1])


def test_fixture_exact_source_transcription():
    theory = load_theory(PATH)
    assert (theory.id, theory.gauge_algebra, theory.gauge_display_name,
            int(theory.number_of_flavours)) == ("su3_nf6_k2_infinite", "A2", "SU(3)", 6)
    assert theory.chern_simons_level == QQ(2)
    assert theory.chern_simons_convention == "absolute value"
    assert theory.simple_factors[0].cartan_name == "D6"
    assert theory.simple_factors[0].display_name == "SO(12)"
    assert theory.source_references[0].equation == "7.3"
    assert theory.pe.original_pe_latex == PE
    assert theory.rational_product.original_rational_product_latex == PRODUCT
    expected = [
        (1, 2, (0, 0, 0, 0, 0, 0), QQ(0)),
        (1, 2, (0, 1, 0, 0, 0, 0), QQ(0)),
        (1, 3, (0, 0, 0, 0, 0, 1), QQ(1)),
        (1, 3, (0, 0, 0, 0, 0, 1), QQ(-1)),
        (1, 4, (0, 0, 0, 1, 0, 0), QQ(0)),
    ]
    assert [signature(item) for item in theory.pe.terms] == expected
    assert len(theory.rational_product.factors) == 5
    assert {x.monomial: -x.power for x in theory.rational_product.factors} == {
        x.monomial: x.coefficient for x in theory.pe.terms}


def test_d6_coroot_conventions_dimensions_and_conjugation():
    d6 = WeylCharacterRing("D6", style="coroots")
    labels = [(0, 0, 0, 0, 0, 0), (1, 0, 0, 0, 0, 0),
              (0, 1, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1),
              (0, 0, 0, 0, 1, 0), (0, 0, 0, 1, 0, 0)]
    assert [d6(label).degree() for label in labels] == [1, 12, 66, 32, 32, 495]
    # D6 half-spinors are distinct and self-dual; terminal-label swapping is
    # therefore not used to alter the source's explicit mu_6 convention.
    assert d6(labels[3]).dual() == d6(labels[3])
    assert d6(labels[4]).dual() == d6(labels[4])
    assert d6(labels[3]) != d6(labels[4])


def test_compact_report_uses_d6_symmetry():
    path = Path("generated/su3_nf6_k2_infinite/order_10/compact_report/compact_results.tex")
    if path.exists():
        report = path.read_text()
        assert r"\mathrm{SO}(12)\times\mathrm{U}(1)_q" in report
        assert r":=[a_1,\ldots,a_{6}]_{D6}" in report
