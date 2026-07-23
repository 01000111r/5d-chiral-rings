"""Transcription checks for the single physical SU(3)+5F fixture."""

from sage.all import QQ, ZZ

from hwg_pipeline import load_theory, render_pe, render_rational_product


PATH = "theories/su3_nf5_k3o2_infinite.yaml"
ORIGINAL_PE = (r"\PE\!\left[(\mu_1\mu_5+1)t^2+(q\mu_2+q^{-1}\mu_4)t^3" "\n"
               r"+\mu_2\mu_4t^4-\mu_2\mu_4t^6\right]")
ORIGINAL_PRODUCT = (r"\frac{1-\mu_2\mu_4t^6}" "\n"
                    r"{(1-t^2)(1-\mu_1\mu_5t^2)(1-q\mu_2t^3)(1-q^{-1}\mu_4t^3)(1-\mu_2\mu_4t^4)}")


def signature(item):
    monomial = item.monomial
    return (int(getattr(item, "coefficient", getattr(item, "power", 0))),
            int(monomial.t_degree),
            tuple(int(x) for x in monomial.representations[0].dynkin_labels),
            monomial.abelian_charges[0][1])


def test_fixture_loads_with_exact_required_metadata():
    theory = load_theory(PATH)
    assert theory.id == "su3_nf5_k3o2_infinite"
    assert theory.chern_simons_level == QQ(3) / 2
    assert theory.gauge_algebra == "A2"
    assert theory.gauge_display_name == "SU(3)"
    assert theory.number_of_flavours == ZZ(5)
    assert theory.chern_simons_convention == "absolute value"
    assert theory.coupling == "infinite"
    assert theory.grading_variable == "t"
    assert theory.simple_factors[0].cartan_name == "A5"
    assert theory.simple_factors[0].display_name == "SU(6)"
    assert [(x.id, x.fugacity) for x in theory.abelian_factors] == [("q", "q")]
    assert theory.source_references[0].equation == "11.3"


def test_exact_six_pe_terms():
    theory = load_theory(PATH)
    assert len(theory.pe.terms) == 6
    assert [signature(x) for x in theory.pe.terms] == [
        (1, 2, (1, 0, 0, 0, 1), QQ(0)),
        (1, 2, (0, 0, 0, 0, 0), QQ(0)),
        (1, 3, (0, 1, 0, 0, 0), QQ(1)),
        (1, 3, (0, 0, 0, 1, 0), QQ(-1)),
        (1, 4, (0, 1, 0, 1, 0), QQ(0)),
        (-1, 6, (0, 1, 0, 1, 0), QQ(0)),
    ]


def test_exact_six_product_factors_and_pe_equivalence():
    theory = load_theory(PATH)
    factors = theory.rational_product.factors
    assert len(factors) == 6
    assert [signature(x) for x in factors] == [
        (1, 6, (0, 1, 0, 1, 0), QQ(0)),
        (-1, 2, (0, 0, 0, 0, 0), QQ(0)),
        (-1, 2, (1, 0, 0, 0, 1), QQ(0)),
        (-1, 3, (0, 1, 0, 0, 0), QQ(1)),
        (-1, 3, (0, 0, 0, 1, 0), QQ(-1)),
        (-1, 4, (0, 1, 0, 1, 0), QQ(0)),
    ]
    # PE[sum c*m] = product (1-m)^(-c), independently of factor order.
    pe_map = {x.monomial: x.coefficient for x in theory.pe.terms}
    product_map = {x.monomial: -x.power for x in factors}
    assert product_map == pe_map
    assert render_pe(theory.pe, theory)
    assert render_rational_product(theory.rational_product, theory)


def test_verbatim_source_strings_are_unchanged():
    theory = load_theory(PATH)
    assert theory.pe.original_pe_latex == ORIGINAL_PE
    assert theory.rational_product.original_rational_product_latex == ORIGINAL_PRODUCT
