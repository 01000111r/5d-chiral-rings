"""Input transcription and independent D10 convention checks."""

from sage.all import QQ, WeylCharacterRing

from hwg_pipeline import expand_pe, load_theory


def test_blind_fixture_matches_supplied_hwg():
    theory = load_theory("theories/su3_nf9_k1o2_infinite.yaml")
    expected = [
        (2, (0, 1, 0, 0, 0, 0, 0, 0, 0, 0)),
        (3, (0, 0, 0, 0, 0, 0, 0, 0, 0, 1)),
        (4, (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
        (4, (0, 0, 0, 1, 0, 0, 0, 0, 0, 0)),
        (5, (0, 0, 0, 0, 0, 0, 0, 0, 0, 1)),
        (6, (0, 0, 0, 0, 0, 1, 0, 0, 0, 0)),
        (8, (0, 0, 0, 0, 0, 0, 0, 1, 0, 0)),
    ]
    assert theory.id == "su3_nf9_k1o2_infinite"
    assert theory.chern_simons_level == QQ(1) / 2
    assert theory.chern_simons_convention == "sign pair ±1/2"
    assert theory.simple_factors[0].cartan_name == "D10"
    assert theory.abelian_factors == ()
    assert len(theory.pe.terms) == 7
    assert [(int(x.monomial.t_degree), tuple(x.monomial.representations[0].dynkin_labels))
            for x in theory.pe.terms] == expected
    assert all(x.coefficient == 1 for x in theory.pe.terms)
    assert expected[-1][0] == 8
    assert all(monomial.t_degree <= 6 for monomial, _ in expand_pe(theory, 6))


def test_sage_d10_coroot_dynkin_label_dimensions():
    ring = WeylCharacterRing("D10", style="coroots")
    fundamental = ring.fundamental_weights()

    def dimension(labels):
        weight = sum((coefficient * fundamental[index + 1]
                      for index, coefficient in enumerate(labels)), ring.space().zero())
        return ring(weight).degree()

    assert dimension((0, 0, 0, 0, 0, 0, 0, 0, 0, 0)) == 1
    assert dimension((1, 0, 0, 0, 0, 0, 0, 0, 0, 0)) == 20
    assert dimension((0, 1, 0, 0, 0, 0, 0, 0, 0, 0)) == 190
    assert dimension((0, 0, 0, 0, 0, 0, 0, 0, 0, 1)) == 512


def test_pe_only_expansion_has_no_manufactured_product():
    theory = load_theory("theories/su3_nf9_k1o2_infinite.yaml")
    assert theory.rational_product is None
    assert expand_pe(theory, 4)
