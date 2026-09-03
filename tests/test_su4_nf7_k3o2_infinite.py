"""Exact transcription and independent PL benchmarks for SU(4)+7F."""

from sage.all import QQ, ZZ

from hwg_pipeline import load_theory
from hwg_pipeline.characters import restore_characters
from hwg_pipeline.expansion import expand_pe
from hwg_pipeline.plethystic import (dimension_refine_virtual,
                                     plethystic_logarithm, unrefine_virtual)


PATH = "theories/su4_nf7_k3o2_infinite.yaml"
ORIGINAL_PE = (r"\PE\!\left[(\mu_1\mu_7+1)t^2+(\mu_2\mu_6+q\mu_3+q^{-1}\mu_5)t^4" "\n"
               r"+\mu_3\mu_5t^6-\mu_3\mu_5t^8\right]")


def signature(item):
    monomial = item.monomial
    return (int(item.coefficient), int(monomial.t_degree),
            tuple(int(x) for x in monomial.representations[0].dynkin_labels),
            monomial.abelian_charges[0][1])


def test_exact_metadata_and_source_string():
    theory = load_theory(PATH)
    assert theory.id == "su4_nf7_k3o2_infinite"
    assert (theory.gauge_algebra, theory.gauge_display_name) == ("A3", "SU(4)")
    assert theory.number_of_flavours == ZZ(7)
    assert theory.chern_simons_level == QQ(3) / 2
    assert theory.chern_simons_convention == "absolute value"
    assert theory.coupling == "infinite"
    assert theory.simple_factors[0].cartan_name == "A7"
    assert theory.simple_factors[0].display_name == "SU(8)"
    assert theory.simple_factors[0].highest_weight_fugacities == tuple(f"mu_{i}" for i in range(1, 8))
    assert [(x.id, x.display_name, x.fugacity) for x in theory.abelian_factors] == [("q", "U(1)", "q")]
    assert theory.source_references[0].equation == "11.3"
    assert theory.pe.original_pe_latex == ORIGINAL_PE
    assert theory.rational_product is None


def test_exact_seven_structured_pe_terms():
    theory = load_theory(PATH)
    assert len(theory.pe.terms) == 7
    assert [signature(x) for x in theory.pe.terms] == [
        (1, 2, (1, 0, 0, 0, 0, 0, 1), QQ(0)),
        (1, 2, (0, 0, 0, 0, 0, 0, 0), QQ(0)),
        (1, 4, (0, 1, 0, 0, 0, 1, 0), QQ(0)),
        (1, 4, (0, 0, 1, 0, 0, 0, 0), QQ(1)),
        (1, 4, (0, 0, 0, 0, 1, 0, 0), QQ(-1)),
        (1, 6, (0, 0, 1, 0, 1, 0, 0), QQ(0)),
        (-1, 8, (0, 0, 1, 0, 1, 0, 0), QQ(0)),
    ]


def test_independent_low_degree_plethystic_logarithm_benchmarks():
    theory = load_theory(PATH)
    pl = plethystic_logarithm(restore_characters(theory, expand_pe(theory, 4)), 4)
    by_degree = {}
    for (degree, charges), content in pl:
        for labels, coefficient in content:
            by_degree.setdefault(int(degree), {})[
                (int(dict(charges)["q"]), labels[0])] = int(coefficient)
    assert by_degree[2] == {
        (0, (0, 0, 0, 0, 0, 0, 0)): 1,
        (0, (1, 0, 0, 0, 0, 0, 1)): 1,
    }
    assert 3 not in by_degree
    assert by_degree[4] == {
        (1, (0, 0, 1, 0, 0, 0, 0)): 1,
        (-1, (0, 0, 0, 0, 1, 0, 0)): 1,
        (0, (0, 0, 0, 0, 0, 0, 0)): -1,
        (0, (1, 0, 0, 0, 0, 0, 1)): -1,
    }
    assert dict(unrefine_virtual(pl)) == {ZZ(2): QQ(64), ZZ(4): QQ(48)}


def test_independent_degree_ten_dimension_benchmark_and_symmetry():
    theory = load_theory(PATH)
    pl = plethystic_logarithm(restore_characters(theory, expand_pe(theory, 10)), 10)
    degree_ten = {int(dict(charges)["q"]): coefficient
                  for (degree, charges), coefficient in dimension_refine_virtual(pl)
                  if degree == 10}
    expected = {-2: ZZ(16576), -1: ZZ(-37024), 0: ZZ(116879),
                1: ZZ(-37024), 2: ZZ(16576)}
    assert degree_ten == expected
    assert all(degree_ten[charge] == degree_ten[-charge] for charge in degree_ten)
    assert dict(unrefine_virtual(pl))[ZZ(10)] == ZZ(75983)
    assert not any(degree <= 9 and degree % 2 for (degree, _), _ in pl)
