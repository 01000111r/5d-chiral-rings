"""Source-transcription and independent D5 convention checks."""
from sage.all import QQ, WeylCharacterRing
from hwg_pipeline import load_theory
from pathlib import Path

PATH = "theories/su3_nf5_k5o2_infinite.yaml"
PE = (r"\PE\!\left[(\mu_2+1)t^2+(q\mu_4+q^{-1}\mu_5)t^3" "\n" r"+\mu_4\mu_5t^4-\mu_4\mu_5t^6\right]")
PRODUCT = (r"\frac{1-\mu_4\mu_5t^6}" "\n" r"{(1-t^2)(1-\mu_2t^2)(1-q\mu_4t^3)(1-q^{-1}\mu_5t^3)(1-\mu_4\mu_5t^4)}")

def signature(item):
    m=item.monomial
    return (int(getattr(item,"coefficient",getattr(item,"power",0))),int(m.t_degree),tuple(map(int,m.representations[0].dynkin_labels)),m.abelian_charges[0][1])

def test_fixture_exact_source_transcription():
    t=load_theory(PATH)
    assert (t.id,t.gauge_algebra,t.gauge_display_name,int(t.number_of_flavours)) == ("su3_nf5_k5o2_infinite","A2","SU(3)",5)
    assert t.chern_simons_level == QQ(5)/2 and t.chern_simons_convention == "absolute value"
    assert t.simple_factors[0].cartan_name == "D5" and t.simple_factors[0].display_name == "SO(10)"
    assert t.source_references[0].equation == "10.3"
    assert t.pe.original_pe_latex == PE and t.rational_product.original_rational_product_latex == PRODUCT
    expected=[(1,2,(0,0,0,0,0),QQ(0)),(1,2,(0,1,0,0,0),QQ(0)),(1,3,(0,0,0,1,0),QQ(1)),(1,3,(0,0,0,0,1),QQ(-1)),(1,4,(0,0,0,1,1),QQ(0)),(-1,6,(0,0,0,1,1),QQ(0))]
    assert [signature(x) for x in t.pe.terms] == expected
    assert {x.monomial:-x.power for x in t.rational_product.factors} == {x.monomial:x.coefficient for x in t.pe.terms}

def test_d5_coroot_conventions_dimensions_and_conjugation():
    D5=WeylCharacterRing("D5",style="coroots")
    labels=[(0,0,0,0,0),(1,0,0,0,0),(0,1,0,0,0),(0,0,0,1,0),(0,0,0,0,1),(0,0,0,1,1)]
    assert [D5(x).degree() for x in labels] == [1,10,45,16,16,210]
    assert D5((0,0,0,1,0)).dual() == D5((0,0,0,0,1))
    assert D5((0,0,0,1,1)).dual() == D5((0,0,0,1,1))

def test_compact_report_uses_d5_symmetry_not_reference_theory_symmetry():
    report = Path("generated/su3_nf5_k5o2_infinite/order_10/compact_report/compact_results.tex").read_text()
    assert r"\mathrm{SO}(10)\times\mathrm{U}(1)_q" in report
    assert r"\mathrm{SU}(6)" not in report
