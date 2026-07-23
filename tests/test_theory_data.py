from copy import deepcopy

import pytest
import yaml
from sage.all import QQ, ZZ

from hwg_pipeline import (dump_theory, load_theory, render_monomial, render_pe,
                          render_pe_exponent, render_rational_product,
                          theory_from_dict, theory_to_dict)


@pytest.fixture
def raw():
    return yaml.safe_load(open("theories/template.yaml", encoding="utf-8"))


def test_template_exact_values_singlets_and_signs(raw):
    theory = theory_from_dict(raw)
    assert theory.chern_simons_level == QQ(3) / 2
    assert theory.pe.terms[1].monomial.abelian_charges[0][1] == QQ(-1)
    assert [x.coefficient for x in theory.pe.terms] == [ZZ(1), ZZ(-1)]
    singlet = deepcopy(raw)
    singlet["pe"]["terms"][0]["monomial"]["representations"]["flavor"] = [0, 0]
    assert theory_from_dict(singlet).pe.terms[0].monomial.representations[0].dynkin_labels == (0, 0)


@pytest.mark.parametrize("location", [
    lambda x: x.update(chern_simons_level=1.5),
    lambda x: x["pe"]["terms"][0]["monomial"]["abelian_charges"].update(topological=1.5),
    lambda x: x["source_references"][0].update(extra=[{"deep": 1.5}]),
])
def test_rejects_floats_anywhere(raw, location):
    location(raw)
    with pytest.raises(ValueError, match="floating-point YAML value.*root"):
        theory_from_dict(raw)


def test_dynkin_label_length_validation(raw):
    raw["pe"]["terms"][0]["monomial"]["representations"]["flavor"] = [1]
    with pytest.raises(ValueError, match="Dynkin-label length 1 does not match rank 2"):
        theory_from_dict(raw)


def test_undeclared_factors_and_fugacities(raw):
    bad = deepcopy(raw)
    bad["pe"]["terms"][0]["monomial"]["representations"]["ghost"] = [0]
    with pytest.raises(ValueError, match="undeclared simple factor.*ghost"):
        theory_from_dict(bad)
    bad = deepcopy(raw)
    bad["pe"]["terms"][0]["monomial"]["abelian_charges"]["ghost"] = 0
    with pytest.raises(ValueError, match="undeclared abelian factor.*ghost"):
        theory_from_dict(bad)
    bad = deepcopy(raw)
    bad["pe"]["terms"][0]["monomial"]["mu_ghost"] = 1
    with pytest.raises(ValueError, match="undeclared highest-weight fugacities.*mu_ghost"):
        theory_from_dict(bad)


def test_other_validation_messages(raw):
    cases = [
        (("simple_factors", 0, "cartan_type"), "AA", "malformed Cartan"),
        (("pe", "terms", 0, "coefficient"), "1/2", "coefficient must be an integer"),
        (("pe", "terms", 0, "monomial", "t_degree"), -1, "t-degree must be nonnegative"),
        (("rational_product", "factors", 0, "power"), 0, "nonzero integer"),
    ]
    for path, value, message in cases:
        item = deepcopy(raw)
        target = item
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        with pytest.raises(ValueError, match=message):
            theory_from_dict(item)


def test_deterministic_rendering_and_product_sign_convention(raw):
    theory = theory_from_dict(raw)
    positive = theory.pe.terms[0].monomial
    assert render_monomial(positive, theory) == r"mu_1 q^{\frac{1}{2}} t^{2}"
    assert render_pe_exponent(theory.pe, theory) == r"mu_1 q^{\frac{1}{2}} t^{2} - mu_2 q^{-1} t^{4}"
    assert render_pe(theory.pe, theory) == r"\operatorname{PE}\!\left[mu_1 q^{\frac{1}{2}} t^{2} - mu_2 q^{-1} t^{4}\right]"
    product = render_rational_product(theory.rational_product, theory)
    assert product == (r"\left(1 - mu_1 q^{\frac{1}{2}} t^{2}\right)^{-1} "
                       r"\left(1 - mu_2 q^{-1} t^{4}\right)")


def test_stable_load_serialize_reload(tmp_path):
    first = load_theory("theories/template.yaml")
    serialized = dump_theory(first)
    path = tmp_path / "copy.yaml"
    path.write_text(serialized, encoding="utf-8")
    second = load_theory(path)
    assert second == first
    assert theory_to_dict(second) == theory_to_dict(first)
    assert dump_theory(second) == serialized
    assert first.pe.original_pe_latex.startswith("mu_1")
    assert first.rational_product.original_rational_product_latex.startswith("(1-")
