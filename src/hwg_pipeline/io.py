"""Safe YAML loading and stable serialization for theory specifications."""

from pathlib import Path
import yaml
from sage.all import QQ, ZZ

from .model import *


def _reject_floats(value, path="root"):
    if isinstance(value, float):
        raise ValueError(f"floating-point YAML value is forbidden at {path}")
    if isinstance(value, dict):
        for key, child in value.items():
            _reject_floats(key, f"{path}.<key>")
            _reject_floats(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _reject_floats(child, f"{path}[{index}]")


def _exact(value, field):
    if isinstance(value, bool) or not isinstance(value, (int, str)):
        raise ValueError(f"{field} must be an integer or rational string")
    try:
        return QQ(value)
    except (TypeError, ValueError):
        raise ValueError(f"invalid exact rational for {field}: {value!r}") from None


def _integer(value, field):
    rational = _exact(value, field)
    if rational.denominator() != 1:
        raise ValueError(f"{field} must be an integer")
    return ZZ(rational)


def _monomial(data):
    allowed = {"t_degree", "representations", "abelian_charges"}
    extra = set(data) - allowed
    if extra:
        raise ValueError(f"undeclared highest-weight fugacities/fields: {', '.join(sorted(extra))}")
    reps = tuple(RepresentationSpec(key, tuple(labels)) for key, labels in data.get("representations", {}).items())
    charges = tuple((key, _exact(value, f"charge {key}")) for key, value in data.get("abelian_charges", {}).items())
    return HighestWeightMonomial(_integer(data["t_degree"], "t_degree"), reps, charges)


def theory_from_dict(data):
    _reject_floats(data)
    simple = tuple(SimpleGroupSpec(x["id"], x["cartan_type"], _integer(x["rank"], "rank"), x["display_name"], tuple(x["highest_weight_fugacities"])) for x in data["simple_factors"])
    abelian = tuple(AbelianFactorSpec(x["id"], x["display_name"], x["fugacity"]) for x in data.get("abelian_factors", []))
    pe_data = data["pe"]
    pe = PlethysticExponentialSpec(tuple(HWGTerm(_integer(x["coefficient"], "HWG term coefficient"), _monomial(x["monomial"])) for x in pe_data["terms"]), pe_data["original_pe_latex"])
    product_data = data.get("rational_product")
    product = None if product_data is None else RationalProductSpec(tuple(RationalProductFactor(_monomial(x["monomial"]), _integer(x["power"], "factor power")) for x in product_data["factors"]), product_data.get("original_rational_product_latex"))
    refs = tuple(SourceReference(x["path"], x["description"], x.get("equation")) for x in data.get("source_references", []))
    return TheorySpec(data["id"], data["title"], bool(data.get("nonphysical", False)), refs, simple, abelian, _exact(data["chern_simons_level"], "chern_simons_level"), pe, product)


def load_theory(path):
    with Path(path).open(encoding="utf-8") as stream:
        data = yaml.safe_load(stream)
    if not isinstance(data, dict):
        raise ValueError("theory YAML root must be a mapping")
    return theory_from_dict(data)


def _rational_text(value):
    return str(QQ(value))


def _monomial_dict(value):
    return {"t_degree": int(value.t_degree), "representations": {x.simple_factor_id: [int(y) for y in x.dynkin_labels] for x in value.representations}, "abelian_charges": {key: _rational_text(x) for key, x in value.abelian_charges}}


def theory_to_dict(value):
    result = {"id": value.id, "title": value.title, "nonphysical": value.nonphysical,
      "source_references": [{k: v for k, v in {"path": x.path, "description": x.description, "equation": x.equation}.items() if v is not None} for x in value.source_references],
      "chern_simons_level": _rational_text(value.chern_simons_level),
      "simple_factors": [{"id": x.id, "cartan_type": x.cartan_type, "rank": int(x.rank), "display_name": x.display_name, "highest_weight_fugacities": list(x.highest_weight_fugacities)} for x in value.simple_factors],
      "abelian_factors": [{"id": x.id, "display_name": x.display_name, "fugacity": x.fugacity} for x in value.abelian_factors],
      "pe": {"original_pe_latex": value.pe.original_pe_latex, "terms": [{"coefficient": int(x.coefficient), "monomial": _monomial_dict(x.monomial)} for x in value.pe.terms]}}
    if value.rational_product:
        result["rational_product"] = {"original_rational_product_latex": value.rational_product.original_rational_product_latex, "factors": [{"power": int(x.power), "monomial": _monomial_dict(x.monomial)} for x in value.rational_product.factors]}
    return result


def dump_theory(value):
    return yaml.safe_dump(theory_to_dict(value), sort_keys=False, allow_unicode=True)
