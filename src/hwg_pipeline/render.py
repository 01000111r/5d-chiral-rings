"""Deterministic normalized LaTeX rendering (not a LaTeX parser)."""

from sage.all import QQ


def _power(name, exponent):
    exponent = QQ(exponent)
    if exponent == 0:
        return ""
    if exponent == 1:
        return name
    text = str(exponent) if exponent.denominator() == 1 else rf"\frac{{{exponent.numerator()}}}{{{exponent.denominator()}}}"
    return rf"{name}^{{{text}}}"


def render_monomial(monomial, theory):
    reps = {x.simple_factor_id: x for x in monomial.representations}
    charges = dict(monomial.abelian_charges)
    pieces = []
    for factor in theory.simple_factors:
        pieces += [_power(name, label) for name, label in zip(factor.highest_weight_fugacities, reps[factor.id].dynkin_labels)]
    pieces += [_power(factor.fugacity, charges[factor.id]) for factor in theory.abelian_factors]
    pieces.append(_power("t", monomial.t_degree))
    return " ".join(x for x in pieces if x) or "1"


def render_pe_exponent(pe, theory):
    ordered = sorted(pe.terms, key=lambda x: (int(x.monomial.t_degree), render_monomial(x.monomial, theory), int(x.coefficient)))
    parts = []
    for term in ordered:
        sign = "-" if term.coefficient < 0 else "+"
        magnitude = abs(term.coefficient)
        body = render_monomial(term.monomial, theory)
        if magnitude != 1:
            body = f"{magnitude} {body}"
        parts.append((sign, body))
    if not parts:
        return "0"
    first_sign, first = parts[0]
    result = ("-" if first_sign == "-" else "") + first
    return result + "".join(f" {sign} {body}" for sign, body in parts[1:])


def render_pe(pe, theory):
    return rf"\operatorname{{PE}}\!\left[{render_pe_exponent(pe, theory)}\right]"


def render_rational_product(product, theory):
    ordered = sorted(product.factors, key=lambda x: (int(x.power), render_monomial(x.monomial, theory)))
    pieces = []
    for factor in ordered:
        piece = rf"\left(1 - {render_monomial(factor.monomial, theory)}\right)"
        if factor.power != 1:
            piece += rf"^{{{factor.power}}}"
        pieces.append(piece)
    return " ".join(pieces) or "1"
