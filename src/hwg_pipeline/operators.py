"""Conservative extraction and first-relation representation analysis."""

from dataclasses import dataclass
from itertools import product

from sage.all import QQ, ZZ

from .sage_backend import adams_decomposition, irrep_dimension, symmetric_power, tensor_product


@dataclass(frozen=True, order=True)
class OperatorRecord:
    t_degree: object
    abelian_charges: tuple
    dynkin_labels: tuple
    signed_multiplicity: object
    representation_dimension: object
    sign: str
    classification: str
    source_theory_id: str
    mixed_degree: bool = False


def _records(pl_series):
    raw = []
    for (degree, charges), content in pl_series:
        for labels, multiplicity in content:
            if not multiplicity:
                continue
            if multiplicity not in ZZ:
                raise ValueError("operator multiplicities must be integral")
            dimension = ZZ(1)
            for factor, label in zip(pl_series.theory.simple_factors, labels):
                dimension *= irrep_dimension(factor, label)
            raw.append((ZZ(degree), tuple(charges), tuple(tuple(ZZ(x) for x in y) for y in labels),
                        ZZ(multiplicity), dimension))
    return raw


def first_negative_degree(pl_series):
    degrees = [d for d, _, _, m, _ in _records(pl_series) if m < 0]
    return min(degrees) if degrees else None


def extract_operator_content(pl_series):
    raw = _records(pl_series)
    negative = min((d for d, _, _, m, _ in raw if m < 0), default=None)
    signs = {}
    for d, _, _, m, _ in raw:
        signs.setdefault(d, set()).add(1 if m > 0 else -1)
    result = []
    for degree, charges, labels, multiplicity, dimension in raw:
        if multiplicity > 0 and (negative is None or degree < negative):
            classification = "low_degree_generator_candidate"
        elif multiplicity < 0 and degree == negative:
            classification = "first_relation_candidate"
        elif multiplicity > 0:
            classification = "higher_positive_correction"
        else:
            classification = "higher_negative_correction"
        result.append(OperatorRecord(degree, charges, labels, multiplicity, dimension,
                                     "positive" if multiplicity > 0 else "negative",
                                     classification, pl_series.theory.id, len(signs[degree]) > 1))
    return tuple(sorted(result))


def candidate_generators(pl_series):
    return tuple(x for x in extract_operator_content(pl_series)
                 if x.classification == "low_degree_generator_candidate")


def first_relation_candidates(pl_series):
    return tuple(x for x in extract_operator_content(pl_series)
                 if x.classification == "first_relation_candidate")


def _combine_factor_decompositions(parts):
    out = {}
    for choice in product(*parts):
        labels = tuple(x[0] for x in choice)
        mult = ZZ.prod(x[1] for x in choice)
        out[labels] = out.get(labels, ZZ(0)) + mult
    return {k: v for k, v in out.items() if v}


def tensor_decomposition(theory, left, right):
    if len(left) != len(theory.simple_factors) or len(right) != len(theory.simple_factors):
        raise ValueError("product data must contain every simple factor")
    return _combine_factor_decompositions([
        tensor_product(f, a, b) for f, a, b in zip(theory.simple_factors, left, right)])


def symmetric_square_decomposition(theory, labels):
    """Sym^2 of an external tensor product, via (chi^2 + psi_2 chi)/2."""
    if len(labels) != len(theory.simple_factors):
        raise ValueError("product data must contain every simple factor")
    square = tensor_decomposition(theory, labels, labels)
    adams = _combine_factor_decompositions([
        adams_decomposition(f, label, 2) for f, label in zip(theory.simple_factors, labels)])
    keys = set(square) | set(adams)
    result = {key: (square.get(key, 0) + adams.get(key, 0)) / 2 for key in keys}
    if any(value not in ZZ or value < 0 for value in result.values()):
        raise ArithmeticError("Sage returned an invalid symmetric-square decomposition")
    return {key: ZZ(value) for key, value in result.items() if value}


def enumerate_quadratic_channels(theory, generators, relations):
    generators, relations = tuple(generators), tuple(relations)
    if any(g.signed_multiplicity <= 0 for g in generators):
        raise ValueError("generator multiplicities must be positive")
    relation_keys = {(r.t_degree, r.abelian_charges, r.dynkin_labels) for r in relations}
    channels = []
    for i, left in enumerate(generators):
        for j in range(i, len(generators)):
            right = generators[j]
            degree = left.t_degree + right.t_degree
            if tuple(k for k, _ in left.abelian_charges) != tuple(k for k, _ in right.abelian_charges):
                raise ValueError("generator charge factors do not agree")
            charges = tuple((k, a + b) for (k, a), (_, b) in zip(left.abelian_charges, right.abelian_charges))
            if not any(d == degree and q == charges for d, q, _ in relation_keys):
                continue
            symmetric = i == j
            if symmetric:
                sym = symmetric_square_decomposition(theory, left.dynkin_labels)
                square = tensor_decomposition(theory, left.dynkin_labels, left.dynkin_labels)
                copies = left.signed_multiplicity
                decomposition = {key: copies * sym.get(key, 0) +
                                  copies * (copies - 1) // 2 * square.get(key, 0)
                                  for key in set(sym) | set(square)}
                decomposition = {key: value for key, value in decomposition.items() if value}
            else:
                decomposition = tensor_decomposition(theory, left.dynkin_labels, right.dynkin_labels)
                decomposition = {key: value * left.signed_multiplicity * right.signed_multiplicity
                                 for key, value in decomposition.items()}
            hits = {labels: decomposition.get(labels, ZZ(0)) for d, q, labels in relation_keys
                    if d == degree and q == charges and decomposition.get(labels, 0)}
            channels.append({"left": left, "right": right,
                             "product_type": "symmetric_square" if symmetric else "tensor",
                             "t_degree": degree, "abelian_charges": charges,
                             "decomposition": tuple(sorted(decomposition.items())),
                             "relation_multiplicities": tuple(sorted(hits.items()))})
    return tuple(channels)
