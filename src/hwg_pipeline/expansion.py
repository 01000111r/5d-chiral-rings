"""Exact sparse expansion of structured highest-weight generating functions."""

from dataclasses import dataclass
from sage.all import ZZ, binomial

from .model import HighestWeightMonomial, RepresentationSpec


def unit_monomial(theory):
    return HighestWeightMonomial(ZZ(0),
        tuple(RepresentationSpec(x.id, (ZZ(0),) * x.rank) for x in theory.simple_factors),
        tuple((x.id, ZZ(0)) for x in theory.abelian_factors))


@dataclass(frozen=True, init=False)
class SparseSeries:
    """An immutable, canonically ordered sparse ZZ-linear combination."""
    terms: tuple
    max_degree: object

    def __init__(self, terms=(), max_degree=0):
        bound = ZZ(max_degree)
        if bound < 0:
            raise ValueError("maximum degree must be nonnegative")
        combined = {}
        source = terms.items() if hasattr(terms, "items") else terms
        for monomial, coefficient in source:
            coefficient = ZZ(coefficient)
            if monomial.t_degree <= bound:
                combined[monomial] = combined.get(monomial, ZZ(0)) + coefficient
        ordered = tuple(sorted(((m, c) for m, c in combined.items() if c), key=lambda item: _key(item[0])))
        object.__setattr__(self, "terms", ordered)
        object.__setattr__(self, "max_degree", bound)

    @classmethod
    def zero(cls, max_degree):
        return cls((), max_degree)

    @classmethod
    def unit(cls, monomial, max_degree):
        return cls(((monomial, ZZ(1)),), max_degree)

    def __iter__(self):
        return iter(self.terms)

    def __len__(self):
        return len(self.terms)

    def coefficient(self, monomial):
        return dict(self.terms).get(monomial, ZZ(0))

    def truncate(self, max_degree):
        return SparseSeries(self.terms, min(ZZ(max_degree), self.max_degree))

    def _combine(self, other, sign):
        bound = min(self.max_degree, other.max_degree)
        return SparseSeries(self.terms + tuple((m, sign * c) for m, c in other), bound)

    def __add__(self, other):
        return self._combine(other, ZZ(1))

    def __sub__(self, other):
        return self._combine(other, ZZ(-1))

    def __neg__(self):
        return self.scalar_mul(-1)

    def scalar_mul(self, scalar):
        scalar = ZZ(scalar)
        return SparseSeries(((m, scalar * c) for m, c in self), self.max_degree)

    def __rmul__(self, scalar):
        return self.scalar_mul(scalar)

    def __mul__(self, other):
        bound = min(self.max_degree, other.max_degree)
        # Products beyond the bound are never inserted (immediate truncation).
        return SparseSeries(((a * b, ca * cb) for a, ca in self for b, cb in other
                             if a.t_degree + b.t_degree <= bound), bound)


def _key(monomial):
    return (monomial.t_degree,
            tuple((x.simple_factor_id, tuple(x.dynkin_labels)) for x in monomial.representations),
            tuple(monomial.abelian_charges))


def _factor(monomial, exponent, max_degree):
    """Expand (1-M)^exponent exactly through max_degree."""
    exponent, max_degree = ZZ(exponent), ZZ(max_degree)
    if monomial.t_degree == 0:
        raise ValueError("degree-zero nonconstant PE monomial has no t-adic expansion")
    limit = max_degree // monomial.t_degree
    if exponent >= 0:
        limit = min(limit, exponent)
        terms = ((monomial ** n, (-1) ** n * binomial(exponent, n)) for n in range(int(limit) + 1))
    else:
        c = -exponent
        terms = ((monomial ** n, binomial(n + c - 1, c - 1)) for n in range(int(limit) + 1))
    return SparseSeries(terms, max_degree)


def expand_pe(theory, max_degree):
    result = SparseSeries.unit(unit_monomial(theory), max_degree)
    for term in theory.pe.terms:
        result = result * _factor(term.monomial, -term.coefficient, max_degree)
    return result


def expand_rational_product(theory, max_degree):
    if theory.rational_product is None:
        raise ValueError("theory has no structured rational-product form")
    result = SparseSeries.unit(unit_monomial(theory), max_degree)
    for factor in theory.rational_product.factors:
        result = result * _factor(factor.monomial, factor.power, max_degree)
    return result


def expand_hwg(theory, max_degree):
    return expand_pe(theory, max_degree)
