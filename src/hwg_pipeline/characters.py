"""Restore irreducible characters from highest-weight monomials."""

from dataclasses import dataclass

from sage.all import QQ, ZZ

from .sage_backend import irrep, irrep_dimension, tensor_product


def _content_key(item):
    return tuple(tuple(x) for x in item)


@dataclass(frozen=True, init=False)
class RepresentationContent:
    """Immutable finite ZZ-linear combination of product-group irreps."""
    group_specs: tuple
    terms: tuple

    def __init__(self, group_specs, terms=()):
        specs = tuple(group_specs)
        combined = {}
        source = terms.items() if hasattr(terms, "items") else terms
        for key, coefficient in source:
            if isinstance(coefficient, float):
                raise ValueError("floating-point multiplicities are forbidden")
            key = _content_key(key)
            if len(key) != len(specs):
                raise ValueError("basis key must contain one label tuple per simple factor")
            # Construction validates labels and eagerly verifies Sage conventions.
            for spec, labels in zip(specs, key):
                irrep(spec, labels)
            coefficient = ZZ(coefficient)
            combined[key] = combined.get(key, ZZ(0)) + coefficient
        object.__setattr__(self, "group_specs", specs)
        object.__setattr__(self, "terms", tuple(sorted((k, v) for k, v in combined.items() if v)))

    @classmethod
    def zero(cls, group_specs):
        return cls(group_specs)

    @classmethod
    def single_irrep(cls, group_specs, labels, multiplicity=1):
        return cls(group_specs, ((tuple(tuple(x) for x in labels), multiplicity),))

    def __iter__(self):
        return iter(self.terms)

    def _combine(self, other, sign):
        if self.group_specs != other.group_specs:
            raise ValueError("representation contents use different simple factors")
        return type(self)(self.group_specs, self.terms + tuple((k, sign * c) for k, c in other))

    def __add__(self, other): return self._combine(other, ZZ(1))
    def __sub__(self, other): return self._combine(other, ZZ(-1))
    def __neg__(self): return self.scalar_mul(-1)
    def scalar_mul(self, scalar):
        if isinstance(scalar, float): raise ValueError("floating-point scalars are forbidden")
        scalar = ZZ(scalar)
        return type(self)(self.group_specs, ((k, scalar * c) for k, c in self))
    def __rmul__(self, scalar): return self.scalar_mul(scalar)

    def total_dimension(self):
        total = ZZ(0)
        for key, multiplicity in self:
            dimension = ZZ(1)
            for spec, labels in zip(self.group_specs, key):
                dimension *= irrep_dimension(spec, labels)
            total += multiplicity * dimension
        return total

    def __mul__(self, other):
        if self.group_specs != other.group_specs:
            raise ValueError("representation contents use different simple factors")
        result = []
        for left, lm in self:
            for right, rm in other:
                partial = [((), ZZ(1))]
                for spec, a, b in zip(self.group_specs, left, right):
                    partial = [(labels + (piece,), mult * piece_mult)
                               for labels, mult in partial
                               for piece, piece_mult in tensor_product(spec, a, b)]
                result.extend((key, lm * rm * mult) for key, mult in partial)
        return type(self)(self.group_specs, result)


@dataclass(frozen=True, init=False)
class CharacterSeries:
    theory: object
    terms: tuple

    def __init__(self, theory, terms=()):
        combined = {}
        source = terms.items() if hasattr(terms, "items") else terms
        for sector, content in source:
            degree, charges = sector
            sector = (ZZ(degree), tuple(sorted((key, QQ(value)) for key, value in charges)))
            combined[sector] = combined.get(sector, RepresentationContent.zero(theory.simple_factors)) + content
        object.__setattr__(self, "theory", theory)
        object.__setattr__(self, "terms", tuple(sorted((k, v) for k, v in combined.items() if v.terms)))

    def __iter__(self): return iter(self.terms)


def restore_characters(theory, hwg_series):
    """Interpret each complete HWG Dynkin tuple as one Sage irrep."""
    restored = []
    factor_ids = tuple(x.id for x in theory.simple_factors)
    for monomial, multiplicity in hwg_series:
        reps = {x.simple_factor_id: x.dynkin_labels for x in monomial.representations}
        labels = tuple(reps[x] for x in factor_ids)
        content = RepresentationContent.single_irrep(theory.simple_factors, labels, multiplicity)
        restored.append(((monomial.t_degree, monomial.abelian_charges), content))
    return CharacterSeries(theory, restored)


def dimension_refine(character_series):
    return tuple((sector, content.total_dimension()) for sector, content in character_series)


def unrefine(character_series):
    coefficients = {}
    for (degree, _charges), dimension in dimension_refine(character_series):
        coefficients[degree] = coefficients.get(degree, ZZ(0)) + dimension
    if not coefficients:
        return ()
    return tuple((ZZ(degree), coefficients.get(ZZ(degree), ZZ(0)))
                 for degree in range(int(max(coefficients)) + 1))
