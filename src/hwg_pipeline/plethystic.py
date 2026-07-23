"""Exact virtual-character arithmetic and truncated plethystic logarithms."""

from dataclasses import dataclass
from math import prod

from sage.all import QQ, ZZ, moebius

from .characters import RepresentationContent
from .sage_backend import adams_decomposition, irrep_dimension, tensor_product


def _key(key):
    return tuple(tuple(ZZ(x) for x in labels) for labels in key)


@dataclass(frozen=True, init=False)
class VirtualRepresentationContent:
    """Immutable QQ-linear combination of product-group irreducibles."""
    group_specs: tuple
    terms: tuple

    def __init__(self, group_specs, terms=()):
        specs, combined = tuple(group_specs), {}
        for key, coefficient in (terms.items() if hasattr(terms, "items") else terms):
            if isinstance(coefficient, float):
                raise ValueError("floating-point coefficients are forbidden")
            key, coefficient = _key(key), QQ(coefficient)
            if len(key) != len(specs):
                raise ValueError("basis key must contain one tuple per simple factor")
            combined[key] = combined.get(key, QQ.zero()) + coefficient
        object.__setattr__(self, "group_specs", specs)
        object.__setattr__(self, "terms", tuple(sorted((k, c) for k, c in combined.items() if c)))

    @classmethod
    def zero(cls, specs): return cls(specs)
    @classmethod
    def trivial(cls, specs): return cls(specs, [((tuple(0 for _ in range(s.rank)) for s in specs), 1)])
    @classmethod
    def from_representation_content(cls, content): return cls(content.group_specs, content.terms)
    @classmethod
    def single_irrep(cls, specs, labels, coefficient=1): return cls(specs, [(labels, coefficient)])
    def __iter__(self): return iter(self.terms)
    def _combine(self, other, sign=1):
        if self.group_specs != other.group_specs: raise ValueError("different simple factors")
        return type(self)(self.group_specs, self.terms + tuple((k, sign*c) for k, c in other))
    def __add__(self, other): return self._combine(other)
    def __sub__(self, other): return self._combine(other, -1)
    def __neg__(self): return self.scalar_mul(-1)
    def scalar_mul(self, scalar):
        if isinstance(scalar, float): raise ValueError("floating-point coefficients are forbidden")
        return type(self)(self.group_specs, ((k, QQ(scalar)*c) for k, c in self))
    def __rmul__(self, scalar): return self.scalar_mul(scalar)
    def __mul__(self, other):
        if self.group_specs != other.group_specs: raise ValueError("different simple factors")
        result = []
        for left, lc in self:
            for right, rc in other:
                pieces = [((), ZZ.one())]
                for spec, a, b in zip(self.group_specs, left, right):
                    pieces = [(labels + (piece,), mult*pm) for labels, mult in pieces
                              for piece, pm in tensor_product(spec, a, b)]
                result.extend((labels, lc*rc*mult) for labels, mult in pieces)
        return type(self)(self.group_specs, result)
    def total_dimension(self):
        return sum((c * prod(irrep_dimension(s, l) for s, l in zip(self.group_specs, key))
                    for key, c in self), QQ.zero())


@dataclass(frozen=True, init=False)
class VirtualCharacterSeries:
    theory: object
    max_degree: object
    terms: tuple
    def __init__(self, theory, terms=(), max_degree=None):
        maximum = ZZ(max_degree if max_degree is not None else max((s[0] for s, _ in terms), default=0))
        combined = {}
        for (degree, charges), content in terms:
            degree = ZZ(degree)
            if degree < 0: raise ValueError("t-degree must be nonnegative")
            if degree > maximum: continue
            sector = degree, tuple(sorted((str(k), QQ(v)) for k, v in charges))
            combined[sector] = combined.get(sector, VirtualRepresentationContent.zero(theory.simple_factors)) + content
        object.__setattr__(self, "theory", theory); object.__setattr__(self, "max_degree", maximum)
        object.__setattr__(self, "terms", tuple(sorted((s, c) for s, c in combined.items() if c.terms)))
    def __iter__(self): return iter(self.terms)
    @classmethod
    def from_character_series(cls, series, max_degree):
        return cls(series.theory, ((s, VirtualRepresentationContent.from_representation_content(c)) for s,c in series), max_degree)
    def _combine(self, other, sign=1):
        return type(self)(self.theory, self.terms + tuple((s, sign*c) for s,c in other), min(self.max_degree, other.max_degree))
    def __add__(self, other): return self._combine(other)
    def __sub__(self, other): return self._combine(other, -1)
    def scalar_mul(self, scalar): return type(self)(self.theory, ((s, scalar*c) for s,c in self), self.max_degree)
    def __rmul__(self, scalar): return self.scalar_mul(scalar)
    def __mul__(self, other):
        terms=[]
        for (d,a), x in self:
            for (e,b), y in other:
                if d + e > min(self.max_degree, other.max_degree):
                    continue
                charges = tuple((k, dict(a)[k]+dict(b)[k]) for k,_ in a)
                terms.append(((d+e, charges), x*y))
        return type(self)(self.theory, terms, min(self.max_degree, other.max_degree))


def _validate_unit(series):
    zero = [(charges, content) for (degree, charges), content in series if degree == 0]
    expected_charges = tuple((a.id, QQ.zero()) for a in series.theory.abelian_factors)
    if len(zero) != 1 or zero[0][0] != expected_charges or zero[0][1] != VirtualRepresentationContent.trivial(series.theory.simple_factors):
        raise ValueError("degree-zero term must be exactly the neutral trivial representation")


def formal_logarithm(series, max_degree=None):
    maximum = ZZ(series.max_degree if max_degree is None else max_degree)
    series = VirtualCharacterSeries(series.theory, series.terms, maximum); _validate_unit(series)
    unit = VirtualCharacterSeries(series.theory, [((0, tuple((a.id, 0) for a in series.theory.abelian_factors)),
             VirtualRepresentationContent.trivial(series.theory.simple_factors))], maximum)
    u, power, result = series-unit, series-unit, VirtualCharacterSeries(series.theory, (), maximum)
    for n in range(1, int(maximum)+1):
        result = result + ((QQ((-1)**(n+1))/n) * power)
        power = power*u
        if not power.terms: break
    return result


def adams_series(series, k, max_degree):
    if isinstance(k, bool) or k not in ZZ or ZZ(k) <= 0: raise ValueError("k must be a positive integer")
    k, maximum, out = ZZ(k), ZZ(max_degree), []
    if k == 1: return VirtualCharacterSeries(series.theory, series.terms, maximum)
    for (degree, charges), content in series:
        if degree*k > maximum: continue
        transformed=[]
        for labels, coefficient in content:
            pieces=[((), coefficient)]
            for spec, factor_labels in zip(content.group_specs, labels):
                pieces=[(key+(piece,), c*m) for key,c in pieces
                        for piece,m in adams_decomposition(spec, factor_labels, k)]
            transformed.extend(pieces)
        out.append(((degree*k, tuple((q, v*k) for q,v in charges)), VirtualRepresentationContent(content.group_specs, transformed)))
    return VirtualCharacterSeries(series.theory, out, maximum)


def plethystic_logarithm(character_series, max_degree):
    source = character_series if isinstance(character_series, VirtualCharacterSeries) else VirtualCharacterSeries.from_character_series(character_series, max_degree)
    result = VirtualCharacterSeries(source.theory, (), max_degree)
    for k in range(1, int(max_degree)+1):
        if moebius(k): result = result + (QQ(moebius(k))/k)*formal_logarithm(adams_series(source,k,max_degree), max_degree)
    return result


def dimension_refine_virtual(series): return tuple((sector, content.total_dimension()) for sector,content in series)
def unrefine_virtual(series):
    out={}
    for (degree,_), value in dimension_refine_virtual(series): out[degree]=out.get(degree,QQ.zero())+value
    return tuple(sorted(out.items()))


def scalar_plethystic_logarithm(coefficients, max_degree):
    """Independent scalar power-series implementation."""
    R = QQ['t']; t=R.gen(); H=sum(QQ(c)*t**ZZ(d) for d,c in coefficients)
    answer=R.zero()
    for k in range(1,int(max_degree)+1):
        if not moebius(k): continue
        u=R(H(t=t**k))-1; power=u
        for n in range(1,int(max_degree)+1):
            answer += (QQ(moebius(k))/k)*(QQ((-1)**(n+1))/n)*power
            power=(power*u).truncate(max_degree+1)
    return tuple((ZZ(d), QQ(answer[d])) for d in range(1,int(max_degree)+1) if answer[d])
