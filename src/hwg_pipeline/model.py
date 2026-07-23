"""Immutable structured data for highest-weight generating functions."""

from dataclasses import dataclass
import re

from sage.all import QQ, ZZ


_CARTAN = re.compile(r"^[A-G]$")


@dataclass(frozen=True)
class SourceReference:
    path: str
    description: str
    equation: str | None = None


@dataclass(frozen=True)
class SimpleGroupSpec:
    id: str
    cartan_type: str
    rank: int
    display_name: str
    highest_weight_fugacities: tuple[str, ...]

    def __post_init__(self):
        if not _CARTAN.fullmatch(self.cartan_type) or self.rank <= 0:
            raise ValueError(f"malformed Cartan type/rank: {self.cartan_type}{self.rank}")
        if len(self.highest_weight_fugacities) != self.rank:
            raise ValueError(f"simple factor {self.id!r}: expected {self.rank} highest-weight fugacities")
        if len(set(self.highest_weight_fugacities)) != self.rank:
            raise ValueError(f"simple factor {self.id!r}: highest-weight fugacities must be unique")

    @property
    def cartan_name(self):
        return f"{self.cartan_type}{self.rank}"


@dataclass(frozen=True)
class AbelianFactorSpec:
    id: str
    display_name: str
    fugacity: str


@dataclass(frozen=True)
class RepresentationSpec:
    simple_factor_id: str
    dynkin_labels: tuple

    def __post_init__(self):
        labels = tuple(ZZ(x) for x in self.dynkin_labels)
        if any(x < 0 for x in labels):
            raise ValueError("Dynkin labels must be nonnegative integers")
        object.__setattr__(self, "dynkin_labels", labels)


@dataclass(frozen=True)
class HighestWeightMonomial:
    t_degree: int
    representations: tuple[RepresentationSpec, ...]
    abelian_charges: tuple[tuple[str, object], ...]

    def __post_init__(self):
        degree = ZZ(self.t_degree)
        if degree < 0:
            raise ValueError("t-degree must be nonnegative")
        object.__setattr__(self, "t_degree", degree)
        object.__setattr__(self, "abelian_charges", tuple((key, QQ(value)) for key, value in self.abelian_charges))


@dataclass(frozen=True)
class HWGTerm:
    coefficient: int
    monomial: HighestWeightMonomial

    def __post_init__(self):
        try:
            coefficient = ZZ(self.coefficient)
        except (TypeError, ValueError):
            raise ValueError("HWG term coefficient must be an integer") from None
        if coefficient != self.coefficient:
            raise ValueError("HWG term coefficient must be an integer")
        object.__setattr__(self, "coefficient", coefficient)


@dataclass(frozen=True)
class PlethysticExponentialSpec:
    terms: tuple[HWGTerm, ...]
    original_pe_latex: str


@dataclass(frozen=True)
class RationalProductFactor:
    monomial: HighestWeightMonomial
    power: int

    def __post_init__(self):
        power = ZZ(self.power)
        if not power:
            raise ValueError("rational-product factor power must be a nonzero integer")
        object.__setattr__(self, "power", power)


@dataclass(frozen=True)
class RationalProductSpec:
    factors: tuple[RationalProductFactor, ...]
    original_rational_product_latex: str | None = None


@dataclass(frozen=True)
class TheorySpec:
    id: str
    title: str
    nonphysical: bool
    source_references: tuple[SourceReference, ...]
    simple_factors: tuple[SimpleGroupSpec, ...]
    abelian_factors: tuple[AbelianFactorSpec, ...]
    chern_simons_level: object
    pe: PlethysticExponentialSpec
    rational_product: RationalProductSpec | None = None

    def __post_init__(self):
        object.__setattr__(self, "chern_simons_level", QQ(self.chern_simons_level))
        simple = {factor.id: factor for factor in self.simple_factors}
        abelian = {factor.id: factor for factor in self.abelian_factors}
        if len(simple) != len(self.simple_factors) or len(abelian) != len(self.abelian_factors):
            raise ValueError("factor identifiers must be unique")
        all_fugacities = [x for f in self.simple_factors for x in f.highest_weight_fugacities]
        all_fugacities += [f.fugacity for f in self.abelian_factors]
        if len(set(all_fugacities)) != len(all_fugacities):
            raise ValueError("fugacity names must be unique")
        monomials = [term.monomial for term in self.pe.terms]
        if self.rational_product:
            monomials += [factor.monomial for factor in self.rational_product.factors]
        for monomial in monomials:
            reps = {rep.simple_factor_id: rep for rep in monomial.representations}
            charges = dict(monomial.abelian_charges)
            unknown_simple = set(reps) - set(simple)
            unknown_abelian = set(charges) - set(abelian)
            if unknown_simple:
                raise ValueError(f"undeclared simple factor(s): {', '.join(sorted(unknown_simple))}")
            if unknown_abelian:
                raise ValueError(f"undeclared abelian factor(s): {', '.join(sorted(unknown_abelian))}")
            if set(reps) != set(simple):
                raise ValueError("every monomial must give Dynkin labels for every simple factor")
            if set(charges) != set(abelian):
                raise ValueError("every monomial must give charges for every abelian factor")
            for factor_id, rep in reps.items():
                if len(rep.dynkin_labels) != simple[factor_id].rank:
                    raise ValueError(f"simple factor {factor_id!r}: Dynkin-label length {len(rep.dynkin_labels)} does not match rank {simple[factor_id].rank}")
