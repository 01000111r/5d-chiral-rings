"""Structured, exact theory data for the HWG pipeline."""

from .io import dump_theory, load_theory, theory_from_dict, theory_to_dict
from .model import (AbelianFactorSpec, HighestWeightMonomial, HWGTerm,
                    PlethysticExponentialSpec, RationalProductFactor,
                    RationalProductSpec, RepresentationSpec, SimpleGroupSpec,
                    SourceReference, TheorySpec)
from .render import render_monomial, render_pe, render_pe_exponent, render_rational_product
from .expansion import (SparseSeries, expand_hwg, expand_pe,
                        expand_rational_product, unit_monomial)

__version__ = "0.1.0"
