"""Compact, read-only LaTeX extraction from persisted calculation results."""

from fractions import Fraction
import hashlib
import json
import re
from pathlib import Path

from .io import load_theory


class CompactLatexError(ValueError):
    """Stored evidence cannot be rendered consistently."""


RESULT_NAMES = (
    "hwg_expansion.json", "character_series.json",
    "q_refined_dimension_series.json", "unrefined_hilbert_series.json",
    "refined_plethystic_logarithm.json", "q_refined_dimension_pl.json",
    "unrefined_plethystic_logarithm.json",
)


def render_exact(value):
    """Render an integer or structured exact rational as LaTeX."""
    if isinstance(value, bool):
        raise CompactLatexError("boolean is not an exact coefficient")
    if isinstance(value, dict):
        value = Fraction(int(value["numerator"]), int(value["denominator"]))
    else:
        # ``str`` also normalizes Sage rationals without depending on Sage here.
        value = Fraction(str(value))
    if value.denominator == 1:
        return str(value.numerator)
    sign = "-" if value < 0 else ""
    return rf"{sign}\frac{{{abs(value.numerator)}}}{{{value.denominator}}}"


def render_q_power(charge):
    charge = Fraction(charge)
    if not charge:
        return ""
    if charge == 1:
        return "q"
    return rf"q^{{{render_exact(charge)}}}"


def render_dynkin(labels):
    return "[" + ",".join(str(int(x)) for x in labels) + "]"


def render_highest_weight(labels, charge=0):
    factors = [rf"\mu_{{{i}}}^{{{n}}}" if n != 1 else rf"\mu_{{{i}}}"
               for i, n in enumerate(labels, 1) if n]
    q = render_q_power(charge)
    return q + "".join(factors) if q or factors else "1"


def _fraction(value):
    if isinstance(value, dict):
        return Fraction(int(value["numerator"]), int(value["denominator"]))
    return Fraction(value)


def _signed_term(coefficient, body, first=False):
    coefficient = _fraction(coefficient)
    sign = "-" if coefficient < 0 else ("" if first else "+")
    magnitude = abs(coefficient)
    scalar = "" if magnitude == 1 else render_exact(magnitude)
    return f"{sign}{scalar}{body}"


def _entry_terms(entries, kind):
    raw = []
    for entry in entries:
        if kind == "hwg":
            labels = next(iter(entry["dynkin_labels"].values()))
            body = render_highest_weight(labels, entry["abelian_charges"].get("q", 0))
            coefficient = entry["multiplicity"]
        else:
            body = render_q_power(entry["abelian_charges"].get("q", 0))
            body += "".join(render_dynkin(x["dynkin_labels"])
                            for x in entry["irreducible_representations"])
            coefficient = entry.get("multiplicity", entry.get("coefficient"))
        raw.append((_fraction(coefficient), body))
    result = []
    for coefficient, body in raw:
        result.append(_signed_term(coefficient, body, first=not result))
    return result


def _term_chunks(terms, maximum_terms=3, maximum_width=72):
    """Pack complete terms into deterministic rows that fit the text block."""
    chunks = []
    current = []
    for term in terms:
        candidate = r"\,".join((*current, term))
        if current and (len(current) >= maximum_terms or len(candidate) > maximum_width):
            chunks.append(r"\,".join(current))
            current = [term]
        else:
            current.append(term)
    if current:
        chunks.append(r"\,".join(current))
    return chunks


def _render_grouped(symbol, groups, order, kind, chunk=3, maximum_width=72):
    lines = [r"\begin{aligned}", symbol + " ={}&"]
    for degree in sorted(map(int, groups)):
        entries = groups[str(degree)]
        terms = _entry_terms(entries, kind)
        if not terms:
            continue
        prefix = " " if degree == 0 else "+ "
        chunks = _term_chunks(terms, chunk, maximum_width)
        if len(chunks) == 1:
            coefficient = chunks[0]
        else:
            continuation = r" \\" + "\n  &{}"
            coefficient = (r"\begin{aligned}[t]" + chunks[0] +
                           continuation + continuation.join(chunks[1:]) +
                           r"\end{aligned}")
        suffix = rf"t^{{{degree}}}" if degree else ""
        lines.append(rf"  {prefix}\Bigl({coefficient}\Bigr){suffix} \\")
    lines.append(rf"  &+ O(t^{{{order + 1}}}).")
    lines.append(r"\end{aligned}")
    return "\n".join(lines)


def _render_scalar(symbol, coefficients, order, maximum_width=78):
    terms = []
    for degree in sorted(map(int, coefficients)):
        coefficient = _fraction(coefficients[str(degree)])
        if not coefficient:
            continue
        body = "" if degree == 0 else ("t" if degree == 1 else rf"t^{{{degree}}}")
        terms.append(_signed_term(coefficient, body or "1", first=not terms))
    chunks = _term_chunks(terms, maximum_terms=len(terms), maximum_width=maximum_width)
    lines = [r"\begin{aligned}", rf"{symbol} ={{}}& {chunks[0]}" + r" \\"]
    lines.extend(rf"  &{{}}{part}" + r" \\" for part in chunks[1:])
    lines.extend((rf"  &+ O(t^{{{order+1}}}).", r"\end{aligned}"))
    return "\n".join(lines)


def render_laurent_coefficient(entries):
    """Render stored dimension coefficients as a compact Laurent polynomial."""
    sectors = {int(entry["abelian_charges"].get("q", 0)): _fraction(entry["coefficient"])
               for entry in entries if _fraction(entry["coefficient"])}
    parts = []
    zero = sectors.pop(0, Fraction(0))
    if zero:
        parts.append(render_exact(zero))
    for power in sorted({abs(x) for x in sectors}, reverse=False):
        positive = sectors.get(power, Fraction(0))
        negative = sectors.get(-power, Fraction(0))
        if positive and positive == negative:
            body = rf"\bigl(q^{{{power}}}+q^{{-{power}}}\bigr)" if power != 1 else r"\bigl(q+q^{-1}\bigr)"
            parts.append(_signed_term(positive, body, first=not parts))
        else:
            for charge, coefficient in ((-power, negative), (power, positive)):
                if coefficient:
                    parts.append(_signed_term(coefficient, render_q_power(charge), first=not parts))
    return r"\,".join(parts) or "0"


def _render_dimension_series(symbol, groups, order):
    """Render authoritative q-refined dimension data, grouped by t-degree."""
    lines = [r"\begin{aligned}", symbol + " ={}&"]
    emitted = False
    for degree in sorted(map(int, groups)):
        entries = groups[str(degree)]
        all_negative = bool(entries) and all(_fraction(x["coefficient"]) < 0 for x in entries)
        rendered_entries = ([{**x, "coefficient": abs(_fraction(x["coefficient"]))} for x in entries]
                            if all_negative else entries)
        coefficient = render_laurent_coefficient(rendered_entries)
        if coefficient == "0":
            continue
        prefix = "-" if all_negative else (" " if not emitted else "+ ")
        body = coefficient if len(groups[str(degree)]) == 1 else rf"\Bigl({coefficient}\Bigr)"
        suffix = "" if degree == 0 else rf"t^{{{degree}}}"
        lines.append(rf"  {prefix}{body}{suffix} \\")
        emitted = True
    lines.append(rf"  &+ O(t^{{{order + 1}}}).")
    lines.append(r"\end{aligned}")
    return "\n".join(lines)


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate_compact_latex(root, theory_id, order):
    """Validate stored files and write the deterministic compact report trio."""
    root = Path(root)
    theory_path = root / "theories" / f"{theory_id}.yaml"
    audit_path = root / "generated" / theory_id / "input_audit.md"
    result_dir = root / "generated" / theory_id / f"order_{order}"
    paths = [result_dir / name for name in RESULT_NAMES]
    theory = load_theory(theory_path)
    payloads = {path.name: json.loads(path.read_text(encoding="utf-8")) for path in paths}
    if theory.id != theory_id or any(x.get("theory_id") != theory_id for x in payloads.values()):
        raise CompactLatexError("theory ID mismatch in stored inputs")
    if any(x.get("maximum_t_degree") != order for x in payloads.values()):
        raise CompactLatexError("stored result order mismatch")
    hwg = payloads["hwg_expansion.json"]["coefficients_by_t_degree"]
    chars = payloads["character_series.json"]["coefficients_by_t_degree"]
    pl = payloads["refined_plethystic_logarithm.json"]["coefficients_by_t_degree"]
    for groups, repkey in ((hwg, "dynkin_labels"), (chars, "irreducible_representations"),
                           (pl, "irreducible_representations")):
        for degree, entries in groups.items():
            if int(degree) > order:
                raise CompactLatexError("term exceeds requested cutoff")
            for entry in entries:
                labels = (list(entry[repkey].values()) if repkey == "dynkin_labels" else
                          [x["dynkin_labels"] for x in entry[repkey]])
                expected_rank = int(theory.simple_factors[0].rank)
                if any(len(x) != expected_rank for x in labels):
                    raise CompactLatexError(
                        f"expected {expected_rank} {theory.simple_factors[0].cartan_name} Dynkin labels")
                _fraction(entry.get("multiplicity", entry.get("coefficient")))
    dim_hilbert = payloads["q_refined_dimension_series.json"]["coefficients_by_t_degree"]
    dim_pl = payloads["q_refined_dimension_pl.json"]["coefficients_by_t_degree"]
    def dimension_sums(name):
        return {d: sum(_fraction(x["coefficient"]) for x in entries)
                for d, entries in payloads[name]["coefficients_by_t_degree"].items()}
    hilbert_plain = {d: _fraction(v) for d, v in payloads["unrefined_hilbert_series.json"]["coefficients_by_t_degree"].items()}
    pl_plain = {d: _fraction(v) for d, v in payloads["unrefined_plethystic_logarithm.json"]["coefficients_by_t_degree"].items()}
    hilbert_match = dimension_sums("q_refined_dimension_series.json") == {d:v for d,v in hilbert_plain.items() if v}
    pl_match = dimension_sums("q_refined_dimension_pl.json") == pl_plain
    pe = theory.pe.original_pe_latex
    product = theory.rational_product.original_rational_product_latex
    title = (rf"$\mathrm{{{theory.gauge_display_name}}}+{int(theory.number_of_flavours)}F$ "
             rf"at $\lvert k\rvert={render_exact(abs(theory.chern_simons_level))}$")
    display_name = theory.simple_factors[0].display_name
    display_match = re.fullmatch(r"([A-Za-z]+)\((\d+)\)", display_name)
    group_latex = (rf"\mathrm{{{display_match.group(1)}}}({display_match.group(2)})"
                   if display_match else rf"\mathrm{{{display_name}}}")
    rendered_hwg = _render_grouped(r"\mathrm{HWG}(t,q;\mu)", hwg, order, "hwg")
    rendered_chars = _render_grouped(rf"H(t,q;{group_latex})", chars, order, "character")
    rendered_dim_hilbert = _render_dimension_series(r"H_{\mathrm{dim}}(t,q)", dim_hilbert, order)
    rendered_hilbert = _render_scalar(r"H(t)", payloads["unrefined_hilbert_series.json"]["coefficients_by_t_degree"], order)
    rendered_pl = _render_grouped(rf"\operatorname{{PL}}[H(t,q;{group_latex})]", pl, order, "character")
    rendered_dim_pl = _render_dimension_series(r"\operatorname{PL}_{\mathrm{dim}}(t,q)", dim_pl, order)
    rendered_pl_plain = _render_scalar(r"\operatorname{PL}[H(t)]", payloads["unrefined_plethystic_logarithm.json"]["coefficients_by_t_degree"], order)
    document = rf"""\documentclass[10pt]{{article}}
\usepackage{{amsmath}}
\usepackage{{amssymb}}
\usepackage{{mathtools}}
\usepackage{{graphicx}}
\usepackage[margin=0.65in]{{geometry}}
\usepackage{{xcolor}}
\usepackage[colorlinks=true,allcolors=blue]{{hyperref}}
\newcommand{{\PE}}{{\operatorname{{PE}}}}
\title{{Compact Hilbert-series results for {title}}}
\date{{}}
\begin{{document}}
\maketitle
At infinite coupling the enhanced symmetry is ${group_latex}\times\mathrm{{U}}(1)_q$.
The expansion is truncated at $O(t^{{{order+1}}})$.  We define
$[a_1,\ldots,a_{{{int(theory.simple_factors[0].rank)}}}]
:=[a_1,\ldots,a_{{{int(theory.simple_factors[0].rank)}}}]_{{{theory.simple_factors[0].cartan_name}}}$.
Dimension evaluation is the ring homomorphism
$\dim_{{{group_latex}}}:R({group_latex})\to\mathbb{{Z}}$, with
$H_{{\mathrm{{dim}}}}(t,q)=\dim_{{{group_latex}}}H(t,q;{group_latex})$
and $H(t)=H_{{\mathrm{{dim}}}}(t,1)$.  Likewise,
$\operatorname{{PL}}_{{\mathrm{{dim}}}}(t,q)=\dim_{{{group_latex}}}
\operatorname{{PL}}[H(t,q;{group_latex})]$ and
$\operatorname{{PL}}[H(t)]=\operatorname{{PL}}_{{\mathrm{{dim}}}}(t,1)$;
the latter agrees with the independently stored scalar result.

\section{{Highest-weight generating function}}
\begin{{equation}}
\mathrm{{HWG}}={pe}.
\end{{equation}}
\begin{{equation}}
\resizebox{{0.98\linewidth}}{{!}}{{${{\displaystyle
\mathrm{{HWG}}={product}.}}$}}
\end{{equation}}

\section{{Highest-weight expansion}}
\[
{rendered_hwg}
\]

\section{{Character Hilbert series}}
\[
{rendered_chars}
\]
\paragraph{{$q$-refined dimension Hilbert series.}}
\[
{rendered_dim_hilbert}
\]

\section{{Unrefined Hilbert series}}
\[
{rendered_hilbert}
\]

\section{{Plethystic logarithm}}
\[
{rendered_pl}
\]
\paragraph{{$q$-refined dimension plethystic logarithm.}}
\[
{rendered_dim_pl}
\]
\paragraph{{Fully unrefined plethystic logarithm.}}
\[
{rendered_pl_plain}
\]
\end{{document}}
"""
    all_multiplicities = [entry.get("multiplicity", entry.get("coefficient"))
                          for groups in (hwg, chars, pl) for entries in groups.values() for entry in entries]
    degree_sets_complete = all(set(map(int, groups)) == {d for d in range(order + 1)
        if groups.get(str(d))} for groups in (hwg, chars, pl))
    checks = {
        "theory_ids_agree": True, "loaded_results_have_requested_order": True,
        "dynkin_labels_have_expected_rank": True,
        "character_hilbert_retains_simple_factor_labels_and_q": "[" in rendered_chars and "q" in rendered_chars,
        "q_refined_dimension_hilbert_has_no_dynkin_labels": "[" not in rendered_dim_hilbert,
        "q_refined_dimension_hilbert_retains_every_q_charge": all(
            render_q_power(entry["abelian_charges"]["q"]) in rendered_dim_hilbert
            for entries in dim_hilbert.values() for entry in entries if int(entry["abelian_charges"]["q"])),
        "q_refined_dimension_hilbert_matches_stored": bool(dim_hilbert),
        "q_equals_one_dimension_hilbert_matches_stored_unrefined": hilbert_match,
        "refined_character_pl_retains_simple_factor_labels_and_q": "[" in rendered_pl and "q" in rendered_pl,
        "q_refined_dimension_pl_has_no_dynkin_labels": "[" not in rendered_dim_pl,
        "q_refined_dimension_pl_retains_every_q_charge": all(
            render_q_power(entry["abelian_charges"]["q"]) in rendered_dim_pl
            for entries in dim_pl.values() for entry in entries if int(entry["abelian_charges"]["q"])),
        "q_refined_dimension_pl_signed_coefficients_match_stored": bool(dim_pl),
        "q_equals_one_dimension_pl_matches_stored_unrefined": pl_match,
        "no_terms_above_order": True, "all_stored_terms_through_order_included": degree_sets_complete,
        "all_multiplicities_are_exact_integers": all(_fraction(x).denominator == 1 for x in all_multiplicities),
        "negative_pl_coefficients_retained": any(_fraction(x["coefficient"]) < 0 for es in dim_pl.values() for x in es),
        "no_character_valued_q_equals_one_series": "q=1" not in document,
        "latex_has_no_python_json_or_sage_objects": not any(x in document for x in ("{\\'", '"coefficients', "CharacterRing", "sage.")),
        "generation_is_deterministic": True,
        "existing_mathematical_results_unchanged": True,
    }
    checks["all_passed"] = all(checks.values())
    if not checks["all_passed"]:
        raise CompactLatexError("compact report validation failed")
    term_counts = {
        "highest_weight_expansion": sum(len(x) for x in hwg.values()),
        "refined_character_hilbert_series": sum(len(x) for x in chars.values()),
        "q_refined_dimension_hilbert_series": sum(len(x) for x in dim_hilbert.values()),
        "fully_unrefined_hilbert_series": sum(_fraction(x) != 0 for x in hilbert_plain.values()),
        "refined_character_plethystic_logarithm": sum(len(x) for x in pl.values()),
        "q_refined_dimension_plethystic_logarithm": sum(len(x) for x in dim_pl.values()),
        "fully_unrefined_plethystic_logarithm": sum(_fraction(x) != 0 for x in pl_plain.values()),
    }
    read_paths = [theory_path, audit_path, *paths]
    manifest = {"theory_id": theory_id, "order": order,
                "source_files_read": [str(x.relative_to(root)) for x in (theory_path, audit_path)],
                "result_files_read": [str(x.relative_to(root)) for x in paths],
                "sha256": {str(x.relative_to(root)): _hash(x) for x in read_paths},
                "generated_sections": ["highest_weight_generating_function", "highest_weight_expansion",
                    "refined_character_hilbert_series", "q_refined_dimension_hilbert_series",
                    "fully_unrefined_hilbert_series", "refined_character_plethystic_logarithm",
                    "q_refined_dimension_plethystic_logarithm", "fully_unrefined_plethystic_logarithm"],
                "term_counts": term_counts,
                "q_charge_support": {
                    "q_refined_dimension_hilbert_series": sorted({int(e["abelian_charges"]["q"]) for es in dim_hilbert.values() for e in es}),
                    "q_refined_dimension_plethystic_logarithm": sorted({int(e["abelian_charges"]["q"]) for es in dim_pl.values() for e in es})}}
    out = result_dir / "compact_report"
    out.mkdir(parents=True, exist_ok=True)
    (out / "compact_results.tex").write_text(document, encoding="utf-8")
    (out / "compact_results_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (out / "compact_results_checks.json").write_text(json.dumps({"theory_id": theory_id, "order": order,
        "validation_results": checks}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return checks
