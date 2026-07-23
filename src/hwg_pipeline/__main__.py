"""Command-line entry point for the exact HWG expansion milestone."""

import argparse
import json
from pathlib import Path

from sage.all import ZZ

from .expansion import expand_hwg, expand_pe, expand_rational_product
from .io import load_theory
from .render import render_monomial
from .characters import dimension_refine, restore_characters, unrefine


def _root():
    current = Path.cwd().resolve()
    for directory in (current, *current.parents):
        if (directory / "theories").is_dir() and (directory / "pyproject.toml").is_file():
            return directory
    raise FileNotFoundError("could not locate project root")


def _groups(series):
    groups = {}
    for monomial, coefficient in series:
        groups.setdefault(str(monomial.t_degree), []).append({
            "dynkin_labels": {x.simple_factor_id: [int(y) for y in x.dynkin_labels]
                               for x in monomial.representations},
            "abelian_charges": {key: str(value) for key, value in monomial.abelian_charges},
            "multiplicity": int(coefficient),
        })
    return groups


def _validations(theory, order, pe, product):
    stability = {}
    for degree in sorted(set((ZZ(0), ZZ(2), ZZ(4), ZZ(6), ZZ(8), ZZ(order)))):
        if degree <= order:
            stability[str(degree)] = (expand_hwg(theory, degree) ==
                                      expand_hwg(theory, degree + 1).truncate(degree))
    unit = next((coefficient for monomial, coefficient in pe if monomial.t_degree == 0), ZZ(0))
    checks = {
        "pe_equals_rational_product": pe == product,
        "constant_coefficient_is_one": unit == 1,
        "no_terms_above_requested_degree": all(m.t_degree <= order for m, _ in pe),
        "all_multiplicities_are_integers": all(c in ZZ for _, c in pe),
        "all_multiplicities_are_nonnegative": all(c >= 0 for _, c in pe),
        "truncation_stability": stability,
    }
    checks["all_passed"] = all(value if isinstance(value, bool) else all(value.values())
                                for value in checks.values())
    return checks


def _write_outputs(theory, order, series, checks, output):
    output.mkdir(parents=True, exist_ok=True)
    payload = {"theory_id": theory.id, "maximum_t_degree": int(order),
               "expansion_route": "plethystic_exponential", "coefficients_by_t_degree": _groups(series),
               "validation_results": checks}
    (output / "hwg_expansion.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    lines = [rf"% Exact highest-weight monomial expansion through $t^{{{order}}}$.",
             rf"\begin{{align*}}", "H(t) = {}"]
    rendered = []
    for monomial, coefficient in series:
        body = render_monomial(monomial, theory)
        rendered.append(("" if not rendered else "+ ") + ("" if coefficient == 1 else f"{coefficient} ") + body)
    lines[-1] += " \\\\\n  ".join(rendered) if rendered else "0"
    lines += [r"\end{align*}", ""]
    (output / "hwg_expansion.tex").write_text("\n".join(lines), encoding="utf-8")
    check_payload = {"theory_id": theory.id, "maximum_t_degree": int(order),
                     "routes_compared": ["plethystic_exponential", "rational_product"],
                     "validation_results": checks}
    (output / "checks.json").write_text(json.dumps(check_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md = [f"# Expansion checks: `{theory.id}` through t^{order}", "",
          "Independent routes: structured PE and structured rational product.", ""]
    for key, value in checks.items():
        if isinstance(value, dict):
            md.append(f"- **{'PASS' if all(value.values()) else 'FAIL'} — {key}**: " +
                      ", ".join(f"D={d}: {'PASS' if passed else 'FAIL'}" for d, passed in value.items()))
        else:
            md.append(f"- **{'PASS' if value else 'FAIL'} — {key}**")
    (output / "checks.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def _charge_dict(charges):
    return {key: str(value) for key, value in charges}


def _character_outputs(theory, order, hwg, output):
    series = restore_characters(theory, hwg)
    refined = dimension_refine(series)
    plain = unrefine(series)
    groups = {}
    latex_groups = {}
    for (degree, charges), content in series:
        entries = groups.setdefault(str(degree), [])
        rendered = latex_groups.setdefault(int(degree), [])
        for labels, multiplicity in content:
            entry = {
                "abelian_charges": _charge_dict(charges),
                "irreducible_representations": [
                    {"cartan_factor_id": factor.id,
                     "dynkin_labels": [int(x) for x in factor_labels]}
                    for factor, factor_labels in zip(theory.simple_factors, labels)],
                "multiplicity": int(multiplicity),
            }
            entries.append(entry)
            reps = " ".join("[" + ",".join(str(x) for x in factor_labels) +
                            rf"]_{{{factor.cartan_type}_{factor.rank}}}"
                            for factor, factor_labels in zip(theory.simple_factors, labels))
            charge = " ".join(_latex_power(key, value) for key, value in charges if value)
            rendered.append(("" if multiplicity == 1 else str(multiplicity), reps, charge))
    character_payload = {"theory_id": theory.id, "maximum_t_degree": int(order),
                         "simple_factors": [{"cartan_factor_id": x.id,
                                             "cartan_type": x.cartan_name}
                                            for x in theory.simple_factors],
                         "coefficients_by_t_degree": groups}
    (output / "character_series.json").write_text(json.dumps(character_payload, indent=2, sort_keys=True) + "\n")
    lines = [rf"% Irreducible-character series through $t^{{{order}}}$.", r"\begin{align*}", "H(t,q) = {}"]
    pieces = []
    for degree in sorted(latex_groups):
        coefficient = " + ".join(" ".join(x for x in entry if x) for entry in latex_groups[degree])
        pieces.append(rf"\left({coefficient}\right)t^{{{degree}}}")
    lines[-1] += " \\\\\n  + ".join(pieces) if pieces else "0"
    lines += [r"\end{align*}", ""]
    (output / "character_series.tex").write_text("\n".join(lines))

    refined_groups = {}
    refined_latex = []
    for (degree, charges), coefficient in refined:
        refined_groups.setdefault(str(degree), []).append({
            "abelian_charges": _charge_dict(charges), "coefficient": int(coefficient)})
        refined_latex.append((int(degree), charges, coefficient))
    refined_payload = {"theory_id": theory.id, "maximum_t_degree": int(order),
                       "coefficients_by_t_degree": refined_groups}
    (output / "q_refined_dimension_series.json").write_text(json.dumps(refined_payload, indent=2, sort_keys=True) + "\n")
    rtex = [rf"% Exact q-refined dimension series through $t^{{{order}}}$.", r"\begin{align*}",
            "H_{\\dim}(t,q) = " + " + ".join(
                f"{coefficient}" + "".join(_latex_power(k, v) for k, v in charges if v) + rf"t^{{{degree}}}"
                for degree, charges, coefficient in refined_latex), r"\end{align*}", ""]
    (output / "q_refined_dimension_series.tex").write_text("\n".join(rtex))

    plain_payload = {"theory_id": theory.id, "maximum_t_degree": int(order),
                     "coefficients_by_t_degree": {str(d): int(c) for d, c in plain}}
    (output / "unrefined_hilbert_series.json").write_text(json.dumps(plain_payload, indent=2, sort_keys=True) + "\n")
    utex = [rf"% Exact unrefined Hilbert series through $t^{{{order}}}$.", r"\begin{align*}",
            "H(t) = " + " + ".join(f"{c}t^{{{d}}}" for d, c in plain), r"\end{align*}", ""]
    (output / "unrefined_hilbert_series.tex").write_text("\n".join(utex))

    checks = {
        "one_irrep_per_hwg_monomial_before_combining": len(hwg) == sum(len(c.terms) for _, c in series),
        "t_degrees_preserved": sorted(set(m.t_degree for m, _ in hwg)) == sorted(set(s[0] for s, _ in series)),
        "q_charges_preserved": sorted(set(m.abelian_charges for m, _ in hwg)) == sorted(set(s[1] for s, _ in series)),
        "representation_multiplicities_nonnegative_integers": all(c in ZZ and c >= 0 for _, content in series for _, c in content),
        "dimensions_nonnegative_integers": all(c in ZZ and c >= 0 for _, c in refined),
        "q_equals_one_matches_unrefined": dict(plain) == {d: sum(c for (sd, _), c in refined if sd == d) for d, _ in plain},
    }
    expected = {0: 1, 1: 0, 2: 36, 3: 30, 4: 630}
    actual = {int(d): int(c) for d, c in plain}
    checks["leading_coefficients"] = all(actual.get(d, 0) == c for d, c in expected.items())
    checks["all_passed"] = all(checks.values())
    check_payload = {"theory_id": theory.id, "maximum_t_degree": int(order), "validation_results": checks}
    (output / "character_checks.json").write_text(json.dumps(check_payload, indent=2, sort_keys=True) + "\n")
    md = [f"# Character checks: `{theory.id}` through t^{order}", ""] + [
        f"- **{'PASS' if value else 'FAIL'} — {key}**" for key, value in checks.items()]
    (output / "character_checks.md").write_text("\n".join(md) + "\n")
    return checks


def _latex_power(name, value):
    if value == 1: return name
    return rf"{name}^{{{value}}}"


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m hwg_pipeline")
    commands = parser.add_subparsers(dest="command", required=True)
    for name, help_text in (("expand", "expand a structured HWG"),
                            ("characters", "restore irreducible characters and dimensions")):
        command = commands.add_parser(name, help=help_text)
        command.add_argument("theory_id")
        command.add_argument("--order", required=True, type=int)
    args = parser.parse_args(argv)
    if args.order < 0:
        parser.error("--order must be nonnegative")
    root = _root()
    theory = load_theory(root / "theories" / f"{args.theory_id}.yaml")
    pe = expand_pe(theory, args.order)
    output = root / "generated" / theory.id / f"order_{args.order}"
    if args.command == "expand":
        product = expand_rational_product(theory, args.order)
        checks = _validations(theory, args.order, pe, product)
        _write_outputs(theory, args.order, pe, checks, output)
    else:
        output.mkdir(parents=True, exist_ok=True)
        checks = _character_outputs(theory, args.order, pe, output)
    if not checks["all_passed"]:
        raise SystemExit("expansion validation failed")


if __name__ == "__main__":
    main()
