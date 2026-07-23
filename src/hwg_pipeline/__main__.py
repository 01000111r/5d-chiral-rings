"""Command-line entry point for the exact HWG expansion milestone."""

import argparse
import json
from pathlib import Path

from sage.all import ZZ

from .expansion import expand_hwg, expand_pe, expand_rational_product
from .io import load_theory
from .render import render_monomial


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


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m hwg_pipeline")
    commands = parser.add_subparsers(dest="command", required=True)
    expand = commands.add_parser("expand", help="expand a structured HWG")
    expand.add_argument("theory_id")
    expand.add_argument("--order", required=True, type=int)
    args = parser.parse_args(argv)
    if args.order < 0:
        parser.error("--order must be nonnegative")
    root = _root()
    theory = load_theory(root / "theories" / f"{args.theory_id}.yaml")
    pe = expand_pe(theory, args.order)
    product = expand_rational_product(theory, args.order)
    checks = _validations(theory, args.order, pe, product)
    _write_outputs(theory, args.order, pe, checks,
                   root / "generated" / theory.id / f"order_{args.order}")
    if not checks["all_passed"]:
        raise SystemExit("expansion validation failed")


if __name__ == "__main__":
    main()
