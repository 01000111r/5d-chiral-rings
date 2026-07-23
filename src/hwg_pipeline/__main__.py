"""Command-line entry point for the exact HWG expansion milestone."""

import argparse
import json
from pathlib import Path

from sage.all import ZZ

from .expansion import expand_hwg, expand_pe, expand_rational_product
from .io import load_theory
from .render import render_monomial
from .characters import dimension_refine, restore_characters, unrefine
from .plethystic import (dimension_refine_virtual, plethystic_logarithm,
                         scalar_plethystic_exponential,
                         scalar_plethystic_logarithm, unrefine_virtual,
                         VirtualCharacterSeries, VirtualRepresentationContent,
                         plethystic_exponential)


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


def _rational(value):
    return int(value) if value.denominator() == 1 else {
        "numerator": int(value.numerator()), "denominator": int(value.denominator())}


def _read_rational(value):
    from sage.all import QQ
    return QQ(value["numerator"], value["denominator"]) if isinstance(value, dict) else QQ(value)


def _load_virtual_pl(theory, path, order):
    payload = json.loads(path.read_text())
    terms = []
    for degree, entries in payload["coefficients_by_t_degree"].items():
        for entry in entries:
            labels = tuple(tuple(rep["dynkin_labels"])
                           for rep in entry["irreducible_representations"])
            content = VirtualRepresentationContent.single_irrep(
                theory.simple_factors, labels, _read_rational(entry["coefficient"]))
            terms.append(((ZZ(degree), tuple(entry["abelian_charges"].items())), content))
    return VirtualCharacterSeries(theory, terms, order)


def _difference(expected, actual):
    left = {(int(d), tuple((k, str(v)) for k, v in q), tuple(tuple(map(int, x)) for x in labels)): c
            for (d, q), content in expected for labels, c in content}
    right = {(int(d), tuple((k, str(v)) for k, v in q), tuple(tuple(map(int, x)) for x in labels)): c
             for (d, q), content in actual for labels, c in content}
    return [{"t_degree": key[0], "abelian_charges": dict(key[1]),
             "dynkin_labels": [list(x) for x in key[2]],
             "expected_multiplicity": _rational(left.get(key, 0)),
             "reconstructed_multiplicity": _rational(right.get(key, 0))}
            for key in sorted(set(left) | set(right)) if left.get(key, 0) != right.get(key, 0)]


def _reconstruction_outputs(theory, order, hwg, output):
    """Load the persisted PL and perform all three exact reconstruction checks."""
    expected = VirtualCharacterSeries.from_character_series(
        restore_characters(theory, hwg), order)
    pl = _load_virtual_pl(theory, output / "refined_plethystic_logarithm.json", order)
    reconstructed = plethystic_exponential(pl, order)
    differences = _difference(expected, reconstructed)
    refined_dim = dimension_refine_virtual(reconstructed)
    expected_dim = dimension_refine_virtual(expected)
    unrefined = unrefine_virtual(reconstructed)
    expected_unrefined = unrefine_virtual(expected)
    scalar_payload = json.loads((output / "unrefined_plethystic_logarithm.json").read_text())
    scalar_pl = tuple((ZZ(d), _read_rational(c)) for d, c in
                      scalar_payload["coefficients_by_t_degree"].items())
    scalar = scalar_plethystic_exponential(scalar_pl, order)

    groups, texgroups = {}, {}
    for (degree, charges), content in reconstructed:
        for labels, coefficient in content:
            groups.setdefault(str(degree), []).append({
                "abelian_charges": _charge_dict(charges),
                "dynkin_labels": [{"cartan_factor_id": factor.id,
                                    "labels": list(map(int, label))}
                                   for factor, label in zip(theory.simple_factors, labels)],
                "multiplicity": _rational(coefficient)})
            rep = " ".join("[" + ",".join(map(str, label)) + "]" for label in labels)
            texgroups.setdefault(int(degree), []).append(
                f"{coefficient}{''.join(_latex_power(k,v) for k,v in charges if v)}{rep}")
    base = {"theory_id": theory.id, "maximum_t_degree": int(order)}
    (output / "reconstructed_character_series.json").write_text(json.dumps(
        {**base, "coefficients_by_t_degree": groups}, indent=2, sort_keys=True) + "\n")
    pieces = ["(" + " + ".join(texgroups[d]) + rf")t^{{{d}}}" for d in sorted(texgroups)]
    (output / "reconstructed_character_series.tex").write_text(
        "% Reconstructed exact character series.\n\\begin{align*}\nH_{rec}=" + " + ".join(pieces) + "\n\\end{align*}\n")

    qgroups = {}
    for (d, q), c in refined_dim:
        qgroups.setdefault(str(d), []).append({"abelian_charges": _charge_dict(q), "coefficient": _rational(c)})
    (output / "reconstructed_q_refined_dimension_series.json").write_text(json.dumps(
        {**base, "coefficients_by_t_degree": qgroups}, indent=2, sort_keys=True) + "\n")
    qtex = " + ".join(f"{c}{''.join(_latex_power(k,v) for k,v in q if v)}t^{{{d}}}" for (d,q),c in refined_dim)
    (output / "reconstructed_q_refined_dimension_series.tex").write_text(
        "% Reconstructed q-refined dimensions.\n\\begin{align*}\nH_{dim,rec}=" + qtex + "\n\\end{align*}\n")
    (output / "reconstructed_unrefined_hilbert_series.json").write_text(json.dumps(
        {**base, "coefficients_by_t_degree": {str(d): _rational(c) for d,c in unrefined}}, indent=2, sort_keys=True) + "\n")
    (output / "reconstructed_unrefined_hilbert_series.tex").write_text(
        "% Reconstructed unrefined series.\n\\begin{align*}\nH_{rec}=" +
        " + ".join(f"{c}t^{{{d}}}" for d,c in unrefined) + "\n\\end{align*}\n")
    (output / "reconstruction_difference.json").write_text(json.dumps(
        {**base, "mismatch_count": len(differences), "mismatches": differences}, indent=2, sort_keys=True) + "\n")
    stable = plethystic_exponential(VirtualCharacterSeries(pl.theory, pl.terms, order + 1), order + 1)
    checks = {"refined_character_series_equal": not differences,
              "q_refined_dimensions_equal": refined_dim == expected_dim,
              "unrefined_dimensions_equal": unrefined == expected_unrefined,
              "independent_scalar_pe_equal": scalar == expected_unrefined,
              "multiplicities_integral": all(c in ZZ for _,x in reconstructed for _,c in x),
              "multiplicities_nonnegative": all(c >= 0 for _,x in reconstructed for _,c in x),
              "difference_is_empty": not differences,
              "order_stability": VirtualCharacterSeries(stable.theory, stable.terms, order) == reconstructed}
    checks["all_passed"] = all(checks.values())
    (output / "reconstruction_checks.json").write_text(json.dumps(
        {**base, "validation_results": checks}, indent=2, sort_keys=True) + "\n")
    (output / "reconstruction_checks.md").write_text("# Reconstruction checks\n\n" +
        "\n".join(f"- **{'PASS' if v else 'FAIL'} — {k}**" for k,v in checks.items()) + "\n")
    return checks


def _pl_outputs(theory, order, hwg, output):
    series = restore_characters(theory, hwg)
    pl = plethystic_logarithm(series, order)
    dim = dimension_refine_virtual(pl); plain = unrefine_virtual(pl)
    groups, texgroups = {}, {}
    for (degree, charges), content in pl:
        for labels, coefficient in content:
            groups.setdefault(str(degree), []).append({
                "abelian_charges": _charge_dict(charges),
                "irreducible_representations": [{"cartan_factor_id": factor.id,
                    "dynkin_labels": [int(x) for x in label]}
                    for factor,label in zip(theory.simple_factors,labels)],
                "coefficient": _rational(coefficient)})
            reps=" ".join("["+",".join(map(str,label))+rf"]_{{{factor.cartan_type}_{factor.rank}}}"
                          for factor,label in zip(theory.simple_factors,labels))
            texgroups.setdefault(int(degree),[]).append((coefficient,charges,reps))
    payload={"theory_id":theory.id,"maximum_t_degree":int(order),"coefficients_by_t_degree":groups}
    (output/"refined_plethystic_logarithm.json").write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    pieces=[]
    for degree in sorted(texgroups):
        body=[]
        for c,charges,reps in texgroups[degree]:
            sign="-" if c<0 else ("+" if body else "")
            magnitude=abs(c); scalar="" if magnitude==1 else str(magnitude)
            body.append(sign+scalar+"".join(_latex_power(k,v) for k,v in charges if v)+reps)
        pieces.append("("+" ".join(body)+rf")t^{{{degree}}}")
    (output/"refined_plethystic_logarithm.tex").write_text("% Exact refined plethystic logarithm.\n\\begin{align*}\nPL[H] = "+" + ".join(pieces)+"\n\\end{align*}\n")
    dgroups={}
    for (degree,charges),coefficient in dim:
        dgroups.setdefault(str(degree),[]).append({"abelian_charges":_charge_dict(charges),"coefficient":_rational(coefficient)})
    (output/"q_refined_dimension_pl.json").write_text(json.dumps({"theory_id":theory.id,"maximum_t_degree":int(order),"coefficients_by_t_degree":dgroups},indent=2,sort_keys=True)+"\n")
    dtex=" + ".join(f"{c}"+"".join(_latex_power(k,v) for k,v in q if v)+rf"t^{{{d}}}" for (d,q),c in dim)
    (output/"q_refined_dimension_pl.tex").write_text("% Exact dimension-refined plethystic logarithm.\n\\begin{align*}\nPL_{\\dim}[H] = "+dtex+"\n\\end{align*}\n")
    upayload={"theory_id":theory.id,"maximum_t_degree":int(order),"coefficients_by_t_degree":{str(d):_rational(c) for d,c in plain}}
    (output/"unrefined_plethystic_logarithm.json").write_text(json.dumps(upayload,indent=2,sort_keys=True)+"\n")
    (output/"unrefined_plethystic_logarithm.tex").write_text("% Exact unrefined plethystic logarithm.\n\\begin{align*}\nPL[H] = "+" + ".join(f"{c}t^{{{d}}}" for d,c in plain)+"\n\\end{align*}\n")
    direct=scalar_plethystic_logarithm(unrefine(series),order)
    expected2={((0,0,0,0,0),):1,((1,0,0,0,1),):1}; actual2={k:int(c) for (d,q),x in pl if d==2 for k,c in x}
    expected3={(1,((0,1,0,0,0),)):1,(-1,((0,0,0,1,0),)):1}
    actual3={(int(dict(q)["q"]),k):int(c) for (d,q),x in pl if d==3 for k,c in x}
    checks={"all_final_coefficients_integral":all(c.denominator()==1 for _,x in pl for _,c in x),
            "negative_coefficients_retained":any(c<0 for _,x in pl for _,c in x),
            "degrees_truncated":all(d<=order for (d,_),_ in pl),
            "degree_2_independent_value":actual2==expected2,"degree_3_independent_value":actual3==expected3,
            "degree_4_independent_value":{k:int(c) for (d,q),x in pl if d==4 for k,c in x}=={((0,0,0,0,0),):-1,((1,0,0,0,1),):-1},
            "direct_scalar_matches_refined_unrefinement":direct==plain}
    checks["all_passed"]=all(checks.values())
    cp={"theory_id":theory.id,"maximum_t_degree":int(order),"direct_scalar_plethystic_logarithm":{str(d):_rational(c) for d,c in direct},"validation_results":checks}
    (output/"plethystic_logarithm_checks.json").write_text(json.dumps(cp,indent=2,sort_keys=True)+"\n")
    (output/"plethystic_logarithm_checks.md").write_text("# Plethystic-logarithm checks\n\n"+"\n".join(f"- **{'PASS' if v else 'FAIL'} — {k}**" for k,v in checks.items())+"\n")
    return checks


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m hwg_pipeline")
    commands = parser.add_subparsers(dest="command", required=True)
    for name, help_text in (("expand", "expand a structured HWG"),
                            ("characters", "restore irreducible characters and dimensions"),
                            ("plethystic-log", "compute the exact refined plethystic logarithm"),
                            ("reconstruct", "reconstruct a character series from its refined PL")):
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
    elif args.command == "characters":
        output.mkdir(parents=True, exist_ok=True)
        checks = _character_outputs(theory, args.order, pe, output)
    elif args.command == "plethystic-log":
        output.mkdir(parents=True, exist_ok=True)
        checks = _pl_outputs(theory, args.order, pe, output)
    else:
        output.mkdir(parents=True, exist_ok=True)
        checks = _reconstruction_outputs(theory, args.order, pe, output)
    if not checks["all_passed"]:
        raise SystemExit("expansion validation failed")


if __name__ == "__main__":
    main()
