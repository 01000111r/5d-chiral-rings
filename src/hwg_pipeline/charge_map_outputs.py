"""Deterministic, theory-independent physical charge-map reports."""
import json
from pathlib import Path
from sage.all import QQ, ZZ, matrix, vector
from .charge_maps import apply_charge_map, apply_charge_map_to_series, rational_json, solve_charge_map


def _write(path, value): path.write_text(json.dumps(value, indent=2, sort_keys=True)+"\n", encoding="utf-8")
def _q(v): return QQ(v["numerator"])/QQ(v["denominator"]) if isinstance(v,dict) else QQ(v)
def _entries(payload):
    field=next(k for k in ("coefficients_by_t_degree","generator_candidates_by_t_degree","relation_candidates_by_t_degree") if k in payload)
    return [x for degree in payload[field].values() for x in degree]
def _fmt(q):
    q=QQ(q); return str(q)
def _texq(q):
    q=QQ(q)
    return str(q) if q.denominator()==1 else ("-" if q<0 else "")+rf"\frac{{{abs(q.numerator())}}}{{{q.denominator()}}}"
def _formula(row,names):
    pieces=[]
    for c,n in zip(row,names):
        if not c: continue
        pieces.append(f"({_fmt(c)}) {n}")
    return " + ".join(pieces) or "0"
def _tex_series(payload,title):
    rows=[]
    for x in _entries(payload):
        m=_q(x.get("signed_multiplicity",x.get("coefficient",x.get("multiplicity"))))
        labels=",".join(map(str,x["child_dynkin_labels"])); B=_q(x["physical_charges"]["B"]); I=_q(x["physical_charges"]["I"])
        rows.append(rf"{_texq(m)}\,\chi_{{[{labels}]}}^{{(B={_texq(B)},I={_texq(I)})}}t^{{{x['t_degree']}}}")
    return "% "+title+"\n\\[\n"+" + ".join(rows)+"\n\\]\n"

def write_charge_map_outputs(spec,theory_id,branching_id,order,source_dir,output):
    solution=solve_charge_map(spec)
    if not solution.diagnostics.unique: raise ValueError("charge map is not unique")
    output.mkdir(parents=True,exist_ok=True)
    inputs={"pl":"branched_refined_plethystic_logarithm.json"}
    if (source_dir/"branched_character_series.json").exists():
        inputs["character"]="branched_character_series.json"
    optional={"generators":"branched_candidate_generators.json","relations":"branched_first_relation_candidates.json"}
    inputs.update({k:v for k,v in optional.items() if (source_dir/v).exists()})
    raw={k:json.loads((source_dir/v).read_text()) for k,v in inputs.items()}
    transformed={k:apply_charge_map_to_series(solution,v) for k,v in raw.items()}
    # An invertible map cannot merge distinct raw charge vectors.  Preserve
    # the compact, parent-level ``raw_provenance`` already carried by each
    # source term instead of duplicating the entire source record.
    for payload in transformed.values():
        for entry in _entries(payload): entry.pop("provenance", None)
    names={"character":"physical_branched_character_series","pl":"physical_branched_refined_plethystic_logarithm","generators":"physical_candidate_generators","relations":"physical_first_relation_candidates"}
    for key,payload in transformed.items():
        _write(output/(names[key]+".json"),payload)
        (output/(names[key]+".tex")).write_text(_tex_series(payload,names[key]),encoding="utf-8")
    validations=[]
    for a in spec.validation_anchors:
        actual=apply_charge_map(solution,a.raw); residual=tuple(x-y for x,y in zip(actual.values,a.physical.values))
        validations.append({"id":a.id,"expected":[rational_json(x) for x in a.physical.values],"actual":[rational_json(x) for x in actual.values],"residual":[rational_json(x) for x in residual],"passed":not any(residual)})
    sectors=sum((_entries(value) for value in transformed.values()),[])
    raw_sectors=sum((_entries(value) for value in raw.values()),[])
    raw_pl=_entries(raw["pl"])
    def anchor_present(anchor):
        return any(x["t_degree"]==anchor.t_degree and tuple(x["child_dynkin_labels"])==anchor.dynkin_labels
                   and all(_q(x["raw_charges"][name])==value for name,value in zip(anchor.raw.names,anchor.raw.values))
                   for x in raw_pl)
    bstep=_q((spec.charge_lattice or {}).get("B_step",1)); istep=_q((spec.charge_lattice or {}).get("I_step",1))
    lattice=all(_q(x["physical_charges"]["B"])/bstep in ZZ and _q(x["physical_charges"]["I"])/istep in ZZ for x in sectors)
    roundtrip=all(tuple(solution.inverse_matrix*vector(QQ,[_q(x["physical_charges"][n]) for n in spec.physical_charge_names]))==tuple(_q(x["raw_charges"][n]) for n in spec.raw_charge_names) for x in sectors)
    pl_entries=_entries(transformed["pl"])
    def pl_key(x): return (x["t_degree"],tuple(x["child_dynkin_labels"]),_q(x["physical_charges"]["B"]),_q(x["physical_charges"]["I"]))
    pl_values={pl_key(x):_q(x.get("coefficient",x.get("signed_multiplicity"))) for x in pl_entries}
    conjugation=all(pl_values.get((d,tuple(reversed(labels)),-B,-I))==m for (d,labels,B,I),m in pl_values.items())
    expected_by_degree={}
    for degree, entries in (spec.physical_pl_benchmarks or {}).items():
        expected_by_degree[int(degree)]={(int(degree),tuple(e["dynkin_labels"]),_q(e["B"]),_q(e["I"])):_q(e["coefficient"]) for e in entries}
    def totals(payload):
        out={}
        for x in _entries(payload):
            m=_q(x.get("signed_multiplicity",x.get("coefficient",x.get("multiplicity"))))
            out[x["t_degree"]]=out.get(x["t_degree"],QQ(0))+m*x["child_representation_dimension"]
        return out
    checks={"solution_unique":True,"matrix_rank_two":solution.matrix.rank()==2,
      "expected_matrix_recovered":spec.expected_matrix is None or solution.matrix==matrix(QQ,spec.expected_matrix),
      "expected_determinant":spec.expected_determinant is None or solution.matrix.det()==spec.expected_determinant,
      "all_defining_residuals_zero":all(not any(x) for x in solution.diagnostics.defining_residuals),"all_validation_anchors_pass":all(x["passed"] for x in validations),
      "all_defining_anchors_exist_in_raw_data":all(anchor_present(a) for a in spec.defining_anchors),
      "all_validation_anchors_exist_in_raw_data":all(anchor_present(a) for a in spec.validation_anchors),
      "configured_charge_lattice":lattice,"exact_raw_physical_raw_round_trip":roundtrip,
      "complete_pl_conjugation_through_order":conjugation,
      "raw_provenance_retained":all(x.get("raw_provenance") for x in sectors),
      "t_degrees_and_dimensions_preserved":all(totals(raw[k])==totals(transformed[k]) for k in raw),
      "complete_mandatory_inputs_translated":len(sectors)==len(raw_sectors)}
    for degree, expected in expected_by_degree.items():
        checks[f"physical_t{degree}_benchmark"]={k:v for k,v in pl_values.items() if k[0]==degree}==expected
    checks["all_passed"]=all(checks.values())
    M,inv=solution.matrix,solution.inverse_matrix
    payload={"charge_map_id":spec.id,"raw_charge_order":list(spec.raw_charge_names),"physical_charge_order":list(spec.physical_charge_names),
      "matrix":[[rational_json(x) for x in row] for row in M.rows()],"inverse_matrix":[[rational_json(x) for x in row] for row in inv.rows()],
      "determinant":rational_json(M.det()),"matrix_rank":M.rank(),"derived_formulas":[f"{n} = {_formula(M.row(i),spec.raw_charge_names)}" for i,n in enumerate(spec.physical_charge_names)],
      "inverse_formulas":[f"{n} = {_formula(inv.row(i),spec.physical_charge_names)}" for i,n in enumerate(spec.raw_charge_names)],"validation_anchors":validations}
    _write(output/"charge_map_solution.json",payload); _write(output/"charge_map_checks.json",{"theory_id":theory_id,"branching_id":branching_id,"validation_results":checks,"validation_anchors":validations})
    md=[f"# Exact charge map: `{spec.id}`","",f"Matrix: `{[[str(x) for x in r] for r in M.rows()]}`.",f"Inverse: `{[[str(x) for x in r] for r in inv.rows()]}`.",f"Rank `{M.rank()}`; determinant `{M.det()}`.",""]+[f"- **{'PASS' if v else 'FAIL'} — {k}**" for k,v in checks.items()]
    (output/"charge_map_solution.md").write_text("\n".join(md)+"\n")
    (output/"charge_map_solution.tex").write_text(r"\[\binom BI="+rf"\begin{{pmatrix}}{_texq(M[0,0])}&{_texq(M[0,1])}\\{_texq(M[1,0])}&{_texq(M[1,1])}\end{{pmatrix}}\binom xq.\]"+"\n")
    (output/"charge_map_checks.md").write_text("# Charge-map checks\n\n"+"\n".join(f"- **{'PASS' if v else 'FAIL'} — {k}**" for k,v in checks.items())+"\n")
    if not checks["all_passed"]: raise ValueError("charge-map checks failed")
    return checks
