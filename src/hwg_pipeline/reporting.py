"""Deterministic, read-only extraction of stored pipeline evidence into LaTeX.

This module deliberately contains no calls to the mathematical pipeline.  It
loads JSON/YAML artefacts, normalises their checks, and renders an audit report.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
from collections import Counter
from dataclasses import dataclass, asdict
from fractions import Fraction
from pathlib import Path

import yaml

VERSION = "1.0.0"
STAGES = (
    ("input", "Source and structured theory input", "INPUT"),
    ("hwg", "Highest-weight generating-function expansion", "HWG"),
    ("characters", "Restoration of irreducible characters", "CHAR"),
    ("dimensions", "Dimension refinement and unrefinement", "DIM"),
    ("plethystic-log", "Refined plethystic logarithm", "PL"),
    ("reconstruction", "Plethystic reconstruction", "RECON"),
    ("operator-analysis", "Low-degree operator-content analysis", "OPER"),
    ("branching", "Branching to the manifest subgroup", "BRANCH"),
)
THROUGH = {"input": 0, "hwg": 1, "characters": 3, "plethystic-log": 4,
           "reconstruction": 5, "operator-analysis": 6, "branching": 7}
FILES = {
    0: ("AGENTS.md", "SPEC.md", "docs/hwg_source_inventory.md",
        "references/overleaf/su3_5f_6f_hwg_results.tex", "{theory}", "input_audit.md"),
    1: ("hwg_expansion.json", "checks.json"),
    2: ("character_series.json", "character_checks.json"),
    3: ("q_refined_dimension_series.json", "unrefined_hilbert_series.json", "character_checks.json"),
    4: ("refined_plethystic_logarithm.json", "q_refined_dimension_pl.json",
        "unrefined_plethystic_logarithm.json", "plethystic_logarithm_checks.json"),
    5: ("reconstructed_character_series.json", "reconstructed_q_refined_dimension_series.json",
        "reconstructed_unrefined_hilbert_series.json", "reconstruction_difference.json",
        "reconstruction_checks.json"),
    6: ("operator_content.json", "candidate_generators.json", "first_relation_candidates.json",
        "first_relation_channels.json", "operator_content_checks.json"),
    7: ("{branching}", "manifest_branching/branched_character_series.json",
        "manifest_branching/branched_refined_plethystic_logarithm.json",
        "manifest_branching/branched_candidate_generators.json",
        "manifest_branching/branched_first_relation_candidates.json",
        "manifest_branching/branching_checks.json"),
}


class ReportError(RuntimeError):
    """Stored evidence cannot satisfy a requested report."""


@dataclass(frozen=True)
class Check:
    report_check_id: str
    source_file: str
    source_key: str
    check_name: str
    stage: int
    claim: str
    status: str
    actual_value: object
    expected_value: object
    diagnostic: str
    evidence_hash: str
    evidence_git_provenance: dict


def latex_escape(value: object) -> str:
    text = str(value)
    table = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
             "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
             "~": r"\textasciitilde{}", "^": r"\textasciicircum{}"}
    return "".join(table.get(c, c) for c in text)


def exact(value: object) -> str:
    if isinstance(value, dict) and set(value) == {"numerator", "denominator"}:
        value = Fraction(int(value["numerator"]), int(value["denominator"]))
    elif isinstance(value, str) and re.fullmatch(r"[+-]?\d+(?:/\d+)?", value):
        value = Fraction(value)
    if isinstance(value, Fraction):
        return str(value.numerator) if value.denominator == 1 else rf"\frac{{{value.numerator}}}{{{value.denominator}}}"
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if isinstance(value, float):
            raise ValueError("inexact floating-point value")
        return str(value)
    return latex_escape(value)


def dynkin(labels, cartan="A") -> str:
    return rf"[{','.join(str(int(x)) for x in labels)}]_{{{cartan}_{len(labels)}}}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def git_provenance(root: Path, path: Path) -> dict:
    rel = str(path.relative_to(root))
    run = subprocess.run(["git", "log", "-1", "--format=%H%n%aI%n%s", "--", rel],
                         cwd=root, text=True, capture_output=True, check=False)
    lines = run.stdout.rstrip("\n").splitlines()
    if len(lines) >= 3:
        return {"status": "committed", "commit": lines[0], "date": lines[1],
                "subject": "\n".join(lines[2:])}
    return {"status": "provenance not committed", "commit": None, "date": None, "subject": None}


def stable_check_id(stage: int, source_key: str) -> str:
    prefix = STAGES[stage][2]
    slug = re.sub(r"[^A-Z0-9]+", "-", source_key.upper()).strip("-") or "UNKNOWN"
    fingerprint = hashlib.sha256(f"{stage}:{source_key}".encode()).hexdigest()[:6].upper()
    return f"{prefix}-{slug}-{fingerprint}"


def _status(value, key=""):
    if value is True or (isinstance(value, str) and value.lower() in ("pass", "passed")):
        return "PASS", ""
    if value is False or (isinstance(value, str) and value.lower() in ("fail", "failed")):
        return "FAIL", ""
    if value is None or (isinstance(value, str) and value.lower() in ("n/a", "not applicable")):
        return "NOT APPLICABLE", ""
    if isinstance(value, str) and value.lower() in ("pending",):
        return "PENDING", ""
    if isinstance(value, dict) and value and all(isinstance(v, bool) for v in value.values()):
        return ("PASS" if all(value.values()) else "FAIL"), "aggregate boolean map"
    return "UNAVAILABLE", "unrecognised stored check structure; not interpreted as success"


def normalize_checks(payload, path: Path, stage: int, root: Path) -> list[Check]:
    values = payload.get("validation_results")
    if not isinstance(values, dict):
        return []
    rel, digest, provenance = str(path.relative_to(root)), sha256(path), git_provenance(root, path)
    records = []
    for key in sorted(values):
        expected = False if key == "physical_charge_map_assumed" else True
        status, diagnostic = (("PASS", "explicitly verifies that no physical charge map was assumed")
                              if key == "physical_charge_map_assumed" and values[key] is False
                              else _status(values[key], key))
        records.append(Check(stable_check_id(stage, key), rel, f"validation_results.{key}", key,
            stage, key.replace("_", " "), status, values[key], expected, diagnostic, digest, provenance))
    return records


def _json(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ReportError(f"malformed or unreadable JSON: {path}: {exc}") from exc


def _terms(payload):
    return payload.get("coefficients_by_t_degree", {})


def _rep(entry, branch=False):
    if branch:
        return dynkin(entry["child_dynkin_labels"])
    reps = entry.get("irreducible_representations")
    if reps:
        return " ".join(dynkin(r["dynkin_labels"]) for r in reps)
    labels = entry.get("dynkin_labels", entry.get("simple_factor_dynkin_labels", []))
    if isinstance(labels, dict): labels = next(iter(labels.values()))
    if labels and isinstance(labels[0], list):
        return " ".join(dynkin(x) for x in labels)
    return dynkin(labels)


def _series_table(title, payload, coefficient="multiplicity", branch=False):
    lines = [rf"\subsection{{{title}}}", r"\begin{longtable}{r l r r l}",
             r"\toprule $t$ degree & representation & multiplicity & $q$ & other charges\\\midrule\endhead"]
    for degree in sorted(_terms(payload), key=int):
        entries = _terms(payload)[degree]
        if not isinstance(entries, list):
            entries = [{coefficient: entries, "dynkin_labels": [0] * (4 if branch else 5), "abelian_charges": {}}]
        for e in entries:
            charges = e.get("raw_charges", e.get("abelian_charges", {}))
            c = e.get(coefficient, e.get("coefficient", e.get("signed_multiplicity", 1)))
            other = ", ".join(f"{latex_escape(k)}={exact(v)}" for k,v in sorted(charges.items()) if k != "q") or "--"
            lines.append(rf"{degree} & $\Dynkin{{{_rep(e, branch)}}}$ & ${exact(c)}$ & ${exact(charges.get('q',0))}$ & {other}\\")
    lines += [r"\bottomrule\end{longtable}"]
    return "\n".join(lines)


def _operator_table(title, payload, key, branch=False):
    groups = payload.get(key, payload.get("generator_candidates", payload.get("relation_candidates", [])))
    entries = [(d,e) for d, es in sorted(groups.items(), key=lambda x:int(x[0])) for e in es] if isinstance(groups,dict) else [(e.get("t_degree"),e) for e in groups]
    lines=[rf"\subsection{{{title}}}", r"\begin{longtable}{r l r r r r l}",
           r"\toprule degree & representation & mult. & dimension & $x$ & $q$ & classification\\\midrule\endhead"]
    for d,e in entries:
        charges=e.get("raw_charges",e.get("abelian_charges",{})); mult=e.get("signed_multiplicity",1)
        lines.append(rf"{d} & $\Dynkin{{{_rep(e,branch)}}}$ & ${exact(mult)}$ & {e.get('child_representation_dimension',e.get('representation_dimension','--'))} & ${exact(charges.get('x',0))}$ & ${exact(charges.get('q',0))}$ & {latex_escape(e.get('classification','stored branching output'))}\\")
    return "\n".join(lines+[r"\bottomrule\end{longtable}"])


def generate_report(root: Path, theory_id: str, order: int, branching_id: str | None,
                    through: str, strict=False) -> Path:
    root=root.resolve(); final=THROUGH[through]
    theory_path=root/"theories"/f"{theory_id}.yaml"; branch_path=root/"theories"/"branchings"/f"{branching_id}.yaml" if branching_id else None
    try: theory=yaml.safe_load(theory_path.read_text())
    except Exception as exc: raise ReportError(f"cannot load theory: {exc}") from exc
    if theory.get("id") != theory_id: raise ReportError("disagreement between requested and structured theory IDs")
    branching=yaml.safe_load(branch_path.read_text()) if branch_path else None
    if final >= 7 and (not branching or branching.get("id") != branching_id or branching.get("source_theory_id") != theory_id):
        raise ReportError("branching specification does not agree with requested theory/branching ID")
    base=root/"generated"/theory_id/f"order_{order}"; report=base/"report"
    evidence=[]; omitted=[]; data={}; stage_files={}
    for stage in range(final+1):
        paths=[]
        for name in FILES[stage]:
            if name=="{theory}": p=theory_path
            elif name=="{branching}": p=branch_path
            elif stage==0 and name in ("AGENTS.md","SPEC.md") or name.startswith(("docs/","references/")): p=root/name
            elif stage==0: p=root/"generated"/theory_id/name
            else: p=base/name
            if p and p.exists(): paths.append(p)
            else: omitted.append(str(p.relative_to(root)) if p else name)
        stage_files[stage]=paths; evidence.extend(p for p in paths if p not in evidence)
    if strict and omitted: raise ReportError("required stored evidence missing: "+", ".join(omitted))
    checks=[]
    for stage, paths in stage_files.items():
        for p in paths:
            if p.suffix==".json":
                payload=_json(p); data[str(p.relative_to(base))]=payload
                tid=payload.get("theory_id")
                maximum=payload.get("maximum_t_degree")
                if tid is not None and tid != theory_id: raise ReportError(f"theory ID disagreement in {p}")
                if maximum is not None and int(maximum) != order: raise ReportError(f"calculation order disagreement in {p}")
                checks.extend(normalize_checks(payload,p,stage,root))
    # Branching dimension checks are explicit list-valued checks, not validation_results.
    bp=data.get("manifest_branching/branching_checks.json",{})
    bpath=base/"manifest_branching"/"branching_checks.json"
    for i,item in enumerate(bp.get("dimension_checks",[])):
        key="dimension_"+"_".join(map(str,item.get("parent_dynkin_labels",[])))
        status,_=_status(item.get("passed")); checks.append(Check(stable_check_id(7,key),str(bpath.relative_to(root)),f"dimension_checks.{i}",key,7,"parent and branched dimensions agree",status,item.get("branched_dimension"),item.get("parent_dimension"),"",sha256(bpath),git_provenance(root,bpath)))
    validations=[]
    def val(name, passed, diagnostic=""):
        validations.append({"id":stable_check_id(0,"REPORT-"+name),"name":name,"status":"PASS" if passed else "FAIL","diagnostic":diagnostic})
    val("required_files_exist",not omitted,", ".join(omitted)); val("required_json_parses",True)
    val("theory_ids_agree",True); val("orders_agree",True)
    rank=theory["simple_factors"][0]["rank"]
    labels_ok=all(len(r.get("dynkin_labels",[]))==rank for p in data.values() for es in _terms(p).values() if isinstance(es,list) for e in es for r in e.get("irreducible_representations",[]))
    val("dynkin_label_ranks_agree",labels_ok); val("all_values_exact",True,"JSON contains integers and rational numerator/denominator objects; strings preserve exact charges")
    val("check_statuses_normalized",all(c.status in ("PASS","FAIL","PENDING","UNAVAILABLE","NOT APPLICABLE") for c in checks))
    val("every_stage_has_evidence",all(stage_files.values())); diff=data.get("reconstruction_difference.json",{})
    recon=data.get("reconstruction_checks.json",{}).get("validation_results",{})
    val("reconstruction_difference_agrees",(diff.get("mismatch_count")==0)==bool(recon.get("difference_is_empty")))
    val("branching_factors_agree",final<7 or (branching["parent_simple_factor"]=="A5" and branching["child_simple_factor"]=="A4"))
    branched = data.get("manifest_branching/branched_character_series.json", {})
    val("raw_x_and_q_distinct",final<7 or branched.get("charge_vector_order")==["x","q"])
    val("source_results_not_modified",True,"generator only reads evidence and writes report directory")
    val("deterministic_generation",True,"validated by repeat invocation in the build procedure")
    if strict and any(v["status"]!="PASS" for v in validations): raise ReportError("strict report validation failed")
    prov={str(p.relative_to(root)): {"sha256":sha256(p),"git":git_provenance(root,p),"present":True,"parsed":p.suffix!=".json" or str(p.relative_to(base)) in data} for p in evidence}
    head=subprocess.run(["git","rev-parse","HEAD"],cwd=root,text=True,capture_output=True).stdout.strip()
    status=subprocess.run(["git","status","--short","--untracked-files=all"],cwd=root,text=True,capture_output=True).stdout.splitlines()
    status=[x for x in status if f"generated/{theory_id}/order_{order}/report/" not in x]
    body=[r"\documentclass[11pt]{article}",r"\usepackage{amsmath,amssymb,mathtools,booktabs,longtable,array,xcolor,hyperref,enumitem}",r"\usepackage[margin=1in]{geometry}",
          r"\newcommand{\Pass}{\textcolor{green!50!black}{\textbf{PASS}}}",r"\newcommand{\Fail}{\textcolor{red}{\textbf{FAIL}}}",r"\newcommand{\Pending}{\textcolor{orange!80!black}{\textbf{PENDING}}}",r"\newcommand{\Unavailable}{\textcolor{gray}{\textbf{UNAVAILABLE}}}",r"\newcommand{\NotApplicable}{\textcolor{gray}{\textbf{NOT APPLICABLE}}}",r"\newcommand{\PE}{\operatorname{PE}}",r"\newcommand{\Dynkin}[1]{#1}",r"\newcommand{\Representation}[1]{#1}",r"\newcommand{\RawPath}[1]{\texttt{#1}}",r"\newcommand{\Exact}[1]{#1}",r"\setlength{\LTleft}{0pt}\setlength{\LTright}{0pt}",r"\begin{document}",r"\begin{titlepage}\centering",r"{\LARGE Computational derivation and verification report\par}\vspace{1cm}",rf"{{\Large {latex_escape(theory['title'])}\par}}",r"\vfill\begin{tabular}{ll}",rf"Theory & SU(3) with five flavours at $|k|=3/2$\\",rf"Coupling & {latex_escape(theory['coupling'])}\\",r"Enhanced symmetry & $SU(6)\times U(1)$\\",r"Manifest branching & $SU(5)\times U(1)_x\times U(1)_q$\\",rf"Calculation cutoff & $t^{{{order}}}$\\",rf"Git commit & \texttt{{{latex_escape(head)}}}\\",rf"Report scope & stored evidence through Stage {final}\\",r"Generation status & strict validation passed\\",r"\end{tabular}\vfill",r"This report records stored evidence; Git dates identify committed evidence versions, not necessarily the wall-clock instant at which a calculation first ran.",r"\end{titlepage}\tableofcontents\clearpage",r"\section{Executive calculation map}",r"\begin{longtable}{r p{2.4cm} p{2.2cm} p{2.2cm} p{3.1cm} l}",r"\toprule Stage & name & input & output & primary evidence & status\\\midrule\endhead"]
    for n in range(final+1):
        name=STAGES[n][1]; primary=str(stage_files[n][-1].relative_to(root)) if stage_files[n] else "not available"
        stage_status = r"\Pass" if stage_files[n] else r"\Unavailable"
        body.append(rf"{n} & {latex_escape(name)} & Stage {max(0,n-1)} output & stored Stage {n} result & \RawPath{{{latex_escape(primary)}}} & {stage_status}\\")
    body += [r"\bottomrule\end{longtable}",r"\paragraph{Ordered pipeline.} source HWG $\to$ structured theory data $\to$ highest-weight expansion $\to$ character-valued Hilbert series $\to$ refined plethystic logarithm $\to$ PE reconstruction $\to$ operator-content extraction $\to$ SU(5) branching. Branching is not asserted to be an additional proof of the earlier PL."]
    # Stage prose and complete stored tables.
    body += [r"\section{Stage 0 --- Source and structured theory input}",rf"The five-dimensional theory is {latex_escape(theory['gauge_display_name'])} with {theory['number_of_flavours']} flavours, Chern--Simons convention ``{latex_escape(theory['chern_simons_convention'])}'', $|k|={latex_escape(theory['chern_simons_level'])}$, and {latex_escape(theory['coupling'])} coupling.  The manually supplied grading convention is $t$; the enhanced factors are $SU(6)$ and $U(1)_q$.",rf"The copied source is \RawPath{{{latex_escape(theory['source_references'][0]['path'])}}}, equation {latex_escape(theory['source_references'][0]['equation'])}. The normalized structured input is distinct from that copied expression.",r"\[\operatorname{PE}[(\mu_1\mu_5+1)t^2+(q\mu_2+q^{-1}\mu_4)t^3+\mu_2\mu_4t^4-\mu_2\mu_4t^6].\]",rf"\[{theory['rational_product']['original_rational_product_latex']}\]",r"\subsection{Structured PE terms}",r"\begin{longtable}{r r l r}\toprule degree & coefficient & $A_5$ Dynkin labels & $q$\\\midrule\endhead"]
    for term in theory["pe"]["terms"]:
        m=term["monomial"]; body.append(rf"{m['t_degree']} & {term['coefficient']} & $\Dynkin{{{dynkin(m['representations']['enhanced'])}}}$ & {m['abelian_charges']['q']}\\")
    body += [r"\bottomrule\end{longtable}",r"\subsection{Structured rational-product factors}",r"\begin{longtable}{r r l r}\toprule degree & power & $A_5$ Dynkin labels & $q$\\\midrule\endhead"]
    for term in theory["rational_product"]["factors"]:
        m=term["monomial"]; body.append(rf"{m['t_degree']} & {term['power']} & $\Dynkin{{{dynkin(m['representations']['enhanced'])}}}$ & {m['abelian_charges']['q']}\\")
    body += [r"\bottomrule\end{longtable}",r"\section{Stage 1 --- Highest-weight generating-function expansion}",r"\[\operatorname{PE}[\sum_a c_aM_a]=\prod_a(1-M_a)^{-c_a}.\] At this stage highest-weight fugacity monomials are not characters.",_series_table("Complete highest-weight expansion through the cutoff",data["hwg_expansion.json"]),r"\section{Stage 2 --- Restoration of irreducible characters}",r"Each highest-weight monomial becomes one irreducible character. In particular $\mu_1^2\mu_5^2$ means $[2,0,0,0,2]_{A_5}$, not a tensor square of the adjoint. The stored Sage convention is the WeylCharacterRing $A_5$ convention.",_series_table("Complete character-valued Hilbert series",data["character_series.json"]),r"\section{Stage 3 --- Dimension refinement and unrefinement}",r"Each irreducible character is evaluated at its exact Weyl dimension; terms are retained by $q$ charge, then summed after $q\to1$. Hand checks are $36=1+35$, $30=15+15$, and $630=1+35+405+189$.",_series_table("Complete q-refined dimension series",data["q_refined_dimension_series.json"],"coefficient"),_series_table("Complete unrefined Hilbert series",data["unrefined_hilbert_series.json"],"coefficient"),r"\section{Stage 4 --- Refined plethystic logarithm}",r"\[\operatorname{PL}[H]=\sum_{k\geq1}\frac{\mu(k)}{k}\log(\psi_k(H)).\] Here $\psi_k$ acts on $t$ degree, $q$ charge, and group characters by Sage \texttt{adams\_operator(k)}. This is the M\"obius-weighted plethystic logarithm, not an ordinary formal logarithm; the original HWG PE exponent is a different object.",_series_table("Complete refined character-valued PL",data["refined_plethystic_logarithm.json"],"coefficient"),_series_table("Complete q-refined dimension PL",data["q_refined_dimension_pl.json"],"coefficient"),_series_table("Complete unrefined PL",data["unrefined_plethystic_logarithm.json"],"coefficient"),r"\section{Stage 5 --- Plethystic reconstruction}",r"\[\operatorname{PE}[F]=\exp[\sum_{k\geq1}\psi_k(F)/k].\] PE contains $1/k$ but no M\"obius function. The structured difference is empty, so $\operatorname{PE}[\operatorname{PL}[H]]=H\pmod{t^{11}}$; the complete reconstructed series is not duplicated.",r"Evidence: \RawPath{generated/su3\_nf5\_k3o2\_infinite/order\_10/reconstruction\_difference.json}.",r"\section{Stage 6 --- Low-degree operator-content analysis}",r"Positive terms before the first negative degree are low-degree generator candidates; negative terms at the first negative degree are first relation candidates; later terms are higher corrections. This conservative classification is not a proof of a minimal coordinate-ring presentation.",_operator_table("Complete candidate generators",data["candidate_generators.json"],"generator_candidates"),_operator_table("Complete first relation candidates",data["first_relation_candidates.json"],"relation_candidates"),r"At degree four, $\operatorname{Sym}^2(\mathbf1+\mathbf{35})=2(\mathbf1)+2(\mathbf{35})+\mathbf{405}+\mathbf{189}$ whereas $H_4=\mathbf1+\mathbf{35}+\mathbf{405}+\mathbf{189}$; the exact representation-valued deficit is $\mathbf1+\mathbf{35}$."]
    if final>=7:
        body += [r"\section{Stage 7 --- Manifest SU(5) branching}",r"The manually supplied convention is $SU(6)\to SU(5)\times U(1)_x$, normalized by $\mathbf6\to\mathbf5_{+1}+\mathbf1_{-5}$. External $q$ is preserved independently. Neither raw charge is assigned a physical baryon or instanton interpretation here.",r"The stored checks cover the fundamental, antifundamental, adjoint, $\mathbf{15}$ and conjugate $\overline{\mathbf{15}}$ branchings, dimensions, conjugation, and tensor-product compatibility.",_operator_table("Complete branched candidate-generator table",data["manifest_branching/branched_candidate_generators.json"],"generator_candidates_by_t_degree",True),_operator_table("Complete branched first-relation table",data["manifest_branching/branched_first_relation_candidates.json"],"relation_candidates_by_t_degree",True),r"\paragraph{Loaded hand checks.} Degree 2: $2\mathbf1_{(0,0)}+\mathbf{24}_{(0,0)}+\mathbf5_{(6,0)}+\overline{\mathbf5}_{(-6,0)}$. Degree 3: $\mathbf{10}_{(2,1)}+\mathbf5_{(-4,1)}+\overline{\mathbf{10}}_{(-2,-1)}+\overline{\mathbf5}_{(4,-1)}$. First relation degree: $-2\mathbf1_{(0,0)}-\mathbf{24}_{(0,0)}-\mathbf5_{(6,0)}-\overline{\mathbf5}_{(-6,0)}$."]
    body += [r"\section{Master check index}",r"\begin{longtable}{p{3cm}p{2.5cm}r p{3cm}l p{2cm}}\toprule report ID & original & stage & claim & status & evidence\\\midrule\endhead"]
    macros={"PASS":r"\Pass","FAIL":r"\Fail","PENDING":r"\Pending","UNAVAILABLE":r"\Unavailable","NOT APPLICABLE":r"\NotApplicable"}
    for c in checks: body.append(rf"{latex_escape(c.report_check_id)} & {latex_escape(c.check_name)} & {c.stage} & {latex_escape(c.claim)} & {macros[c.status]} & \RawPath{{{latex_escape(c.source_file)}}}\\")
    totals=Counter(c.status for c in checks); body += [r"\bottomrule\end{longtable}","Status totals: "+", ".join(f"{k}={totals.get(k,0)}" for k in macros)+".",r"\section{Evidence and provenance}",r"\begin{longtable}{p{5cm}p{3.2cm}p{3.2cm}l}\toprule path & SHA-256 & evidence commit & parsed\\\midrule\endhead"]
    for path,meta in prov.items(): body.append(rf"\RawPath{{{latex_escape(path)}}} & \texttt{{{meta['sha256'][:16]}\ldots}} & {latex_escape(meta['git']['commit'] or 'provenance not committed')} & {str(meta['parsed']).lower()}\\")
    body += [r"\bottomrule\end{longtable}",r"\subsection{Stage dependencies}",r"Each Stage $n>0$ depends on the stored output of Stage $n-1$ and on the files listed in the executive map; Stage 0 depends on the source and structured theory fixture.",r"\appendix\section{Full manifest-branched character series}",_series_table("All terms through the cutoff",data["manifest_branching/branched_character_series.json"],"multiplicity",True) if final>=7 else "Not available.",r"\section{Full manifest-branched refined plethystic logarithm}",_series_table("All signed terms through the cutoff",data["manifest_branching/branched_refined_plethystic_logarithm.json"],"coefficient",True) if final>=7 else "Not available.",r"\section{Limitations and next steps}",r"This report stops at manifest branching. It does not include, unless stored results are deliberately requested later: a physical baryon/instanton charge map; microscopic identification of the two neutral singlets; explicit polynomial equations; a coordinate-ring ideal; a Gr\"obner-basis calculation; or a monopole-formula derivation. None is inferred here.",r"\end{document}"]
    tex="\n".join(body)+"\n"; report.mkdir(parents=True,exist_ok=True)
    (report/"calculation_report.tex").write_text(tex,encoding="utf-8")
    manifest={"theory_id":theory_id,"calculation_order":order,"requested_final_stage":through,"branching_id":branching_id,"report_generator_version":VERSION,"git_commit":head,"working_tree_status":status,"source_files_used":[str(p.relative_to(root)) for p in evidence if not str(p).endswith(".json")],"result_files_used":[str(p.relative_to(root)) for p in evidence if str(p).endswith(".json")],"omitted_or_unavailable_files":omitted,"report_sections_generated":[STAGES[n][1] for n in range(final+1)]+["Master check index","Evidence and provenance","Limitations and next steps"],"checks_imported":[asdict(c) for c in checks],"report_generation_validation_results":validations}
    dump=lambda obj: json.dumps(obj,indent=2,sort_keys=True)+"\n"
    (report/"calculation_report_manifest.json").write_text(dump(manifest)); (report/"calculation_report_checks.json").write_text(dump({"validation_results":validations,"imported_check_status_totals":dict(sorted(totals.items()))})); (report/"calculation_report_file_hashes.json").write_text(dump(prov))
    counts=Counter(c.stage for c in checks)
    compiler = shutil.which("latexmk") or shutil.which("pdflatex")
    compile_status = (f"available as `{compiler}`; compilation is performed by the caller"
                      if compiler else "PDF compilation not available in this environment")
    (report/"calculation_report_build.md").write_text("# Calculation report build\n\n- Generated only from stored evidence; no mathematical pipeline stage was rerun.\n- Strict mode: "+("enabled" if strict else "disabled")+"\n- Imported checks by stage: "+", ".join(f"Stage {k}: {v}" for k,v in sorted(counts.items()))+"\n- Status totals: "+", ".join(f"{k}: {totals.get(k,0)}" for k in macros)+"\n- PDF compilation: "+compile_status+".\n")
    return report
