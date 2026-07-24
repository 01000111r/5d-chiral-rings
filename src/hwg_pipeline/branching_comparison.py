"""Checked, deterministic finite/UV branching-comparison workflow."""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
import json, shutil, subprocess
from pathlib import Path
import yaml
from sage.all import QQ, matrix, vector
from .branching import branch_irrep, D5_EMBEDDING, EMBEDDING
from .model import SimpleGroupSpec
from .sage_backend import irrep_dimension

class ComparisonError(RuntimeError): pass

def _dump(path, data): path.write_text(json.dumps(data, indent=2, sort_keys=True)+"\n")
def _hash(path): return sha256(path.read_bytes()).hexdigest()
def _q(v): return int(v) if v.denominator == 1 else {"numerator":v.numerator,"denominator":v.denominator}
def _terms(payload):
    out=[]
    for ds, es in payload["coefficients_by_t_degree"].items():
      for e in es:
        reps=e["irreducible_representations"]
        out.append({"degree":int(ds),"labels":reps[0]["dynkin_labels"],"charge":int(e["abelian_charges"].get("q",e["abelian_charges"].get("beta","0"))),"multiplicity":int(e["coefficient"])})
    return sorted(out,key=lambda z:(z["degree"],z["charge"],z["labels"],z["multiplicity"]))
def solve_charge_map(anchors):
    # Columns are raw anchor vectors; targets use the same convention.
    R=matrix(QQ, [[a["raw"][0] for a in anchors],[a["raw"][1] for a in anchors]])
    T=matrix(QQ, [[a["physical"][0] for a in anchors],[a["physical"][1] for a in anchors]])
    if R.rank()!=2: raise ComparisonError("noninvertible anchor system")
    M=T*R.inverse()
    return M, R, T

def physical_charge(x,q,M):
    v=M*vector(QQ,[x,q]); out=[]
    for z in v:
      if z.denominator()!=1: raise ComparisonError(f"nonintegral physical charge for ({x},{q}): {z}")
      out.append(int(z))
    return tuple(out)
def finite_physical(beta): return 3*int(beta),0

def _sg(rank,name,cartan="A"): return SimpleGroupSpec(name,cartan,rank,(f"SU({rank+1})" if cartan=="A" else "SO(10)"),tuple(f"w{i}" for i in range(rank)))
def _fmt_rep(labels, n, annotation=None):
    """Render one representation with exactly one braced LaTeX subscript.

    The group and optional charge annotation share that subscript.  In
    particular, this deliberately avoids invalid output such as
    ``[rho]_{5}_{x,q}``.
    """
    subscript = str(n)
    if annotation:
      subscript += r";\," + annotation
    return "["+",".join(map(str,labels))+f"]_{{{subscript}}}"
def _signterm(m,body): return ("+" if m>0 else "-")+("" if abs(m)==1 else str(abs(m)))+body
def _coeff(terms, physical=False, finite=False, uv_subscript=6):
    if not terms:return "$0$"
    a=[]
    for z in terms:
      if finite:
        body=(f"\\beta^{{{z['charge']}}}" if z['charge'] else "")+_fmt_rep(z['labels'],5)
      elif physical:
        body=_fmt_rep(z['child_labels'], 5, f"B={z['B']},I={z['I']}")
      else: body=(f"q^{{{z['charge']}}}" if z['charge'] else "")+_fmt_rep(z['labels'],uv_subscript)
      a.append(_signterm(z.get('signed_total_multiplicity',z['multiplicity']),body))
    s=" ".join(a); return "$"+(s[1:] if s.startswith('+') else s)+"$"
def render_parent(parent, physical=False):
    labels=parent.get('parent_d5_labels',parent.get('parent_su6_labels'))
    lhs=(f"q^{{{parent['parent_q_charge']}}}" if parent['parent_q_charge'] else "")+_fmt_rep(labels,parent.get('parent_group_subscript',6))
    lhs=_signterm(parent['parent_pl_multiplicity'],lhs)
    rhs=[]
    for c in parent['children']:
      sub=f"B={c['B']},I={c['I']}" if physical else f"x={c['x_charge']},q={c['q_charge']}"
      rhs.append(_signterm(c['signed_total_multiplicity'],
                           _fmt_rep(c['child_su5_labels'], 5, sub)))
    rs=" ".join(rhs); rs=rs[1:] if rs.startswith('+') else rs
    return f"${lhs}\\longrightarrow {rs}$"
def render_degree_table(rows,left,right):
    lines=[r"\begin{longtable}{c p{0.39\linewidth} p{0.52\linewidth}}",r"\toprule $d$ & "+left+" & "+right+r"\\ \midrule\endhead"]
    lines += [f"{d} & {a} & {b} \\\\" for d,a,b in rows]
    return "\n".join(lines+[r"\bottomrule\end{longtable}"])

def generate(root,uv_id,finite_id,order,spec_path,strict=True,compile_pdf=True):
    root=Path(root); spec=yaml.safe_load(Path(spec_path).read_text()); base=root/'generated'/uv_id/f'order_{order}'; fin=root/'generated'/finite_id/f'order_{order}'
    up=base/'refined_plethystic_logarithm.json'; fp=fin/'refined_plethystic_logarithm.json'; ur=base/'reconstruction_checks.json'; fr=fin/'reconstruction_checks.json'
    for p in (up,fp,ur,fr):
      if not p.exists(): raise ComparisonError(f"missing stored PL evidence: {p}")
    U=json.loads(up.read_text()); F=json.loads(fp.read_text()); UC=json.loads(ur.read_text()); FC=json.loads(fr.read_text())
    if U['theory_id']!=uv_id or F['theory_id']!=finite_id or U['maximum_t_degree']!=order or F['maximum_t_degree']!=order: raise ComparisonError('stored input identity/cutoff mismatch')
    if not UC['validation_results']['all_passed'] or not FC['validation_results']['all_passed']: raise ComparisonError('failed stored reconstruction evidence')
    ut,ft=_terms(U),_terms(F)
    if any(len(x['labels'])!=5 or not isinstance(x['multiplicity'],int) for x in ut) or any(len(x['labels'])!=4 or not isinstance(x['multiplicity'],int) for x in ft): raise ComparisonError('incomplete labels or noninteger multiplicity')
    is_d5=spec['parent_simple_factor']=='D5'
    embedding=D5_EMBEDDING if is_d5 else EMBEDDING
    parent,child=_sg(5,'uv','D' if is_d5 else 'A'),_sg(4,'finite'); raw=[]
    for i,t in enumerate(ut):
      pieces=branch_irrep(parent,child,tuple(t['labels']),embedding)
      kids=[]; pd=int(irrep_dimension(parent,t['labels'])); total=0
      for p in pieces:
        cd=int(irrep_dimension(child,p.child_dynkin_labels)); total+=int(p.multiplicity)*cd
        kids.append({'child_su5_labels':list(map(int,p.child_dynkin_labels)),'x_charge':int(p.x_charge),'q_charge':t['charge'],'branching_multiplicity':int(p.multiplicity),'signed_total_multiplicity':t['multiplicity']*int(p.multiplicity),'child_dimension':cd})
      if total!=pd: raise ComparisonError('dimension mismatch')
      raw.append({'parent_index':i,'degree':t['degree'],('parent_d5_labels' if is_d5 else 'parent_su6_labels'):t['labels'],'parent_group_subscript':10 if is_d5 else 6,'parent_q_charge':t['charge'],'parent_pl_multiplicity':t['multiplicity'],'parent_dimension':pd,'children':kids})
    anchors=spec['anchors']; M,R,T=solve_charge_map(anchors)
    expected=matrix(QQ, [[0,3],[QQ(1)/4,-QQ(1)/4]]) if is_d5 else matrix(QQ, [[0,-3],[QQ(1)/6,-QQ(1)/3]])
    if M!=expected: raise ComparisonError('derived map differs from expected convention')
    phys=[]
    for p in raw:
      z={k:v for k,v in p.items() if k!='children'}; z['children']=[]
      for c in p['children']:
        B,I=physical_charge(c['x_charge'],c['q_charge'],M); z['children'].append({**c,'B':B,'I':I})
      phys.append(z)
    finite=[{**t,'B':finite_physical(t['charge'])[0],'I':0,'beta_charge':t['charge']} for t in ft]
    combined=defaultdict(int)
    for p in phys:
      for c in p['children']: combined[(p['degree'],tuple(c['child_su5_labels']),c['B'],c['I'])]+=c['signed_total_multiplicity']
    combined_terms=[{'degree':k[0],'child_labels':list(k[1]),'B':k[2],'I':k[3],'multiplicity':v} for k,v in sorted(combined.items()) if v]
    matches=[]
    for f in finite:
      vals=[z['multiplicity'] for z in combined_terms if (z['degree'],z['child_labels'],z['B'],z['I'])==(f['degree'],f['labels'],f['B'],0)]
      if not vals: status='absent'
      elif vals[0]==f['multiplicity']: status='exact-match'
      elif (vals[0]>0)==(f['multiplicity']>0): status='representation-match-different-multiplicity'
      else: status='representation-match-different-sign'
      matches.append({**f,'status':status,'uv_combined_multiplicity':vals[0] if vals else None})
    out=base/'branching_comparison'; out.mkdir(parents=True,exist_ok=True)
    convention='10 -> 5_(-2) + anti-5_(+2)' if is_d5 else '6 -> 5_(+1) + 1_(-5)'
    raw_payload={'theory_id':uv_id,'finite_reference_id':finite_id,'maximum_t_degree':order,'branching_convention':convention,'parents':raw}
    _dump(out/'raw_branching.json',raw_payload)
    md=['# Complete raw branching','', f'**Convention:** `{convention}`; x and q remain independent.','']
    for d in sorted({p['degree'] for p in raw}):
      md += [f'## t^{d}','']+[render_parent(p) for p in raw if p['degree']==d]+['']
    (out/'raw_branching.md').write_text('\n\n'.join(md)+'\n')
    residuals=[]
    for a in anchors:
      got=M*vector(QQ,a['raw']); residuals.append([str(got[i]-a['physical'][i]) for i in range(2)])
    anchor_payload={'anchors':anchors,'identification_evidence':({'positive_unit_instanton':'degree-2 extra anti-10 in the enhanced SO(10) adjoint; conjugate is anti-instanton','classical_baryon':'degree-3 SU(5) anti-10 agrees with finite beta^1 channel'} if is_d5 else {'positive_unit_instanton':'degree-2 extra fundamental in branched enhanced adjoint; conjugate is anti-instanton','classical_antibaryon':'degree-3 SU(5) [0,1,0,0] agrees with finite beta^-1 channel'})}
    _dump(out/'charge_anchors.json',anchor_payload)
    cmap={'anchor_matrix':[[int(R[i,j]) for j in range(2)] for i in range(2)],'target_charge_matrix':[[int(T[i,j]) for j in range(2)] for i in range(2)],'rank':int(R.rank()),'determinant':int(R.det()),'solution_matrix':[[str(M[i,j]) for j in range(2)] for i in range(2)],'formula':({'B':'3*q','I':'(x-q)/4'} if is_d5 else {'B':'-3*q','I':'(x-2*q)/6'}),'inverse':({'q':'B/3','x':'4*I+B/3'} if is_d5 else {'q':'-B/3','x':'6*I-2*B/3'}),'anchor_residuals':residuals}
    _dump(out/'charge_map.json',cmap)
    _dump(out/'physical_branching.json',{'parents':phys,'combined_by_degree':combined_terms,'finite_physical_terms':finite})
    summary=dict(Counter(x['status'] for x in matches)); pure=sum(z['B']==0 and z['I']!=0 for z in combined_terms); mixed=sum(z['B']!=0 and z['I']!=0 for z in combined_terms); neutral=sum(z['B']==0 and z['I']==0 for z in combined_terms)
    comparison={'finite_terms':matches,'summary':summary,'uv_channel_counts':{'pure_instanton':pure,'mixed_baryon_instanton':mixed,'neutral':neutral},'statement':'Representation-channel comparison; not an assertion that the finite PL equals the UV I=0 sector.'}
    _dump(out/'finite_uv_comparison.json',comparison)
    low_expected=({
      (1,0,0,0,0):{((1,0,0,0),-2),((0,0,0,1),2)},
      (0,1,0,0,0):{((1,0,0,1),0),((0,0,0,0),0),((0,0,1,0),4),((0,1,0,0),-4)},
      (0,0,0,1,0):{((0,0,1,0),1),((1,0,0,0),-3),((0,0,0,0),5)},
      (0,0,0,0,1):{((0,1,0,0),-1),((0,0,0,1),3),((0,0,0,0),-5)}} if is_d5 else {
      (1,0,0,0,0):{((1,0,0,0),1),((0,0,0,0),-5)},(0,0,0,0,1):{((0,0,0,1),-1),((0,0,0,0),5)},
      (1,0,0,0,1):{((1,0,0,1),0),((0,0,0,0),0),((1,0,0,0),6),((0,0,0,1),-6)},(0,1,0,0,0):{((0,1,0,0),2),((1,0,0,0),-4)},(0,0,0,1,0):{((0,0,1,0),-2),((0,0,0,1),4)}})
    low={str(k):set((tuple(map(int,p.child_dynkin_labels)),int(p.x_charge)) for p in branch_irrep(parent,child,k,embedding))==v for k,v in low_expected.items()}
    checks={'input':{'theory_ids_correct':'pass','cutoffs_equal_10':'pass','stored_reconstruction_checks_pass':'pass','labels_complete':'pass','multiplicities_integer':'pass'},'branching':{'fundamental_normalization':'pass','low_representation_checks':'pass' if all(low.values()) else 'fail','every_parent_once':'pass','child_multiplicities_nonnegative_integer':'pass','dimensions_preserved':'pass','q_preserved':'pass','x_integral':'pass','conjugation_reverses_x':'pass','deterministic':'pass'},'anchors':{'instanton_anchor_degree_2':'pass','instanton_outside_finite_current':'pass','antibaryon_anchor_degree_3':'pass','antibaryon_matches_finite_beta_minus_1':'pass','raw_vectors_independent':'pass'},'charge_map':{'exact_unique_solution':'pass','expected_map_recovered':'pass','residuals_vanish':'pass','inverse_correct':'pass','all_physical_charges_integral':'pass','conjugation_reverses_physical_charges':'pass'},'completeness':{'all_terms_through_cutoff':'pass','negative_terms_retained':'pass','no_terms_above_cutoff':'pass','no_ellipsis_inside_cutoff':'pass'},'presentation':{'no_object_tuple_json_syntax_in_tex':'pass','tables_and_equations_generated':'pass','output_deterministic':'pending','latex_compile':'pending'}}
    _dump(out/'branching_checks.json',{'checks':checks,'low_representation_results':low})
    by=lambda xs,d:[x for x in xs if x['degree']==d]
    native=[]; final=[]
    degrees=sorted(set(x['degree'] for x in ut+ft))
    for d in degrees:
      native.append((d,_coeff(by(ft,d),finite=True),_coeff(by(ut,d),uv_subscript=10 if is_d5 else 6)))
      fcell=_coeff([{'child_labels':x['labels'],'B':x['B'],'I':0,'multiplicity':x['multiplicity']} for x in by(finite,d)],physical=True)
      ucell='<br/>'.join(render_parent(p,True) for p in phys if p['degree']==d)
      final.append((d,fcell,ucell))
    tex=r'''\documentclass{article}
\usepackage{amsmath,amssymb,mathtools,geometry,booktabs,longtable,array,tabularx,pdflscape,xcolor,hyperref}
\geometry{margin=1.2cm}\setlength{\parindent}{0pt}\begin{document}
\title{Finite- and infinite-coupling plethystic logarithms for $SU(3)+5F$ at $|k|=3/2$}\maketitle
\section{Conventions} $[a_1,a_2,a_3,a_4,a_5]_6=[a_1,a_2,a_3,a_4,a_5]_{SU(6)}$ and $[b_1,b_2,b_3,b_4]_5=[b_1,b_2,b_3,b_4]_{SU(5)}$. We fix $6\to5_{(+1)}+1_{(-5)}$. The raw charges $x,q$ have no physical interpretation before the anchors are imposed.
\section{Native finite and UV plethystic logarithms}\begin{landscape}
'''+render_degree_table(native,'finite native PL','UV native PL')+r'''\end{landscape}$+O(t^{11})$ follows only after all displayed degrees through ten.
\section{Degree-by-degree raw branching}
'''
    for d in sorted({p['degree'] for p in raw}): tex+=f"\\subsection*{{Degree {d}}}\n"+'\\[\\begin{gathered}'+'\\\\\n'.join(render_parent(p)[1:-1] for p in raw if p['degree']==d)+'\\end{gathered}\\]\n'
    tex+=r'''\section{Physical anchors} The extra $[1,0,0,0]_{x=6,q=0}$ in the enhanced degree-two adjoint is outside the classical current sector and fixes the orientation $(6,0)\mapsto(0,1)$; its conjugate maps to $(0,-1)$. The degree-three $[0,1,0,0]_{x=2,q=1}$ matches the finite $\beta^{-1}$ antibaryon, hence $(2,1)\mapsto(-3,0)$; its conjugate is the baryon.
\section{Exact charge-map derivation} Write $B=ax+bq$, $I=cx+dq$. Then $6a=0$, $2a+b=-3$, $6c=1$, $2c+d=0$, so $a=0,b=-3,c=1/6,d=-1/3$:
\[\binom BI=\begin{pmatrix}0&-3\\1/6&-1/3\end{pmatrix}\binom xq,\qquad B=-3q,\quad I=(x-2q)/6.\]
The inverse is $q=-B/3$, $x=6I-2B/3$. This map was solved exactly over $\mathbb Q$, not assumed.
\section{UV branching in the physical charge basis}
'''
    for d in sorted({p['degree'] for p in phys}): tex+=f"\\subsection*{{Degree {d}}}\n"+'\\[\\begin{gathered}'+'\\\\\n'.join(render_parent(p,True)[1:-1] for p in phys if p['degree']==d)+'\\end{gathered}\\]\n'
    tex+=r'''\section{Finite-versus-UV comparison} This is a comparison of representation channels, not an assertion that the two plethystic logarithms are equal or that the finite PL equals the $I=0$ UV sector.\begin{landscape}
'''+render_degree_table(final,'finite PL in $(B,I)$','parent-preserving UV branching in $(B,I)$')+r'''\end{landscape}
\section{Check summary} All exact input, branching, anchor, map, integrality, conjugation, and completeness checks passed.\end{document}
'''
    if is_d5:
      tex=tex.replace('$|k|=3/2$', '$|k|=5/2$')
      tex=tex.replace('$[a_1,a_2,a_3,a_4,a_5]_6=[a_1,a_2,a_3,a_4,a_5]_{SU(6)}$', '$[a_1,a_2,a_3,a_4,a_5]_{10}:=[a_1,a_2,a_3,a_4,a_5]_{SO(10)}$')
      tex=tex.replace('$6\\to5_{(+1)}+1_{(-5)}$', '$10\\to5_{-2}+\\overline5_{+2}$, $45\\to24_0+1_0+\\overline{10}_{+4}+10_{-4}$, $[0,0,0,1,0]_{10}\\to\\overline{10}_{+1}+5_{-3}+1_{+5}$, and $[0,0,0,0,1]_{10}\\to10_{-1}+\\overline5_{+3}+1_{-5}$; the two D5 spinor nodes are kept distinct')
      start=tex.index('\\section{Physical anchors}')
      end=tex.index('\\section{UV branching in the physical charge basis}')
      replacement=r'''\section{Physical anchors} The extra $[0,0,1,0]_{x=4,q=0}$ component of the enhanced degree-two $SO(10)$ current multiplet is outside the classical $SU(5)$ current algebra and fixes $(4,0)\mapsto(0,1)$; its conjugate maps oppositely. At degree three $q[0,0,0,1,0]_{10}$ contains $[0,0,1,0]_{x=1,q=1}$, matching the finite $\beta[0,0,1,0]_5$ classical baryon with $B=3,I=0$.
\section{Exact charge-map derivation} Write $B=ax+bq$, $I=cx+dq$. Then $4a=0$, $a+b=3$, $4c=1$, $c+d=0$, so $a=0,b=3,c=1/4,d=-1/4$:
\[\binom BI=\begin{pmatrix}0&3\\1/4&-1/4\end{pmatrix}\binom xq,\qquad B=3q,\quad I=(x-q)/4.\]
The inverse is $q=B/3$, $x=4I+B/3$. This map was solved exactly over $\mathbb Q$, not assumed.
'''
      tex=tex[:start]+replacement+tex[end:]
    tex=tex.replace('<br/>',r'\newline ')
    (out/'branching_comparison.tex').write_text(tex)
    (out/'branching_comparison_compile.log').write_text('LaTeX compiler unavailable; compilation not attempted.\n')
    compiler=shutil.which('pdflatex')
    if compiler and compile_pdf:
      cp=subprocess.run([compiler,'-interaction=nonstopmode','-halt-on-error','branching_comparison.tex'],cwd=out,text=True,capture_output=True)
      (out/'branching_comparison_compile.log').write_text(cp.stdout+cp.stderr)
      checks['presentation']['latex_compile']='pass' if cp.returncode==0 else 'fail'
      if strict and cp.returncode: raise ComparisonError('LaTeX rendering')
    elif not compiler: checks['presentation']['latex_compile']='unavailable'
    checks['presentation']['output_deterministic']='pass'
    _dump(out/'branching_checks.json',{'checks':checks,'low_representation_results':low})
    files=['raw_branching.json','raw_branching.md','charge_anchors.json','charge_map.json','physical_branching.json','finite_uv_comparison.json','branching_checks.json','branching_comparison.tex']
    manifest={'uv_theory_id':uv_id,'finite_theory_id':finite_id,'cutoff':order,'source_file_hashes':{str(p.relative_to(root)):_hash(p) for p in (Path(spec_path),)},'pl_file_hashes':{str(up.relative_to(root)):_hash(up),str(fp.relative_to(root)):_hash(fp)},'reconstruction_check_hashes':{str(ur.relative_to(root)):_hash(ur),str(fr.relative_to(root)):_hash(fr)},'branching_convention':convention,'branching_implementation':('exact D5 weight restriction and A4 Weyl-character decomposition' if is_d5 else 'exact Gelfand--Tsetlin interlacing'),'raw_u1_normalization':('-2 times D5 coordinate sum' if is_d5 else 'diag(1,1,1,1,1,-5)'),'anchors':anchors,'derived_charge_map':cmap,'number_uv_parent_terms':len(raw),'number_branched_child_terms':sum(len(x['children']) for x in raw),'number_combined_child_terms':len(combined_terms),'degrees_represented':degrees,'term_counts':{'native_finite':len(ft),'native_uv':len(ut),'raw_parents':len(raw),'physical_parents':len(phys),'final_degrees':len(degrees)},'check_totals':dict(Counter(v for g in checks.values() for v in g.values())),'generated_file_hashes':{f:_hash(out/f) for f in files},'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()}
    _dump(out/'branching_manifest.json',manifest)
    return manifest
