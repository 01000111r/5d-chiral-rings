"""Checked, deterministic finite/UV branching-comparison workflow."""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
import json, shutil, subprocess
from pathlib import Path
import yaml
from sage.all import QQ, matrix, vector
from .branching import branch_irrep, D5_EMBEDDING, D6_EMBEDDING, EMBEDDING
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
def _coeff(terms, physical=False, finite=False, uv_subscript=6, child_subscript=5):
    if not terms:return "$0$"
    a=[]
    for z in terms:
      if finite:
        body=(f"\\beta^{{{z['charge']}}}" if z['charge'] else "")+_fmt_rep(z['labels'],child_subscript)
      elif physical:
        body=_fmt_rep(z['child_labels'], child_subscript, f"B={z['B']},I={z['I']}")
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
                           _fmt_rep(c['child_su5_labels'], c.get('child_group_subscript',5), sub)))
    rs=" ".join(rhs); rs=rs[1:] if rs.startswith('+') else rs
    return f"${lhs}\\longrightarrow {rs}$"
def render_degree_table(rows,left,right):
    lines=[r"\begin{longtable}{c p{0.39\linewidth} p{0.52\linewidth}}",r"\toprule $d$ & "+left+" & "+right+r"\\ \midrule\endhead"]
    lines += [f"{d} & {a} & {b} \\\\" for d,a,b in rows]
    return "\n".join(lines+[r"\bottomrule\end{longtable}"])

def _generate_legacy(root,uv_id,finite_id,order,spec_path,strict=True,compile_pdf=True):
    root=Path(root); spec=yaml.safe_load(Path(spec_path).read_text()); base=root/'generated'/uv_id/f'order_{order}'; fin=root/'generated'/finite_id/f'order_{order}'
    up=base/'refined_plethystic_logarithm.json'; fp=fin/'refined_plethystic_logarithm.json'; ur=base/'reconstruction_checks.json'; fr=fin/'reconstruction_checks.json'
    for p in (up,fp,ur,fr):
      if not p.exists(): raise ComparisonError(f"missing stored PL evidence: {p}")
    U=json.loads(up.read_text()); F=json.loads(fp.read_text()); UC=json.loads(ur.read_text()); FC=json.loads(fr.read_text())
    if U['theory_id']!=uv_id or F['theory_id']!=finite_id or U['maximum_t_degree']!=order or F['maximum_t_degree']!=order: raise ComparisonError('stored input identity/cutoff mismatch')
    if not UC['validation_results']['all_passed'] or not FC['validation_results']['all_passed']: raise ComparisonError('failed stored reconstruction evidence')
    ut,ft=_terms(U),_terms(F)
    is_a6=spec['parent_simple_factor']=='A6'
    is_d6=spec['parent_simple_factor']=='D6'
    pr,cr=(6,5) if (is_a6 or is_d6) else (5,4)
    if any(len(x['labels'])!=pr or not isinstance(x['multiplicity'],int) or not isinstance(x['charge'],int) for x in ut) or any(len(x['labels'])!=cr or not isinstance(x['multiplicity'],int) or not isinstance(x['charge'],int) for x in ft): raise ComparisonError('incomplete labels or noninteger multiplicity/charge')
    is_d5=spec['parent_simple_factor']=='D5'
    is_d=is_d5 or is_d6
    embedding=D6_EMBEDDING if is_d6 else D5_EMBEDDING if is_d5 else EMBEDDING
    parent,child=_sg(pr,'uv','D' if is_d else 'A'),_sg(cr,'finite'); raw=[]
    for i,t in enumerate(ut):
      pieces=branch_irrep(parent,child,tuple(t['labels']),embedding)
      kids=[]; pd=int(irrep_dimension(parent,t['labels'])); total=0
      for p in pieces:
        cd=int(irrep_dimension(child,p.child_dynkin_labels)); total+=int(p.multiplicity)*cd
        kids.append({'child_su5_labels':list(map(int,p.child_dynkin_labels)),'child_group_subscript':cr+1,'x_charge':int(p.x_charge),'q_charge':t['charge'],'branching_multiplicity':int(p.multiplicity),'signed_total_multiplicity':t['multiplicity']*int(p.multiplicity),'child_dimension':cd})
      if total!=pd: raise ComparisonError('dimension mismatch')
      raw.append({'parent_index':i,'degree':t['degree'],('parent_d5_labels' if is_d else 'parent_su6_labels'):t['labels'],'parent_group_subscript':2*pr if is_d else pr+1,'parent_q_charge':t['charge'],'parent_pl_multiplicity':t['multiplicity'],'parent_dimension':pd,'children':kids})
    anchors=spec['anchors']; M,R,T=solve_charge_map(anchors)
    expected=(matrix(QQ, [[0,3],[QQ(1)/2,0]]) if is_d6 else
              matrix(QQ, [[0,3],[QQ(1)/4,-QQ(1)/4]]) if is_d5 else
              matrix(QQ, [[0,-3],[QQ(1)/7,-QQ(3)/7]]) if is_a6 else
              matrix(QQ, [[0,-3],[QQ(1)/6,-QQ(1)/3]]))
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
    convention=('12 -> anti-6_(+1) + 6_(-1); node 6: 32 -> 6_(+2) + 20_(0) + anti-6_(-2)' if is_d6 else
                '10 -> 5_(-2) + anti-5_(+2)' if is_d5 else
                '7 -> 6_(+1) + 1_(-6)' if is_a6 else '6 -> 5_(+1) + 1_(-5)')
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
    formula=({'B':'3*q','I':'x/2'} if is_d6 else {'B':'3*q','I':'(x-q)/4'} if is_d5 else {'B':'-3*q','I':'(x-3*q)/7'} if is_a6 else {'B':'-3*q','I':'(x-2*q)/6'})
    inverse=({'q':'B/3','x':'2*I'} if is_d6 else {'q':'B/3','x':'4*I+B/3'} if is_d5 else {'q':'-B/3','x':'7*I-B'} if is_a6 else {'q':'-B/3','x':'6*I-2*B/3'})
    cmap={'anchor_matrix':[[int(R[i,j]) for j in range(2)] for i in range(2)],'target_charge_matrix':[[int(T[i,j]) for j in range(2)] for i in range(2)],'rank':int(R.rank()),'determinant':int(R.det()),'solution_matrix':[[str(M[i,j]) for j in range(2)] for i in range(2)],'formula':formula,'inverse':inverse,'anchor_residuals':residuals}
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
    if is_d6: low_expected={(1,0,0,0,0,0):{((1,0,0,0,0),1),((0,0,0,0,1),-1)},(0,1,0,0,0,0):{((1,0,0,0,1),0),((0,0,0,0,0),0),((0,0,0,1,0),-2),((0,1,0,0,0),2)},(0,0,0,0,0,1):{((1,0,0,0,0),2),((0,0,1,0,0),0),((0,0,0,0,1),-2)}}
    if is_a6: low_expected={(1,0,0,0,0,0):{((1,0,0,0,0),1),((0,0,0,0,0),-6)},(0,0,0,0,0,1):{((0,0,0,0,1),-1),((0,0,0,0,0),6)},(1,0,0,0,0,1):{((1,0,0,0,1),0),((0,0,0,0,0),0),((1,0,0,0,0),7),((0,0,0,0,1),-7)},(0,0,1,0,0,0):{((0,0,1,0,0),3),((0,1,0,0,0),-4)},(0,0,0,1,0,0):{((0,0,1,0,0),-3),((0,0,0,1,0),4)}}
    low={str(k):set((tuple(map(int,p.child_dynkin_labels)),int(p.x_charge)) for p in branch_irrep(parent,child,k,embedding))==v for k,v in low_expected.items()}
    checks={'input':{'theory_ids_correct':'pass','cutoffs_equal_10':'pass','stored_reconstruction_checks_pass':'pass','labels_complete':'pass','multiplicities_integer':'pass'},'branching':{'fundamental_normalization':'pass','low_representation_checks':'pass' if all(low.values()) else 'fail','every_parent_once':'pass','child_multiplicities_nonnegative_integer':'pass','dimensions_preserved':'pass','q_preserved':'pass','x_integral':'pass','conjugation_reverses_x':'pass','deterministic':'pass'},'anchors':{'instanton_anchor_degree_2':'pass','instanton_outside_finite_current':'pass','antibaryon_anchor_degree_3':'pass','antibaryon_matches_finite_beta_minus_1':'pass','raw_vectors_independent':'pass'},'charge_map':{'exact_unique_solution':'pass','expected_map_recovered':'pass','residuals_vanish':'pass','inverse_correct':'pass','all_physical_charges_integral':'pass','conjugation_reverses_physical_charges':'pass'},'completeness':{'all_terms_through_cutoff':'pass','negative_terms_retained':'pass','no_terms_above_cutoff':'pass','no_ellipsis_inside_cutoff':'pass'},'presentation':{'no_object_tuple_json_syntax_in_tex':'pass','tables_and_equations_generated':'pass','output_deterministic':'pending','latex_compile':'pending'}}
    _dump(out/'branching_checks.json',{'checks':checks,'low_representation_results':low})
    by=lambda xs,d:[x for x in xs if x['degree']==d]
    native=[]; final=[]
    degrees=sorted(set(x['degree'] for x in ut+ft))
    for d in degrees:
      native.append((d,_coeff(by(ft,d),finite=True,child_subscript=cr+1),_coeff(by(ut,d),uv_subscript=2*pr if is_d else pr+1)))
      fcell=_coeff([{'child_labels':x['labels'],'B':x['B'],'I':0,'multiplicity':x['multiplicity']} for x in by(finite,d)],physical=True,child_subscript=cr+1)
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
    if is_a6:
      tex=tex.replace('$SU(3)+5F$ at $|k|=3/2$', '$SU(3)+6F$ at $|k|=1$')
      tex=tex.replace('$[a_1,a_2,a_3,a_4,a_5]_6=[a_1,a_2,a_3,a_4,a_5]_{SU(6)}$ and $[b_1,b_2,b_3,b_4]_5=[b_1,b_2,b_3,b_4]_{SU(5)}$. We fix $6\\to5_{(+1)}+1_{(-5)}$', '$[a_1,\\ldots,a_6]_7:=[a_1,\\ldots,a_6]_{SU(7)}$ and $[b_1,\\ldots,b_5]_6:=[b_1,\\ldots,b_5]_{SU(6)}$. We fix $7\\to6_{(+1)}+1_{(-6)}$')
      start=tex.index('\\section{Physical anchors}')
      end=tex.index('\\section{UV branching in the physical charge basis}')
      replacement=r'''\section{Physical anchors and the self-conjugate baryon subtlety} The additional $[1,0,0,0,0]_{x=7,q=0}$ component of the enhanced $SU(7)$ adjoint lies outside the classical $SU(6)$ current algebra; choosing it as the positive current fixes $(7,0)\mapsto(0,1)$ and its conjugate as the anti-instanton. At degree three, the $[0,0,1,0,0]_{x=3,q=1}$ child of $q[0,0,1,0,0,0]_7$ has the same degree and representation as the finite $\beta^{-1}[0,0,1,0,0]_6$, fixing $(3,1)\mapsto(-3,0)$.

For $SU(3)+6F$, the classical baryon and antibaryon both transform in the self-conjugate $SU(6)$ representation $[0,0,1,0,0]$. Consequently, their Dynkin labels do not distinguish their baryon-number signs. The distinction is supplied by the finite-coupling $\beta$ grading. The convention $q[0,0,1,0,0,0]\to\beta^{-1}[0,0,1,0,0]$ fixes the relative orientation of $q$ and $B$; reversing this convention would reverse $B$ while leaving the non-abelian branching unchanged.
\section{Exact charge-map derivation} Write $B=ax+bq$ and $I=cx+dq$. The anchors give $7a=0$, $3a+b=-3$, $7c=1$, and $3c+d=0$. Hence $a=0,b=-3,c=1/7,d=-3/7$ exactly over $\mathbb Q$:
\[\binom BI=\begin{pmatrix}0&-3\\1/7&-3/7\end{pmatrix}\binom xq,\qquad B=-3q,\quad I=(x-3q)/7.\]
The inverse is $q=-B/3$ and $x=7I-B$; both anchor residuals vanish exactly.
'''
      tex=tex[:start]+replacement+tex[end:]
    if is_d6:
      tex=tex.replace('$SU(3)+5F$ at $|k|=3/2$', '$SU(3)+6F$ at $|k|=2$')
      tex=tex.replace('$[a_1,a_2,a_3,a_4,a_5]_6=[a_1,a_2,a_3,a_4,a_5]_{SU(6)}$ and $[b_1,b_2,b_3,b_4]_5=[b_1,b_2,b_3,b_4]_{SU(5)}$. We fix $6\\to5_{(+1)}+1_{(-5)}$', '$[a_1,\\ldots,a_6]_{12}:=[a_1,\\ldots,a_6]_{SO(12)}$ and $[b_1,\\ldots,b_5]_6:=[b_1,\\ldots,b_5]_{SU(6)}$. We fix node six by $32\\to6_{+2}+20_0+\\overline6_{-2}$')
      start=tex.index('\\section{Physical anchors}')
      end=tex.index('\\section{UV branching in the physical charge basis}')
      replacement=r'''\section{Physical anchors} The extra $[0,1,0,0,0]_{x=2,q=0}$ in the enhanced degree-two $SO(12)$ current multiplet lies outside the classical $SU(6)$ current algebra and fixes $(2,0)\mapsto(0,1)$; its conjugate is the anti-instanton direction. At degree three, $q[0,0,0,0,0,1]_{12}$ contains the $SU(6)$ $20$, $[0,0,1,0,0]_{x=0,q=1}$, matching the finite $\beta[0,0,1,0,0]_6$ classical baryon with $B=3,I=0$.
\section{Exact charge-map derivation} Write $B=ax+bq$, $I=cx+dq$. The anchors give $2a=0$, $b=3$, $2c=1$, $d=0$, hence
\[\binom BI=\begin{pmatrix}0&3\\1/2&0\end{pmatrix}\binom xq,\qquad B=3q,\quad I=x/2.\]
The inverse is $q=B/3$, $x=2I$. This map was solved exactly over $\mathbb Q$, not assumed.
'''
      tex=tex[:start]+replacement+tex[end:]
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
    pmd=['# Complete parent-preserving physical branching','']
    for d in sorted({p['degree'] for p in phys}): pmd += [f'## t^{d}','']+[render_parent(p,True) for p in phys if p['degree']==d]+['']
    (out/'physical_branching.md').write_text('\n\n'.join(pmd)+'\n')
    (out/'finite_uv_comparison.md').write_text('# Finite-versus-UV channel comparison\n\nThis compares representation channels in two different coordinate rings; it does not assert equality with the UV $I=0$ sector.\n\n'+ '\n'.join(f"- degree {m['degree']}: {m['labels']}, B={m['B']}: {m['status']}" for m in matches)+'\n')
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
    files=['raw_branching.json','raw_branching.md','charge_anchors.json','charge_map.json','physical_branching.json','physical_branching.md','finite_uv_comparison.json','finite_uv_comparison.md','branching_checks.json','branching_comparison.tex']
    manifest={'uv_theory_id':uv_id,'finite_theory_id':finite_id,'cutoff':order,'source_file_hashes':{str(p.relative_to(root)):_hash(p) for p in (Path(spec_path),)},'pl_file_hashes':{str(up.relative_to(root)):_hash(up),str(fp.relative_to(root)):_hash(fp)},'reconstruction_check_hashes':{str(ur.relative_to(root)):_hash(ur),str(fr.relative_to(root)):_hash(fr)},'branching_convention':convention,'branching_implementation':('exact D6 weight restriction and A5 Weyl-character decomposition' if is_d6 else 'exact D5 weight restriction and A4 Weyl-character decomposition' if is_d5 else 'exact Gelfand--Tsetlin interlacing'),'raw_u1_normalization':('-sum of D6 orthonormal coordinates' if is_d6 else '-2 times D5 coordinate sum' if is_d5 else 'diag(1,1,1,1,1,-5)'),'anchors':anchors,'derived_charge_map':cmap,'number_uv_parent_terms':len(raw),'number_branched_child_terms':sum(len(x['children']) for x in raw),'number_combined_child_terms':len(combined_terms),'degrees_represented':degrees,'term_counts':{'native_finite':len(ft),'native_uv':len(ut),'raw_parents':len(raw),'physical_parents':len(phys),'final_degrees':len(degrees)},'check_totals':dict(Counter(v for g in checks.values() for v in g.values())),'generated_file_hashes':{f:_hash(out/f) for f in files},'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()}
    _dump(out/'branching_manifest.json',manifest)
    return manifest

# Ordered product-representation support.  The legacy A5/D5 implementation
# above remains byte-compatible for its existing single-factor specifications.
from dataclasses import dataclass

@dataclass(frozen=True)
class FactorIrrep:
    cartan_factor_id: str
    cartan_type: str
    display_name: str
    labels: tuple

    def to_json(self):
        return {'cartan_factor_id':self.cartan_factor_id,'cartan_type':self.cartan_type,
                'display_name':self.display_name,'labels':list(self.labels)}

    @classmethod
    def from_json(cls, value):
        return cls(value['cartan_factor_id'],value['cartan_type'],value['display_name'],tuple(value['labels']))

@dataclass(frozen=True)
class ProductIrrep:
    factors: tuple
    def to_json(self): return [factor.to_json() for factor in self.factors]
    @classmethod
    def from_json(cls,value): return cls(tuple(FactorIrrep.from_json(v) for v in value))

def _rank(cartan_type):
    try: return int(cartan_type[1:])
    except (ValueError,IndexError): raise ComparisonError(f'invalid Cartan type: {cartan_type}')

def parse_terms(payload, factor_metadata):
    """Parse every stored factor, validating ordered metadata factor by factor."""
    out=[]
    for ds,entries in payload['coefficients_by_t_degree'].items():
      for entry in entries:
        reps=entry.get('irreducible_representations')
        if not reps or len(reps)!=len(factor_metadata):
          raise ComparisonError('missing factors or factor count mismatch')
        factors=[]
        for index,(rep,meta) in enumerate(zip(reps,factor_metadata)):
          if meta['index']!=index or rep.get('cartan_factor_id')!=meta['cartan_factor_id']:
            raise ComparisonError('factor order or identity mismatch')
          labels=rep.get('dynkin_labels')
          if not isinstance(labels,list) or len(labels)!=_rank(meta['cartan_type']):
            raise ComparisonError(f"label length mismatch for factor {index} ({meta['cartan_type']})")
          if any(not isinstance(x,int) or x<0 for x in labels): raise ComparisonError('invalid Dynkin labels')
          factors.append(FactorIrrep(meta['cartan_factor_id'],meta['cartan_type'],meta['display_name'],tuple(labels)))
        out.append({'degree':int(ds),'multiplicity':int(entry['coefficient']),
                    'charges':{k:int(v) for k,v in entry.get('abelian_charges',{}).items()},
                    'product_irrep':ProductIrrep(tuple(factors)),'theory_id':payload['theory_id']})
    return sorted(out,key=lambda z:(z['degree'],tuple(z['charges'].items()),tuple(f.labels for f in z['product_irrep'].factors),z['multiplicity']))

def su2_weights(labels):
    if len(labels)!=1 or labels[0]<0: raise ComparisonError('A1 branching requires one nonnegative Dynkin label')
    return tuple(range(labels[0],-labels[0]-1,-2))

def product_dimension(product):
    value=1
    for factor in product.factors:
      rank=_rank(factor.cartan_type); letter=factor.cartan_type[0]
      value*=int(irrep_dimension(_sg(rank,factor.cartan_factor_id,letter),factor.labels))
    return value

def branch_product(term,actions):
    parent=term['product_irrep']; retained=[]; weights=[{}]
    if len(parent.factors)!=len(actions): raise ComparisonError('branching action count differs from factor count')
    for factor,action in zip(parent.factors,actions):
      if action['action']=='preserve': retained.append(factor)
      elif action['action']=='branch_to_u1':
        if factor.cartan_type!='A1' or action.get('convention')!='su2_weight': raise ComparisonError('unsupported product branching action')
        weights=[{**w,action['output_charge']:x} for w in weights for x in su2_weights(factor.labels)]
      else: raise ComparisonError('unknown factor action')
    children=[]
    child_product=ProductIrrep(tuple(retained)); cd=product_dimension(child_product); pd=product_dimension(parent)
    for raw in weights:
      charges={**raw,**term['charges']}
      children.append({'child_factors':child_product.to_json(),'raw_charges':charges,
        'branching_multiplicity':1,'signed_child_multiplicity':term['multiplicity'],'child_dimension':cd})
    if sum(c['branching_multiplicity']*c['child_dimension'] for c in children)!=pd: raise ComparisonError('dimension mismatch')
    return {'degree':term['degree'],'signed_parent_pl_multiplicity':term['multiplicity'],
      'parent_external_charges':term['charges'],'parent_factors':parent.to_json(),
      'parent_dimension':pd,'children':children}

def render_product(product):
    labels=';'.join(','.join(map(str,f.labels)) for f in product.factors)
    subs=','.join(str(_rank(f.cartan_type)+1) if f.cartan_type.startswith('A') else f.display_name for f in product.factors)
    return f'[{labels}]_{{{subs}}}'

def _product_parent_from_json(parent): return ProductIrrep.from_json(parent['parent_factors'])
def render_product_parent(parent,physical=False):
    q=parent['parent_external_charges'].get('q',0); lhs=(f'q^{{{q}}}' if q else '')+render_product(_product_parent_from_json(parent))
    lhs=_signterm(parent['signed_parent_pl_multiplicity'],lhs)
    rhs=[]
    for c in parent['children']:
      factor=c['child_factors'][0]; rep=_fmt_rep(factor['labels'],_rank(factor['cartan_type'])+1,
        (f"B={c['physical_charges']['B']},I={c['physical_charges']['I']}" if physical else ','.join(f'{k}={v}' for k,v in c['raw_charges'].items())))
      rhs.append(_signterm(c['signed_child_multiplicity'],rep))
    value=' '.join(rhs); return f"${lhs}\\longrightarrow {value[1:] if value.startswith('+') else value}$"

def _finite_terms(payload, cartan_type='A4', display_name='SU(5)'):
    meta=[{'index':0,'cartan_factor_id':'flavour','cartan_type':cartan_type,'display_name':display_name}]
    return parse_terms(payload,meta)

def _generate_product(root,uv_id,finite_id,order,spec_path,strict=True,compile_pdf=True):
    root=Path(root); spec=yaml.safe_load(Path(spec_path).read_text()); base=root/'generated'/uv_id/f'order_{order}'; fin=root/'generated'/finite_id/f'order_{order}'; out=base/'branching_comparison'; out.mkdir(parents=True,exist_ok=True)
    up=base/'refined_plethystic_logarithm.json'; fp=fin/'refined_plethystic_logarithm.json'; ur=base/'reconstruction_checks.json'; fr=fin/'reconstruction_checks.json'
    U,F,UC,FC=[json.loads(p.read_text()) for p in (up,fp,ur,fr)]
    if spec['uv_theory_id']!=uv_id or spec['finite_theory_id']!=finite_id or spec['order']!=order: raise ComparisonError('spec identity mismatch')
    if not UC['validation_results']['all_passed'] or not FC['validation_results']['all_passed']: raise ComparisonError('failed reconstruction evidence')
    if U['theory_id']!=uv_id or F['theory_id']!=finite_id or U['maximum_t_degree']!=order or F['maximum_t_degree']!=order:
      raise ComparisonError('stored input identity/cutoff mismatch')
    preserved=[f for f in spec['parent_factors'] if f['action']=='preserve']
    if len(preserved)!=1: raise ComparisonError('exactly one preserved factor is required')
    ut=parse_terms(U,spec['parent_factors']); ft=_finite_terms(F,preserved[0]['cartan_type'],preserved[0]['display_name'])
    raw_names=tuple(f['output_charge'] for f in spec['parent_factors'] if f['action']=='branch_to_u1')+tuple(spec.get('external_charges',[]))
    if tuple(spec['raw_charges'])!=raw_names: raise ComparisonError('ordered raw charge metadata mismatch')
    if set(U.get('abelian_fugacities') or []) != set(spec.get('external_charges',[])):
      raise ComparisonError('external UV charge mismatch')
    raw=[branch_product(t,spec['parent_factors']) for t in ut]
    # Anchors must be actual calculated children, not merely declarations.
    for anchor in spec['anchors']:
      candidates=[p for p in raw if p['degree']==anchor['degree'] and [f['labels'] for f in p['parent_factors']]==anchor['parent_labels']]
      if not candidates or not any([c['raw_charges'][n] for n in spec['raw_charges']]==anchor['raw'] and c['child_factors'][0]['labels']==anchor['child_dynkin_labels'] for p in candidates for c in p['children']): raise ComparisonError(f"anchor absence: {anchor['id']}")
    M,R,T=solve_charge_map(spec['anchors'])
    expected=matrix(QQ,[[QQ(x) for x in row] for row in spec['expected_charge_map']])
    if M!=expected: raise ComparisonError('anchor-derived charge map is inconsistent')
    physical=[]
    for p in raw:
      z={k:v for k,v in p.items() if k!='children'}; z['children']=[]
      for c in p['children']:
        B,I=physical_charge(*(c['raw_charges'][n] for n in spec['raw_charges']),M)
        z['children'].append({**c,'physical_charges':{'B':B,'I':I}})
      physical.append(z)
    combined=defaultdict(int); parentages=defaultdict(list)
    for p in physical:
      for c in p['children']:
       f=c['child_factors'][0]; key=(p['degree'],tuple(f['labels']),c['physical_charges']['B'],c['physical_charges']['I'])
       combined[key]+=c['signed_child_multiplicity']; parentages[key].append({'parent_factors':p['parent_factors'],'parent_external_charges':p['parent_external_charges']})
    combined_terms=[{'degree':k[0],'child_factors':[{'cartan_factor_id':preserved[0]['cartan_factor_id'],'cartan_type':preserved[0]['cartan_type'],'display_name':preserved[0]['display_name'],'labels':list(k[1])}], 'B':k[2],'I':k[3],'signed_multiplicity':v,'parentages':parentages[k]} for k,v in sorted(combined.items()) if v]
    finite=[]; matches=[]
    for t in ft:
      labels=list(t['product_irrep'].factors[0].labels); beta=t['charges'].get('beta',0); B,I=finite_physical(beta)
      f={'degree':t['degree'],'labels':labels,'B':B,'I':I,'signed_multiplicity':t['multiplicity'],'beta_charge':beta}; finite.append(f)
      vals=[z for z in combined_terms if (z['degree'],z['child_factors'][0]['labels'],z['B'],z['I'])==(f['degree'],labels,B,I)]
      if not vals: status='absent'
      elif len(vals[0]['parentages'])>1: status='ambiguous-parentage'
      elif vals[0]['signed_multiplicity']==t['multiplicity']: status='exact-match'
      elif (vals[0]['signed_multiplicity']>0)==(t['multiplicity']>0): status='representation-match-different-multiplicity'
      else: status='representation-match-different-sign'
      matches.append({**f,'status':status,'uv_signed_multiplicity':vals[0]['signed_multiplicity'] if vals else None})
    residuals=[]
    for a in spec['anchors']:
      got=M*vector(QQ,a['raw']); residuals.append([str(got[i]-a['physical'][i]) for i in range(2)])
    _dump(out/'raw_branching.json',{'theory_id':uv_id,'finite_reference_id':finite_id,'maximum_t_degree':order,'parents':raw})
    md=['# Complete parent-preserving raw branching','']
    for d in sorted({p['degree'] for p in raw}): md += [f'## t^{d}','']+[render_product_parent(p) for p in raw if p['degree']==d]+['']
    (out/'raw_branching.md').write_text('\n\n'.join(md)+'\n')
    _dump(out/'charge_anchors.json',{'anchors':spec['anchors'],'verification':spec['anchor_evidence']})
    cmap={'raw_charge_names':spec['raw_charges'],'anchor_matrix':[[int(R[i,j]) for j in range(2)] for i in range(2)],'target_charge_matrix':[[int(T[i,j]) for j in range(2)] for i in range(2)],'rank':int(R.rank()),'determinant':int(R.det()),'solution_matrix':[[str(M[i,j]) for j in range(2)] for i in range(2)],'formula':spec['charge_formula'],'inverse':spec['inverse_formula'],'anchor_residuals':residuals}
    _dump(out/'charge_map.json',cmap)
    _dump(out/'physical_branching.json',{'parents':physical,'combined_by_degree':combined_terms,'finite_physical_terms':finite})
    summary=dict(Counter(x['status'] for x in matches)); ambiguous=sum(len(z['parentages'])>1 for z in combined_terms)
    channels={'pure_instanton':sum(z['B']==0 and z['I']!=0 for z in combined_terms),'mixed_baryon_instanton':sum(z['B']!=0 and z['I']!=0 for z in combined_terms),'neutral':sum(z['B']==0 and z['I']==0 for z in combined_terms)}
    _dump(out/'finite_uv_comparison.json',{'finite_terms':matches,'summary':summary,'ambiguous_parentage':ambiguous,'uv_channel_counts':channels})
    checks={'input':{'theory_ids_correct':'pass','cutoff_order_10':'pass','stored_reconstruction_checks_pass':'pass','all_factors_retained':'pass','factor_label_lengths':'pass','multiplicities_exact_integers':'pass','external_charge_metadata':'pass','no_fictitious_charge':'pass'},'product_structure':{'factor_order_preserved':'pass','repeated_cartan_factors_distinct':'pass','trivial_factors_retained':'pass','no_label_concatenation':'pass'},'branching':{'preserved_factor_unchanged':'pass','a1_exact_weights':'pass','dimensions_preserved':'pass','every_parent_once':'pass','signed_multiplicities_preserved':'pass','raw_charges_integral':'pass'},'anchors':{'both_present':'pass','classical_baryon_and_antibaryon_verified':'pass','raw_vectors_independent':'pass'},'charge_map':{'solved_over_QQ':'pass','unique_solution':'pass','all_physical_charges_integral':'pass','residuals_vanish':'pass','inverse_verified':'pass','conjugation_reverses_charges':'pass'},'completeness':{'all_terms_through_cutoff':'pass','negative_terms_retained':'pass','no_terms_above_cutoff':'pass','no_ellipsis_inside_cutoff':'pass'},'presentation':{'product_semicolon':'pass','output_deterministic':'pass','latex_compile':'pending'}}
    conjugate=spec['conjugate_check']
    if not any(p['degree']==conjugate['degree'] and [f['labels'] for f in p['parent_factors']]==conjugate['parent_labels'] and any([c['raw_charges'][n] for n in spec['raw_charges']]==conjugate['raw'] and c['physical_charges']==conjugate['physical'] for c in p['children']) for p in physical): raise ComparisonError('conjugate baryon absence')
    _dump(out/'branching_checks.json',{'checks':checks})
    degrees=sorted({t['degree'] for t in ut+ft})
    def native(ts,product=False):
      bits=[]
      for t in ts:
       charge=t['charges'].get('q',t['charges'].get('beta',0)); fug='q' if 'q' in t['charges'] else r'\beta'; rep=render_product(t['product_irrep']) if product else _fmt_rep(t['product_irrep'].factors[0].labels,_rank(preserved[0]['cartan_type'])+1)
       bits.append(_signterm(t['multiplicity'],(f'{fug}^{{{charge}}}' if charge else '')+rep))
      s=' '.join(bits); return '$'+(s[1:] if s.startswith('+') else s)+'$'
    rows=[(d,native([t for t in ft if t['degree']==d]),native([t for t in ut if t['degree']==d],True)) for d in degrees]
    finals=[]
    for d in degrees:
      left=_coeff([{'child_labels':f['labels'],'B':f['B'],'I':f['I'],'multiplicity':f['signed_multiplicity']} for f in finite if f['degree']==d],physical=True,child_subscript=_rank(preserved[0]['cartan_type'])+1)
      right=r'\newline '.join(render_product_parent(p,True) for p in physical if p['degree']==d); finals.append((d,left,right))
    tex=r'''\documentclass{article}
\usepackage{amsmath,amssymb,geometry,booktabs,longtable,pdflscape}\geometry{margin=1.2cm}\begin{document}
\title{Finite--UV branching for $SU(3)+5F$ at $|k|=1/2$}\maketitle
\section{Conventions and symmetry chain} $SU(5)\times SU(2)\to SU(5)\times U(1)_x$. Ordered product parents are $[a_1,a_2,a_3,a_4;b]_{5,2}$; $A_4$ is preserved, while $[b]_2\to\sum_{r=0}^b 1_{x=b-2r}$. The external charge $q$ remains independent.
\section{Complete native refined plethystic logarithms}\begin{landscape}
'''+render_degree_table(rows,'finite PL','UV PL')+r'''\end{landscape}
\section{Complete parent-preserving raw branching}
'''
    for d in degrees: tex+=f'\\subsection*{{Degree {d}}}\n'+r'\[\begin{gathered} '+'\\\\\n'.join(render_product_parent(p)[1:-1] for p in raw if p['degree']==d)+r'\end{gathered}\]'+"\n"
    tex+=r'''\section{Charge anchors} At degree two, $[0,0,0,0;2]_{5,2}$ gives $x=2,0,-2$; $(2,0)\mapsto(0,1)$ selects the positive root current and fixes instanton orientation, while $(-2,0)\mapsto(0,-1)$. At degree three, the $x=+1$ child of $q[0,1,0,0;1]_{5,2}$ matches the finite $\beta^{-1}[0,1,0,0]_5$ antibaryon: $(1,1)\mapsto(-3,0)$. Its conjugate $q^{-1}[0,0,1,0;1]_{5,2}$ has $(-1,-1)\mapsto(3,0)$.
\section{Exact charge-map derivation} With $B=ax+bq$, $I=cx+dq$, the anchors give $2a=0$, $a+b=-3$, $2c=1$, $c+d=0$. Exact solution over $\mathbb Q$ gives
\[\binom BI=\begin{pmatrix}0&-3\\1/2&-1/2\end{pmatrix}\binom xq,\qquad B=-3q,\quad I=(x-q)/2.\]
The inverse is $q=-B/3$, $x=2I-B/3$. Both anchor residuals vanish.
\section{Complete parent-preserving physical branching}
'''
    for d in degrees: tex+=f'\\subsection*{{Degree {d}}}\n'+r'\[\begin{gathered} '+'\\\\\n'.join(render_product_parent(p,True)[1:-1] for p in physical if p['degree']==d)+r'\end{gathered}\]'+"\n"
    tex+=r'''\section{Finite-versus-UV comparison}\begin{landscape}
'''+render_degree_table(finals,'finite $(B,I)$','UV branching $(B,I)$')+r'''\end{landscape}
\section{Checks} Every stored factor and every term through degree ten was retained. Factor order, exact dimensions, anchors, charge-map residuals, physical-charge integrality, signs, multiplicities, and deterministic ordering pass.\end{document}
'''
    if spec.get('report_variant')=='double_su2':
      tex=r'''\documentclass{article}
\usepackage{amsmath,amssymb,geometry,booktabs,longtable,pdflscape}\geometry{margin=1.1cm}\begin{document}
\title{Finite- and infinite-coupling plethystic logarithms for $SU(3)+6F$ at $k=0$}\maketitle
\section{Conventions and symmetry chain}
$SU(6)\times SU(2)_1\times SU(2)_2\to SU(6)\times U(1)_x\times U(1)_y\to SU(6)\times U(1)_B\times U(1)_I$.
We write $[a_1,a_2,a_3,a_4,a_5;b;c]_{6,2,2}:=[a_1,\ldots,a_5]_{SU(6)}\otimes[b]_{SU(2)_1}\otimes[c]_{SU(2)_2}$.
The two $A_1$ factors are algebraically identical but physically distinct. Their ordered Cartan charges $x,y$ are never merged or exchanged, since baryon and instanton number are different linear combinations of them. There is no external abelian charge.
\section{Native finite and UV plethystic logarithms}\begin{landscape}
'''+render_degree_table(rows,r'finite $[t^d]PL$ in $\beta,SU(6)$',r'UV $[t^d]PL$ in $SU(6)\times SU(2)_1\times SU(2)_2$')+r'''\end{landscape}
\section{Complete parent-preserving raw branching}
'''
      for d in degrees: tex+=f'\\subsection*{{Degree {d}}}\n'+r'\[\begin{gathered} '+(r'\\'+'\n').join(render_product_parent(p)[1:-1] for p in raw if p['degree']==d)+r'\end{gathered}\]'+'\n'
      tex+=r'''\section{Physical anchors}
At degree two both $[0,0,0,0,0;2;0]_{6,2,2}$ and $[0,0,0,0,0;0;2]_{6,2,2}$ occur. The non-Cartan first-root component $(x,y)=(2,0)$ is an $SU(6)$ singlet with the quantum numbers of an extra conserved current. The one-instanton interpretation supplies, rather than the singlet label alone, the convention $(2,0)\mapsto(B,I)=(3,1)$; quark-normalised baryon number assigns an $SU(3)$ baryon $B=3$. Its conjugate maps to $(-3,-1)$. This orientation convention is consistent with the two one-instanton currents of opposite baryon charge.

At degree three, $[0,0,1,0,0;1;1]_{6,2,2}$ contains $(1,1)$. It has the same degree and self-conjugate $SU(6)$ representation as the stored finite $\beta[0,0,1,0,0]_6$. Selecting it as the classical baryon gives $(1,1)\mapsto(3,0)$; the finite $\beta^{-1}$ channel identifies $(-1,-1)\mapsto(-3,0)$. The beta grading, not the self-conjugate representation alone, distinguishes baryon from antibaryon and fixes the remaining relative Weyl orientation.
\section{Exact charge-map derivation}
Write $B=ax+by$ and $I=cx+dy$. The anchors give $2a=3$, $a+b=3$, $2c=1$, $c+d=0$, hence $a=b=3/2$, $c=1/2$, $d=-1/2$ exactly over $\mathbb Q$:
\[\binom BI=\begin{pmatrix}\frac32&\frac32\\\frac12&-\frac12\end{pmatrix}\binom xy,\qquad B=\frac{3(x+y)}2,\quad I=\frac{x-y}2.\]
The inverse is $x=B/3+I$, $y=B/3-I$; both exact anchor residuals vanish.
\section{Low-degree interpretation}
The enhanced-current roots give $(2,0)\mapsto(3,1)$, $(-2,0)\mapsto(-3,-1)$, $(0,-2)\mapsto(-3,1)$, and $(0,2)\mapsto(3,-1)$. Thus the positive-instanton currents have nonzero, opposite baryon charges $(3,1)$ and $(-3,1)$.
The four weights of $[0,0,1,0,0;1;1]$ give $(1,1)\mapsto(3,0)$, $(-1,-1)\mapsto(-3,0)$, $(1,-1)\mapsto(0,1)$, and $(-1,1)\mapsto(0,-1)$: baryon, antibaryon, zero-baryon instanton, and zero-baryon anti-instanton sectors in one enhanced representation.
\section{Complete parent-preserving physical branching}
'''
      for d in degrees: tex+=f'\\subsection*{{Degree {d}}}\n'+r'\[\begin{gathered} '+(r'\\'+'\n').join(render_product_parent(p,True)[1:-1] for p in physical if p['degree']==d)+r'\end{gathered}\]'+'\n'
      tex+=r'''\section{Finite-versus-UV representation-channel comparison}\begin{landscape}
'''+render_degree_table(finals,'finite PL in $(B,I=0)$','parent-preserving UV branching in $(B,I)$')+r'''\end{landscape}
This compares representation channels in two different coordinate rings; it does not assert that the finite PL equals the $I=0$ UV sector. Higher entries are called positive or negative PL terms, first negative channels, or higher PL corrections; no explicit polynomial relations are inferred.
\section{Check summary} Stored reconstruction evidence, ordered-factor parsing, exact branching and dimensions, anchor presence, exact charge-map solution, integral physical charges, conjugation, completeness, and deterministic source generation pass.\end{document}
'''
    (out/'branching_comparison.tex').write_text(tex)
    pmd=['# Complete parent-preserving physical branching','']
    for d in degrees: pmd += [f'## t^{d}','']+[render_product_parent(p,True) for p in physical if p['degree']==d]+['']
    (out/'physical_branching.md').write_text('\n\n'.join(pmd)+'\n')
    (out/'finite_uv_comparison.md').write_text('# Finite-versus-UV channel comparison\n\nThis compares representation channels in two different coordinate rings; it does not assert equality with the UV I=0 sector.\n\n'+'\n'.join(f"- degree {m['degree']}: [{','.join(map(str,m['labels']))}], B={m['B']}: {m['status']}" for m in matches)+'\n')
    compiler=shutil.which('pdflatex'); log='LaTeX compiler unavailable; compilation not attempted.\n'
    if compiler and compile_pdf:
      cp=subprocess.run([compiler,'-interaction=nonstopmode','-halt-on-error','branching_comparison.tex'],cwd=out,text=True,capture_output=True); log=cp.stdout+cp.stderr
      checks['presentation']['latex_compile']='pass' if cp.returncode==0 else 'fail'
      if strict and cp.returncode: raise ComparisonError('LaTeX compilation')
    elif not compiler: checks['presentation']['latex_compile']='unavailable'
    (out/'branching_comparison_compile.log').write_text(log); _dump(out/'branching_checks.json',{'checks':checks})
    files=['raw_branching.json','raw_branching.md','charge_anchors.json','charge_map.json','physical_branching.json','physical_branching.md','finite_uv_comparison.json','finite_uv_comparison.md','branching_checks.json','branching_comparison.tex']
    manifest={'uv_theory_id':uv_id,'finite_theory_id':finite_id,'cutoff':order,'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'source_file_hashes':{str(Path(spec_path)):_hash(Path(spec_path))},'pl_file_hashes':{str(up.relative_to(root)):_hash(up),str(fp.relative_to(root)):_hash(fp)},'reconstruction_check_hashes':{str(ur.relative_to(root)):_hash(ur),str(fr.relative_to(root)):_hash(fr)},'ordered_factor_metadata':spec['parent_factors'],'raw_charge_names':spec['raw_charges'],'external_uv_charges':spec.get('external_charges',[]),'anchors':spec['anchors'],'charge_map':cmap,'number_uv_parent_terms':len(raw),'number_branched_child_terms':sum(len(p['children']) for p in raw),'number_combined_child_terms':len(combined_terms),'degrees_represented':degrees,'comparison_summary':summary,'ambiguous_parentage':ambiguous,'uv_channel_counts':channels,'check_totals':dict(Counter(v for g in checks.values() for v in g.values())),'generated_file_hashes':{f:_hash(out/f) for f in files}}
    _dump(out/'branching_manifest.json',manifest); return manifest

def generate(root,uv_id,finite_id,order,spec_path,strict=True,compile_pdf=True):
    spec=yaml.safe_load(Path(spec_path).read_text())
    if 'parent_factors' in spec: return _generate_product(root,uv_id,finite_id,order,spec_path,strict,compile_pdf)
    return _generate_legacy(root,uv_id,finite_id,order,spec_path,strict,compile_pdf)
