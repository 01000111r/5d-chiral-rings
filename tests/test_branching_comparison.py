"""Focused exact tests for the reusable comparison workflow."""
from pathlib import Path
import hashlib
import pytest
from sage.all import QQ, matrix
from hwg_pipeline.branching import branch_irrep, D5_EMBEDDING
from hwg_pipeline.branching_comparison import (ComparisonError, finite_physical,
    physical_charge, render_degree_table, render_parent, solve_charge_map,
    _fmt_rep)
from hwg_pipeline.model import SimpleGroupSpec
from hwg_pipeline.sage_backend import irrep_dimension

def g(n): return SimpleGroupSpec(str(n),'A',n,f'SU({n+1})',tuple(map(str,range(n))))
def pieces(l): return {(tuple(map(int,p.child_dynkin_labels)),int(p.x_charge)):int(p.multiplicity) for p in branch_irrep(g(5),g(4),l)}
def test_fundamental(): assert pieces((1,0,0,0,0))=={((1,0,0,0),1):1,((0,0,0,0),-5):1}
def test_adjoint(): assert pieces((1,0,0,0,1))=={((1,0,0,1),0):1,((0,0,0,0),0):1,((1,0,0,0),6):1,((0,0,0,1),-6):1}
def test_antisymmetric(): assert pieces((0,1,0,0,0))=={((0,1,0,0),2):1,((1,0,0,0),-4):1}
def test_conjugate(): assert pieces((0,0,0,1,0))=={((0,0,1,0),-2):1,((0,0,0,1),4):1}
def test_dimension_preservation():
 p=branch_irrep(g(5),g(4),(2,1,0,0,1)); assert irrep_dimension(g(5),(2,1,0,0,1))==sum(x.multiplicity*irrep_dimension(g(4),x.child_dynkin_labels) for x in p)
def test_external_q_and_negative_multiplicity_preserved():
 q=-2;m=-3; assert all(q==-2 and m*x.multiplicity<0 for x in branch_irrep(g(5),g(4),(0,1,0,0,0)))
def anchors(): return [{'raw':[6,0],'physical':[0,1]},{'raw':[2,1],'physical':[-3,0]}]
def test_exact_anchor_solution():
 M,R,T=solve_charge_map(anchors()); assert M==matrix(QQ, [[0,-3],[QQ(1)/6,-QQ(1)/3]]) and M*R==T
def test_allowed_sublattice_integrality(): assert physical_charge(2,1,solve_charge_map(anchors())[0])==(-3,0)
def test_nonintegral_rejected():
 with pytest.raises(ComparisonError): physical_charge(1,0,solve_charge_map(anchors())[0])
def test_finite_beta_conversion(): assert finite_physical(-1)==(-3,0)
def parent(): return {'parent_q_charge':1,'parent_su6_labels':[0,1,0,0,0],'parent_pl_multiplicity':-2,'children':[{'child_su5_labels':[0,1,0,0],'x_charge':2,'q_charge':1,'signed_total_multiplicity':-2,'B':-3,'I':0}]}
def test_parent_preserving_rendering():
 s=render_parent(parent(),True)
 assert '\\longrightarrow' in s and '[0,1,0,0]' in s and '-2' in s
 assert ']_{5;\\,B=-3,I=0}' in s and '}_{' not in s
def test_charge_annotation_uses_one_braced_subscript():
 assert _fmt_rep((1,0,0,0),5,'x=6,q=0') == '[1,0,0,0]_{5;\\,x=6,q=0}'
def test_side_by_side_table_rendering():
 s=render_degree_table([(2,'a','b')],'finite','UV'); assert 'longtable' in s and 'finite & UV' in s
def test_complete_cutoff_inclusion():
 s=render_degree_table([(d,str(d),str(d)) for d in range(11)],'f','u'); assert all(f'{d} &' in s for d in range(11))
def test_deterministic_render_generation():
 rows=[(2,'a','b'),(10,'c','d')]; assert render_degree_table(rows,'f','u')==render_degree_table(rows,'f','u')
def test_branching_deterministic(): assert branch_irrep(g(5),g(4),(2,1,0,1,0))==branch_irrep(g(5),g(4),(2,1,0,1,0))

def a6pieces(labels):
 return {(tuple(map(int,p.child_dynkin_labels)),int(p.x_charge)):int(p.multiplicity) for p in branch_irrep(g(6),g(5),labels)}
def test_a6_fundamental_and_antifundamental():
 assert a6pieces((1,0,0,0,0,0))=={((1,0,0,0,0),1):1,((0,0,0,0,0),-6):1}
 assert a6pieces((0,0,0,0,0,1))=={((0,0,0,0,1),-1):1,((0,0,0,0,0),6):1}
def test_a6_adjoint_and_antisymmetrics():
 assert a6pieces((1,0,0,0,0,1))=={((1,0,0,0,1),0):1,((0,0,0,0,0),0):1,((1,0,0,0,0),7):1,((0,0,0,0,1),-7):1}
 assert a6pieces((0,0,1,0,0,0))=={((0,0,1,0,0),3):1,((0,1,0,0,0),-4):1}
 assert a6pieces((0,0,0,1,0,0))=={((0,0,1,0,0),-3):1,((0,0,0,1,0),4):1}
def test_a6_twenty_self_conjugate_and_exact_map():
 assert tuple(reversed((0,0,1,0,0)))==(0,0,1,0,0)
 M,R,T=solve_charge_map([{'raw':[7,0],'physical':[0,1]},{'raw':[3,1],'physical':[-3,0]}])
 assert M==matrix(QQ,[[0,-3],[QQ(1)/7,-QQ(3)/7]]) and M*R==T

def d5pieces(labels):
 return {(tuple(map(int,p.child_dynkin_labels)),int(p.x_charge)):int(p.multiplicity)
         for p in branch_irrep(SimpleGroupSpec('d','D',5,'SO(10)',tuple('abcde')),g(4),labels,D5_EMBEDDING)}
def test_d5_vector_convention(): assert d5pieces((1,0,0,0,0))=={((1,0,0,0),-2):1,((0,0,0,1),2):1}
def test_d5_adjoint_convention(): assert d5pieces((0,1,0,0,0))=={((1,0,0,1),0):1,((0,0,0,0),0):1,((0,0,1,0),4):1,((0,1,0,0),-4):1}
def test_d5_spinor_nodes_remain_distinct():
 assert d5pieces((0,0,0,1,0))=={((0,0,1,0),1):1,((1,0,0,0),-3):1,((0,0,0,0),5):1}
 assert d5pieces((0,0,0,0,1))=={((0,1,0,0),-1):1,((0,0,0,1),3):1,((0,0,0,0),-5):1}

from hwg_pipeline.branching_comparison import (FactorIrrep, ProductIrrep,
    branch_product, parse_terms, product_dimension, render_product,
    render_product_parent, su2_weights)

def product_meta(): return [
 {'index':0,'cartan_factor_id':'enhanced','cartan_type':'A4','display_name':'SU(5)','action':'preserve'},
 {'index':1,'cartan_factor_id':'su2','cartan_type':'A1','display_name':'SU(2)','action':'branch_to_u1','output_charge':'x','convention':'su2_weight'}]
def stored(reps,coefficient=-2,q=1): return {'theory_id':'synthetic','coefficients_by_t_degree':{'3':[{'abelian_charges':{'q':str(q)},'coefficient':coefficient,'irreducible_representations':reps}]}}
def reps(): return [{'cartan_factor_id':'enhanced','dynkin_labels':[0,1,0,0]},{'cartan_factor_id':'su2','dynkin_labels':[1]}]
def test_parse_single_factor_stored_term():
 t=parse_terms(stored([{'cartan_factor_id':'flavour','dynkin_labels':[1,0,0,0]}]),[{'index':0,'cartan_factor_id':'flavour','cartan_type':'A4','display_name':'SU(5)'}])[0]; assert len(t['product_irrep'].factors)==1
def test_parse_two_factor_and_order_preserved():
 t=parse_terms(stored(reps()),product_meta())[0]; assert [f.cartan_factor_id for f in t['product_irrep'].factors]==['enhanced','su2'] and [f.labels for f in t['product_irrep'].factors]==[(0,1,0,0),(1,)]
@pytest.mark.parametrize('bad',[[{'cartan_factor_id':'enhanced','dynkin_labels':[0,1,0,0,1]},{'cartan_factor_id':'su2','dynkin_labels':[]}],[{'cartan_factor_id':'enhanced','dynkin_labels':[0,1,0,0,1]}],[]])
def test_reject_bad_per_factor_lengths_missing_or_concatenated(bad):
 with pytest.raises(ComparisonError): parse_terms(stored(bad),product_meta())
@pytest.mark.parametrize(('label','expected'),[(0,(0,)),(1,(1,-1)),(2,(2,0,-2)),(3,(3,1,-1,-3))])
def test_exact_a1_weights(label,expected): assert su2_weights((label,))==expected
def product_term(b=2,m=-2,q=1):
 return {'degree':3,'multiplicity':m,'charges':{'q':q},'product_irrep':ProductIrrep((FactorIrrep('enhanced','A4','SU(5)',(0,1,0,0)),FactorIrrep('su2','A1','SU(2)',(b,))))}
def test_product_branch_preserves_a4_q_sign_and_dimensions():
 p=branch_product(product_term(),product_meta()); assert all(c['child_factors'][0]['labels']==[0,1,0,0] and c['raw_charges']['q']==1 and c['signed_child_multiplicity']==-2 for c in p['children']); assert p['parent_dimension']==sum(c['child_dimension'] for c in p['children'])
def test_product_dimension(): assert product_dimension(product_term()['product_irrep'])==product_dimension(ProductIrrep((product_term()['product_irrep'].factors[0],)))*3
def test_product_json_roundtrip():
 p=product_term()['product_irrep']; assert ProductIrrep.from_json(p.to_json())==p
def test_product_latex_semicolon_and_single_factor_compatibility():
 assert render_product(product_term()['product_irrep'])=='[0,1,0,0;2]_{5,2}'
 assert render_product(ProductIrrep((product_term()['product_irrep'].factors[0],)))=='[0,1,0,0]_{5}'
def test_product_parent_physical_rendering():
 p=branch_product(product_term(1),product_meta()); p['children']=[{**c,'physical_charges':{'B':-3,'I':0}} for c in p['children']]; assert '[0,1,0,0;1]_{5,2}' in render_product_parent(p,True) and 'B=-3,I=0' in render_product_parent(p,True)
def test_product_deterministic(): assert branch_product(product_term(),product_meta())==branch_product(product_term(),product_meta())
def test_product_matching_key_discards_parent_a1():
 a=branch_product(product_term(1),product_meta())['children'][0]; b=branch_product(product_term(3),product_meta())['children'][0]; assert a['child_factors']==b['child_factors']

def double_meta(): return [
 {'index':0,'cartan_factor_id':'a5','cartan_type':'A5','display_name':'SU(6)','action':'preserve'},
 {'index':1,'cartan_factor_id':'a1_1','cartan_type':'A1','display_name':'SU(2)_1','action':'branch_to_u1','output_charge':'x','convention':'su2_weight'},
 {'index':2,'cartan_factor_id':'a1_2','cartan_type':'A1','display_name':'SU(2)_2','action':'branch_to_u1','output_charge':'y','convention':'su2_weight'}]
def double_term(b=1,c=1,m=-2):
 return {'degree':3,'multiplicity':m,'charges':{},'product_irrep':ProductIrrep((FactorIrrep('a5','A5','SU(6)',(0,0,1,0,0)),FactorIrrep('a1_1','A1','SU(2)_1',(b,)),FactorIrrep('a1_2','A1','SU(2)_2',(c,))))}
def test_three_factor_parse_preserves_repeated_a1_identity_order_and_trivial_factors():
 payload=stored([{'cartan_factor_id':'a5','dynkin_labels':[0,0,1,0,0]},{'cartan_factor_id':'a1_1','dynkin_labels':[1]},{'cartan_factor_id':'a1_2','dynkin_labels':[0]}],q=0)
 payload['coefficients_by_t_degree']['3'][0]['abelian_charges']={}
 t=parse_terms(payload,double_meta())[0]
 assert [(f.cartan_factor_id,f.labels) for f in t['product_irrep'].factors]==[('a5',(0,0,1,0,0)),('a1_1',(1,)),('a1_2',(0,))]
def test_repeated_factor_exchange_is_rejected():
 payload=stored([{'cartan_factor_id':'a5','dynkin_labels':[0,0,1,0,0]},{'cartan_factor_id':'a1_2','dynkin_labels':[1]},{'cartan_factor_id':'a1_1','dynkin_labels':[0]}])
 with pytest.raises(ComparisonError): parse_terms(payload,double_meta())
@pytest.mark.parametrize(('b','c'),[(1,1),(2,1),(1,2),(2,2)])
def test_double_a1_branching_uses_distinct_arbitrary_charge_names(b,c):
 p=branch_product(double_term(b,c),double_meta())
 assert [(z['raw_charges']['x'],z['raw_charges']['y']) for z in p['children']]==[(x,y) for x in su2_weights((b,)) for y in su2_weights((c,))]
 assert len(p['children'])==(b+1)*(c+1) and p['parent_dimension']==sum(z['child_dimension'] for z in p['children'])
 assert all(z['signed_child_multiplicity']==-2 and z['child_factors'][0]['labels']==[0,0,1,0,0] for z in p['children'])
def test_three_factor_json_latex_and_no_external_q():
 p=branch_product(double_term(),double_meta())
 assert ProductIrrep.from_json(double_term()['product_irrep'].to_json())==double_term()['product_irrep']
 assert render_product(double_term()['product_irrep'])=='[0,0,1,0,0;1;1]_{6,2,2}'
 assert 'q' not in str(p) and 'q' not in render_product_parent(p)
def test_double_a1_exact_charge_map_and_integrality_rejection():
 M,R,T=solve_charge_map([{'raw':[2,0],'physical':[3,1]},{'raw':[1,1],'physical':[3,0]}])
 assert M==matrix(QQ,[[QQ(3)/2,QQ(3)/2],[QQ(1)/2,-QQ(1)/2]]) and M*R==T
 assert physical_charge(1,-1,M)==(0,1)
 with pytest.raises(ComparisonError): physical_charge(1,0,M)
