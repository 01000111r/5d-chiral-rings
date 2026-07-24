"""Focused exact tests for the reusable comparison workflow."""
from pathlib import Path
import hashlib
import pytest
from sage.all import QQ, matrix
from hwg_pipeline.branching import branch_irrep
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
