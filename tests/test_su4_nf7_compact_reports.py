from pathlib import Path
ROOT=Path(__file__).parents[1]
def test_finite_and_uv_compact_conventions():
 finite=(ROOT/'generated/su4_nf7_finite/order_10/compact_report/compact_results.tex').read_text()
 uv=(ROOT/'generated/su4_nf7_k3o2_infinite/order_10/compact_report/compact_results.tex').read_text()
 assert 'N_c=4,N_f=7' in finite and 'B=4B_\\beta' in finite
 assert '\\mu_4\\mu_3t^8-\\mu_4\\mu_3t^8=0' in finite and '\\beta' in finite
 assert 'q' in uv and 'Rational-product form' in uv and 'N/A' in uv
 assert not any(x in finite+uv for x in ('CharacterRing','sage.','"coefficients'))
