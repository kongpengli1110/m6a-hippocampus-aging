"""Synthetic-only mathematical tests; never load a real expression matrix or run actual H2."""

# Release-only safety/path adapter; see provenance/adaptation_diffs.
from pathlib import Path as _ReleasePath
import sys as _release_sys
_release_home = next(p for p in _ReleasePath(__file__).resolve().parents if (p / 'release_support.py').is_file())
_release_sys.path.insert(0, str(_release_home))
from release_support import gate, release_root, resolve_source
gate(__file__, __name__)
import unittest,sys,math,itertools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
import pandas as pd
from scipy import sparse
from statsmodels.stats.multitest import multipletests
from h2_core import *


def fixture():
    # 4 fixed candidate slots, 8 controls, two occupied bins plus optional empty bins.
    ids=['C1','C2','C3','C4']+[f'N{i}' for i in range(8)]
    bg=pd.DataFrame({'HGNC_ID':ids,'ensembl_id':['E'+x for x in ids],
        'expression_bin':[1,1,2,2]+[1]*4+[2]*4,'detection_bin':[1]*12,
        'noncandidate_pool':[False]*4+[True]*8})
    m=pd.DataFrame({'HGNC_ID':ids[:4],'weight':[1.,4.,2.,8.]})
    union=set(ids[:4]);p=build_sampling_plan(bg,m,('expression_bin','detection_bin'),union,min_count=1)
    return bg,m,union,p

class ReferenceTests(unittest.TestCase):
    def test_full_feature_denominator_equal_donors_not_count_weighted(self):
        # First two rows are same donor/class and must sum; third is a second donor with very different depth.
        s=pd.DataFrame({'donor':['a','a','b','a','b','a','b'],
          'Class':['G','G','G','I','I','N','N'],'n_nuclei':[10,10,20,20,20,20,20]})
        X=sparse.csr_matrix([[1,0,4],[1,0,4],[9,0,1],[0,3,7],[0,4,6],[0,0,5],[0,0,8]])
        a=aggregate_reference(X,s,['G','I','N'],min_nuclei=20,min_donors=2)
        self.assertEqual(a['n_by_class'],{'G':2,'I':2,'N':2})
        self.assertAlmostEqual(a['means'][0,0],550000.) # mean(2/10,9/10); not candidate-only denominator
        self.assertAlmostEqual(a['means'][0,2],450000.)
        self.assertEqual(a['manifest'].n_nuclei.sum(),120)
        # Different depths make pooled-cell/library means demonstrably different.
        X2=sparse.csr_matrix([[2,0,8],[90,0,10],[0,3,7],[0,4,6],[0,0,5],[0,0,8]])
        s2=s.drop(index=1).reset_index(drop=True);s2.loc[0,'n_nuclei']=20
        aa=aggregate_reference(X2,s2,['G','I','N'],20,2)
        self.assertAlmostEqual(aa['means'][0,0],550000.)
        self.assertNotAlmostEqual(aa['means'][0,0],92/110*1e6)
    def test_filter_after_sum_and_min_donors(self):
        s=pd.DataFrame({'donor':['a','a','b'],'Class':['G']*3,'n_nuclei':[10,10,19]})
        a=aggregate_reference([[1,9],[2,8],[5,5]],s,['G'],20,1)
        self.assertEqual(a['n_by_class']['G'],1)
        with self.assertRaises(NotTestable):aggregate_reference([[1,9],[2,8],[5,5]],s,['G'],20,2)
    def test_zero_library_not_sample(self):
        s=pd.DataFrame({'donor':['a','b'],'Class':['G']*2,'n_nuclei':[20]*2})
        a=aggregate_reference([[0,0],[1,1]],s,['G'],20,1)
        self.assertEqual(a['n_by_class']['G'],1)
    def test_share_sum_and_zero_gene(self):
        shares,ok=specificity_shares([[2,0,1],[1,0,1],[1,0,1]])
        np.testing.assert_array_equal(ok,[True,False,True]);self.assertTrue(np.isnan(shares[:,1]).all())
        np.testing.assert_allclose(shares[:,ok].sum(axis=0),[1,1]);self.assertEqual(shares[0,0],.5)
        self.assertAlmostEqual(shares[0,2],1/3)
    def test_bad_means_fail(self):
        for x in ([[1,-1]],[[1,float('nan')]]):
            with self.assertRaises(ValueError):specificity_shares(x)

class MatchingTests(unittest.TestCase):
    def test_no_candidate_controls_and_no_duplicates_every_draw(self):
        bg,m,u,p=fixture();rng=np.random.default_rng(30)
        for _ in range(200):
            ix=draw_positions(p,rng)
            self.assertEqual(len(set(ix)),p.n_slots)
            self.assertFalse(set(bg.iloc[ix].HGNC_ID)&u)
            for slot,j in enumerate(ix):self.assertEqual(tuple(bg.loc[j,list(p.bin_columns)]),p.slot_bins[slot])
    def test_candidate_union_not_only_current_set(self):
        bg,m,u,p=fixture();bg.loc[4,'noncandidate_pool']=False;u.add('N0')
        pp=build_sampling_plan(bg,m,p.bin_columns,u,min_count=1)
        for _,_,pool in pp.bins:self.assertNotIn(4,pool)
    def test_inconsistent_pool_fails(self):
        bg,m,u,p=fixture();bg.loc[0,'noncandidate_pool']=True
        with self.assertRaises(ValueError):build_sampling_plan(bg,m,p.bin_columns,u,min_count=1)
    def test_weights_and_bin_slots_preserved(self):
        bg,m,u,p=fixture();values=np.linspace(.05,.9,len(bg));null,audit=generate_null(p,values,np.random.default_rng(4),40,(0,1,17,39))
        self.assertEqual(sorted(p.weights),sorted(m.weight));self.assertEqual(p.weight_sum,m.weight.sum())
        for b in [0,1,17,39]:
            a=[r for r in audit if r['draw_0based']==b]
            self.assertEqual(Counter(r['bin'] for r in a),Counter(':'.join(map(str,x)) for x in p.slot_bins))
            self.assertEqual([r['slot_weight'] for r in a],list(p.weights))
            independent=sum(r['slot_weight']*r['control_Glut_share'] for r in a)/sum(r['slot_weight'] for r in a)
            self.assertAlmostEqual(independent,null[b],places=14)
    def test_random_slot_assignment_not_fixed_by_weight_or_gene(self):
        bg,m,u,p=fixture();rng=np.random.default_rng(9);seen=set()
        for _ in range(100):seen.add(tuple(draw_positions(p,rng)[:2]))
        # Both assignments to two differently weighted slots are observed in this deterministic fixture.
        self.assertTrue(any((b,a) in seen for a,b in seen));self.assertGreater(len(seen),8)
    def test_empty_unoccupied_bins_ok(self):
        bg,m,u,p=fixture();self.assertEqual(len(p.bins),2)
        # No bin 3 rows exist and no candidate requests it; there is no forced draw from empty bins.
        generate_null(p,np.linspace(.1,.9,len(bg)),np.random.default_rng(1),5,())
    def test_empty_requested_bin_is_not_testable(self):
        bg,m,u,p=fixture();bg.loc[0,'expression_bin']=99
        with self.assertRaisesRegex(NotTestable,'insufficient_pool'):build_sampling_plan(bg,m,p.bin_columns,u,min_count=1)
    def test_insufficient_pool_no_replacement_fallback(self):
        bg,m,u,p=fixture();bg=bg.drop(index=[4,5,6]).reset_index(drop=True)
        with self.assertRaisesRegex(NotTestable,'insufficient_pool'):build_sampling_plan(bg,m,p.bin_columns,u,min_count=1)
    def test_coverage_failure_not_rescue(self):
        bg,m,u,p=fixture();bg=bg[bg.HGNC_ID!='C4']
        with self.assertRaisesRegex(NotTestable,'coverage'):build_sampling_plan(bg,m,p.bin_columns,u,min_count=1)
    def test_three_dimensional_length_bins_are_preserved(self):
        bg,m,u,p=fixture()
        bg['length_bin']=[1,2,1,2,1,1,2,2,1,1,2,2]
        pp=build_sampling_plan(bg,m,('expression_bin','detection_bin','length_bin'),u,min_count=1)
        for _ in range(50):
            ix=draw_positions(pp,np.random.default_rng(_))
            for slot,j in enumerate(ix):self.assertEqual(tuple(bg.loc[j,list(pp.bin_columns)]),pp.slot_bins[slot])
    def test_unweighted_and_scaled_weights(self):
        self.assertAlmostEqual(weighted_score([.1,.3,.8],[1,1,1]),.4)
        self.assertAlmostEqual(weighted_score([.1,.3,.8],[2,4,6]),weighted_score([.1,.3,.8],[20,40,60]))
    def test_seed_reproducible_and_substreams_separate(self):
        bg,m,u,p=fixture();v=np.linspace(.01,.99,len(bg))
        c1=np.random.SeedSequence(20261012).spawn(16);c2=np.random.SeedSequence(20261012).spawn(16)
        a,_=generate_null(p,v,np.random.default_rng(c1[0]),40,());b,_=generate_null(p,v,np.random.default_rng(c2[0]),40,());d,_=generate_null(p,v,np.random.default_rng(c2[8]),40,())
        np.testing.assert_array_equal(a,b);self.assertFalse(np.array_equal(a,d))
        self.assertEqual(len({tuple(c.generate_state(4)) for c in c1}),16)
    def test_tiny_enumerated_null_matches_allowed_support(self):
        # Exactly one candidate slot, two legal controls; synthetic RNG outputs must lie in enumerated support.
        bg=pd.DataFrame({'HGNC_ID':['C','A','B'],'ensembl_id':['EC','EA','EB'],'expression_bin':[1]*3,'noncandidate_pool':[False,True,True]})
        p=build_sampling_plan(bg,pd.DataFrame({'HGNC_ID':['C'],'weight':[7.]}),('expression_bin',),{'C'},1)
        n,_=generate_null(p,np.array([.4,.2,.8]),np.random.default_rng(2),100,())
        np.testing.assert_allclose(np.unique(n),[.2,.8])

class InferenceFormulaTests(unittest.TestCase):
    def test_ties_in_both_tails(self):
        s=empirical_summary(2.,[1,2,2,3]);self.assertEqual(s['p_lower'],4/5);self.assertEqual(s['p_upper'],4/5);self.assertEqual(s['raw_p'],1.)
    def test_extreme_two_tailed_plus_one_not_zero(self):
        for obs in [-1,10]:self.assertEqual(empirical_summary(obs,[1,2,3,4])['raw_p'],2/5)
    def test_null_reference_band_not_effect_ci(self):
        s=empirical_summary(.8,[.1,.2,.3,.4]);self.assertAlmostEqual(s['Delta'],.55)
        self.assertEqual(s['null_median'],.25);self.assertAlmostEqual(s['null_reference_lower_025'],.1075)
        self.assertIn('NOT_effect_CI',s['null_reference_interval_type'])
    def test_constant_null(self):
        self.assertEqual(empirical_summary(.5,[.5]*9)['raw_p'],1)
    def test_bh_matches_independent_library(self):
        p=np.array([.01,.06,.02,.6,1,.003,.8,.02]);q=bh_fixed(p,8)
        np.testing.assert_allclose(q,multipletests(p,method='fdr_bh')[1],rtol=0,atol=1e-15)
    def test_bh_keeps_missing_position_and_fixed_family(self):
        p=np.array([.001,np.nan,.02,.3,.4,.5,.8,1]);q=bh_fixed(p,8)
        np.testing.assert_allclose(q,multipletests(np.nan_to_num(p,nan=1),method='fdr_bh')[1])
        self.assertEqual(q[1],1)
        with self.assertRaises(ValueError):bh_fixed(p[:7],8)
    def test_invalid_p_fails(self):
        with self.assertRaises(ValueError):bh_fixed([-.1]+[1]*7,8)

if __name__=='__main__':unittest.main(verbosity=2)
