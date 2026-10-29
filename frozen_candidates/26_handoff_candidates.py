"""Task 5: freeze mouse candidate list v1 for human-validation handoff. Rules are fixed here BEFORE any human data are seen.
Core list (<=50): mouse genes satisfying ALL of
  E1 MeRIP: hypomethylated with age (PNAS Dataset S17) in >=2 of CA1/CA3/DG
  E2 cell type: Allen-HPF excitatory specificity (CA1+CA3+DG+SUB share) >=0.5
then ranked by evidence score = #regions hypo (2-3) + 1[in m6A-enriched excitatory module with kME>=0.6] + 1[METTL3 target set (cKO HIP down / Cb m6A-loss) or YTHDF1 CLIP]
 + 1[Tegowski cortex Glut DM-RNA] + 1[Allen excitatory pseudobulk age effect padj<0.1 in any of CA1/CA3/DG/SUB]; ties broken by max kME; take top 50 with score>=4.
Human orthologs: MGI HOM_MouseHumanSequence.rpt (one-to-one only flagged; many-to-many kept with flag)."""
import pandas as pd,numpy as np
D='./'; T=D+'analysis/tables/'; H=D+'handoff/'
gs=pd.read_csv(T+'genesets_m6a_aging.csv'); G={k:set(v) for k,v in gs.groupby('set').gene}
rg=pd.read_csv(T+'genesets_regulator_targets.csv'); RG={k:set(v) for k,v in rg.groupby('set').gene}
spec=pd.read_csv(T+'allen_HPF_specificity_EWCE.csv',index_col=0)
mm=pd.read_csv(T+'mod_gene_membership_kME.csv'); hub=pd.read_csv(T+'mod_m6A_enriched_hubs.csv'); ag=pd.read_csv(T+'mod_age_effect_donor.csv').set_index('module')
R=pd.read_csv(T+'allen_HPF_pseudobulk_DESeq2_age.csv'); R=R[R.celltype.isin(['CA1','CA3','DG','SUB'])]
hom=pd.read_csv(D+'data/ref/HOM_MouseHumanSequence.rpt',sep='\t',dtype=str)
mo=hom[hom['NCBI Taxon ID']=='10090'][['DB Class Key','Symbol']]; hu=hom[hom['NCBI Taxon ID']=='9606'][['DB Class Key','Symbol','HGNC ID','EntrezGene ID']]
om=mo.merge(hu,on='DB Class Key',suffixes=('_mouse','_human'))
nmo=mo.groupby('DB Class Key').size(); nhu=hu.groupby('DB Class Key').size()
om['ortholog_type']=[('1:1' if nmo[k]==1 and nhu[k]==1 else f'{nmo[k]}:{nhu[k]}') for k in om['DB Class Key']]
# synonym fallback: old symbols (e.g. Rnf165 -> Ark2c) via MGI MRK_List2.rpt (current Gene markers)
mk=pd.read_csv(D+'data/ref/MRK_List2.rpt',sep='\t',dtype=str,usecols=['Marker Symbol','Status','Marker Type','Marker Synonyms (pipe-separated)'])
mk=mk[(mk.Status=='O')&(mk['Marker Type']=='Gene')].dropna(subset=['Marker Synonyms (pipe-separated)'])
SYN={}
for cur,syn in zip(mk['Marker Symbol'],mk['Marker Synonyms (pipe-separated)']):
    for x in syn.split('|'): SYN.setdefault(x,cur)
O=om.groupby('Symbol_mouse').agg(human_symbol=('Symbol_human',lambda x:';'.join(sorted(set(x)))),HGNC_ID=('HGNC ID',lambda x:';'.join(sorted(set(x.dropna())))),human_Entrez=('EntrezGene ID',lambda x:';'.join(sorted(set(x.dropna())))),ortholog_type=('ortholog_type','first'))
nreg=pd.Series({g:sum(g in G[f'A1_{r}_hypo'] for r in ['CA1','CA3','DG']) for g in set().union(*[G[f'A1_{r}_hypo'] for r in ['CA1','CA3','DG']])})
regions=pd.Series({g:';'.join(r for r in ['CA1','CA3','DG'] if g in G[f'A1_{r}_hypo']) for g in nreg.index})
exc=spec[['CA1','CA3','DG','SUB']].sum(axis=1); top=spec[['CA1','CA3','DG','SUB']].idxmax(axis=1)
m6amods=set(hub.module)
mmo=mm[mm.module.isin(m6amods)].sort_values('kME_own',ascending=False)
best=mmo.groupby('gene').agg(modules=('module',lambda x:';'.join(x)),max_kME=('kME_own','max'))
mett=RG['METTL3cKO_HIP_down']|RG['METTL3KO_Cb_m6Aloss']; y1=RG['YTHDF1_CLIP_HIP']; a3=G['A3_ctx_Glut_hyper']|G['A3_ctx_Glut_hypo']
ae=R.sort_values('pvalue').groupby('gene').first()
cand=[g for g in nreg.index if nreg[g]>=2 and exc.get(g,0)>=0.5]
rows=[]
for g in cand:
    km=best.max_kME.get(g,np.nan); inmod=(km>=0.6) if not pd.isna(km) else False
    agep=R[(R.gene==g)]; sig=agep[agep.padj<0.1]
    e=dict(mouse_symbol=g,MeRIP_hypo_regions=regions[g],n_regions_hypo=int(nreg[g]),excitatory_specificity=round(exc[g],3),top_celltype=top[g],
           in_m6A_module=bool(inmod),modules=best.modules.get(g,''),max_kME=round(km,3) if not pd.isna(km) else np.nan,
           METTL3_target=g in mett,YTHDF1_CLIP=g in y1,Tegowski_ctx_Glut_DM=('hyper' if g in G['A3_ctx_Glut_hyper'] else 'hypo' if g in G['A3_ctx_Glut_hypo'] else ''),
           Allen_age_log2FC_bestType=round(ae.log2FoldChange.get(g,np.nan),3) if g in ae.index else np.nan,Allen_age_bestType=ae.celltype.get(g,''),Allen_age_padj_best=ae.padj.get(g,np.nan) if g in ae.index else np.nan,
           Allen_age_sig_types=';'.join(sig.celltype))
    e['evidence_score']=e['n_regions_hypo']+int(e['in_m6A_module'])+int(e['METTL3_target'] or e['YTHDF1_CLIP'])+int(e['Tegowski_ctx_Glut_DM']!='')+int(len(sig)>0)
    rows.append(e)
C=pd.DataFrame(rows).sort_values(['evidence_score','max_kME'],ascending=[False,False])
C['mgi_current_symbol']=[g if g in O.index else SYN.get(g,g) for g in C.mouse_symbol]; C=C.join(O,on='mgi_current_symbol')
C['tier']=np.select([C.evidence_score>=5,C.evidence_score>=4],['Tier1','Tier2'],'Tier3')
core=C[C.evidence_score>=4].head(50).copy()
core['predicted_direction_human_aging_expression']='down in aged hippocampus/cortex excitatory neurons (or bulk), relative to young'
core['predicted_direction_m6A_if_measurable']='lower m6A with age (hypomethylation), conditional evidence only'
print('candidates passing E1&E2:',len(C),' core (score>=4):',len(core),' tier1:',(core.tier=='Tier1').sum())
print(core[['mouse_symbol','human_symbol','ortholog_type','n_regions_hypo','excitatory_specificity','modules','max_kME','METTL3_target','YTHDF1_CLIP','Tegowski_ctx_Glut_DM','Allen_age_sig_types','evidence_score']].to_string())
# module gene sets
ms=[]
for m in m6amods:
    g=mm[mm.module==m].sort_values('kME_own',ascending=False); a=ag.loc[m]
    for _,r in g.iterrows(): ms.append(dict(module=m,celltype=m.split('-')[0],mouse_symbol=r.gene,kME=round(r.kME_own,3),module_age_beta_SD=round(a.beta_aged,3),module_age_p=a.p))
MS=pd.DataFrame(ms); MS['mgi_current_symbol']=[g if g in O.index else SYN.get(g,g) for g in MS.mouse_symbol]; MS=MS.join(O,on='mgi_current_symbol')
# program sets
ps=[]
for nm,s in [('METTL3cKO_HIP_down (GSE113796, Nex-CreERT2 CA1/CA3 basal, padj<0.05, log2FC<0)',RG['METTL3cKO_HIP_down']),('METTL3KO_Cb_m6Aloss (RM2Target GSE100528 MeRIP, top500 peaks lost)',RG['METTL3KO_Cb_m6Aloss']),
             ('A1_core87_hypo (PNAS Dataset S17, hypo in CA1&CA3&DG)',G['A1_CA1_hypo']&G['A1_CA3_hypo']&G['A1_DG_hypo'])]:
    for g in sorted(s): ps.append(dict(program=nm.split(' ')[0],definition=nm,mouse_symbol=g))
PS=pd.DataFrame(ps); PS['mgi_current_symbol']=[g if g in O.index else SYN.get(g,g) for g in PS.mouse_symbol]; PS=PS.join(O,on='mgi_current_symbol')
readme=pd.DataFrame({'item':[
 'Version / freeze','Purpose','E1 MeRIP evidence','E2 cell-type evidence','Module evidence','METTL3/YTHDF1 program evidence','scDART evidence','Allen age evidence','Evidence score / tiers','Core list','Ortholog mapping',
 'Pre-specified human hypotheses','H1','H2','H3','H4','What NOT to claim','Sheets'],
 'description':[
 'mouse_candidates_v1, frozen 2026-10-12 (phase 2b). No modification allowed after human data are inspected; any later change -> new version file with changelog.',
 'Candidate genes/sets from the mouse storyline: aging m6A hypomethylation in hippocampus attributable to excitatory neurons + weakening METTL3-dependent program; for human transcriptomic/genetic tests (performed by assistant 2).',
 'Hypomethylated with age (3 vs 16 mo) in >=2 of CA1/CA3/DG: GSE198526 / PNAS 2023 Dataset S17 (FC>=1.2, padj<=0.05, authors\' consistent-peak calls). Tissue-level MeRIP; not cell-specific.',
 'Allen Zeng-Aging-Mouse-10Xv3 HPF (unbiased libraries, 8 adult vs 15 aged donors): sum of specificity (mean CPM share) in CA1+CA3+DG+SUB glutamatergic types >=0.5.',
 'Signed WGCNA-like networks on metacells of Allen HPF CA1/DG/CA3/SUB neurons (scripts/24); modules enriched (BH q<0.01) for A1 hypo or METTL3 target sets; member with kME>=0.6 counts.',
 'METTL3cKO_HIP_down (GSE113796) or METTL3 KO cerebellum m6A-loss peaks (RM2Target/GSE100528) or YTHDF1 CLIP hippocampus (GSE106607).',
 'Tegowski 2024 Supp Table 10 cortex glutamatergic aging DM-RNA (either direction). Conditional/orthogonal only (2 reps/age; age confounded with kit).',
 'Allen excitatory pseudobulk DESeq2 aged-vs-adult padj<0.1 in any of CA1/CA3/DG/SUB.',
 'score = n_regions_hypo(2-3) + module + program + scDART + Allen_age (max 7). Tier1 >=5, Tier2 =4. Core = top 50 with score>=4, ranked by score then kME.',
 'See sheet core_candidates.',
 'MGI HOM_MouseHumanSequence.rpt (downloaded 2026-10-12, md5 89dc473ae85b11204fec5318bbf9ed81); outdated mouse symbols mapped to current MGI symbol via MRK_List2.rpt synonyms (md5 30e1088b64f45d2fd45c6befb3867e36). ortholog_type 1:1 recommended for tests; others flagged; blank = no MGI human ortholog.',
 'All directions stated as human-equivalent of mouse findings; one-sided tests allowed only in stated direction.',
 'H1 (expression, aging): core-list score (mean z) is LOWER in older vs younger human hippocampus/cortex, preferentially in excitatory-neuron pseudobulk (snRNA) rather than glia.',
 'H2 (cell-type specificity): core list and module sets are enriched among human excitatory-neuron-specific genes (EWCE-like with expression-matched background).',
 'H3 (program): METTL3cKO_HIP_down ortholog set shows lower activity with age in human excitatory neurons; direction of effect only, no m6A claim.',
 'H4 (genetics, exploratory): no directional prior for cognition GWAS/eQTL colocalisation; report as exploratory; mouse cognition association was NULL (GSE302188 within-old q>=0.39; GSE13799 non-replicating).',
 'Do not describe any human expression result as validating m6A changes. Mouse sc/bulk RNA supports expression/cell state only. Bulk MeRIP is tissue-level.',
 'README; core_candidates; all_E1E2_candidates; module_genesets; program_genesets']})
with pd.ExcelWriter(H+'mouse_candidates_v1.xlsx',engine='xlsxwriter') as w:
    readme.to_excel(w,sheet_name='README',index=False); core.to_excel(w,sheet_name='core_candidates',index=False); C.to_excel(w,sheet_name='all_E1E2_candidates',index=False); MS.to_excel(w,sheet_name='module_genesets',index=False); PS.to_excel(w,sheet_name='program_genesets',index=False)
core.to_csv(H+'mouse_candidates_v1_core.csv',index=False); MS.to_csv(H+'mouse_candidates_v1_modules.csv',index=False); PS.to_csv(H+'mouse_candidates_v1_programs.csv',index=False); readme.to_csv(H+'mouse_candidates_v1_README.csv',index=False)
