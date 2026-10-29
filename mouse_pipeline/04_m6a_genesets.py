"""Build m6A gene sets (symbols, mouse) from authors' processed files. No re-alignment.
Outputs analysis/tables/genesets_m6a.csv (set, gene) and analysis/tables/merip_backgrounds.csv"""
import pandas as pd,numpy as np,gzip
O='analysis/tables/'
genes=pd.read_csv('data/allen/genes.csv'); e2s=dict(zip(genes.gene_identifier,genes.gene_symbol))
S=[]
# A1 PNAS Dataset S17 (FC>=1.2, padj<=0.05, consistent peaks; authors' curated)
a1=pd.read_csv(O+'A1_PNAS_DatasetS17_diffmethyl_genes_long.csv')
for (r,d),g in a1.groupby(['region','direction']): S+= [(f'A1_{r}_{d}',x) for x in g.gene.unique()]
# backgrounds: genes with any MeTPeak peak in young or old (GEO BEDs, ENSEMBL ids)
bg=[]
for r in ['CA1','CA3','DG','ACC']:
    ids=set()
    for t in ['young','old']:
        d=pd.read_csv(f'data/A1_GSE198526/GSE198526_{r}{t}_peak.bed.gz',sep='\t'); ids|=set(d.iloc[:,3])
    d=pd.read_csv(f'data/A1_GSE198526/GSE198526_{r}_sig_diff_peak.bed.gz',sep='\t'); ids|=set(d.iloc[:,3])
    sy={e2s.get(i) for i in ids}-{None}
    bg+=[(f'A1_{r}',x) for x in sy]
# A2 Huang Aging Cell 2023 whole hippocampus 3 vs 20 mo; m6alogFC = aged vs young (checked: neg rows have higher Young IP counts)
a2=pd.read_excel('data/A2_GSE217884/GSE217882_Ageing_MeRIPseq.xlsx')
a2=a2[a2.m6aFDR<0.05]
hyp=set(a2[a2.m6alogFC>0].genename.dropna()); hpo=set(a2[a2.m6alogFC<0].genename.dropna())
S+=[('A2_HIP_hyper',x) for x in hyp-hpo]+[('A2_HIP_hypo',x) for x in hpo-hyp]
# A3 Tegowski ST10 cortex cell-type DM-RNAs (aged vs young), conditional discovery only
a3=pd.read_csv(O+'A3_Tegowski_ST10_celltype_aging_DMRNAs.csv')
for (c,d),g in a3.groupby(['celltype','direction']):
    if len(g)>=10: S+=[(f'A3_ctx_{c}_{d.replace("_in_aged","")}',x) for x in g.rna.unique()]
# A3 site-level GEO file (2,055 sites, glutamatergic neurons) -> genes
s=pd.read_csv('data/A3_GSE268159/GSE268159_Ageddiffmethylsites.txt.gz',sep='\t')
S+=[('A3_ctx_GlutSites_hyper',x) for x in s[s.YoungL2FC<0].V14.unique()]+[('A3_ctx_GlutSites_hypo',x) for x in s[s.YoungL2FC>0].V14.unique()]
gs=pd.DataFrame(S,columns=['set','gene']).drop_duplicates()
gs.to_csv(O+'genesets_m6a_aging.csv',index=False)
pd.DataFrame(bg,columns=['region','gene']).to_csv(O+'merip_backgrounds.csv',index=False)
print(gs.groupby('set').size().to_string()); print(pd.DataFrame(bg,columns=['region','gene']).groupby('region').size())
