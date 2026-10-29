"""Regulator-dependent target gene sets (mouse brain only). Sources: GSE113796 authors' DESeq2 (Nex-CreERT2 Mettl3/Fto cKO, CA1/CA3, basal),
GSE106607 Ythdf1 high-confidence CLIP peaks (hippocampus), RM2Target v2 mouse m6A perturbation tables (brain entries)."""
import pandas as pd,numpy as np
D='./'; T=D+'analysis/tables/'
S=[]
a=pd.read_csv(D+'data/A9/analysis.csv',index_col=0)
for k,nm in [('Mettl3_genotypeeffect_in_basal','METTL3cKO_HIP'),('Fto_genotypeeffect_in_basal','FTOcKO_HIP'),('Mettl3_genotypeeffect_in_fearcondtioned','METTL3cKO_HIP_FC'),('Fto_genotypeeffect_in_fearcondtioned','FTOcKO_HIP_FC')]:
    p=a['padj_'+k]; l=a['log2FoldChange_'+k]
    S+=[(f'{nm}_up',g,'GSE113796') for g in a.Gene[(p<0.05)&(l>0)]]+[(f'{nm}_down',g,'GSE113796') for g in a.Gene[(p<0.05)&(l<0)]]
c=pd.read_excel(D+'data/A9/GSE106607_A_Ythdf1_CLIP_seq_peak_list.xlsx',sheet_name='High-confidence Y1-CLIP peaks')
S+=[('YTHDF1_CLIP_HIP',g,'GSE106607') for g in c.Symbol.dropna().unique()]
e=pd.read_csv(D+'data/rm2target/Perturbation_Expr_Mus_musculus_m6A.csv',usecols=['wer','target_gene','cell_line','gse_id','padj','log2fc'],low_memory=False)
m=pd.read_csv(D+'data/rm2target/Perturbation_Mod-RIP_mm10_Mus_musculus_m6A.csv',usecols=['wer','target_gene','cell_line','gse_id','padj','log2fc'])
def top(df,n=500):
    df=df[df.padj<0.05].sort_values('padj'); return df.drop_duplicates('target_gene').head(n).target_gene
sel=[('expr',e,'Mettl3','Cortex','GSE154992','METTL3KO_Cortex'),('expr',e,'Mettl14','Oligodendrocytes','GSE124243','METTL14mut_Oligo'),
     ('expr',e,'Fmr1','Hippocampus','GSE94559','FMR1KO_HIP'),('expr',e,'Ythdf2','Cerebellar-GCs','GSE153688','YTHDF2KD_CGN'),
     ('mod',m,'Mettl3','Cerebellum','GSE100528','METTL3KO_Cb_m6Aloss'),('mod',m,'Virma','Forebrain','GSE279609','VIRMAKO_FB_m6Aloss')]
for kind,df,w,cl,g,nm in sel:
    x=df[(df.wer==w)&(df.cell_line==cl)&(df.gse_id==g)]
    if kind=='mod': x=x[x.log2fc<0]; S+=[(nm,t,f'RM2Target:{g}') for t in top(x)]
    else:
        S+=[(nm+'_up',t,f'RM2Target:{g}') for t in top(x[x.log2fc>0],300)]+[(nm+'_down',t,f'RM2Target:{g}') for t in top(x[x.log2fc<0],300)]
R=pd.DataFrame(S,columns=['set','gene','source']).drop_duplicates(['set','gene'])
R.to_csv(T+'genesets_regulator_targets.csv',index=False); print(R.groupby(['set','source']).size().to_string())
