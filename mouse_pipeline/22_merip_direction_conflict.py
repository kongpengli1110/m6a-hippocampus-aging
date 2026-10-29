"""Task 2: why GSE198526 (A1; mostly hypo) and GSE217882 (A2; mostly hyper) disagree. Uses only authors' processed files + ENA run metadata."""
import pandas as pd,numpy as np,pyranges as pr,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
D='./'; T=D+'analysis/tables/'; F=D+'analysis/figures/'
genes=pd.read_csv(D+'data/allen/genes.csv'); e2s=dict(zip(genes.gene_identifier,genes.gene_symbol))
# ---- 1. design / depth table
e1=pd.read_csv(D+'data/A_sra/PRJNA815761_ena.tsv',sep='\t'); e2=pd.read_csv(D+'data/A_sra/PRJNA901175_ena.tsv',sep='\t')
e1['region']=e1.sample_title.str.split('_').str[0]; e1['age']=e1.sample_title.str.extract(r'_(young|old)')[0]; e1['lib']=np.where(e1.sample_title.str.contains('meRIP'),'IP','input')
e2['region']='wholeHIP'; e2['age']=np.where(e2.sample_title.str.startswith('Young'),'young','old'); e2['lib']=np.where(e2.sample_title.str.contains('IP'),'IP','input')
E=pd.concat([e1,e2]); E['M_reads']=E.read_count/1e6
dep=E.groupby(['region','age','lib']).M_reads.mean().unstack(['age','lib']).round(1)
dep['IP_old/young']=(dep[('old','IP')]/dep[('young','IP')]).round(2); dep['input_old/young']=(dep[('old','input')]/dep[('young','input')]).round(2)
dep['IP/input_young']=(dep[('young','IP')]/dep[('young','input')]).round(2); dep['IP/input_old']=(dep[('old','IP')]/dep[('old','input')]).round(2)
# peaks per age
pk={}
for r in ['CA1','CA3','DG','ACC']:
    pk[r]={a:len(pd.read_csv(D+f'data/A1_GSE198526/GSE198526_{r}{a}_peak.bed.gz',sep='\t')) for a in ['young','old']}
pk['wholeHIP']={'young':8084,'old':7527}   # Huang 2023 text
dep['peaks_young']=[pk[r]['young'] for r in dep.index]; dep['peaks_old']=[pk[r]['old'] for r in dep.index]; dep['peaks_old/young']=(dep.peaks_old/dep.peaks_young).round(2)
s17=pd.read_csv(T+'A1_PNAS_DatasetS17_diffmethyl_genes_long.csv'); c=s17.groupby(['region','direction']).size().unstack()
pct={r:c.loc[r,'hypo']/c.loc[r].sum()*100 for r in c.index}; pct['wholeHIP']=102/528*100
dep['pct_hypo_genes']=[round(pct[r],1) for r in dep.index]
dep.columns=['_'.join(x) if isinstance(x,tuple) and x[1] else (x[0] if isinstance(x,tuple) else x) for x in dep.columns]
dep.to_csv(T+'2b_conflict_depth_peaks_direction.csv'); print(dep.to_string())
rho=stats.spearmanr(dep['peaks_old/young'],dep.pct_hypo_genes); rho2=stats.spearmanr(dep['IP_old/young'],dep.pct_hypo_genes)
print('spearman peaks old/young vs %hypo',rho,' IP depth ratio vs %hypo',rho2)
# ---- 2. design differences (from methods text; recorded verbatim-summarised)
design=pd.DataFrame([
 ['Age contrast','3 vs 16 mo (C57BL/6, sex not stated in GEO)','3 vs 20 mo (male C57BL/6)'],
 ['Tissue','microdissected CA1 / CA3 / DG (+ACC); 10 mice/age pooled into 3 replicates','whole hippocampus, 3 individual mice/age'],
 ['RNA input','rRNA-depleted total RNA (NEBNext rRNA removal), 10 ug/IP, fragmented ~80 nt','poly(A)-selected RNA (oligo-dT x2), ~1.4 ug, fragmented ~100 nt'],
 ['Antibody','anti-m6A, Synaptic Systems (10 ug)','anti-m6A, Synaptic Systems 202003'],
 ['Library','TruSeq Stranded Total RNA; HiSeq 2000 SE50','SMARTer Stranded Total RNA-Seq Pico v2; HiSeq 4000 PE150'],
 ['Peak calling','MeTPeak 1.0.0 (STAR 2.5.2b, mm10)','not specified in text; 200-nt windows (fixed-width bins) in Table S1'],
 ['Differential test','ExomePeak 2.16 (IP vs input ratio), FC>=1.2, padj<=0.05, consistent peaks only','edgeR-like on IP window counts (columns m6alogFC/LR/FDR); FC>1.5, FDR<0.05'],
 ['Reported result','hypo >> hyper (70-94% hypo genes)','hyper >> hypo (426 vs 102 genes)']],columns=['aspect','A1_GSE198526_PNAS2023','A2_GSE217882_AgingCell2023'])
design.to_csv(T+'2b_conflict_design_comparison.csv',index=False)
# ---- 3. A2: is the reported m6A logFC input-normalised? compare with IP-only logFC and RNA logFC
a2=pd.read_excel(D+'data/A2_GSE217884/GSE217882_Ageing_MeRIPseq.xlsx'); ip=a2[['IP_Y1','IP_Y2','IP_Y3','IP_A1','IP_A2','IP_A3']]
ipn=ip/ip.sum()*1e6; a2['ipOnly_lfc']=np.log2(ipn[['IP_A1','IP_A2','IP_A3']].mean(axis=1)+1)-np.log2(ipn[['IP_Y1','IP_Y2','IP_Y3']].mean(axis=1)+1)
rna=pd.read_excel(D+'data/A2_GSE217884/GSE217883_Ageing_RNAseq.xlsx').drop_duplicates('genename').set_index('genename')
a2['rna_lfc']=a2.genename.map(rna.rnalogFC)
r_ip=np.corrcoef(a2.ipOnly_lfc,a2.m6alogFC)[0,1]; sr=stats.spearmanr(a2.m6alogFC,a2.rna_lfc,nan_policy='omit')
a2['m6a_minus_rna']=a2.m6alogFC-a2.rna_lfc
print('A2 corr(reported, IP-only lfc)=%.3f ; spearman(m6A lfc, RNA lfc)=%.2f p=%.1e'%(r_ip,sr.statistic,sr.pvalue))
print('A2 hyper fraction raw',(a2.m6alogFC>0).mean().round(3),' after subtracting RNA lfc',(a2.m6a_minus_rna.dropna()>0).mean().round(3))
# ---- 4. threshold dependence
thr=[]
for f in [0.05,0.01,0.001]:
    for l in [np.log2(1.5),1,1.5]:
        s=a2[(a2.m6aFDR<f)&(a2.m6alogFC.abs()>=l)]; thr.append(dict(dataset='A2',level='peak',FDR=f,min_abs_log2FC=round(l,2),n=len(s),frac_hyper=(s.m6alogFC>0).mean()))
# region annotation dependence (A2)
reg=a2.groupby('type').m6alogFC.agg(n='size',frac_hyper=lambda x:(x>0).mean()).reset_index(); reg.to_csv(T+'2b_conflict_A2_by_peak_region.csv',index=False); print(reg)
# A1: peak-level direction inferred from presence in young-only vs old-only MeTPeak peaks, across p thresholds of the GEO diff BED
def rd(f):
    d=pd.read_csv(f,sep='\t'); d.columns=['Chromosome','Start','End','gene','score','Strand']+list(d.columns[6:]); return d
A1dir=[]
for r in ['CA1','CA3','DG','ACC']:
    df=rd(D+f'data/A1_GSE198526/GSE198526_{r}_sig_diff_peak.bed.gz'); df['pid']=range(len(df))
    Y=pr.PyRanges(rd(D+f'data/A1_GSE198526/GSE198526_{r}young_peak.bed.gz')[['Chromosome','Start','End']]); O=pr.PyRanges(rd(D+f'data/A1_GSE198526/GSE198526_{r}old_peak.bed.gz')[['Chromosome','Start','End']])
    P=pr.PyRanges(df[['Chromosome','Start','End','pid']]); iy=set(P.overlap(Y).df.pid); io=set(P.overlap(O).df.pid)
    df['dir']=np.select([df.pid.isin(iy)&~df.pid.isin(io),df.pid.isin(io)&~df.pid.isin(iy)],['young_only(hypo)','old_only(hyper)'],'both/none')
    for f in [0.05,0.01,1e-3,1e-5]:
        s=df[df.score<=f]; n_h=(s.dir=='young_only(hypo)').sum(); n_y=(s.dir=='old_only(hyper)').sum()
        thr.append(dict(dataset=f'A1_{r}',level='peak (presence-inferred)',FDR=f,min_abs_log2FC=np.nan,n=n_h+n_y,frac_hyper=n_y/max(n_h+n_y,1)))
    A1dir.append(df.assign(region=r))
TH=pd.DataFrame(thr); TH.to_csv(T+'2b_conflict_threshold_dependence.csv',index=False); print(TH.round(3).to_string())
# ---- 5. peak-level coordinate overlap: A2 DM peaks vs A1 DM peaks (any HIP region)
A1=pd.concat(A1dir); A1=A1[A1.region!='ACC']
a2g=a2.assign(Chromosome='chr'+a2.Chr.astype(str),pid2=range(len(a2)))
ov=pr.PyRanges(a2g[['Chromosome','Start','End','pid2']]).join(pr.PyRanges(A1[['Chromosome','Start','End','dir','region']])).df
ov=ov.merge(a2g[['pid2','m6alogFC','genename']],on='pid2')
ct=pd.crosstab(np.where(ov.m6alogFC>0,'A2 hyper','A2 hypo'),[ov.region,ov.dir]); ct.to_csv(T+'2b_conflict_peak_overlap_A2_vs_A1.csv'); print(ct)
print('A2 peaks overlapping any A1 HIP DM peak:',ov.pid2.nunique(),'/',len(a2))
# ---- 6. gene-level overlap
G={k:set(v) for k,v in s17.groupby(['region','direction']).gene}
a2h=set(a2[a2.m6alogFC>0].genename)-set(a2[a2.m6alogFC<0].genename); a2o=set(a2[a2.m6alogFC<0].genename)-a2h
gl=[]
for r in ['CA1','CA3','DG']:
    for dname in ['hypo','hyper']:
        gl.append(dict(A1=f'{r}_{dname}',n=len(G[(r,dname)]),with_A2_hyper=len(G[(r,dname)]&a2h),with_A2_hypo=len(G[(r,dname)]&a2o)))
GL=pd.DataFrame(gl); GL.to_csv(T+'2b_conflict_gene_overlap.csv',index=False); print(GL)
# ---- 7. cell-type of A2 hyper genes vs A1 hypo (Allen specificity): neuron vs glia share
spec=pd.read_csv(T+'allen_HPF_specificity_EWCE.csv',index_col=0)
neu=['CA1','CA3','DG','SUB','ENT_Glut','GABA']; gli=['Astro','Oligo','OPC','Microglia']
def share(gs): s=spec.reindex(list(gs)).dropna(); return pd.Series({'neuronal':s[neu].sum(axis=1).mean(),'glial':s[gli].sum(axis=1).mean(),'n':len(s)})
SH=pd.DataFrame({'A1_CA1_hypo':share(G[('CA1','hypo')]),'A1_DG_hypo':share(G[('DG','hypo')]),'A1_DG_hyper':share(G[('DG','hyper')]),'A2_hyper':share(a2h),'A2_hypo':share(a2o)}).T
SH.to_csv(T+'2b_conflict_celltype_share.csv'); print(SH.round(3))
# ---- figure
fig,ax=plt.subplots(2,3,figsize=(16,9))
x=dep.index; ax[0,0].scatter(dep['peaks_old/young'],dep.pct_hypo_genes,s=60)
for i,r in enumerate(x): ax[0,0].annotate(r,(dep['peaks_old/young'].iloc[i],dep.pct_hypo_genes.iloc[i]))
ax[0,0].set_xlabel('peaks called old / young'); ax[0,0].set_ylabel('% hypo among DM genes'); ax[0,0].set_title(f'Peak-calling asymmetry vs direction\n(Spearman rho={rho.statistic:.2f}, n=5 contrasts)')
w=0.35; xi=np.arange(len(x)); ax[0,1].bar(xi-w/2,dep['IP_old/young'],w,label='IP reads old/young'); ax[0,1].bar(xi+w/2,dep['input_old/young'],w,label='input reads old/young'); ax[0,1].axhline(1,color='k',lw=.5)
ax[0,1].set_xticks(xi); ax[0,1].set_xticklabels(x); ax[0,1].legend(fontsize=8); ax[0,1].set_title('Sequencing depth asymmetry (ENA run metadata)')
ax[0,2].scatter(a2.rna_lfc,a2.m6alogFC,s=6,alpha=.5); ax[0,2].axhline(0,color='k',lw=.5); ax[0,2].axvline(0,color='k',lw=.5)
ax[0,2].set_xlabel('A2 RNA-seq log2FC (aged/young)'); ax[0,2].set_ylabel('A2 reported m6A log2FC'); ax[0,2].set_title(f'A2: m6A change tracks RNA (rho={sr.statistic:.2f})\nreported FC ~ IP-only FC (r={r_ip:.3f})')
th=TH.copy(); th['lab']=th.dataset+' FDR'+th.FDR.astype(str)+np.where(th.min_abs_log2FC.notna(),' |lfc|>'+th.min_abs_log2FC.round(2).astype(str),'')
for ds,g in th.groupby('dataset'):
    ax[1,0].plot(range(len(g)),g.frac_hyper,marker='o',label=ds)
ax[1,0].axhline(.5,color='grey',ls='--'); ax[1,0].set_ylabel('fraction hyper (peaks)'); ax[1,0].set_xlabel('increasingly stringent threshold'); ax[1,0].legend(fontsize=7); ax[1,0].set_title('Direction vs threshold (A2 stable; A1 CA3/DG flip to hyper at strict p)')
GLp=GL.set_index('A1')[['with_A2_hyper','with_A2_hypo']]; GLp.plot.bar(ax=ax[1,1],color=['#d1495b','#3b75af']); ax[1,1].set_title('Gene-level overlap A1 vs A2'); ax[1,1].tick_params(axis='x',rotation=30)
SH[['neuronal','glial']].plot.bar(ax=ax[1,2],color=['#2b8cbe','#e6550d']); ax[1,2].set_title('Mean Allen-HPF specificity share'); ax[1,2].tick_params(axis='x',rotation=30)
plt.tight_layout(); plt.savefig(F+'P2b_Fig3_MeRIP_direction_conflict.png',dpi=300)
