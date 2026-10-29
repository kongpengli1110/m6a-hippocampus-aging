"""Verify sample-ID correspondence between GEO GSE302188 brain counts and Zenodo 15849341 VST matrix (same animals?)
by expression-profile correlation. Output: analysis/tables/B11_sample_matching.csv"""
import pandas as pd,numpy as np
g=pd.read_csv('data/B11_GSE302188/GSE302188_brain_RNA_gene_counts.csv.gz',index_col=0)
z=pd.read_csv('data/B11_GSE302188/zenodo/varianceStabilizedRNASeq_Brain.tsv',sep='\t',index_col=0)
print('GEO',g.shape,'Zenodo VST',z.shape)
common=g.index.intersection(z.index)
gl=np.log2(g.loc[common]/g.sum()*1e6+1)
keep=(gl.mean(1)>3)
gl=gl[keep]; zz=z.loc[gl.index]
# rank correlation on centered values (remove gene means) to emphasise sample-specific signal
gc=gl.sub(gl.mean(1),axis=0); zc=zz.sub(zz.mean(1),axis=0)
C=pd.DataFrame(np.corrcoef(gc.T.values,zc.T.values)[:gc.shape[1],gc.shape[1]:],index=gc.columns,columns=zc.columns)
best=C.idxmax(1); bestv=C.max(1); second=C.apply(lambda r:r.nlargest(2).iloc[1],axis=1)
rev=C.idxmax(0)
res=pd.DataFrame({'geo_id':C.index,'zenodo_best':best.values,'r_best':bestv.values,'r_second':second.values})
res['reciprocal']=[rev[b]==a for a,b in zip(res.geo_id,res.zenodo_best)]
res['geo_age']=res.geo_id.str.split('.').str[0]; res['zen_age']=res.zenodo_best.str.split('_').str[0]
res.to_csv('analysis/tables/B11_sample_matching.csv',index=False)
print(res.sort_values('r_best').head(10).to_string())
print('reciprocal',res.reciprocal.sum(),'/',len(res),'unique zen',res.zenodo_best.nunique())
print(pd.crosstab(res.geo_age,res.zen_age))
# nominal-name agreement after renaming
def ren(x):
    a,n=x.split('.N'); return {'3M':'2M','28M':'30M'}.get(a,a)+'_'+n
res['nominal']=res.geo_id.map(ren); print('nominal==best',(res.nominal==res.zenodo_best).sum())
print(res[res.nominal!=res.zenodo_best][['geo_id','nominal','zenodo_best','r_best','r_second']].to_string())
