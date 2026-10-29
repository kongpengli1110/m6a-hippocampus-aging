import os
import pandas as pd,numpy as np
D=os.environ.get('PROJECT_DIR','./')
MAJOR={'016 CA1-ProS Glut':'CA1','017 CA3 Glut':'CA3','037 DG Glut':'DG','023 SUB-ProS Glut':'SUB','036 HPF CR Glut':'CR',
 '319 Astro-TE NN':'Astro','327 Oligo NN':'Oligo','326 OPC NN':'OPC','334 Microglia NN':'Microglia','333 Endo NN':'Endo',
 '331 Peri NN':'Mural','332 SMC NN':'Mural','330 VLMC NN':'VLMC','335 BAM NN':'BAM','038 DG-PIR Ex IMN':'DG_IMN'}
def major(s):
    if s in MAJOR: return MAJOR[s]
    if s.endswith('Gaba'): return 'GABA'
    if 'ENT' in s and 'Glut' in s: return 'ENT_Glut'
    if s.endswith('Glut'): return 'OtherGlut'
    return 'Other'
REGS={'Writer':['Mettl3','Mettl14','Wtap','Virma','Rbm15','Rbm15b','Zc3h13','Cbll1','Mettl16','Mettl5','Zcchc4'],
 'Eraser':['Fto','Alkbh5'],
 'Reader':['Ythdf1','Ythdf2','Ythdf3','Ythdc1','Ythdc2','Igf2bp1','Igf2bp2','Igf2bp3','Hnrnpa2b1','Hnrnpc','Rbmx','Prrc2a','Fmr1','Elavl1','Eif3a']}
REG2CLS={g:k for k,v in REGS.items() for g in v}
def load_pb(regions=('HIP','PAR','ENT')):
    """donor x sampling x subclass pseudobulk -> unbiased only, aggregated to donor x major type"""
    genes=pd.read_csv(D+'data/allen/genes.csv')
    Xs=[];Ms=[]
    for r in regions:
        X=pd.read_parquet(D+f'data/allen/pseudobulk_{r}_donor_sampling_subclass.parquet'); M=pd.read_csv(D+f'data/allen/pseudobulk_{r}_meta.csv',index_col=0)
        M['region']=r; keep=M.population_sampling=='unbiased'; X=X[keep.values]; M=M[keep]
        M['major']=M.subclass.map(major); M['key']=M.donor_label+'|'+M.major
        Xg=X.groupby(M.key.values).sum(); Mg=M.groupby('key').agg(donor=('donor_label','first'),major=('major','first'),region=('region','first'),sex=('donor_sex','first'),age=('age','first'),n_cells=('n_cells','sum'),donor_line=('donor_label',lambda s:s.iloc[0].split('-')[0]))
        Xs.append(Xg);Ms.append(Mg.loc[Xg.index])
    X=pd.concat(Xs);M=pd.concat(Ms)
    return X,M,genes
