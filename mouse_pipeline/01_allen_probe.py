import fsspec,h5py,time,numpy as np,pandas as pd
t=time.time()
url='https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/expression_matrices/Zeng-Aging-Mouse-10Xv3/20241130/Zeng-Aging-Mouse-10Xv3-raw.h5ad'
f=fsspec.open(url,'rb',block_size=2**23).open()
h=h5py.File(f,'r')
o=h['obs']
for k in o: print(k, type(o[k]), dict(o[k].attrs), flush=True)
def rd(node):
    if isinstance(node,h5py.Group):
        c=node['categories'][:].astype(str); return c[node['codes'][:]]
    return node[:].astype(str)
cl=rd(o['cell_label']); print(len(cl), cl[:3], time.time()-t, flush=True)
np.save('./data/allen/h5ad_cell_label.npy', cl)
ip=h['X']['indptr'][:]; np.save('./data/allen/h5ad_indptr.npy',ip); print('indptr',time.time()-t,flush=True)
v=h['var']; print({k:dict(v[k].attrs) for k in v})
g=pd.DataFrame({'gene_identifier':rd(v['gene_identifier']),'gene_symbol':rd(v['gene_symbol'])}); g.to_csv('./data/allen/genes.csv',index=False); print(g.head(), time.time()-t)
