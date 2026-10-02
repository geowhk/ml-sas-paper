from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch,PathPatch
from matplotlib.path import Path as MPath
import sys
ROOT=Path(__file__).resolve().parents[1]/'reference/figure_inputs'
E=ROOT/'dgp'
if len(sys.argv) not in (2,3): raise SystemExit('Usage: python scripts/render_dgp_figures.py NEW_OUTPUT_DIR')
if len(sys.argv)==3: ROOT=Path(sys.argv[2]);E=ROOT/'dgp'
A=Path(sys.argv[1]); A.mkdir(parents=True,exist_ok=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'pdf.fonttype':42})
v=pd.read_csv(E/'hex_vertices.csv');polys=[g[['X','Y']].to_numpy() for _,g in v.groupby('id')]
def boundary(t):
 edges={}
 for p in np.array(polys)[t]:
  for a,b in zip(p[:-1],p[1:]):
   key=tuple(sorted([tuple(np.round(a,9)),tuple(np.round(b,9))]));edges[key]=edges.get(key,0)+1
 return [list(k) for k,n in edges.items() if n==1]
def panel(fig,ax,values,title,truth=None,labels=None,colors=None):
 if labels is not None:
  cmap=ListedColormap(colors);norm=BoundaryNorm(np.arange(len(labels)+1)-.5,len(labels))
 else:
  lim=max(abs(np.asarray(values).min()),abs(np.asarray(values).max()));cmap='RdBu_r';norm=plt.Normalize(-lim,lim)
 p=PolyCollection(polys,array=np.array(values),cmap=cmap,norm=norm,edgecolors='white',linewidths=.12);ax.add_collection(p);ax.autoscale();ax.set_aspect('equal');ax.axis('off');ax.set_title(title)
 if truth is not None:ax.add_collection(LineCollection(boundary(truth),colors='#333333',linewidths=.75))
 if labels is None:fig.colorbar(p,ax=ax,orientation='horizontal',fraction=.04,pad=.025)
 else:ax.legend(handles=[Patch(facecolor=colors[i],label=l) for i,l in enumerate(labels)],loc='lower center',bbox_to_anchor=(.5,-.15),ncol=2 if len(labels)>2 else len(labels),fontsize=7,frameon=False)
def save(fig,name):
 fig.savefig(A/(name+'.png'),dpi=350,bbox_inches='tight');fig.savefig(A/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
prod=ROOT/'02_revision-simulation/runs/revision_v2/production'
f=pd.read_csv(E/'fixed_dgp1_attributes.csv');q=pd.read_csv(ROOT/'02_revision-simulation/runs/revision_v2/corrections/dgp1_local_geary_greater_perm99999_20260719/tables/dgp1_univariate_lisa_results_corrected.csv');m=pd.read_csv(prod/'dgp1_section4_v2_ntree100_perm99999_20260718_133406/tables/dgp1_univariate_mlsas_local_results.csv')
fig,axes=plt.subplots(3,2,figsize=(7.1,8.6),layout='constrained');ax=axes.ravel()
panel(fig,ax[0],f.x,r'(a) Input $X$',f.true_cluster)
panel(fig,ax[1],f.cluster_type.map({'background':0,'HH':1,'LL':2,'MM':3}),'(b) Planted regions',f.true_cluster,['Background','High','Low','Moderate'],['#eeeeee','#d55e00','#0072b2','#009e73'])
flags=[q.moran_flag,q.geary_flag,q.gi_flag,m.cluster_flag_strict];counts=[]
for i,(fl,title) in enumerate(zip(flags,[r"Local Moran’s $I$",r"Local Geary’s $c$",r"Getis–Ord $G_i^*$","Path ML-SAS"])):
 t=fl!='Not significant';counts.append(int(t.sum()));panel(fig,ax[i+2],t.astype(int),f'({chr(99+i)}) {title}',f.true_cluster,['Not significant','Significant'],['#eeeeee','#d55e00'])
assert counts==[61,83,66,69],counts
save(fig,'Figure_05_dgp1_inputs_results')
for d,tag,dirname,methodfiles in [(3,'dgp3.1','dgp3.1_section4_v2_ntree100_perm99999_20260718_152439',['continuous_multigeary','continuous_mlsas']),(4,'dgp4.3','dgp4.3_section4_v2_ntree100_perm99999_20260718_155829',['continuous_multigeary','continuous_mlsas','mixed_mlsas'])]:
 f=pd.read_csv(E/f'fixed_dgp{d}_attributes.csv');fig,axes=plt.subplots(2,3,figsize=(8.0,5.7),layout='constrained');fig.set_constrained_layout_pads(h_pad=0.10,w_pad=0.04,hspace=0.10,wspace=0.03);ax=axes.ravel()
 for i,col in enumerate(['x1','x2','x3' if d==3 else 'cat1']):
  if col=='cat1':panel(fig,ax[i],f[col].map({'A':0,'B':1,'C':2}),f'({chr(97+i)}) CAT1',f.true_cluster,['A','B','C'],['#e69f00','#56b4e9','#009e73'])
  else:panel(fig,ax[i],f[col],f'({chr(97+i)}) $X_{col[-1]}$',f.true_cluster)
 counts=[]
 for i,meth in enumerate(methodfiles):
  g=pd.read_csv(prod/dirname/'tables'/f'{tag}_{meth}_local_results.csv');t=(g.cluster_flag if 'multigeary' in meth else g.cluster_flag_strict)!='Not significant';counts.append(int(t.sum()))
  title='Multivariate Local Geary' if 'multigeary' in meth else ('ML-SAS continuous' if meth.startswith('continuous') else 'ML-SAS mixed')
  panel(fig,ax[i+3],t.astype(int),f'({chr(100+i)}) {title}',f.true_cluster,['Not significant','Significant'],['#eeeeee','#d55e00'])
 if d==3:panel(fig,ax[5],f.true_cluster.astype(int),'(f) Planted region',None,['Background','Target'],['#eeeeee','#d55e00'])
 assert counts==({3:[33,31],4:[42,19,46]}[d]),counts
 save(fig,f'Figure_0{d+3}_dgp{d}_inputs_results')

f=pd.read_csv(E/'fixed_dgp2_attributes.csv');fig,axes=plt.subplots(2,2,figsize=(7.1,6.4),layout='constrained')
for i,(ax,col) in enumerate(zip(axes.ravel(),['x1','x2','x3','cat1'])):
 if col=='cat1':panel(fig,ax,f[col].map({'C1':0,'C2':1,'C3':2}),f'({chr(97+i)}) CAT1',None,['C1','C2','C3'],['#e69f00','#56b4e9','#009e73'])
 else:panel(fig,ax,f[col],f'({chr(97+i)}) $X_{col[-1]}$')
save(fig,'Figure_S1_dgp2_inputs')
