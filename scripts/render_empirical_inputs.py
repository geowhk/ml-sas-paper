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
if len(sys.argv)!=3: raise SystemExit('Usage: python scripts/render_empirical_inputs.py INPUT_GEOJSON NEW_OUTPUT_DIR')
input_geojson=Path(sys.argv[1]); A=Path(sys.argv[2]); A.mkdir(parents=True,exist_ok=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10,'pdf.fonttype':42})
def save(fig,name):
 fig.savefig(A/(name+'.png'),dpi=350,bbox_inches='tight');fig.savefig(A/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
# Empirical maps: five quantile classes for continuous variables; original units retained.
j=json.loads(input_geojson.read_text());features=j['features'];assert len(features)==337
fig,axes=plt.subplots(2,2,figsize=(7.3,8.8),layout='constrained')
config=[('population','(a) Population (persons)'),('dwelling_value','(b) Dwelling value (thousand euro)'),('employment_rate','(c) Employment rate (%)'),('urbanity','(d) Ordered urbanity')]
meta={}
for ax,(col,title) in zip(axes.ravel(),config):
 vals=[f['properties'][col] for f in features]
 if col=='urbanity':
  labs=['extremely urbanised','strongly urbanised','moderately urbanised','hardly urbanised','not urbanised'];cats=np.array([labs.index(x) for x in vals]);colors=plt.colormaps['viridis'](np.linspace(.05,.9,5));labels=['Extremely','Strongly','Moderately','Hardly','Not urbanised']
 else:
  edges=np.unique(np.quantile(vals,np.linspace(0,1,6)));cats=np.digitize(vals,edges[1:-1],right=True);colors=plt.colormaps['YlGnBu'](np.linspace(.15,.95,len(edges)-1));labels=[f'{a:,.0f}–{b:,.0f}' for a,b in zip(edges[:-1],edges[1:])];meta[col]={'breaks':edges.tolist(),'unit':title}
 for f,cat in zip(features,cats):
  geom=f['geometry'];parts=geom['coordinates'] if geom['type']=='MultiPolygon' else [geom['coordinates']]
  for rings in parts:
   paths=[]
   for ring in rings:
    vs=np.asarray(ring)[:,:2];codes=np.full(len(vs),MPath.LINETO);codes[0]=MPath.MOVETO;codes[-1]=MPath.CLOSEPOLY;paths.append(MPath(vs,codes))
   ax.add_patch(PathPatch(MPath.make_compound_path(*paths),facecolor=colors[cat],edgecolor='white',linewidth=.15))
 ax.autoscale();ax.set_aspect('equal');ax.axis('off');ax.set_title(title)
 ax.legend(handles=[Patch(facecolor=colors[i],label=l) for i,l in enumerate(labels)],loc='upper center',bbox_to_anchor=(.5,-.015),fontsize=7,frameon=False,ncol=2)
(A/'empirical_map_breaks.json').write_text(json.dumps(meta,indent=2))
save(fig,'Figure_08_empirical_inputs')
