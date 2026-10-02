"""Appendix C diagnostics from local results; v6 F1, score orientation applied once."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
if len(sys.argv)!=3: raise SystemExit('Usage: analyze_appendix.py LOCAL_EVIDENCE_CSV NEW_OUTPUT_DIR')
P=Path(sys.argv[2]);P.mkdir(parents=True,exist_ok=False)
a=pd.read_csv(sys.argv[1],low_memory=False)
assert len(a)==240000 and not a.duplicated(['dgp','replication_id','method','id']).any()
def auc(s,t):
    s=np.asarray(s);t=np.asarray(t,dtype=bool)
    if not t.any() or t.all():return np.nan
    ranks=pd.Series(s).rank().to_numpy(); n=t.sum()
    return (ranks[t].sum()-n*(n+1)/2)/(n*(len(t)-n))
def bh(p):
    p=np.asarray(p);o=np.argsort(p);v=p[o]*len(p)/np.arange(1,len(p)+1)
    q=np.empty(len(p));q[o]=np.minimum(1,np.minimum.accumulate(v[::-1])[::-1]);return q
assert auc([1,2,3,4],[0,0,1,1])==1 and auc([1,1],[0,1])==.5
rows=[];strata=[];pairauc=[]
for (d,r,m),g in a.groupby(['dgp','replication_id','method']):
    t=g.truth.to_numpy();det=g.detected.to_numpy();s=g.statistic.to_numpy() * (-1 if m=='local_geary' else 1)
    assert len(g)==400 and t.sum()==60
    assert np.allclose(bh(g.p_raw),g.p_adjusted,atol=1e-12)
    if m!='local_geary':assert np.array_equal(det,(g.p_adjusted<=.05)&(g.statistic>g.permutation_mean))
    tp=(det&t).sum();fp=(det&~t).sum();pr=tp/det.sum() if det.sum() else np.nan;rec=tp/60
    row=dict(dgp=d,replication_id=r,method=m,true_positive=tp,false_positive=fp,discoveries=det.sum(),recall=rec,background_fpr=fp/340,precision=pr,f1=2*tp/(2*tp+fp+(60-tp)),rank_auc=auc(s,t),p_rank_auc=auc(-g.p_raw,t),target_raw_p05=((g.p_raw<=.05)&g.truth).sum()/60,target_median_p=g.loc[g.truth,'p_raw'].median(),bh_critical_p=g.loc[g.p_adjusted<=.05,'p_raw'].max(),top60_recall=t[np.argsort(-s,kind='stable')[:60]].sum()/60)
    if m!='local_geary':
        row.update(centered_auc=auc(g.statistic-g.permutation_mean,t),target_mean_statistic=g.loc[g.truth,'statistic'].mean(),target_mean_expectation=g.loc[g.truth,'permutation_mean'].mean(),background_mean_statistic=g.loc[~g.truth,'statistic'].mean(),background_mean_expectation=g.loc[~g.truth,'permutation_mean'].mean())
    rows.append(row)
    for target in [True,False]:
        cats=['all','A','B','C'] if d==4 else ['all']
        for cat in cats:
            h=g[(g.truth==target)&((g.cat1==cat) if cat!='all' else True)]
            if len(h):strata.append(dict(dgp=d,replication_id=r,method=m,target=target,category=cat,n=len(h),detections=h.detected.sum(),statistic=h.statistic.mean(),expectation=h.permutation_mean.mean(),median_p=h.p_raw.median()))
    if d==4:
        for tc in ['A','non-A']:
            for bc in ['A','non-A']:
                x=s[t & ((g.cat1=='A').to_numpy() if tc=='A' else (g.cat1!='A').to_numpy())]
                y=s[~t & ((g.cat1=='A').to_numpy() if bc=='A' else (g.cat1!='A').to_numpy())]
                if len(x) and len(y):pairauc.append(dict(replication_id=r,method=m,target_category=tc,background_category=bc,pairs=len(x)*len(y),auc=((x[:,None]>y).sum()+.5*(x[:,None]==y).sum())/(len(x)*len(y))))
r=pd.DataFrame(rows);r.to_csv(P/'recomputed_metrics.csv',index=False)
r.groupby(['dgp','method']).mean(numeric_only=True).drop(columns='replication_id').to_csv(P/'diagnostic_summary.csv')
s=pd.DataFrame(strata);s.to_csv(P/'category_diagnostics_by_replication.csv',index=False)
ss=s.groupby(['dgp','method','target','category']).agg(mean_n=('n','mean'),mean_detections=('detections','mean'),mean_statistic=('statistic','mean'),mean_expectation=('expectation','mean'),mean_median_p=('median_p','mean'),total_n=('n','sum'),total_detections=('detections','sum'));ss['pooled_detection_rate']=ss.total_detections/ss.total_n;ss.to_csv(P/'category_diagnostics_summary.csv')
pa=pd.DataFrame(pairauc);pa.to_csv(P/'category_pair_auc_by_replication.csv',index=False)
pa.groupby(['method','target_category','background_category']).agg(mean_auc=('auc','mean'),mean_pairs=('pairs','mean')).to_csv(P/'category_pair_auc_summary.csv')
x=a[(a.dgp==4)&(a.method=='path_continuous')];y=a[(a.dgp==4)&(a.method=='path_mixed')]
w=x.merge(y,on=['replication_id','id'],suffixes=('_continuous','_mixed'))
w['transition']=np.select([~w.detected_continuous&w.detected_mixed,w.detected_continuous&~w.detected_mixed,w.detected_continuous&w.detected_mixed],['gain','loss','both'],default='neither')
w.groupby(['truth_continuous','cat1_continuous','transition']).size().rename('cells_across_100_datasets').to_csv(P/'dgp4_detection_transitions.csv')
# Post-hoc p-value tails and category contribution to whole-sample AUC.
tails=[]
for (method,rep),g in a[a.dgp==3].groupby(['method','replication_id']):
    for cut in [.0001,.0005,.001,.002,.003,.005,.01,.05]:
        tails.append(dict(method=method,replication_id=rep,p_cut=cut,target_rate=(g.loc[g.truth,'p_raw']<=cut).mean(),background_rate=(g.loc[~g.truth,'p_raw']<=cut).mean()))
pd.DataFrame(tails).groupby(['method','p_cut']).mean().drop(columns='replication_id').to_csv(P/'dgp3_pvalue_tail_rates.csv')
pa['auc_contribution']=pa['auc']*pa['pairs']/(60*340)
pa.groupby(['method','target_category']).agg(mean_auc_contribution=('auc_contribution','sum')).div(100).to_csv(P/'category_auc_contributions.csv')

print('PASS: Appendix C diagnostics generated from 240,000 cells; BH checks passed')
