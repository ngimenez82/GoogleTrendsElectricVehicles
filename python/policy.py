"""Subsidy-deadline indicators from a policy calendar (CSV: iso2, effective_date, type, ...).
pre  = last month in which buyers could still register under the old scheme
       (the month before the effective date, or the same month if the change takes effect after the 15th);
post = first month fully under the new scheme."""
import sys, json; sys.path.insert(0,'.')
from common import *
import warnings; warnings.filterwarnings('ignore')
import pandas as pd, numpy as np
CAL = sys.argv[1] if len(sys.argv)>1 else '../policy/policy_calendar_provisional.csv'
OUT = sys.argv[2] if len(sys.argv)>2 else '.'
def add_policy(d, cal):
    d['pol_pre']=0.0; d['pol_post']=0.0
    for _,row in cal.iterrows():
        e=pd.Timestamp(row.effective_date)
        pre = (e - pd.offsets.MonthBegin(1)) if e.day<=15 else e.replace(day=1)
        post = pre + pd.offsets.MonthBegin(1)
        d.loc[(d.iso2==row.iso2)&(d.tm==pre),'pol_pre']=1.0
        d.loc[(d.iso2==row.iso2)&(d.tm==post),'pol_post']=1.0
    return d
if __name__=='__main__':
    d=load(); cal=pd.read_csv(CAL); d=add_policy(d,cal)
    d['cm']=d.iso2+'_'+d.tm.dt.month.astype(str)
    out={'calendar':cal.to_dict('records'),'n_pre':int(d.pol_pre.sum()),'n_post':int(d.pol_post.sum())}
    for lab,X,fe in [('FE',['L_dlnbev']+ATT+MAC+MDUM,'iso2'),('FE_cm',['L_dlnbev']+ATT+MAC,'cm')]:
        r0=fe_ols(d,'d_ln_bev',X,fe=fe); r1=fe_ols(d,'d_ln_bev',X+['pol_pre','pol_post'],fe=fe)
        out[lab]=dict(sum_without=sum_test(r0,ATT), sum_with=sum_test(r1,ATT), joint_with=wald(r1,ATT),
                      pre=(float(r1['b']['pol_pre']),float(np.sqrt(r1['V'].loc['pol_pre','pol_pre']))),
                      post=(float(r1['b']['pol_post']),float(np.sqrt(r1['V'].loc['pol_post','pol_post']))), N=int(r1['N']))
        print(lab, out[lab])
    json.dump(out,open(f'{OUT}/out_policy.json','w'),indent=1,default=float)
