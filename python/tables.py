import sys, json; sys.path.insert(0,'.')
from common import *
import warnings; warnings.filterwarnings('ignore')
d = load(); g=d.groupby('iso2')
d['l_gasdies']=np.log(d.gasdies_ktons.where(d.gasdies_ktons>0)); d['L_dgasdies']=g.l_gasdies.diff().groupby(d.iso2).shift(1)
d['l_n2o']=np.log(d.n2o_gg.where(d.n2o_gg>0)); d['L_dn2o']=g.l_n2o.diff().groupby(d.iso2).shift(1)
lo,hi=d.d_ln_bev.quantile([.01,.99]); d['yw']=d.d_ln_bev.clip(lo,hi); d['L_yw']=g.yw.shift(1)
LEAD=[f'F{k}_att' for k in range(1,7)]
def col(r, lag='L_dlnbev', att=ATT, kind='FE', extra=None):
    o={'coef':{}, 'N':int(r['N']), 'G':int(r['G'])}
    for v in r['b'].index:
        if v.startswith('mdum'): continue
        b=float(r['b'][v]); se=float(np.sqrt(r['V'].loc[v,v])); p=2*(1-stats.t.cdf(abs(b/se), r['df']))
        o['coef']['L_dlnbev' if v==lag else v]=(b,se,p)
    o['r2']= r.get('r2w', r.get('r2mg'))
    if att and all(a in r['b'].index for a in att):
        est,se,p=sum_test(r,att); F,pj=wald(r,att); rho=r['b'][lag]
        idx={v:i for i,v in enumerate(r['b'].index)}; gg=np.zeros(len(idx))
        for v in att: gg[idx[v]]=1/(1-rho)
        gg[idx[lag]]=est/(1-rho)**2; lrse=float(np.sqrt(gg@r['V'].values@gg)); lr=est/(1-rho)
        o.update(sum=(est,se,p), joint=(F,pj), LR=(lr,lrse,2*(1-stats.t.cdf(abs(lr/lrse),r['df']))))
    return o
T={}
y='d_ln_bev'
T['T2']=[col(fe_ols(d,y,['L_dlnbev']+MAC+MDUM),att=None),
         col(fe_ols(d,y,['L_dlnbev']+ATT+MAC+MDUM)),
         col(cce_mg(d,y,['L_dlnbev']+MAC,['d_ln_bev']),att=None),
         col(cce_mg(d,y,['L_dlnbev']+ATT+MAC,['d_ln_bev','d_att']))]
T['T3']=[T['T2'][3],
         col(cce_mg(d,y,['L_dlnbev']+ATT+MAC+['L_dgasdies'],['d_ln_bev','d_att'])),
         col(cce_mg(d,y,['L_dlnbev']+ATT+MAC+['L_dn2o'],['d_ln_bev','d_att']))]
T['A1']=[col(fe_ols(d,'yw',['L_yw']+MAC+MDUM),lag='L_yw',att=None),
         col(fe_ols(d,'yw',['L_yw']+ATT+MAC+MDUM),lag='L_yw'),
         col(cce_mg(d,'yw',['L_yw']+MAC,['yw']),lag='L_yw',att=None),
         col(cce_mg(d,'yw',['L_yw']+ATT+MAC,['yw','d_att']),lag='L_yw')]
ra=fe_ols(d,y,['L_dlnbev']+LEAD+ATT+MAC+MDUM); a2=col(ra)
a2['joint_leads']=wald(ra,LEAD); a2['sum_leads']=sum_test(ra,LEAD); T['A2']=[a2]
# Descriptives on FE col (2) estimation sample
s=fe_ols(d,y,['L_dlnbev']+ATT+MAC+MDUM)['sample']; idx=s.index
dd=d.loc[idx]
desc={}
for v,lab in [('ev_bev','New BEV registrations (monthly)'),('d_ln_bev','Δln(new BEV registrations)'),('lr_elec_vs_base','ln(relative EV searches), LR'),('d_att','Δln(relative EV searches)'),('unem_rate','Unemployment rate (%)'),('cpi_all','CPI inflation, YoY (%)'),('cci','Consumer confidence index'),('d_rel_price','Δ relative energy price (diesel vs electricity)')]:
    x=dd[v].astype(float); desc[lab]=(int(x.notna().sum()),float(x.mean()),float(x.std()),float(x.min()),float(x.max()))
T['T1']=desc; T['T1_countries']=sorted(dd.iso2.unique()); T['T1_period']=(str(dd.tm.min().date()),str(dd.tm.max().date()))
ins=json.load(open('out_insample.json')); T['ins']=ins
T['oos']=json.load(open('out_oos.json')); T['oosv']=json.load(open('out_oos_variants.json'))
json.dump(T,open('out_tables.json','w'),indent=1,default=float)
for k in ['T2','T3','A1']:
    for c in T[k]: print(k, c['N'], c['G'], round(c['r2'],3), {kk:tuple(round(x,3) for x in v) for kk,v in c.items() if kk in ('sum','joint','LR')})
print('A2', T['A2'][0]['joint'], T['A2'][0]['joint_leads'], T['A2'][0]['N'])
print(desc, T['T1_period'])
