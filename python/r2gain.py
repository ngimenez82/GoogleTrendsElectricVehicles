import sys, json; sys.path.insert(0,'.')
from common import *
import warnings; warnings.filterwarnings('ignore')
OUT = sys.argv[1] if len(sys.argv)>1 else '.'
d=load(); d['cm']=d.iso2+'_'+d.tm.dt.month.astype(str)
out={}
for fe,lab in [('iso2','country FE + month dummies'),('cm','country x month FE')]:
    full=['L_dlnbev']+ATT+MAC+([] if fe=='cm' else MDUM)
    base=['L_dlnbev']+MAC+([] if fe=='cm' else MDUM)
    s=d[['d_ln_bev','iso2','cm']+full].dropna()
    r1=fe_ols(s,'d_ln_bev',base,fe=fe); r2=fe_ols(s,'d_ln_bev',full,fe=fe)
    out[fe]=dict(label=lab,N=int(r2['N']),r2_base=r2['r2w'] if False else r1['r2w'],r2_full=r2['r2w'],gain=r2['r2w']-r1['r2w'],
                 sum=sum_test(r2,ATT),joint=wald(r2,ATT),rho=float(r2['b']['L_dlnbev']))
    print(fe, out[fe])
# keep the full column for Table 2 (country x month FE)
r=fe_ols(d,'d_ln_bev',['L_dlnbev']+ATT+MAC,fe='cm')
o={'coef':{v:(float(r['b'][v]),float(np.sqrt(r['V'].loc[v,v])),float(2*(1-stats.t.cdf(abs(r['b'][v]/np.sqrt(r['V'].loc[v,v])),r['df'])))) for v in r['b'].index},'N':int(r['N']),'G':int(r['G']),'r2':r['r2w']}
est,se,p=sum_test(r,ATT); F,pj=wald(r,ATT); rho=r['b']['L_dlnbev']
idx={v:i for i,v in enumerate(r['b'].index)}; gg=np.zeros(len(idx))
for v in ATT: gg[idx[v]]=1/(1-rho)
gg[idx['L_dlnbev']]=est/(1-rho)**2; lrse=float(np.sqrt(gg@r['V'].values@gg)); lr=est/(1-rho)
o.update(sum=(est,se,p),joint=(F,pj),LR=(lr,lrse,float(2*(1-stats.t.cdf(abs(lr/lrse),r['df'])))))
out['T2_cm']=o
json.dump(out,open(f'{OUT}/out_r2gain.json','w'),indent=1,default=float)
