import sys, json; sys.path.insert(0,'.')
from common import *
import warnings; warnings.filterwarnings('ignore')
d = load()
rng = np.random.default_rng(12345)

def wild_boot_sum(d, y, X, names, B=1999):
    """Wild cluster bootstrap-t (Webb 6-point weights), restricted under H0: sum=0."""
    r = fe_ols(d,y,X); s = r['sample']
    est,se,_ = sum_test(r,names); t0 = est/se
    cols=[y]+X; dm = s[cols]-s.groupby('iso2')[cols].transform('mean')
    Xm = dm[X].values; Y = dm[y].values
    R = np.array([1.0 if v in names else 0.0 for v in X])
    # restricted estimate: impose R b = 0
    XtX_inv = np.linalg.pinv(Xm.T@Xm); bu = XtX_inv@Xm.T@Y
    br = bu - XtX_inv@R*(R@bu)/(R@XtX_inv@R)
    er = Y - Xm@br; fit = Xm@br
    groups = pd.factorize(s['iso2'])[0]; G = groups.max()+1
    w6 = np.array([-np.sqrt(1.5),-1,-np.sqrt(.5),np.sqrt(.5),1,np.sqrt(1.5)])
    N,K = len(Y), len(X)+1; adj = G/(G-1)*(N-1)/(N-K)
    tb = np.empty(B)
    for b in range(B):
        w = rng.choice(w6, G)[groups]
        Yb = fit + er*w; bb = XtX_inv@Xm.T@Yb; eb = Yb-Xm@bb
        meat = np.zeros((K-1,K-1))
        for g in range(G):
            idx = groups==g; sc = Xm[idx].T@eb[idx]; meat += np.outer(sc,sc)
        Vb = XtX_inv@meat@XtX_inv*adj
        tb[b] = (R@bb)/np.sqrt(R@Vb@R)
    return float(np.mean(np.abs(tb)>=abs(t0)))

out = {}
def summarize(name, r, att=ATT, lag='L_dlnbev'):
    est,se,p = sum_test(r,att); F,pj = wald(r,att)
    rho = r['b'][lag]
    # long-run multiplier on the level of log registrations: sum(beta)/(1-rho)
    g = np.zeros(len(r['b'])); idx={v:i for i,v in enumerate(r['b'].index)}
    lr = est/(1-rho)
    # delta method
    for v in att: g[idx[v]] = 1/(1-rho)
    g[idx[lag]] = est/(1-rho)**2
    lr_se = np.sqrt(g@r['V'].values@g)
    o = dict(N=int(r['N']), G=int(r['G']), sum=est, sum_se=se, p_sum=p, F_joint=F, p_joint=pj,
             rho=rho, rho_se=float(np.sqrt(r['V'].loc[lag,lag])), LR=lr, LR_se=lr_se,
             coefs={v: (float(r['b'][v]), float(np.sqrt(r['V'].loc[v,v]))) for v in r['b'].index if not v.startswith('mdum')})
    if 'r2w' in r: o['r2']=r['r2w']
    if 'r2mg' in r: o['r2']=r['r2mg']
    out[name]=o
    print(f"{name:28s} N={o['N']} G={o['G']} sum={est:.3f}({se:.3f}) p_sum={p:.3f} p_joint={pj:.4f} rho={rho:.3f} LR={lr:.3f}({lr_se:.3f})")
    return o

y='d_ln_bev'; base=['L_dlnbev']
r = fe_ols(d,y,base+MAC+MDUM); summarize('FE_macro', r, att=[]) if False else None
out['FE_macro'] = dict(N=int(r['N']),G=int(r['G']),r2=r['r2w'],coefs={v:(float(r['b'][v]),float(np.sqrt(r['V'].loc[v,v]))) for v in r['b'].index if not v.startswith('mdum')})
r = fe_ols(d,y,base+ATT+MAC+MDUM); summarize('FE_att', r)
out['FE_att']['p_wildboot_sum'] = wild_boot_sum(d,y,base+ATT+MAC+MDUM,ATT); print(' wild bootstrap p(sum)=',out['FE_att']['p_wildboot_sum'])
c = cce_mg(d,y,base+MAC,['d_ln_bev'])
out['CCE_macro']=dict(N=int(c['N']),G=int(c['G']),r2=c['r2mg'],coefs={v:(float(c['b'][v]),float(np.sqrt(c['V'].loc[v,v]))) for v in c['b'].index})
c = cce_mg(d,y,base+ATT+MAC,['d_ln_bev','d_att']); summarize('CCE_att', c)
c4 = cce_mg(d,y,base+ATT+MAC,['d_ln_bev','d_att'],cr_lags=4); summarize('CCE_att_crlags4', c4)
cj = cce_jackknife(d,y,base+ATT+MAC,['d_ln_bev','d_att']); summarize('CCE_att_jackknife', cj)

# winsorized DV (correct implementation)
lo,hi = d[y].quantile([.01,.99]); d['yw']=d[y].clip(lo,hi); d['L_yw']=d.groupby('iso2').yw.shift(1)
summarize('FE_att_winsor', fe_ols(d,'yw',['L_yw']+ATT+MAC+MDUM), lag='L_yw')
summarize('CCE_att_winsor', cce_mg(d,'yw',['L_yw']+ATT+MAC,['yw','d_att']), lag='L_yw')

# Leads placebo
LEAD=[f'F{k}_att' for k in range(1,7)]
ra = fe_ols(d,y,base+LEAD+ATT+MAC+MDUM)
Fl,pl = wald(ra,LEAD); Fg,pg = wald(ra,ATT)
out['FE_leads'] = dict(N=int(ra['N']), p_joint_leads=pl, p_joint_lags=pg, sum_leads=sum_test(ra,LEAD), sum_lags=sum_test(ra,ATT),
     coefs={v:(float(ra['b'][v]),float(np.sqrt(ra['V'].loc[v,v]))) for v in LEAD+ATT+['L_dlnbev']})
print('leads p=',pl,' lags p=',pg)

# Placebo attention: other mobility log-ratios (public transport, bike) vs same base
for nm,var in [('public','lr_public_vs_base'),('bike','lr_bike_vs_base')]:
    d['dd_'+nm] = d.groupby('iso2')[var].diff()
    P=[f'L{k}_{nm}' for k in range(1,7)]
    for k in range(1,7): d[P[k-1]] = d.groupby('iso2')['dd_'+nm].shift(k)
    rp = fe_ols(d,y,base+P+MAC+MDUM); est,se,p = sum_test(rp,P); F,pj=wald(rp,P)
    out['placebo_'+nm]=dict(N=int(rp['N']),sum=est,sum_se=se,p_sum=p,p_joint=pj)
    print('placebo',nm,est,se,p,pj)
    rb = fe_ols(d,y,base+ATT+P+MAC+MDUM); e2=sum_test(rb,ATT); e3=sum_test(rb,P)
    out['placebo_'+nm+'_horse']=dict(att=e2,plac=e3); print('  horse race att',e2,'placebo',e3)

# Market-share DV: growth of BEV share in new registrations
d['ln_share'] = np.log(d.ev_bev/d.registered_cars); d['d_ln_share']=d.groupby('iso2').ln_share.diff(); d['L_dshare']=d.groupby('iso2').d_ln_share.shift(1)
summarize('FE_att_share', fe_ols(d,'d_ln_share',['L_dshare']+ATT+MAC+MDUM), lag='L_dshare')
summarize('CCE_att_share', cce_mg(d,'d_ln_share',['L_dshare']+ATT+MAC,['d_ln_share','d_att']), lag='L_dshare')

# Without CCI -> expands sample (adds non-EU countries)
MAC2=['L_unem_rate','L_cpi_all','L_d_rel_price']
r=fe_ols(d,y,base+ATT+MAC2+MDUM); summarize('FE_att_noCCI', r); out['FE_att_noCCI']['countries']=sorted(r['sample'].iso2.unique())
# Year FE (common shocks, simple alternative to CCE)
YD=[f'ydum{k}' for k in range(2,9)]
summarize('FE_att_yearFE', fe_ols(d,y,base+ATT+MAC+MDUM+YD))
# time FE (country + month-year FE): within-time variation only
d['tt']=d['t']; tdum=pd.get_dummies(d.t,prefix='T',drop_first=True,dtype=float); d=pd.concat([d,tdum],axis=1)
summarize('FE_att_twoway', fe_ols(d,y,base+ATT+MAC+list(tdum.columns)))

json.dump(out, open('out_insample.json','w'), indent=1, default=float)
