import pandas as pd, numpy as np
from scipy import stats

def load(path='../dataset.dta'):
    d = pd.read_stata(path)
    d['iso2'] = d['iso2'].astype(str)
    d = d.sort_values(['iso2','tm']).reset_index(drop=True)
    d['t'] = (d.tm.dt.year-2018)*12 + d.tm.dt.month - 1
    g = d.groupby('iso2')
    def L(v,k=1): return g[v].shift(k)
    d['L_dlnbev'] = L('d_ln_bev')
    for k in range(1,7):
        d[f'L{k}_att'] = L('d_att',k)
        d[f'F{k}_att'] = g['d_att'].shift(-k)
    d["L0_att"] = d["d_att"]; d = d.copy()
    for v in ['unem_rate','cci','cpi_all','d_rel_price']:
        d['L_'+v] = L(v)
    return d

MDUM = [f'mdum{m}' for m in range(2,13)]
ATT = [f'L{k}_att' for k in range(1,7)]
MAC = ['L_unem_rate','L_cci','L_cpi_all','L_d_rel_price']

def fe_ols(d, y, X, cluster='iso2', fe='iso2'):
    """Within estimator with country FE, Stata-style CR1 cluster SE (xtreg, fe vce(cluster))."""
    s = d[[y]+X+[cluster]+([fe] if fe!=cluster else [])].dropna().copy()
    cols=[y]+X
    dm = s[cols] - s.groupby(fe)[cols].transform('mean')
    Y = dm[y].values; Xm = dm[X].values
    XtX_inv = np.linalg.pinv(Xm.T@Xm)
    b = XtX_inv@Xm.T@Y
    e = Y - Xm@b
    G = s[cluster].nunique(); N=len(s); K=len(X)+1
    meat = np.zeros((len(X),len(X)))
    for _,idx in s.groupby(cluster).indices.items():
        sc = Xm[idx].T@e[idx]; meat += np.outer(sc,sc)
    V = XtX_inv@meat@XtX_inv * G/(G-1)*(N-1)/(N-K)
    r2w = 1 - (e@e)/((Y-Y.mean())@(Y-Y.mean()))
    return dict(b=pd.Series(b,X), V=pd.DataFrame(V,X,X), N=N, G=G, r2w=r2w, resid=e, sample=s, df=G-1)

def wald(res, names, R=None):
    b = res['b'][names].values; V = res['V'].loc[names,names].values
    if R is None:  # joint
        W = b@np.linalg.pinv(V)@b; q=len(names)
        F = W/q; p = 1-stats.f.cdf(F, q, res['df'])
        return F,p
    R=np.asarray(R); est=R@b; se=np.sqrt(R@V@R)
    tt=est/se; p=2*(1-stats.t.cdf(abs(tt),res['df']))
    return est,se,p

def sum_test(res, names):
    return wald(res, names, R=np.ones(len(names)))

def cce_mg(d, y, X, csa_vars, cr_lags=2, extra_country=MDUM, min_obs=None):
    """Dynamic CCE mean-group (Chudik-Pesaran 2015), mimicking xtdcce2 defaults:
    country-specific OLS with intercept, cross-sectional averages of csa_vars (lags 0..cr_lags)."""
    d = d.copy()
    csa_cols=[]
    for v in csa_vars:
        m = d.groupby('t')[v].transform('mean')
        d['csa_'+v] = m
    d = d.sort_values(['iso2','t'])
    for v in csa_vars:
        for l in range(cr_lags+1):
            c=f'csa_{v}_l{l}'; d[c]=d.groupby('iso2')['csa_'+v].shift(l); csa_cols.append(c)
    allX = X+extra_country+csa_cols
    s = d[[y,'iso2','t']+allX].dropna()
    # xtdcce2 builds CSA on estimation sample; recompute on that sample
    s2 = d.loc[s.index].copy()
    for v in csa_vars:
        sub = d.loc[d[[y]+X].dropna().index]
        m = sub.groupby('t')[v].mean()
        d['csa_'+v] = d['t'].map(m)
    for v in csa_vars:
        for l in range(cr_lags+1):
            c=f'csa_{v}_l{l}'; d[c]=d.groupby('iso2')['csa_'+v].shift(l)
    s = d[[y,'iso2','t']+allX].dropna()
    coefs=[]; r2=[]; Ns=0; ctry=[]
    for c,gdf in s.groupby('iso2'):
        Z = np.column_stack([np.ones(len(gdf)), gdf[allX].values])
        if len(gdf) <= Z.shape[1]+2: continue
        bb,_,_,_ = np.linalg.lstsq(Z, gdf[y].values, rcond=None)
        e = gdf[y].values - Z@bb
        r2.append(1-(e@e)/(((gdf[y]-gdf[y].mean())**2).sum()))
        coefs.append(bb[1:1+len(X)]); Ns+=len(gdf); ctry.append(c)
    B = np.array(coefs); Nc=len(B)
    bmg = B.mean(0); V = np.cov(B.T, ddof=1)/Nc
    return dict(b=pd.Series(bmg,X), V=pd.DataFrame(V,X,X), N=Ns, G=Nc, df=Nc-1,
                B=pd.DataFrame(B,index=ctry,columns=X), r2mg=np.mean(r2), sample=s)

def cce_jackknife(d, y, X, csa_vars, cr_lags=2):
    full = cce_mg(d,y,X,csa_vars,cr_lags)
    tmid = d['t'].median()
    a = cce_mg(d[d.t<=tmid],y,X,csa_vars,cr_lags); b = cce_mg(d[d.t>tmid],y,X,csa_vars,cr_lags)
    bj = 2*full['b'] - 0.5*(a['b']+b['b'])
    out = dict(full); out['b']=bj
    return out

def stars(p): return '***' if p<0.01 else '**' if p<0.05 else '*' if p<0.1 else ''
