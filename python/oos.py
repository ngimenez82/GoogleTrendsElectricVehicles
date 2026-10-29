import sys, json; sys.path.insert(0,'.')
from common import *
import warnings; warnings.filterwarnings('ignore')
import statsmodels.api as sm
d = load()
d['ln_bev'] = d['ln_bev'].astype(float)
g = d.groupby('iso2')
# direct h-step targets: cumulative log growth from t-1 to t+h-1
for h in [1,3,6]:
    d[f'y_h{h}'] = g['ln_bev'].shift(-(h-1)) - g['ln_bev'].shift(1)
SPECS = {
 'M0_AR':     ['L_dlnbev']+MDUM,
 'M1_macro':  ['L_dlnbev']+MAC+MDUM,
 'M2_GT':     ['L_dlnbev']+ATT+MAC+MDUM,
 'M3_nowcast':['L_dlnbev','L0_att']+ATT+MAC+MDUM,
}
def fit_predict(train, test, y, X):
    tr = train[[y,'iso2']+X].dropna(); te = test[['iso2',y]+X].dropna()
    mu_y = tr.groupby('iso2')[y].mean(); mu_X = tr.groupby('iso2')[X].mean()
    Xm = tr[X]-mu_X.loc[tr.iso2].values; Ym = tr[y]-mu_y.loc[tr.iso2].values
    b = np.linalg.lstsq(Xm.values, Ym.values, rcond=None)[0]
    te = te[te.iso2.isin(mu_y.index)]
    a = mu_y.loc[te.iso2].values - mu_X.loc[te.iso2].values@b
    return te.assign(pred = a + te[X].values@b)

def run(h, start=48, cols_needed=None):
    # forecast origin: target month t, info up to t-1 (for nowcast also GT of month t)
    res=[]
    common_ix = d[[f'y_h{h}','iso2']+SPECS['M3_nowcast']].dropna().index
    for t in range(start, d.t.max()-h+2):
        # training: targets fully observed before forecast origin -> t' + h - 1 <= t-1
        train = d[(d.t <= t-h)].loc[lambda x: x.index.isin(common_ix)]
        test = d[(d.t==t)].loc[lambda x: x.index.isin(common_ix)]
        if len(test)==0: continue
        row = test[['iso2','t',f'y_h{h}']].rename(columns={f'y_h{h}':'y'}).copy()
        for m,X in SPECS.items():
            p = fit_predict(train,test,f'y_h{h}',X).set_index('iso2')['pred']
            row[m] = row.iso2.map(p).values
        row['RW'] = 0.0  # no-change in log level (zero growth) benchmark
        res.append(row)
    return pd.concat(res)

def nw_t(x, lags=None):
    x = np.asarray(x); T=len(x); lags = lags if lags is not None else int(np.floor(4*(T/100)**(2/9)))
    m = sm.OLS(x, np.ones(T)).fit(cov_type='HAC', cov_kwds={'maxlags':max(lags,1)})
    return float(m.params[0]), float(m.tvalues[0]), float(m.pvalues[0])

out={}
for h in [1,3,6]:
    r = run(h)
    e = {m: r['y']-r[m] for m in list(SPECS)+['RW']}
    rmse = {m: float(np.sqrt((e[m]**2).mean())) for m in e}
    o=dict(h=h, n=len(r), T=r.t.nunique(), G=r.iso2.nunique(), first=str(d.tm[d.t==r.t.min()].iloc[0].date()), last=str(d.tm[d.t==r.t.max()].iloc[0].date()), rmse=rmse)
    for alt,ben in [('M2_GT','M1_macro'),('M3_nowcast','M1_macro'),('M2_GT','M0_AR'),('M3_nowcast','M0_AR')]:
        # Clark-West adjusted loss differential, averaged across countries each month (panel), HAC over time
        fb, fa = r[ben], r[alt]
        cw = (r['y']-fb)**2 - ((r['y']-fa)**2 - (fb-fa)**2)
        dm = (r['y']-fb)**2 - (r['y']-fa)**2
        cwt = r.assign(cw=cw).groupby('t').cw.mean(); dmt = r.assign(dm=dm).groupby('t').dm.mean()
        _,tcw,_ = nw_t(cwt.values, lags=h)  ; _,tdm,pdm = nw_t(dmt.values, lags=h)
        from scipy import stats as st
        pcw = 1-st.norm.cdf(tcw)  # one-sided
        # per-country share with lower MSE
        cs = r.assign(la=(r['y']-fa)**2, lb=(r['y']-fb)**2).groupby('iso2')[['la','lb']].mean()
        o[f'{alt}_vs_{ben}'] = dict(ratio=rmse[alt]/rmse[ben], CW_t=tcw, CW_p=pcw, DM_t=tdm, DM_p=pdm, share_countries_better=float((cs.la<cs.lb).mean()))
    out[f'h{h}']=o
    print(json.dumps(o, indent=1, default=float))
    if h==1: r.to_csv('oos_h1_forecasts.csv', index=False)
    if h==3: r.to_csv('oos_h3_forecasts.csv', index=False)
json.dump(out, open('out_oos.json','w'), indent=1, default=float)
