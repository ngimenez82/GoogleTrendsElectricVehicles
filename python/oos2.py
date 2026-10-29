"""Extended pseudo out-of-sample evaluation for the Journal of Forecasting version."""
import sys, json; sys.path.insert(0,'.')
from common import *
import warnings; warnings.filterwarnings('ignore')
import statsmodels.api as sm
from scipy import stats as st
OUTDIR = sys.argv[1] if len(sys.argv)>1 else '.'
d = load(); g = d.groupby('iso2')
d['ln_bev'] = d.ln_bev.astype(float)
d['cm'] = d.iso2 + '_' + d.tm.dt.month.astype(str)
# hybrid seasonality: common month dummies, except Ireland (semi-annual registration calendar) which gets its own month effects
d['ck'] = np.where(d.iso2 == 'IE', d['cm'], d['iso2'])
HD = [f'hd{m}' for m in range(2, 13)]
for m in range(2, 13): d[f'hd{m}'] = d[f'mdum{m}'] * (d.iso2 != 'IE')
d['att_ma3'] = d[['L1_att','L2_att','L3_att']].sum(1, min_count=3)/3
d['SN'] = g['d_ln_bev'].shift(12)            # seasonal naive forecast of monthly growth
d['ln_cars'] = np.log(d.registered_cars); d['d_cars'] = g.ln_cars.diff(); d['L_dcars'] = g.d_cars.shift(1)
for h in [3]:
    d[f'y_h{h}'] = g['ln_bev'].shift(-(h-1)) - g['ln_bev'].shift(1)
    d[f'SN_h{h}'] = g['ln_bev'].shift(12-(h-1)) - g['ln_bev'].shift(13)
B0 = ['L_dlnbev']+MDUM
BH = ['L_dlnbev']+HD
SPEC = {'M0c': (B0,'iso2'), 'GT3c': (B0+['att_ma3'],'iso2'),
        'M0h': (BH,'ck'), 'GT3h': (BH+['att_ma3'],'ck'), 'GT6h': (BH+ATT,'ck'), 'GT3h+now': (BH+['L0_att','att_ma3'],'ck'), 'M0h+macro': (BH+MAC,'ck'),
        'M0': (B0,'cm'), 'M0+macro': (B0+MAC,'cm'), 'GT3': (B0+['att_ma3'],'cm'), 'GT3+now': (B0+['L0_att','att_ma3'],'cm'), 'GT6': (B0+ATT,'cm')}
SPEC_CARS = {'M0h': (BH,'ck'), 'GT3h': (BH+['att_ma3'],'ck'), 'M0h+cars': (BH+['L_dcars'],'ck'), 'GT3h+cars': (BH+['L_dcars','att_ma3'],'ck')}
POLICY = {('DE','2022-12'),('DE','2023-01'),('DE','2023-08'),('DE','2023-09'),('DE','2023-12'),('DE','2024-01'),('SE','2022-11'),('SE','2022-12'),('SE','2023-01')}

def fit_pooled(tr, te, y, X, key='iso2'):
    """Pooled slopes; intercepts by `key` (country, or country x calendar month = country-specific seasonality)."""
    X=[v for v in X if not (key=='cm' and v.startswith('mdum'))]
    tr=tr[[y,key]+X].dropna(); te=te[[key]+X].dropna(); te=te[te[key].isin(tr[key].unique())]
    my=tr.groupby(key)[y].mean(); mX=tr.groupby(key)[X].mean()
    b=np.linalg.lstsq((tr[X]-mX.loc[tr[key]].values).values,(tr[y]-my.loc[tr[key]].values).values,rcond=None)[0]
    return pd.Series(my.loc[te[key]].values-mX.loc[te[key]].values@b+te[X].values@b, index=te.index)
def fit_het(tr, te, y, X):
    out={}
    for c,te_c in te.groupby('iso2'):
        tc=tr[tr.iso2==c][[y]+X].dropna(); tec=te_c[X].dropna()
        if len(tc)<len(X)+8 or len(tec)==0: continue
        Z=np.column_stack([np.ones(len(tc)),tc[X].values]); b=np.linalg.lstsq(Z,tc[y].values,rcond=None)[0]
        out.update(dict(zip(tec.index, np.column_stack([np.ones(len(tec)),tec[X].values])@b)))
    return pd.Series(out)

def run(y='d_ln_bev', start=48, window=None, specs=SPEC, het=False, need=None):
    need = need or sorted(set(sum([v[0] for v in specs.values()],[])))
    ix = d[[y,'iso2','SN']+need].dropna().index
    rows=[]
    for t in range(start, d.t.max()+1):
        lo = -1 if window is None else t-1-window
        tr=d[(d.t<=t-1)&(d.t>lo)&d.index.isin(ix)]; te=d[(d.t==t)&d.index.isin(ix)]
        if len(te)==0: continue
        r=te[['iso2','t',y,'SN']].rename(columns={y:'y'}).copy(); r['ZERO']=0.0
        for m,(X,key) in specs.items():
            r[m]=fit_pooled(tr,te,y,X,key)
        if het:
            r['M0_het']=fit_het(tr,te,y,B0); r['GT3_het']=fit_het(tr,te,y,B0+['att_ma3'])
        rows.append(r)
    r=pd.concat(rows)
    if 'GT3' in r: r['COMB']=(r['M0']+r['GT3'])/2
    if 'GT3h' in r: r['COMBh']=(r['M0h']+r['GT3h'])/2
    r['tm']=d.loc[r.index,'tm'].dt.strftime('%Y-%m'); r['policy']=[(a,b) in POLICY for a,b in zip(r.iso2,r.tm)]
    return r.dropna()

def nw(x,l=1):
    m=sm.OLS(np.asarray(x,float),np.ones(len(x))).fit(cov_type='HAC',cov_kwds={'maxlags':l}); return float(m.tvalues[0])
def compare(r, alt, ben, nested=True, lag=1):
    e_a=r.y-r[alt]; e_b=r.y-r[ben]
    o=dict(rmse_ratio=float(np.sqrt((e_a**2).mean()/(e_b**2).mean())), mae_ratio=float(e_a.abs().mean()/e_b.abs().mean()))
    dm=r.assign(z=e_b**2-e_a**2).groupby('t').z.mean(); o['DM_t']=nw(dm,lag); o['DM_p']=float(2*(1-st.norm.cdf(abs(o['DM_t']))))
    if nested:
        cw=r.assign(z=e_b**2-(e_a**2-(r[ben]-r[alt])**2)).groupby('t').z.mean(); o['CW_t']=nw(cw,lag); o['CW_p']=float(1-st.norm.cdf(o['CW_t']))
    return o
def pt_test(y,f):
    """Pesaran-Timmermann (1992) directional accuracy test."""
    y=np.sign(y); f=np.sign(f); n=len(y); P=np.mean(y==f); py=np.mean(y>0); pf=np.mean(f>0)
    Ps=py*pf+(1-py)*(1-pf); vP=Ps*(1-Ps)/n
    vPs=((2*py-1)**2*pf*(1-pf)+(2*pf-1)**2*py*(1-py))/n+4*py*pf*(1-py)*(1-pf)/n**2
    S=(P-Ps)/np.sqrt(vP-vPs); return float(P), float(1-st.norm.cdf(S))
def absolute(r, models):
    return {m: dict(rmse=float(np.sqrt(((r.y-r[m])**2).mean())), mae=float((r.y-r[m]).abs().mean()), hit=pt_test(r.y.values, r[m].values)) for m in models}

out={}
r = run(het=True)
r.to_csv(f'{OUTDIR}/oos2_h1_forecasts.csv', index=False)
models=['ZERO','SN','M0c','GT3c','M0','M0+macro','GT3','GT3+now','GT6','COMB','M0_het','GT3_het','M0h','GT3h','GT6h','GT3h+now','M0h+macro','COMBh']
out['n']=len(r); out['T']=int(r.t.nunique()); out['G']=int(r.iso2.nunique()); out['n_policy']=int(r.policy.sum())
out['abs']=absolute(r, models)
cmp={}
for alt,ben,nest in [('ZERO','M0',False),('SN','M0',False),('M0c','M0',False),('M0_het','M0',False),('GT3c','M0c',True),('M0+macro','M0',True),('GT3','M0',True),('GT3+now','M0',True),('GT6','M0',True),('COMB','M0',False),('GT3_het','M0_het',True),
                     ('M0c','M0h',False),('M0','M0h',False),('ZERO','M0h',False),('SN','M0h',False),('M0_het','M0h',False),('M0h+macro','M0h',True),
                     ('GT3h','M0h',True),('GT6h','M0h',True),('GT3h+now','M0h',True),('COMBh','M0h',False)]:
    cmp[f'{alt}|{ben}']=compare(r,alt,ben,nest)
out['cmp']=cmp
# policy episodes vs rest
ep={}
for lab,sub in [('policy',r[r.policy]),('rest',r[~r.policy]),('all',r)]:
    ep[lab]=dict(n=len(sub), rmse_M0=float(np.sqrt(((sub.y-sub.M0h)**2).mean())), rmse_GT3=float(np.sqrt(((sub.y-sub.GT3h)**2).mean())),
                 msfe_share=float(((sub.y-sub.M0h)**2).sum()/((r.y-r.M0h)**2).sum()))
out['episodes']=ep
out['rest_cmp']=compare(r[~r.policy],'GT3h','M0h',True)
# Dec/Jan months for all countries (year-end pull-forward)
dj=r.tm.str[-2:].isin(['12','01'])
out['decjan']=dict(share_msfe=float(((r.y-r.M0h)**2)[dj].sum()/((r.y-r.M0h)**2).sum()), share_obs=float(dj.mean()),
                   ratio_decjan=float(np.sqrt(((r.y-r.GT3h)**2)[dj].mean()/((r.y-r.M0h)**2)[dj].mean())),
                   ratio_other=float(np.sqrt(((r.y-r.GT3h)**2)[~dj].mean()/((r.y-r.M0h)**2)[~dj].mean())))
# CSSED for figure
cs=r.assign(z=(r.y-r.M0h)**2-(r.y-r.GT3h)**2).groupby('tm').z.sum().cumsum()
out['cssed']=cs.to_dict()
rc=run(specs=SPEC_CARS); out['cars']=dict(n=len(rc), cars=compare(rc,'M0h+cars','M0h',True), gt=compare(rc,'GT3h+cars','M0h+cars',True))
# rolling window
rr=run(window=36, specs={'M0h':(BH,'ck'),'GT3h':(BH+['att_ma3'],'ck')})
out['rolling36']=compare(rr,'GT3h','M0h',True); out['rolling36']['n']=len(rr)
# alternative starts
for s,lab in [(36,'start2021'),(60,'start2023')]:
    rs=run(start=s, specs={'M0h':(BH,'ck'),'GT3h':(BH+['att_ma3'],'ck')}); out[lab]=compare(rs,'GT3h','M0h',True); out[lab]['n']=len(rs)
# h=3 direct
r3=run(y='y_h3', specs={'M0h':(BH,'ck'),'GT3h':(BH+['att_ma3'],'ck'),'GT6h':(BH+ATT,'ck')})
out['h3']={'GT3':compare(r3,'GT3h','M0h',True,lag=3),'GT6':compare(r3,'GT6h','M0h',True,lag=3),'n':len(r3)}
# BEV share target
d['ln_share']=np.log(d.ev_bev/d.registered_cars); d['d_share']=g.ln_share.diff(); d['L_dshare']=g.d_share.shift(1)
S0=['L_dshare']+MDUM
d_keep=d.copy()
rsh=run(y='d_share', specs={'M0h':(['L_dshare']+HD,'ck'),'GT3h':(['L_dshare']+HD+['att_ma3'],'ck')})
out['share']=compare(rsh,'GT3h','M0h',True); out['share']['n']=len(rsh)
json.dump(out, open(f'{OUTDIR}/out_oos2.json','w'), indent=1, default=float)
print(json.dumps({k:v for k,v in out.items() if k!='cssed'}, indent=1, default=float))
# Excluding Ireland from the evaluation
rx = r[r.iso2 != 'IE']
out['exIE'] = dict(n=len(rx), gt_M0h=compare(rx,'GT3h','M0h',True), gt_M0=compare(rx,'GT3','M0',True), M0_M0h=compare(rx,'M0','M0h',False))
json.dump(out, open(f'{OUTDIR}/out_oos2.json','w'), indent=1, default=float)
print('exIE', out['exIE'])
