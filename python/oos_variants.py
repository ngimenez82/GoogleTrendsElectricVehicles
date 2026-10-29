import sys, json; sys.path.insert(0,'.')
from common import *
import warnings; warnings.filterwarnings('ignore')
import statsmodels.api as sm
from scipy import stats as st
d = load(); g=d.groupby('iso2')
d['ln_bev']=d.ln_bev.astype(float)
d['att_ma3'] = d[['L1_att','L2_att','L3_att']].sum(1, min_count=3)
d['att_ma6'] = d[ATT].sum(1, min_count=6)
d['lvl_att_L1'] = g['lr_elec_vs_base'].shift(1)
d['ln_share']=np.log(d.ev_bev/d.registered_cars); d['d_share']=g.ln_share.diff(); d['L_dshare']=d.groupby('iso2').d_share.shift(1)
def fitpred(tr,te,y,X):
    tr=tr[[y,'iso2']+X].dropna(); te=te[['iso2',y]+X].dropna(); te=te[te.iso2.isin(tr.iso2.unique())]
    my=tr.groupby('iso2')[y].mean(); mX=tr.groupby('iso2')[X].mean()
    b=np.linalg.lstsq((tr[X]-mX.loc[tr.iso2].values).values,(tr[y]-my.loc[tr.iso2].values).values,rcond=None)[0]
    return pd.Series(my.loc[te.iso2].values-mX.loc[te.iso2].values@b+te[X].values@b, index=te.index)
def nw(x,l): m=sm.OLS(np.asarray(x),np.ones(len(x))).fit(cov_type='HAC',cov_kwds={'maxlags':l}); return m.tvalues[0]
def evaluate(name, y, ben, alt, start=48, lagdv=None):
    ix = d[[y,'iso2']+list(set(ben+alt))].dropna().index
    rows=[]
    for t in range(start, d.t.max()+1):
        tr=d[(d.t<=t-1)&d.index.isin(ix)]; te=d[(d.t==t)&d.index.isin(ix)]
        if len(te)==0: continue
        rows.append(pd.DataFrame({'t':t,'iso2':te.iso2,'y':te[y],'fb':fitpred(tr,te,y,ben),'fa':fitpred(tr,te,y,alt)}))
    r=pd.concat(rows).dropna()
    rb=np.sqrt(((r.y-r.fb)**2).mean()); ra=np.sqrt(((r.y-r.fa)**2).mean())
    cw=((r.y-r.fb)**2-((r.y-r.fa)**2-(r.fb-r.fa)**2)); cwt=r.assign(c=cw).groupby('t').c.mean()
    dmt=r.assign(c=(r.y-r.fb)**2-(r.y-r.fa)**2).groupby('t').c.mean()
    tcw=nw(cwt,1); tdm=nw(dmt,1)
    o=dict(n=len(r),T=r.t.nunique(),ratio=ra/rb,CW_t=tcw,CW_p=1-st.norm.cdf(tcw),DM_t=tdm,DM_p=2*(1-st.norm.cdf(abs(tdm))))
    print(f'{name:40s} n={o["n"]} ratio={o["ratio"]:.4f} CW p={o["CW_p"]:.3f} DM t={tdm:.2f}'); return o
out={}
B=['L_dlnbev']+MAC+MDUM; B0=['L_dlnbev']+MDUM
out['ma3_vs_macro']=evaluate('MA(1-3) attention vs macro','d_ln_bev',B,B+['att_ma3'])
out['ma6_vs_macro']=evaluate('Sum(1-6) attention vs macro','d_ln_bev',B,B+['att_ma6'])
out['ma3_vs_AR']=evaluate('MA(1-3) attention vs AR','d_ln_bev',B0,B0+['att_ma3'])
out['now_vs_AR']=evaluate('L0 + MA3 vs AR (nowcast)','d_ln_bev',B0,B0+['L0_att','att_ma3'])
out['lvl_vs_AR']=evaluate('GT level L1 vs AR','d_ln_bev',B0,B0+['lvl_att_L1'])
out['ma3_start2023']=evaluate('MA3 vs AR, start 2023m1','d_ln_bev',B0,B0+['att_ma3'],start=60)
out['ma3_start2021']=evaluate('MA3 vs AR, start 2021m1','d_ln_bev',B0,B0+['att_ma3'],start=36)
S0=['L_dshare']+MDUM
out['share_ma3']=evaluate('Share target: MA3 vs AR','d_share',S0,S0+['att_ma3'])
out['share_now']=evaluate('Share target: L0+MA3 vs AR','d_share',S0,S0+['L0_att','att_ma3'])
json.dump(out,open('out_oos_variants.json','w'),indent=1,default=float)
