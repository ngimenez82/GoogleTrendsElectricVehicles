import sys, json; import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
IN = sys.argv[1] if len(sys.argv)>1 else '.'; OUT = sys.argv[2] if len(sys.argv)>2 else '../figures'
r=pd.read_csv(f'{IN}/oos2_h1_forecasts.csv')
plt.rcParams.update({'font.family':'serif','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(7.2,3.5))
for alt,ben,ls,lab in [('GT3h','M0h','-','GT3 (average of attention lags 1–3) vs. benchmark M0'),('GT6h','M0h','--','GT6 (six attention lags) vs. benchmark M0')]:
    cs=r.assign(z=(r.y-r[ben])**2-(r.y-r[alt])**2).groupby('tm').z.sum().cumsum()
    ax.plot(pd.to_datetime(cs.index),cs.values,color='#1f3b73',ls=ls,lw=1.6,label=lab)
ax.axhline(0,color='#888888',lw=0.8)
for dt in ['2022-12-01','2023-08-01','2023-12-01']: ax.axvline(pd.Timestamp(dt),color='#bbbbbb',ls=':',lw=0.8)
ax.set_ylabel('Cumulative sum of squared\nforecast error differences'); ax.grid(axis='y',color='#e5e5e5',lw=0.6)
ax.legend(frameon=False,fontsize=8,loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2); fig.tight_layout(); fig.savefig(f'{OUT}/figure2.png',dpi=300)
