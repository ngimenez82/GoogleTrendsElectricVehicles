import sys; sys.path.insert(0,'.')
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.dates as mdates
d=load()
plt.rcParams.update({'font.family':'serif','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
eps={'DE':[('2022-12-01','A'),('2023-08-01','B'),('2023-12-01','C')],'SE':[('2022-11-01','D')]}
names={'DE':'Germany','SE':'Sweden'}
fig,ax=plt.subplots(2,2,figsize=(7.2,4.6),sharex=True)
for j,c in enumerate(['DE','SE']):
    x=d[(d.iso2==c)&(d.tm>='2021-07-01')&(d.tm<='2024-06-01')]
    ax[0,j].plot(x.tm,x.ev_bev/1000,color='#1f3b73',lw=1.6)
    ax[1,j].plot(x.tm,x.lr_elec_vs_base,color='#1f3b73',lw=1.6)
    ax[0,j].set_title(names[c],fontsize=10)
    for k in range(2):
        for dt,lab in eps[c]:
            ax[k,j].axvline(pd.Timestamp(dt),color='#888888',ls='--',lw=0.8)
        ax[k,j].grid(axis='y',color='#e5e5e5',lw=0.6)
    for dt,lab in eps[c]:
        ax[0,j].annotate(lab,(pd.Timestamp(dt),ax[0,j].get_ylim()[1]*0.97),fontsize=8,ha='right',fontweight='bold',va='top',color='#444444',xytext=(-3,0),textcoords='offset points')
ax[0,0].set_ylabel('New BEV registrations\n(thousands)'); ax[1,0].set_ylabel('Relative EV search\nintensity (log-ratio)')
for a in ax[1]: a.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m')); a.tick_params(axis='x',rotation=45)
fig.tight_layout(); fig.savefig('../figures/figure1.png',dpi=300)
