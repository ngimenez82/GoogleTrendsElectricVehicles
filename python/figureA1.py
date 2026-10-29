"""Figure A1: relative EV search intensity and new BEV registrations, by country (standardized, 3-month moving averages)."""
import sys, json; sys.path.insert(0,'.')
from common import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.dates as mdates
OUT = sys.argv[1] if len(sys.argv)>1 else '../figures'; JOUT = sys.argv[2] if len(sys.argv)>2 else '.'
d = load()
C = ['AT','BE','CZ','DE','DK','EE','ES','FI','FR','GR','HU','IE','IT','LT','LU','LV','NL','PL','PT','SE','SI','SK']
NAMES = dict(AT='Austria',BE='Belgium',CZ='Czech Republic',DE='Germany',DK='Denmark',EE='Estonia',ES='Spain',FI='Finland',FR='France',GR='Greece',HU='Hungary',IE='Ireland',IT='Italy',LT='Lithuania',LU='Luxembourg',LV='Latvia',NL='Netherlands',PL='Poland',PT='Portugal',SE='Sweden',SI='Slovenia',SK='Slovakia')
plt.rcParams.update({'font.family':'serif','font.size':7,'axes.spines.top':False,'axes.spines.right':False})
fig, axs = plt.subplots(6, 4, figsize=(7.2, 9.0), sharex=True)
corr = {}
for k, c in enumerate(C):
    ax = axs.flat[k]; x = d[d.iso2==c].sort_values('tm')
    z = lambda s: (s - s.mean())/s.std()
    reg = z(x.ln_bev.astype(float).rolling(3, center=True).mean()); att = z(x.lr_elec_vs_base.rolling(3, center=True).mean())
    ax.plot(x.tm, reg, color='#1f3b73', lw=1.1, label='ln(new BEV registrations)')
    ax.plot(x.tm, att, color='#b85c00', lw=1.1, ls='--', label='Relative EV search intensity')
    ax.set_title(NAMES[c], fontsize=7.5, pad=2); ax.grid(axis='y', color='#e5e5e5', lw=0.5); ax.tick_params(length=2)
    # within-country correlations: levels and 12-month changes
    yy = x.ln_bev.astype(float); aa = x.lr_elec_vs_base
    corr[c] = dict(levels=float(np.corrcoef(yy, aa)[0,1]),
                   yoy=float(pd.concat([yy.diff(12), aa.diff(12)], axis=1).dropna().corr().iloc[0,1]))
for ax in axs.flat[len(C):]: ax.axis('off')
for ax in axs[-1]: ax.xaxis.set_major_locator(mdates.YearLocator(2)); ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
for ax in axs[-2][2:]: ax.xaxis.set_tick_params(labelbottom=True); ax.xaxis.set_major_locator(mdates.YearLocator(2)); ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
h, l = axs.flat[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower right', bbox_to_anchor=(0.97, 0.06), frameon=False, fontsize=8)
fig.tight_layout(h_pad=0.6, w_pad=0.6); fig.savefig(f'{OUT}/figureA1.png', dpi=300)
lv = np.array([v['levels'] for v in corr.values()]); yo = np.array([v['yoy'] for v in corr.values()])
summary = dict(by_country=corr, median_levels=float(np.median(lv)), median_yoy=float(np.median(yo)), n_pos_levels=int((lv>0).sum()), n_pos_yoy=int((yo>0).sum()))
json.dump(summary, open(f'{JOUT}/out_figA1.json','w'), indent=1)
print({k:v for k,v in summary.items() if k!='by_country'})
print({c:(round(v['levels'],2),round(v['yoy'],2)) for c,v in corr.items()})
