# Figure 1: study flow chart. Counts are taken from sql/04_flowchart_counts.sql.
# Contains aggregate counts only (no patient-level data).
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib import font_manager as fm
fam = 'Liberation Sans' if any('Liberation Sans' in f.name for f in fm.fontManager.ttflist) else 'DejaVu Sans'
plt.rcParams['font.family'] = fam
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'

W, H = 10.65, 11.9
Y0 = 1.0
fig = plt.figure(figsize=(W*0.72, H*0.72))
ax = fig.add_axes([0,0,1,1]); ax.set_xlim(0,W); ax.set_ylim(Y0, Y0+H); ax.axis('off')
LW = 1.0; FS = 10.5; LH = 0.30   # line height in units

def box(cx, cy, w, h, fc='white', ec='black', lw=LW):
    ax.add_patch(FancyBboxPatch((cx-w/2, cy-h/2), w, h,
        boxstyle='round,pad=0,rounding_size=0.12', fc=fc, ec=ec, lw=lw))

def main(cx, cy, w, lines, fc='white', bold_first=False):
    h = LH*len(lines) + 0.36
    box(cx, cy, w, h, fc=fc)
    y0 = cy + LH*(len(lines)-1)/2
    for i, t in enumerate(lines):
        ax.text(cx, y0 - i*LH, t, ha='center', va='center', fontsize=FS,
                fontweight='bold' if (bold_first and i == 0) else 'normal')
    return h

def excl(xl, cy, w, header, bullets):
    n = 1 + len(bullets)
    h = LH*n + 0.34
    box(xl+w/2, cy, w, h, fc='#f4f4f4', ec='#7a7a7a')
    y0 = cy + LH*(n-1)/2
    ax.text(xl+0.2, y0, header, ha='left', va='center', fontsize=FS, fontweight='bold')
    for i, b in enumerate(bullets):
        ax.text(xl+0.2, y0-(i+1)*LH, '\u2022  ' + b, ha='left', va='center', fontsize=FS)
    return h

def arrow(x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(arrowstyle='-|>,head_length=0.45,head_width=0.22',
                        color='black', lw=LW, shrinkA=0, shrinkB=0))

XM, MW = 3.4, 3.7
XE, EW = 5.65, 4.85

ys = {'b0': 12.3, 'b1': 10.3, 'b2': 8.0, 'b3': 5.8, 'b4': 3.7}
hb = {}
hb['b0'] = main(XM, ys['b0'], MW, ['ICU stays in MIMIC-IV v3.1', '(n = 94,458)'])
hb['b1'] = main(XM, ys['b1'], MW, ['First ICU stay of adult patients', '(n = 65,366)'])
hb['b2'] = main(XM, ys['b2'], MW, ['Eligible ICU stays', '(n = 29,915)'])
hb['b3'] = main(XM, ys['b3'], MW, ['Serial acid-base data available', '(n = 7,699)'])
hb['b4'] = main(XM, ys['b4'], MW, ['Study cohort', '(n = 6,173)'], fc='#e6eef8', bold_first=True)

order = ['b0','b1','b2','b3','b4']
for a, b in zip(order[:-1], order[1:]):
    arrow(XM, ys[a]-hb[a]/2, XM, ys[b]+hb[b]/2)

ex = [
    ((ys['b0']+ys['b1'])/2, 'Excluded (n = 29,092)', ['Repeat ICU stays']),
    ((ys['b1']+ys['b2'])/2, 'Excluded (n = 35,451)',
        ['ICU length of stay < 48 h (n = 34,224)', 'Died within 48 h (n = 148)',
         'ESKD or chronic dialysis (n = 1,079)']),
    ((ys['b2']+ys['b3'])/2, 'Excluded (n = 22,216)',
        ['No paired BE and lactate, 0-24 h (n = 11,487)',
         'No paired BE and lactate, 24-48 h (n = 10,729)']),
    ((ys['b3']+ys['b4'])/2, 'Excluded (n = 1,526)',
        ['No creatinine at ICU admission (n = 202)', 'KRT started within 48 h (n = 432)',
         'KDIGO stage 2-3 AKI within 48 h (n = 892)']),
]
for cy, hd, bl in ex:
    excl(XE, cy, EW, hd, bl)
    arrow(XM, cy, XE, cy)

# outcome split
ysplit = 2.75
yo = 1.75
OW = 3.0
xl, xr = XM-1.6, XM+1.6
ax.plot([XM, XM], [ys['b4']-hb['b4']/2, ysplit], color='black', lw=LW)
ax.plot([xl, xr], [ysplit, ysplit], color='black', lw=LW)
ho = LH*3 + 0.36
for x in (xl, xr):
    arrow(x, ysplit, x, yo+ho/2)
main(xl, yo, OW, ['KDIGO stage 2-3 AKI or KRT', '(48 h to day 7)', 'n = 629 (10.2%)'])
main(xr, yo, OW, ['No primary outcome', '(48 h to day 7)', 'n = 5,544 (89.8%)'])

for ext, kw in [('png', dict(dpi=600)), ('pdf', {}), ('svg', {}), ('tiff', dict(dpi=600, pil_kwargs={'compression':'tiff_lzw'}))]:
    fig.savefig(f'../figures/Figure1_flowchart.{ext}', facecolor='white', **kw)

print(fam)
