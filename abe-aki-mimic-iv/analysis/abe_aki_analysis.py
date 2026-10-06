# =====================================================================
# abe_aki_analysis.py
#
# Early trajectory of alactic base excess (ABE) and acute kidney injury
# in critically ill adults: a retrospective cohort study using MIMIC-IV.
#
# Statistical analysis pipeline. Run after sql/01_cohort.sql and
# sql/02_covariates.sql have created the tables `abe.kohort` and
# `abe.kovaryat` in your own Google Cloud project.
#
# Designed for Google Colaboratory. Only aggregate results are printed;
# no patient-level data are written to disk.
#
# Variable names are in Turkish (as used during the analysis); an
# English glossary is provided in README.md.
# =====================================================================

PROJECT_ID = 'YOUR_PROJECT_ID'      # <-- replace with your Google Cloud project ID
N_IMPUTATIONS = 20

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import chi2, chi2_contingency, fisher_exact, norm
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from matplotlib import gridspec
from google.cloud import bigquery

warnings.filterwarnings('ignore')

try:  # authenticate when running in Google Colab
    from google.colab import auth
    auth.authenticate_user()
except ImportError:
    pass

client = bigquery.Client(project=PROJECT_ID)
SQL_DIR = Path('sql')                     # path to the sql/ folder of this repository
FIG_DIR = Path('figures'); FIG_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------
df = client.query(f"""
SELECT k.*, v.* EXCEPT(stay_id)
FROM `{PROJECT_ID}.abe.kohort` k
JOIN `{PROJECT_ID}.abe.kovaryat` v ON k.stay_id = v.stay_id
""").to_dataframe()

ek_sql = (SQL_DIR / '03_additional_covariates.sql').read_text().replace('YOUR_PROJECT_ID', PROJECT_ID)
df = df.merge(client.query(ek_sql).to_dataframe(), on='stay_id', how='left')

# Fluid balance winsorised at the 1st and 99th percentiles, in litres
lo_fb, hi_fb = df['sivi_dengesi_48s'].quantile([0.01, 0.99])
df['sivi_dengesi_48s_w'] = df['sivi_dengesi_48s'].clip(lo_fb, hi_fb) / 1000

print('Cohort size:', len(df), '| events:', int(df['sonlanim_aki23_rrt'].sum()))


# ---------------------------------------------------------------------
# 2. Descriptive tables (Tables 1 and 2)
# ---------------------------------------------------------------------
from tableone import TableOne  # pip install tableone

t1_cols = ['yas', 'cinsiyet', 'kabul_tipi', 'yb_tipi',
           'abe_ort_0_24', 'abe_min_0_24', 'abe_ort_24_48', 'delta_abe',
           'laktat_ort_0_24', 'be_ort_0_24', 'kre_bazal', 'kre_ilk', 'bazal_kaynak',
           'vazopressor_48s', 'imv_48s', 'sivi_dengesi_48s_w', 'idrar_ml_kg_saat',
           'map_min_24s', 'trombosit_min', 'bilirubin_max', 'bun_max', 'klor_max',
           'bikarbonat_min', 'hemoglobin_min', 'lokosit_max',
           'vankomisin_iv', 'aminoglikozid', 'pip_tazo', 'loop_diuretik', 'nsaid',
           'kbh', 'diyabet', 'kky', 'siroz', 'sepsis_icd', 'hospital_expire_flag']
t1_cat = ['cinsiyet', 'kabul_tipi', 'yb_tipi', 'bazal_kaynak', 'vazopressor_48s', 'imv_48s',
          'vankomisin_iv', 'aminoglikozid', 'pip_tazo', 'loop_diuretik', 'nsaid',
          'kbh', 'diyabet', 'kky', 'siroz', 'sepsis_icd', 'hospital_expire_flag']
t1 = TableOne(df, columns=t1_cols, categorical=t1_cat,
              nonnormal=[c for c in t1_cols if c not in t1_cat],
              groupby='sonlanim_aki23_rrt', pval=True, missing=True)
print(t1.tabulate(tablefmt='github'))

print('\nAdditional binary variables (Fisher exact test when expected counts < 5):')
y0 = df['sonlanim_aki23_rrt'] == 0
for v in ['bikarbonat_inf', 'eritrosit_tx', 'taze_plazma_tx', 'asetazolamid']:
    a, b = df.loc[y0, v], df.loc[~y0, v]
    tab = [[(a == 0).sum(), (a == 1).sum()], [(b == 0).sum(), (b == 1).sum()]]
    use_fisher = (chi2_contingency(tab)[3] < 5).any()
    p = fisher_exact(tab)[1] if use_fisher else chi2_contingency(tab)[1]
    print(f"{v}: overall {int(df[v].sum())} ({df[v].mean()*100:.1f}) | "
          f"no outcome {int(a.sum())} ({a.mean()*100:.1f}) | outcome {int(b.sum())} ({b.mean()*100:.1f}) | P = {p:.4f}")


# ---------------------------------------------------------------------
# 3. Variable preparation
# ---------------------------------------------------------------------
def yb_grup(x):
    if 'CVICU' in x: return 'CVICU'
    if 'CCU' in x: return 'CCU'
    if 'TSICU' in x: return 'TSICU'
    if x.startswith('Neuro'): return 'Noro'
    if x.startswith('Medical/Surgical'): return 'MICU_SICU'
    if x.startswith('Medical Intensive'): return 'MICU'
    if x.startswith('Surgical Intensive'): return 'SICU'
    return 'Diger'


def kabul_grup(x):
    if x in ('ELECTIVE', 'SURGICAL SAME DAY ADMISSION'): return 'Planli'
    if x == 'URGENT': return 'Urgent'
    if x == 'OBSERVATION ADMIT': return 'Observasyon'
    return 'Acil'


d0 = df.copy()
d0['yb_grup'] = d0['yb_tipi'].apply(yb_grup)
d0['kabul_grup'] = d0['kabul_tipi'].apply(kabul_grup)
d0['erkek'] = (d0['cinsiyet'] == 'M').astype(int)
d0['bazal_onceki'] = (d0['bazal_kaynak'] == 'onceki_7gun_min').astype(int)
d0['log_kre_bazal'] = np.log(d0['kre_bazal'].clip(lower=0.2))
d0['log_laktat'] = np.log(d0['laktat_ort_0_24'].clip(lower=0.3))
d0['log_bilirubin'] = np.log(d0['bilirubin_max'].clip(lower=0.1))
d0['log_bun'] = np.log(d0['bun_max'].clip(lower=1))
d0['idrar'] = d0['idrar_ml_kg_saat'].clip(0, 5)
d0['y'] = d0['sonlanim_aki23_rrt'].astype(int)
d0['cvicu'] = (d0['yb_grup'] == 'CVICU').astype(int)
d0['olum_7g'] = (d0['deathtime'].notna() &
                 (d0['deathtime'] <= d0['intime'] + pd.Timedelta(days=7))).astype(int)
d0['y_bilesik'] = ((d0['y'] == 1) | (d0['olum_7g'] == 1)).astype(int)
d0['delta_kat'] = pd.cut(d0['delta_abe'], [-np.inf, -2, 0, 2, np.inf],
                         labels=['dec_ge2', 'dec_lt2', 'rise_0_2', 'rise_gt2']).astype(str)

IMP_COLS = ['yas', 'log_kre_bazal', 'log_laktat', 'abe_ort_0_24', 'delta_abe',
            'sivi_dengesi_48s_w', 'idrar', 'map_min_24s', 'trombosit_min',
            'log_bilirubin', 'log_bun', 'hemoglobin_min', 'y']

COV = ("yas + erkek + C(kabul_grup) + C(yb_grup) + log_kre_bazal + bazal_onceki"
       " + log_laktat + vazopressor_48s + imv_48s + sivi_dengesi_48s_w + idrar"
       " + map_min_24s + trombosit_min + log_bilirubin + log_bun + hemoglobin_min"
       " + vankomisin_iv + pip_tazo + loop_diuretik + bikarbonat_inf + eritrosit_tx"
       " + kbh + kky + siroz + sepsis_icd")


def drop_term(cov, term):
    return cov.replace(f" + {term}", "").replace(f"{term} + ", "")


def fit(formula, data):
    try:
        return smf.logit(formula, data).fit(disp=0, method='newton', maxiter=100)
    except Exception:
        return smf.logit(formula, data).fit(disp=0, method='bfgs', maxiter=2000)


def fmt(q, se, p, dig=3):
    ps = '<0.001' if p < 0.001 else f'{p:.3f}'
    return f"{np.exp(q):.{dig}f} ({np.exp(q - 1.96*se):.{dig}f}-{np.exp(q + 1.96*se):.{dig}f}), P = {ps}"


# ---------------------------------------------------------------------
# 4. Multiple imputation (20 datasets) and Rubin's rules
# ---------------------------------------------------------------------
imps = []
for i in range(N_IMPUTATIONS):
    t = d0.copy()
    t[IMP_COLS] = IterativeImputer(random_state=i, sample_posterior=True,
                                   max_iter=20).fit_transform(t[IMP_COLS])
    t['y'] = d0['y']
    imps.append(t)


def pool(formula, terms, subset=None):
    b = {k: [] for k in terms}
    w = {k: [] for k in terms}
    for t in imps:
        x = t if subset is None else t[subset(t)]
        m = fit(formula, x)
        for k in terms:
            b[k].append(m.params[k]); w[k].append(m.bse[k] ** 2)
    out = {}
    for k in terms:
        bb, ww = np.array(b[k]), np.array(w[k])
        q = bb.mean()
        se = np.sqrt(ww.mean() + (1 + 1 / N_IMPUTATIONS) * bb.var(ddof=1))
        out[k] = (q, se, 2 * norm.sf(abs(q / se)))
    return out


# ---------------------------------------------------------------------
# 5. Primary analysis (Table 3)
# ---------------------------------------------------------------------
print('\n=== Primary analysis ===')
m_unadj = smf.logit('y ~ delta_abe', d0).fit(disp=0)
print('Unadjusted, per 1 mmol/L increase:',
      fmt(m_unadj.params['delta_abe'], m_unadj.bse['delta_abe'], m_unadj.pvalues['delta_abe']))
r = pool(f"y ~ abe_ort_0_24 + delta_abe + {COV}", ['delta_abe', 'abe_ort_0_24'])
print('Adjusted (MI), per 1 mmol/L increase:', fmt(*r['delta_abe']))
q, se, p = r['delta_abe']
print('Adjusted (MI), per 1 mmol/L decrease:', fmt(-q, se, p))
print('Day-1 ABE in the same model:', fmt(*r['abe_ort_0_24']))

print('\n=== Categorical delta ABE (MI; reference 0 < dABE <= 2) ===')
ref = "C(delta_kat, Treatment(reference='rise_0_2'))"
terms = [f"{ref}[T.{k}]" for k in ['dec_ge2', 'dec_lt2', 'rise_gt2']]
r = pool(f"y ~ abe_ort_0_24 + {ref} + {COV}", terms)
for k in ['dec_ge2', 'dec_lt2', 'rise_0_2', 'rise_gt2']:
    sub = d0[d0['delta_kat'] == k]
    est = 'reference' if k == 'rise_0_2' else fmt(*r[f"{ref}[T.{k}]"], dig=2)
    print(f"{k}: n = {len(sub)}, events = {int(sub['y'].sum())} ({sub['y'].mean()*100:.1f}%), OR {est}")


# ---------------------------------------------------------------------
# 6. Incremental value and trajectory vs level (Supplementary Table S1)
#    Single deterministic imputation
# ---------------------------------------------------------------------
d1 = d0.copy()
d1[IMP_COLS] = IterativeImputer(random_state=1, max_iter=20).fit_transform(d1[IMP_COLS])
d1['y'] = d0['y']


def lrt(big, small):
    stat = 2 * (big.llf - small.llf)
    df_ = big.df_model - small.df_model
    return f"LR chi2 = {stat:.2f}, df = {df_:.0f}, P = {chi2.sf(stat, df_):.4f}"


m0 = fit(f"y ~ {COV}", d1)
m1 = fit(f"y ~ abe_ort_0_24 + {COV}", d1)
m2 = fit(f"y ~ abe_ort_0_24 + delta_abe + {COV}", d1)
m_d2 = fit(f"y ~ abe_ort_24_48 + {COV}", d1)
m_both = fit(f"y ~ abe_ort_0_24 + abe_ort_24_48 + {COV}", d1)
print('\n=== Incremental value (single imputation) ===')
print('Model 1 vs Model 0 (add day-1 ABE):', lrt(m1, m0))
print('Model 2 vs Model 1 (add delta ABE):', lrt(m2, m1))
for name, m in [('Model 0', m0), ('Model 1', m1), ('Model 2', m2),
                ('Day-2 ABE only', m_d2), ('Day-1 + day-2 ABE', m_both)]:
    print(f"AUROC {name}: {roc_auc_score(d1['y'], m.predict()):.4f}")
for label, m, v in [('Day-2 ABE only', m_d2, 'abe_ort_24_48'),
                    ('Both days, day-1', m_both, 'abe_ort_0_24'),
                    ('Both days, day-2', m_both, 'abe_ort_24_48')]:
    print(f"{label}: OR {fmt(m.params[v], m.bse[v], m.pvalues[v])}")
print('Added value of day-1 ABE given day-2 ABE:', lrt(m_both, m_d2))


# ---------------------------------------------------------------------
# 7. Non-linearity (natural cubic regression spline, 4 df, centred)
# ---------------------------------------------------------------------
f_lin = f"y ~ abe_ort_0_24 + delta_abe + {COV}"
f_spl = f"y ~ abe_ort_0_24 + cr(delta_abe, df=4, constraints='center') + {COV}"
lr_stats, df_diff = [], None
for t in imps:
    y1, X1 = patsy.dmatrices(f_lin, t, return_type='dataframe')
    y2, X2 = patsy.dmatrices(f_spl, t, return_type='dataframe')
    a1 = sm.Logit(y1, X1).fit(disp=0, method='newton', maxiter=100)
    a2 = sm.Logit(y2, X2).fit(disp=0, method='newton', maxiter=100)
    lr_stats.append(2 * (a2.llf - a1.llf))
    df_diff = np.linalg.matrix_rank(X2.values) - np.linalg.matrix_rank(X1.values)
lr_stats = np.array(lr_stats)
print(f"\n=== Non-linearity: median LR chi2 {np.median(lr_stats):.2f} "
      f"(range {lr_stats.min():.2f}-{lr_stats.max():.2f}), df = {df_diff}, "
      f"P = {chi2.sf(np.median(lr_stats), df_diff):.4f}")


# ---------------------------------------------------------------------
# 8. Sensitivity analyses (Table 4; MI)
# ---------------------------------------------------------------------
print('\n=== Sensitivity analyses (OR per 1 mmol/L increase in delta ABE) ===')
scenarios = [
    ('Primary analysis',                         None,                               COV,                              'y'),
    ('Without sodium bicarbonate',               lambda t: t['bikarbonat_inf'] == 0, drop_term(COV, 'bikarbonat_inf'), 'y'),
    ('Without loop diuretics',                   lambda t: t['loop_diuretik'] == 0,  drop_term(COV, 'loop_diuretik'),  'y'),
    ('Arterial samples only',                    lambda t: t['n_venoz'] == 0,        COV,                              'y'),
    ('Baseline creatinine obtained before ICU',  lambda t: t['bazal_onceki'] == 1,   drop_term(COV, 'bazal_onceki'),   'y'),
    ('Urine output removed from the model',      None,                               drop_term(COV, 'idrar'),          'y'),
    ('Composite outcome incl. death within 7 d', None,                               COV,                              'y_bilesik'),
]
for name, sub, cov, outcome in scenarios:
    x = d0 if sub is None else d0[sub(d0)]
    r = pool(f"{outcome} ~ abe_ort_0_24 + delta_abe + {cov}", ['delta_abe'], sub)
    print(f"{name}: n = {len(x)}, events = {int(x[outcome].sum())}, OR {fmt(*r['delta_abe'])}")


# ---------------------------------------------------------------------
# 9. Subgroup analyses (Supplementary Table S2, Figure 3; MI)
# ---------------------------------------------------------------------
print('\n=== Subgroups ===')
subgroups = [
    ('Sepsis',                          'sepsis_icd',      'sepsis_icd',      ['No', 'Yes']),
    ('Chronic kidney disease',          'kbh',             'kbh',             ['No', 'Yes']),
    ('Invasive mechanical ventilation', 'imv_48s',         'imv_48s',         ['No', 'Yes']),
    ('Vasopressors',                    'vazopressor_48s', 'vazopressor_48s', ['No', 'Yes']),
    ('Loop diuretics',                  'loop_diuretik',   'loop_diuretik',   ['No', 'Yes']),
    ('ICU type',                        'cvicu',           'C(yb_grup)',      ['Other ICUs', 'Cardiovascular ICU']),
]
overall = pool(f"y ~ abe_ort_0_24 + delta_abe + {COV}", ['delta_abe'])['delta_abe']
fp_rows = []
for name, var, term, labels in subgroups:
    cov = drop_term(COV, term)
    p_int = pool(f"y ~ abe_ort_0_24 + delta_abe * {var} + {cov}",
                 [f'delta_abe:{var}'])[f'delta_abe:{var}'][2]
    fp_rows.append(dict(label=name, header=True))
    for lev in [0, 1]:
        x = d0[d0[var] == lev]
        q, se, p = pool(f"y ~ abe_ort_0_24 + delta_abe + {cov}", ['delta_abe'],
                        lambda t, var=var, lev=lev: t[var] == lev)['delta_abe']
        fp_rows.append(dict(label=labels[lev], header=False, n=len(x), ev=int(x['y'].sum()),
                            q=q, se=se, pint=(p_int if lev == 0 else None)))
        print(f"{name}: {labels[lev]} | n = {len(x)} | events = {int(x['y'].sum())} | OR {fmt(q, se, p)}"
              + (f" | P interaction = {p_int:.3f}" if lev == 0 else ''))


# ---------------------------------------------------------------------
# 10. Figures 2 and 3
# ---------------------------------------------------------------------
plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                     'font.size': 10, 'axes.linewidth': 0.8, 'pdf.fonttype': 42})


def save(fig, name):
    fig.savefig(FIG_DIR / f'{name}.png', dpi=600, bbox_inches='tight', facecolor='white')
    fig.savefig(FIG_DIR / f'{name}.pdf', bbox_inches='tight', facecolor='white')
    fig.savefig(FIG_DIR / f'{name}.tiff', dpi=600, bbox_inches='tight', facecolor='white',
                pil_kwargs={'compression': 'tiff_lzw'})


# Figure 2: adjusted OR relative to delta ABE = 0, pooled over imputations
lo, hi = d0['delta_abe'].quantile([0.02, 0.98])
grid = np.linspace(lo, hi, 200)
E, V, risks = [], [], []
for t in imps:
    yy, XX = patsy.dmatrices(f_spl, t, return_type='dataframe')
    di = XX.design_info
    m = sm.Logit(yy, XX).fit(disp=0, method='newton', maxiter=100)
    base = t.iloc[[0] * len(grid)].copy().reset_index(drop=True)
    refr = base.copy()
    base['delta_abe'] = grid
    refr['delta_abe'] = 0.0
    D = (np.asarray(patsy.build_design_matrices([di], base)[0])
         - np.asarray(patsy.build_design_matrices([di], refr)[0]))
    E.append(D @ m.params.values)
    V.append(np.einsum('ij,jk,ik->i', D, m.cov_params().values, D))
    rp = []
    for v in [-4, 0, 7]:  # g-computation: average predicted risk
        tt = t.copy(); tt['delta_abe'] = v
        Xt = np.asarray(patsy.build_design_matrices([di], tt)[0])
        rp.append((1 / (1 + np.exp(-(Xt @ m.params.values)))).mean() * 100)
    risks.append(rp)
E, V = np.array(E), np.array(V)
q = E.mean(0)
se = np.sqrt(V.mean(0) + (1 + 1 / N_IMPUTATIONS) * E.var(0, ddof=1))
print('\nAverage adjusted risk (%) at dABE -4 / 0 / +7:', np.round(np.mean(risks, 0), 1))
for v in [-4, -2, 2, 4]:
    i = np.argmin(abs(grid - v))
    print(f"dABE {v:+d}: OR {np.exp(q[i]):.2f} ({np.exp(q[i]-1.96*se[i]):.2f}-{np.exp(q[i]+1.96*se[i]):.2f})")

fig = plt.figure(figsize=(6.5, 5.2))
gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1], hspace=0.06)
ax = fig.add_subplot(gs[0]); axh = fig.add_subplot(gs[1], sharex=ax)
ax.fill_between(grid, np.exp(q - 1.96 * se), np.exp(q + 1.96 * se), color='#9ecae1', alpha=0.55, lw=0)
ax.plot(grid, np.exp(q), color='#08519c', lw=2)
ax.axhline(1, color='black', lw=0.8, ls='--')
ax.axvline(0, color='gray', lw=0.6, ls=':')
ax.set_yscale('log')
ax.set_yticks([0.3, 0.4, 0.5, 0.7, 1, 1.5, 2])
ax.set_yticklabels(['0.3', '0.4', '0.5', '0.7', '1.0', '1.5', '2.0'])
ax.minorticks_off()
ax.set_ylabel('Adjusted odds ratio (95% CI)')
ax.tick_params(labelbottom=False)
for s in ['top', 'right']:
    ax.spines[s].set_visible(False)
ax.text(0.98, 0.95, 'Reference: ΔABE = 0 mmol/L', transform=ax.transAxes,
        ha='right', va='top', fontsize=9, color='dimgray')
inside = d0['delta_abe'][(d0['delta_abe'] >= lo) & (d0['delta_abe'] <= hi)]
axh.hist(inside, bins=50, range=(lo, hi), color='#bdbdbd', edgecolor='white', lw=0.3)
axh.set_ylabel('Patients, n')
axh.set_xlabel('ΔABE (mmol/L), mean ABE 24-48 h minus mean ABE 0-24 h')
for s in ['top', 'right']:
    axh.spines[s].set_visible(False)
save(fig, 'Figure2_dABE_OR_curve')
plt.close(fig)

# Figure 3: forest plot
rows = [dict(label='Overall', header=False, overall=True, n=len(d0), ev=int(d0['y'].sum()),
             q=overall[0], se=overall[1], pint=None)] + fp_rows
nR = len(rows)
fig = plt.figure(figsize=(10, 0.36 * nR + 1.2))
gs = gridspec.GridSpec(1, 3, width_ratios=[1.55, 1.15, 1.2], wspace=0.02)
axL, axP, axR = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
ys = np.arange(nR)[::-1]
for a in (axL, axR):
    a.set_xlim(0, 1); a.set_ylim(-0.7, nR + 0.3); a.axis('off')
axP.set_ylim(-0.7, nR + 0.3)
H = nR - 0.1
axL.text(0.00, H, 'Subgroup', weight='bold', va='bottom')
axL.text(0.70, H, 'n', weight='bold', ha='right', va='bottom')
axL.text(0.95, H, 'Events', weight='bold', ha='right', va='bottom')
axR.text(0.05, H, 'OR (95% CI)', weight='bold', va='bottom')
axR.text(0.98, H, 'P for\ninteraction', weight='bold', ha='right', va='bottom', linespacing=1.0)
for yv, r in zip(ys, rows):
    if r['header']:
        axL.text(0.00, yv, r['label'], weight='bold', va='center')
        continue
    ind = 0.0 if r.get('overall') else 0.04
    axL.text(ind, yv, r['label'], va='center', weight='bold' if r.get('overall') else 'normal')
    axL.text(0.70, yv, f"{r['n']:,}", ha='right', va='center')
    axL.text(0.95, yv, f"{r['ev']:,}", ha='right', va='center')
    o = np.exp(r['q'])
    l, u = np.exp(r['q'] - 1.96 * r['se']), np.exp(r['q'] + 1.96 * r['se'])
    axP.plot([l, u], [yv, yv], color='black', lw=1)
    axP.plot(o, yv, marker='D' if r.get('overall') else 's', color='black',
             ms=6 if r.get('overall') else 5)
    axR.text(0.05, yv, f"{o:.2f} ({l:.2f}-{u:.2f})", va='center')
    if r['pint'] is not None:
        axR.text(0.98, yv - 0.5, '<0.001' if r['pint'] < 0.001 else f"{r['pint']:.3f}",
                 ha='right', va='center')
axP.axvline(1, color='black', lw=0.8, ls='--')
axP.axvline(np.exp(overall[0]), color='#08519c', lw=0.8, ls=':')
axP.set_xlim(0.75, 1.15)
axP.set_xticks([0.8, 0.9, 1.0, 1.1])
axP.set_yticks([])
for s in ['top', 'right', 'left']:
    axP.spines[s].set_visible(False)
axP.set_xlabel('Adjusted OR per 1 mmol/L increase in ΔABE')
save(fig, 'Figure3_forest_subgroups')
plt.close(fig)

# Package versions (for the Methods section)
import importlib.metadata as md
import platform
print('\nPython', platform.python_version())
for pkg in ['pandas', 'numpy', 'statsmodels', 'scikit-learn', 'patsy', 'scipy',
            'tableone', 'matplotlib', 'google-cloud-bigquery']:
    try:
        print(f'{pkg}=={md.version(pkg)}')
    except md.PackageNotFoundError:
        print(f'{pkg}: not installed')

print('\nDone. Figures saved to', FIG_DIR.resolve())
