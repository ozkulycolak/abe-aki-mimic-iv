# =====================================================================
# abe_aki_analysis.py  (version 2.2)
#
# Early changes in alactic base excess (ABE) and subsequent acute kidney
# injury in critically ill adults: a retrospective cohort study using
# MIMIC-IV.
#
# Run AFTER the SQL scripts 01, 02 and 05 have created the tables
# `abe.kohort`, `abe.kovaryat` and `abe.sofa_d1` in your own Google
# Cloud project (see README.md). Designed for Google Colaboratory:
#     !pip install -q tableone lifelines
#     %run analysis/abe_aki_analysis.py
# Only aggregate results are printed and written to results_v2_1.txt;
# no patient-level data are written to disk.
# Variable names are partly Turkish; see the glossary in README.md.
# =====================================================================

import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import statsmodels.formula.api as smf
from scipy.stats import chi2_contingency, fisher_exact, norm
from scipy.stats import f as f_dist
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from matplotlib import gridspec

warnings.filterwarnings('ignore')

PROJECT_ID = os.environ.get('ABE_PROJECT_ID', 'YOUR_PROJECT_ID')   # <-- your Google Cloud project
N_IMP = int(os.environ.get('ABE_N_IMP', 40))      # number of imputed datasets
N_BOOT = int(os.environ.get('ABE_N_BOOT', 200))   # bootstrap replicates
SEED = 2026

try:
    BASE = Path(__file__).resolve().parent.parent
except NameError:
    BASE = Path('.')
SQL_DIR = BASE / 'sql'
FIG_DIR = BASE / 'figures'
FIG_DIR.mkdir(exist_ok=True)
OUT = open(BASE / 'results_v2_2.txt', 'w')


def P(*args):
    s = ' '.join(str(a) for a in args)
    print(s)
    OUT.write(s + '\n')
    OUT.flush()


def H(title):
    P('\n' + '=' * 78 + '\n' + title + '\n' + '=' * 78)


if 'client' not in globals():
    try:
        from google.colab import auth
        auth.authenticate_user()
    except ImportError:
        pass
    from google.cloud import bigquery
    client = bigquery.Client(project=PROJECT_ID)


def q(sql):
    return client.query(sql).to_dataframe()


def read_sql(name):
    return (SQL_DIR / name).read_text().replace('YOUR_PROJECT_ID', PROJECT_ID)


def to_dt(s):
    if str(s.dtype) in ('dbdate', 'object'):
        return pd.to_datetime(s.astype(str), errors='coerce')
    return pd.to_datetime(s, errors='coerce')


def med_iqr(x, d=1):
    x = pd.Series(x).dropna()
    return f"{x.median():.{d}f} ({x.quantile(.25):.{d}f} to {x.quantile(.75):.{d}f})"


def pfmt(p):
    return '<0.001' if p < 0.001 else f'{p:.3f}'


# =====================================================================
# 1. DATA
# =====================================================================
df = q(f"""
SELECT k.*, v.* EXCEPT(stay_id)
FROM `{PROJECT_ID}.abe.kohort` k
JOIN `{PROJECT_ID}.abe.kovaryat` v ON k.stay_id = v.stay_id
""")
for f_ in ['03_additional_covariates.sql', '07_extra_covariates.sql']:
    df = df.merge(q(read_sql(f_)), on='stay_id', how='left')
sofa = q(f"""
SELECT stay_id, pf_min, gcs_min,
       sofa_resp + sofa_coag + sofa_liver + sofa_cv + sofa_cns AS sofa_nonrenal
FROM `{PROJECT_ID}.abe.sofa_d1`
""")
df = df.merge(sofa, on='stay_id', how='left')

for c in ['intime', 'outtime', 'dischtime', 'deathtime', 'aki_zaman', 'aki_kre_zaman',
          'rrt_ilk', 'kre_bazal_zaman']:
    if c in df.columns:
        df[c] = to_dt(df[c])
df['dod_dt'] = to_dt(df['dod']) + pd.Timedelta(hours=23, minutes=59)   # date-only deaths: end of day
df['olum_zaman'] = df['deathtime'].fillna(df['dod_dt'])
df['olum_7g'] = (df['olum_zaman'] <= df['intime'] + pd.Timedelta(days=7)).astype(int)

lo_fb, hi_fb = df['sivi_dengesi_48s'].quantile([0.01, 0.99])
df['sivi_dengesi_48s_w'] = df['sivi_dengesi_48s'].clip(lo_fb, hi_fb) / 1000
df['salin_L'] = df['salin_ml'].clip(0, df['salin_ml'].quantile(.99)) / 1000
df['dengeli_L'] = df['dengeli_ml'].clip(0, df['dengeli_ml'].quantile(.99)) / 1000
df['kre_oran_48'] = df['kre_max_0_48'] / df['kre_bazal']
df['evre1_48s'] = df['evre1_48s_kdigo'].astype(int)             # KDIGO, rolling 48-h / 7-day references
df['evre1_sabit'] = ((df['kre_max_0_48'] >= 1.5 * df['kre_bazal']) |
                     (df['kre_max_0_48'] - df['kre_bazal'] >= 0.3)).astype(int)   # version-2 fixed baseline
df['n_gaz'] = df['n_0_24'] + df['n_24_48']
df['delta_be'] = df['be_ort_24_48'] - df['be_ort_0_24']
df['delta_lak'] = df['laktat_ort_24_48'] - df['laktat_ort_0_24']
df['delta_siki'] = df['abe_ort_24_48'] - df['abe_ort_0_24_siki']
df['delta_ilk_son'] = df['abe_son_g2'] - df['abe_ilk_g1']
df['delta_klor'] = df['klor_g2'] - df['klor_g1']
df['delta_na'] = df['sodyum_g2'] - df['sodyum_g1']
df['kre_bazal_saat'] = (df['intime'] - df['kre_bazal_zaman']).dt.total_seconds() / 3600
df['y'] = df['sonlanim_aki23_rrt'].astype(int)
df['y_evre3'] = df['sonlanim_evre3'].astype(int)
df['y_sabit'] = df['sonlanim_sabit'].astype(int)
# Single chronological rule for all analyses: an outcome time-stamped after the recorded
# time of death is not an outcome (the patient is classified as dying without the outcome).
olum_once = df['olum_zaman'].notna() & df['aki_zaman'].notna() & (df['aki_zaman'] > df['olum_zaman'])
N_POST = int(olum_once.sum())
df.loc[olum_once, ['y', 'y_evre3', 'kre_sonlanim', 'krt_sonlanim']] = 0
df.loc[olum_once, 'aki_zaman'] = pd.NaT
krt_win = (df['rrt_ilk'] > df['intime'] + pd.Timedelta(hours=48)) & (df['rrt_ilk'] <= df['intime'] + pd.Timedelta(days=7))
sabit_t = pd.concat([df['aki_kre_zaman'], df['rrt_ilk'].where(krt_win)], axis=1).min(axis=1)
sabit_post = (df['y_sabit'] == 1) & df['olum_zaman'].notna() & sabit_t.notna() & (sabit_t > df['olum_zaman'])
df.loc[sabit_post, 'y_sabit'] = 0
df['y_bilesik'] = ((df['y'] == 1) | (df['olum_7g'] == 1)).astype(int)


def alt_outcome(base):
    """Outcome and early-AKI flag recomputed with an alternative baseline creatinine."""
    late = df['kre_max_48s_7g']
    early = df['kre_max_0_48']
    y_ = ((df['krt_sonlanim'] == 1) | (late >= 2 * base) | ((late >= 4.0) & (late >= base + 0.3))).astype(int)
    excl = (early >= 2 * base) | ((early >= 4.0) & (early >= base + 0.3))
    return y_, excl


df['kre_bazal_a'] = df['kre_son_onceki'].fillna(df['kre_ilk'])
df['y_bazal_a'], df['dis_bazal_a'] = alt_outcome(df['kre_bazal_a'])
df['kre_bazal_b'] = df['kre_uzak_medyan']
df['y_bazal_b'], df['dis_bazal_b'] = alt_outcome(df['kre_bazal_b'])
df.loc[olum_once, ['y_bazal_a', 'y_bazal_b']] = 0

N = len(df)
H('1. COHORT')
P(f"Cohort size: {N} | primary events: {int(df['y'].sum())} ({df['y'].mean()*100:.1f}%)")
P(f"Outcome events time-stamped after the recorded time of death, reclassified as no outcome in all analyses: {N_POST} "
  f"(fixed-baseline definition: {int(sabit_post.sum())})")


# =====================================================================
# 2. DESCRIPTIVE STATISTICS
# =====================================================================
H('2. DESCRIPTIVE STATISTICS (Tables 1 and 2)')
cols = ['yas', 'cinsiyet', 'kabul_tipi', 'yb_tipi', 'kalp_cerrahisi', 'cerrahi_servis',
        'kbh', 'diyabet', 'kky', 'siroz', 'sepsis_icd',
        'kre_bazal', 'bazal_kaynak', 'kre_ilk',
        'abe_ort_0_24', 'abe_min_0_24', 'abe_ort_24_48', 'delta_abe', 'laktat_ort_0_24',
        'be_ort_0_24', 'n_0_24', 'n_24_48',
        'sofa_nonrenal', 'kre_oran_48', 'evre1_48s',
        'vazopressor_48s', 'imv_48s', 'map_min_24s', 'sivi_dengesi_48s_w', 'idrar_ml_kg_saat',
        'salin_L', 'dengeli_L', 'klor_g1', 'delta_klor',
        'trombosit_min', 'bilirubin_max', 'bun_max', 'bikarbonat_min', 'hemoglobin_min',
        'vankomisin_iv', 'pip_tazo', 'aminoglikozid', 'nsaid', 'loop_diuretik',
        'bikarbonat_inf', 'eritrosit_tx', 'taze_plazma_tx', 'asetazolamid',
        'olum_7g', 'hospital_expire_flag']
cat = ['cinsiyet', 'kabul_tipi', 'yb_tipi', 'kalp_cerrahisi', 'cerrahi_servis', 'kbh', 'diyabet',
       'kky', 'siroz', 'sepsis_icd', 'bazal_kaynak', 'evre1_48s', 'vazopressor_48s', 'imv_48s',
       'vankomisin_iv', 'pip_tazo', 'aminoglikozid', 'nsaid', 'loop_diuretik', 'bikarbonat_inf',
       'eritrosit_tx', 'taze_plazma_tx', 'asetazolamid', 'olum_7g', 'hospital_expire_flag']
try:
    from tableone import TableOne
    t1 = TableOne(df, columns=cols, categorical=cat, nonnormal=[c for c in cols if c not in cat],
                  groupby='y', pval=True, smd=True, missing=True)
    P(t1.tabulate(tablefmt='github'))
except Exception as e:
    P('TableOne not available or failed:', e)

for v in ['bikarbonat_inf', 'eritrosit_tx', 'taze_plazma_tx', 'asetazolamid', 'aminoglikozid', 'nsaid']:
    a, b = df.loc[df['y'] == 0, v], df.loc[df['y'] == 1, v]
    tab = [[(a == 0).sum(), (a == 1).sum()], [(b == 0).sum(), (b == 1).sum()]]
    test = 'Fisher' if (chi2_contingency(tab)[3] < 5).any() else 'chi2'
    p = fisher_exact(tab)[1] if test == 'Fisher' else chi2_contingency(tab)[1]
    P(f"{v}: {int(df[v].sum())} ({df[v].mean()*100:.1f}) | no outcome {int(a.sum())} ({a.mean()*100:.1f}) | "
      f"outcome {int(b.sum())} ({b.mean()*100:.1f}) | P = {pfmt(p)} [{test}]")

P('\n--- Exposure measurement ---')
P('Paired measurements, day 1, median (IQR):', med_iqr(df['n_0_24'], 0), '| day 2:', med_iqr(df['n_24_48'], 0))
P('At least 2 measurements on both days:', int(((df['n_0_24'] >= 2) & (df['n_24_48'] >= 2)).sum()))
P('Day-1 measurement within the strict 0-24 h window:', int((df['n_0_24_siki'] >= 1).sum()))

P('\n--- Baseline creatinine ---')
pre = df['bazal_kaynak'] == 'onceki_7gun_min'
P(f"Pre-ICU baseline: {int(pre.sum())} ({pre.mean()*100:.1f}%); hours before ICU admission, median (IQR): "
  f"{med_iqr(df.loc[pre, 'kre_bazal_saat'])}")
P(f"Remote baseline (median 8-365 days before ICU) available: {int(df['kre_uzak_medyan'].notna().sum())}")

P('\n--- Outcome ascertainment (48 h to day 7) ---')
P(f"At least one creatinine after 48 h: {int((df['n_kre_48s_7g'] >= 1).sum())} "
  f"({(df['n_kre_48s_7g'] >= 1).mean()*100:.1f}%); measurements per patient, median (IQR): {med_iqr(df['n_kre_48s_7g'], 0)}")
d7 = df['intime'] + pd.Timedelta(days=7)
P(f"Died by day 7: {int(df['olum_7g'].sum())} | discharged alive from hospital by day 7: "
  f"{int(((df['dischtime'] <= d7) & (df['olum_7g'] == 0)).sum())}")
P(f"Events by component: creatinine only {int(((df['kre_sonlanim'] == 1) & (df['krt_sonlanim'] == 0)).sum())}, "
  f"KRT only {int(((df['kre_sonlanim'] == 0) & (df['krt_sonlanim'] == 1)).sum())}, "
  f"both {int(((df['kre_sonlanim'] == 1) & (df['krt_sonlanim'] == 1)).sum())}")
ev_day = np.floor((df.loc[df['y'] == 1, 'aki_zaman'] - df.loc[df['y'] == 1, 'intime']).dt.total_seconds() / 86400) + 1
P('Events by ICU day:', ev_day.value_counts().sort_index().astype(int).to_dict())
P(f"Version-2 fixed-baseline definition in this cohort: {int(df['y_sabit'].sum())} events; both definitions "
  f"{int(((df['y'] == 1) & (df['y_sabit'] == 1)).sum())}, KDIGO rolling only {int(((df['y'] == 1) & (df['y_sabit'] == 0)).sum())}, "
  f"fixed only {int(((df['y'] == 0) & (df['y_sabit'] == 1)).sum())}")
P(f"Events meeting the acute absolute criterion (>= 4.0 mg/dL with acute rise): {int(((df['y'] == 1) & (df['mutlak_akut_flag'] == 1)).sum())}")
P('\n--- Phenotype by outcome definition (final cohort) ---')
grp = np.select([(df['y'] == 1) & (df['y_sabit'] == 1), (df['y'] == 1) & (df['y_sabit'] == 0),
                 (df['y'] == 0) & (df['y_sabit'] == 1)], ['both', 'rolling_only', 'fixed_only'], 'neither')
for g_ in ['both', 'rolling_only', 'fixed_only', 'neither']:
    s_ = df[grp == g_]
    P(f"{g_}: n = {len(s_)}; baseline creatinine {med_iqr(s_['kre_bazal'], 2)}; admission creatinine {med_iqr(s_['kre_ilk'], 2)}; "
      f"peak creatinine 48 h-day 7 {med_iqr(s_['kre_max_48s_7g'], 2)}; peak/baseline {med_iqr(s_['kre_max_48s_7g'] / s_['kre_bazal'], 2)}; "
      f"fluid balance (L) {med_iqr(s_['sivi_dengesi_48s_w'])}; loop diuretics {s_['loop_diuretik'].mean()*100:.1f}%; "
      f"KRT {s_['krt_sonlanim'].mean()*100:.1f}%; dABE {med_iqr(s_['delta_abe'])}; in-hospital mortality {s_['hospital_expire_flag'].mean()*100:.1f}%")
P('\n--- KDIGO definition reconciliation (abe.kohort_genis) ---')
rec = q(read_sql('09_kdigo_reconciliation.sql'))
for c_ in rec.columns:
    P(f"{c_}: {int(rec[c_].iloc[0])}")
P(f"Secondary outcome (KDIGO stage 3 or KRT): {int(df['y_evre3'].sum())} ({df['y_evre3'].mean()*100:.1f}%)")
for g in [0, 1]:
    s_ = df[df['y'] == g]
    P(f"KDIGO stage 1 criteria by 48 h, outcome={g}: {int(s_['evre1_48s'].sum())} ({s_['evre1_48s'].mean()*100:.1f}%)")
P(f"PaO2/FiO2 not available: {int(df['pf_min'].isna().sum())} | GCS not assessable: {int(df['gcs_min'].isna().sum())}")


# =====================================================================
# 3. ANALYSIS DATASET AND MULTIPLE IMPUTATION
# =====================================================================
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
d0['bazal_onceki'] = pre.astype(int)
d0['cvicu'] = (d0['yb_grup'] == 'CVICU').astype(int)
d0['log_kre_bazal'] = np.log(d0['kre_bazal'].clip(lower=0.2))
d0['log_kre_bazal_a'] = np.log(d0['kre_bazal_a'].clip(lower=0.2))
d0['log_kre_bazal_b'] = np.log(d0['kre_bazal_b'].clip(lower=0.2))
d0['log_laktat'] = np.log(d0['laktat_ort_0_24'].clip(lower=0.3))
d0['log_bilirubin'] = np.log(d0['bilirubin_max'].clip(lower=0.1))
d0['log_bun'] = np.log(d0['bun_max'].clip(lower=1))
d0['log_albumin'] = np.log(d0['albumin_min'].clip(lower=0.5))
d0['idrar'] = d0['idrar_ml_kg_saat'].clip(0, 5)
d0['log_kre_oran_48'] = np.log(d0['kre_oran_48'])
d0['log_n_gaz'] = np.log(d0['n_gaz'])
d0['klor'] = d0['klor_max']
d0['delta_kat'] = pd.cut(d0['delta_abe'], [-np.inf, -2, 0, 2, np.inf],
                         labels=['dec_ge2', 'dec_lt2', 'rise_0_2', 'rise_gt2']).astype(str)

# Restricted cubic spline (Harrell) for dABE, knots at the 5th, 35th, 65th and 95th percentiles
KNOTS = np.percentile(d0['delta_abe'], [5, 35, 65, 95])


def rcs_basis(x, k=KNOTS):
    x = np.asarray(x, dtype=float)
    norm_ = (k[-1] - k[0]) ** 2
    cols = []
    for j in range(len(k) - 2):
        b = (np.clip(x - k[j], 0, None) ** 3
             - np.clip(x - k[-2], 0, None) ** 3 * (k[-1] - k[j]) / (k[-1] - k[-2])
             + np.clip(x - k[-1], 0, None) ** 3 * (k[-2] - k[j]) / (k[-1] - k[-2])) / norm_
        cols.append(b)
    return np.column_stack(cols)


def set_rcs(t):
    B = rcs_basis(t['delta_abe'])
    t['rcs_n1'], t['rcs_n2'] = B[:, 0], B[:, 1]
    return t


d0 = set_rcs(d0)

IMP_COLS = ['yas', 'log_kre_bazal', 'log_laktat', 'abe_ort_0_24', 'delta_abe', 'sivi_dengesi_48s_w',
            'idrar', 'map_min_24s', 'trombosit_min', 'log_bilirubin', 'log_bun', 'hemoglobin_min',
            'klor', 'delta_klor', 'delta_na', 'sofa_nonrenal', 'log_kre_oran_48', 'y']

COV = ("yas + erkek + C(kabul_grup) + C(yb_grup) + log_kre_bazal + bazal_onceki"
       " + log_laktat + vazopressor_48s + imv_48s + sivi_dengesi_48s_w + idrar"
       " + map_min_24s + trombosit_min + log_bilirubin + log_bun + hemoglobin_min"
       " + vankomisin_iv + pip_tazo + loop_diuretik + bikarbonat_inf + eritrosit_tx"
       " + kbh + diyabet + kky + siroz + sepsis_icd")
PARS = "yas + erkek + log_kre_bazal + log_laktat + sofa_nonrenal"


def drop_term(c, term):
    return c.replace(f" + {term}", "").replace(f"{term} + ", "")


# All analysis-model variables are predictors in the imputation model; only variables with missing values are imputed.
BIN_COLS = ['erkek', 'bazal_onceki', 'vazopressor_48s', 'imv_48s', 'vankomisin_iv', 'pip_tazo', 'loop_diuretik',
            'bikarbonat_inf', 'eritrosit_tx', 'kbh', 'diyabet', 'kky', 'siroz', 'sepsis_icd', 'kalp_cerrahisi']
DUM = pd.get_dummies(d0[['kabul_grup', 'yb_grup']], drop_first=True).astype(float)
imps = []
for i in range(N_IMP):
    t = d0.copy()
    X_imp = pd.concat([t[IMP_COLS], t[BIN_COLS].astype(float), DUM], axis=1)
    X_done = IterativeImputer(random_state=i, sample_posterior=True, max_iter=20).fit_transform(X_imp)
    t[IMP_COLS] = X_done[:, :len(IMP_COLS)]
    t['y'] = d0['y']
    imps.append(t)
H('3. MULTIPLE IMPUTATION')
P(f"{N_IMP} imputed datasets (chained equations, Bayesian ridge, posterior sampling, 20 iterations); "
  f"predictors: {len(IMP_COLS)} continuous/outcome variables, {len(BIN_COLS)} binary covariates and {DUM.shape[1]} category indicators.")
P('Missing values before imputation:', {c: int(d0[c].isna().sum()) for c in IMP_COLS if d0[c].isna().sum() > 0})
mis = d0['log_bilirubin'].isna()
P('Bilirubin (log), observed mean vs mean of imputed values:',
  f"{d0['log_bilirubin'].mean():.3f} vs {np.mean([t.loc[mis, 'log_bilirubin'].mean() for t in imps]):.3f}")


def fit(formula, data):
    try:
        return smf.logit(formula, data).fit(disp=0, method='newton', maxiter=100)
    except Exception:
        return smf.logit(formula, data).fit(disp=0, method='bfgs', maxiter=3000)


def pool(formula, terms, subset=None, data=None):
    """Rubin's rules for single coefficients."""
    data = imps if data is None else data
    b = {k: [] for k in terms}; w = {k: [] for k in terms}
    for t in data:
        x = t if subset is None else t[subset(t)]
        m = fit(formula, x)
        for k in terms:
            b[k].append(m.params[k]); w[k].append(m.bse[k] ** 2)
    out = {}
    for k in terms:
        bb, ww = np.array(b[k]), np.array(w[k])
        qb = bb.mean(); se = np.sqrt(ww.mean() + (1 + 1 / len(data)) * bb.var(ddof=1))
        out[k] = (qb, se, 2 * norm.sf(abs(qb / se)))
    return out


def pool_d1(formula, terms, subset=None, data=None):
    """Multivariate Wald test pooled across imputations (D1 statistic)."""
    data = imps if data is None else data
    Q, U = [], []
    for t in data:
        x = t if subset is None else t[subset(t)]
        m = fit(formula, x)
        Q.append(m.params[terms].values); U.append(m.cov_params().loc[terms, terms].values)
    Q, U = np.array(Q), np.array(U)
    m_, k = len(Q), len(terms)
    qb, Ub = Q.mean(0), U.mean(0)
    B = np.atleast_2d(np.cov(Q.T, ddof=1))
    T = Ub + (1 + 1 / m_) * B
    r = max((1 + 1 / m_) * np.trace(B @ np.linalg.inv(Ub)) / k, 1e-8)
    D1 = float(qb @ np.linalg.inv(Ub) @ qb) / (k * (1 + r))
    t_ = k * (m_ - 1)
    nu = 4 + (t_ - 4) * (1 + (1 - 2 / t_) / r) ** 2 if t_ > 4 else t_ * (1 + 1 / k) * (1 + 1 / r) ** 2 / 2
    return dict(D1=D1, k=k, nu=nu, p=f_dist.sf(D1, k, nu), q=qb, T=T)


def fmt(q_, se, p, dig=3):
    return f"{np.exp(q_):.{dig}f} ({np.exp(q_-1.96*se):.{dig}f}-{np.exp(q_+1.96*se):.{dig}f}), P = {pfmt(p)}"


# =====================================================================
# 4. PRIMARY ANALYSIS (Table 3)
# =====================================================================
H('4. PRIMARY ANALYSIS (Table 3)')
F_MAIN = f"y ~ abe_ort_0_24 + delta_abe + {COV}"
m_un = smf.logit("y ~ delta_abe", d0).fit(disp=0)
P('Unadjusted, per 1 mmol/L increase:', fmt(m_un.params['delta_abe'], m_un.bse['delta_abe'], m_un.pvalues['delta_abe']))
prim = pool(F_MAIN, ['delta_abe', 'abe_ort_0_24'])
P('Adjusted, per 1 mmol/L increase:', fmt(*prim['delta_abe']))
q_, se_, p_ = prim['delta_abe']
P('Adjusted, per 1 mmol/L decrease (reciprocal of the above):', fmt(-q_, se_, p_))
P('Day-1 ABE in the same model:', fmt(*prim['abe_ort_0_24']))
REF = "C(delta_kat, Treatment(reference='rise_0_2'))"
catr = pool(f"y ~ abe_ort_0_24 + {REF} + {COV}", [f"{REF}[T.{k}]" for k in ['dec_ge2', 'dec_lt2', 'rise_gt2']])
for k in ['dec_ge2', 'dec_lt2', 'rise_0_2', 'rise_gt2']:
    s_ = d0[d0['delta_kat'] == k]
    est = 'reference' if k == 'rise_0_2' else fmt(*catr[f"{REF}[T.{k}]"], dig=2)
    P(f"Category {k}: n = {len(s_)}, events = {int(s_['y'].sum())} ({s_['y'].mean()*100:.1f}%), OR {est}")


# =====================================================================
# 5. INCREMENTAL VALUE (pooled across all imputations)
# =====================================================================
def midrank(x):
    J = np.argsort(x); Z = x[J]; n = len(x); T = np.zeros(n); i = 0
    while i < n:
        j = i
        while j < n and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(n); T2[J] = T
    return T2


def delong(y, p1, p2):
    """Paired DeLong test: AUCs, difference (p2 - p1) and its variance."""
    y = np.asarray(y)
    order = np.argsort(-y, kind='mergesort'); m = int(y.sum()); n = len(y) - m
    preds = np.vstack((p1, p2))[:, order]
    tx = np.array([midrank(r[:m]) for r in preds]); ty = np.array([midrank(r[m:]) for r in preds])
    tz = np.array([midrank(r) for r in preds])
    aucs = tz[:, :m].sum(1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n; v10 = 1.0 - (tz[:, m:] - ty) / m
    S = np.cov(v01) / m + np.cov(v10) / n
    return aucs, aucs[1] - aucs[0], S[0, 0] + S[1, 1] - 2 * S[0, 1]


def compare_models(f_small, f_big, label):
    a1, a2, dd, vv, b1, b2 = [], [], [], [], [], []
    for t in imps:
        ms, mb = fit(f_small, t), fit(f_big, t)
        ps, pb = ms.predict(t), mb.predict(t)
        aucs, d_, v_ = delong(t['y'].values, ps.values, pb.values)
        a1.append(aucs[0]); a2.append(aucs[1]); dd.append(d_); vv.append(v_)
        b1.append(np.mean((t['y'] - ps) ** 2)); b2.append(np.mean((t['y'] - pb) ** 2))
    dd, vv = np.array(dd), np.array(vv)
    qb = dd.mean(); se = np.sqrt(vv.mean() + (1 + 1 / len(dd)) * dd.var(ddof=1))
    P(f"{label}: AUROC {np.mean(a1):.3f} -> {np.mean(a2):.3f}; difference {qb:.4f} "
      f"(95% CI {qb-1.96*se:.5f} to {qb+1.96*se:.5f}), P = {pfmt(2*norm.sf(abs(qb/se)))}; "
      f"Brier {np.mean(b1):.4f} -> {np.mean(b2):.4f}")


def optimism(formula, t, B, rng):
    """Bootstrap optimism-corrected AUROC and calibration slope (Harrell)."""
    m = fit(formula, t); p_app = m.predict(t); auc_app = roc_auc_score(t['y'], p_app)
    opt_auc, slopes = [], []
    for _ in range(B):
        tb = t.iloc[rng.integers(0, len(t), len(t))]
        try:
            mb = fit(formula, tb)
        except Exception:
            continue
        pb_b = mb.predict(tb); pb_o = mb.predict(t).clip(1e-6, 1 - 1e-6)
        opt_auc.append(roc_auc_score(tb['y'], pb_b) - roc_auc_score(t['y'], pb_o))
        lp = np.log(pb_o / (1 - pb_o))
        slopes.append(smf.logit('y ~ lp', pd.DataFrame({'y': t['y'].values, 'lp': lp.values})).fit(disp=0).params['lp'])
    return auc_app, auc_app - np.mean(opt_auc), np.mean(slopes)


H('5. INCREMENTAL VALUE (Supplementary Table S1)')
F_M0 = f"y ~ {COV}"
F_M1 = f"y ~ abe_ort_0_24 + {COV}"
F_C = f"y ~ abe_ort_0_24 + {COV} + log_kre_oran_48"
F_CD = f"y ~ abe_ort_0_24 + delta_abe + {COV} + log_kre_oran_48"
P('Wald test, day-1 ABE added to covariates:', fmt(*pool(F_M1, ['abe_ort_0_24'])['abe_ort_0_24']))
P('Wald test, dABE added to covariates + day-1 ABE:', fmt(*prim['delta_abe']))
P('Wald test, dABE added to covariates + day-1 ABE + 48-h creatinine change:',
  fmt(*pool(F_CD, ['delta_abe'])['delta_abe']))
compare_models(F_M0, F_M1, 'Covariates -> + day-1 ABE')
compare_models(F_M1, F_MAIN, 'Covariates + day-1 ABE -> + dABE')
compare_models(F_C, F_CD, 'Covariates + day-1 ABE + 48-h creatinine change -> + dABE')
rng = np.random.default_rng(SEED)
for lab, f_ in [('Creatinine model', F_C), ('Creatinine model + dABE', F_CD)]:
    app, corr, slope = optimism(f_, imps[0], max(N_BOOT // 2, 20), rng)
    P(f"{lab} (bootstrap, imputation 1): apparent AUROC {app:.3f}, optimism-corrected {corr:.3f}, "
      f"calibration slope {slope:.3f}")
lvl = pool(f"y ~ abe_ort_24_48 + {COV}", ['abe_ort_24_48'])
both = pool(f"y ~ abe_ort_0_24 + abe_ort_24_48 + {COV}", ['abe_ort_0_24', 'abe_ort_24_48'])
P('Day-2 ABE alone:', fmt(*lvl['abe_ort_24_48']))
P('Both days, day-1 ABE:', fmt(*both['abe_ort_0_24']), '| day-2 ABE:', fmt(*both['abe_ort_24_48']))


# =====================================================================
# 6. NON-LINEARITY AND ABSOLUTE RISKS (Figure 2)
# =====================================================================
H('6. NON-LINEARITY AND ABSOLUTE RISKS (Figure 2)')
F_SPL = f"y ~ abe_ort_0_24 + delta_abe + rcs_n1 + rcs_n2 + {COV}"
P('Knots (mmol/L):', np.round(KNOTS, 2).tolist())
nl = pool_d1(F_SPL, ['rcs_n1', 'rcs_n2'])
P(f"Non-linearity (pooled D1 Wald test): F = {nl['D1']:.2f}, df = {nl['k']}, {min(nl['nu'], 99999):.0f}; P = {pfmt(nl['p'])}")
cv = pool_d1(F_SPL, ['delta_abe', 'rcs_n1', 'rcs_n2'])
P(f"Overall association, spline model (pooled D1): F = {cv['D1']:.2f}, P = {pfmt(cv['p'])}")


def basis3(x):
    x = np.atleast_1d(np.asarray(x, dtype=float))
    return np.column_stack([x, rcs_basis(x)])


grid_lo, grid_hi = np.percentile(d0['delta_abe'], [2, 98])
grid = np.linspace(grid_lo, grid_hi, 200)
Dg = basis3(grid) - basis3(0.0)
lor = Dg @ cv['q']; lse = np.sqrt(np.einsum('ij,jk,ik->i', Dg, cv['T'], Dg))
for v in [-4, -2, 2, 4]:
    D_ = (basis3(v) - basis3(0.0))[0]
    e_ = float(D_ @ cv['q']); s_ = float(np.sqrt(D_ @ cv['T'] @ D_))
    P(f"dABE {v:+d} vs 0: OR {np.exp(e_):.2f} ({np.exp(e_-1.96*s_):.2f}-{np.exp(e_+1.96*s_):.2f})")


def gcomp(m, t, vals):
    out = []
    for v in vals:
        tt = t.copy(); tt['delta_abe'] = v; tt = set_rcs(tt)
        out.append(m.predict(tt).mean() * 100)
    return out


VALS = [-4, 0, 7]
risk_pool = np.mean([gcomp(fit(F_SPL, t), t, VALS) for t in imps], axis=0)
boot = []
t0 = imps[0]
for b_ in range(N_BOOT):
    tb = t0.iloc[rng.integers(0, len(t0), len(t0))]
    try:
        boot.append(gcomp(fit(F_SPL, tb), tb, VALS))
    except Exception:
        pass
boot = np.array(boot)
for i, v in enumerate(VALS):
    lo_, hi_ = np.percentile(boot[:, i], [2.5, 97.5])
    P(f"Average adjusted risk at dABE {v:+d}: {risk_pool[i]:.1f}% (95% CI {lo_:.1f}-{hi_:.1f})")
for i, v in [(0, -4), (2, 7)]:
    dif = boot[:, i] - boot[:, 1]
    P(f"Risk difference dABE {v:+d} vs 0: {risk_pool[i]-risk_pool[1]:.1f} percentage points "
      f"(95% CI {np.percentile(dif, 2.5):.1f} to {np.percentile(dif, 97.5):.1f})")
P(f"Patients below -4 mmol/L: {int((d0['delta_abe'] < -4).sum())}; above +7 mmol/L: {int((d0['delta_abe'] > 7).sum())}")


# =====================================================================
# 7. SENSITIVITY ANALYSES (Table 4)
# =====================================================================
H('7. SENSITIVITY ANALYSES (Table 4; OR per 1 mmol/L increase in dABE)')
ALL_ADD = "sofa_nonrenal + log_kre_oran_48 + klor + delta_klor + salin_L + dengeli_L + log_n_gaz"
COV_NL = (COV.replace("yas +", "cr(yas, df=3, constraints='center') +", 1)
             .replace("log_kre_bazal +", "cr(log_kre_bazal, df=3, constraints='center') +", 1)
             .replace("log_laktat +", "cr(log_laktat, df=3, constraints='center') +", 1)
             .replace("map_min_24s +", "cr(map_min_24s, df=3, constraints='center') +", 1)
             .replace("sivi_dengesi_48s_w +", "cr(sivi_dengesi_48s_w, df=3, constraints='center') +", 1))
STAGE_A = ("yas + erkek + C(kabul_grup) + C(yb_grup) + log_kre_bazal + bazal_onceki"
           " + kbh + diyabet + kky + siroz + sepsis_icd")
STAGE_B = STAGE_A + " + log_laktat + map_min_24s + trombosit_min + log_bilirubin + log_bun + hemoglobin_min + sofa_nonrenal"
COV_A = COV.replace('log_kre_bazal', 'log_kre_bazal_a')
COV_B = drop_term(COV.replace('log_kre_bazal', 'log_kre_bazal_b'), 'bazal_onceki')

scen = [
    ('Primary analysis', F_MAIN, None, 'y', 'delta_abe'),
    ('[Adjustment] + non-renal SOFA', F_MAIN + " + sofa_nonrenal", None, 'y', 'delta_abe'),
    ('[Adjustment] Parsimonious model', f"y ~ abe_ort_0_24 + delta_abe + {PARS}", None, 'y', 'delta_abe'),
    ('[Adjustment] + creatinine change 0-48 h', F_CD, None, 'y', 'delta_abe'),
    ('[Adjustment] + day-1 chloride', F_MAIN + " + klor", None, 'y', 'delta_abe'),
    ('[Adjustment] + change in chloride', F_MAIN + " + delta_klor", None, 'y', 'delta_abe'),
    ('[Adjustment] + change in sodium', F_MAIN + " + delta_na", None, 'y', 'delta_abe'),
    ('[Adjustment] + saline and balanced crystalloid volumes', F_MAIN + " + salin_L + dengeli_L", None, 'y', 'delta_abe'),
    ('[Adjustment] + number of blood gases', F_MAIN + " + log_n_gaz", None, 'y', 'delta_abe'),
    ('[Adjustment] + all additional variables', F_MAIN + " + " + ALL_ADD, None, 'y', 'delta_abe'),
    ('[Staged] A: baseline characteristics', f"y ~ abe_ort_0_24 + delta_abe + {STAGE_A}", None, 'y', 'delta_abe'),
    ('[Staged] B: A + day-1 severity', f"y ~ abe_ort_0_24 + delta_abe + {STAGE_B}", None, 'y', 'delta_abe'),
    ('[Model] Urine output removed', f"y ~ abe_ort_0_24 + delta_abe + {drop_term(COV, 'idrar')}", None, 'y', 'delta_abe'),
    ('[Model] Non-linear continuous covariates', f"y ~ abe_ort_0_24 + delta_abe + {COV_NL}", None, 'y', 'delta_abe'),
    ('[Exposure] Strict day-1 window (0-24 h)', f"y ~ abe_ort_0_24_siki + delta_siki + {COV}",
     lambda t: t['n_0_24_siki'] >= 1, 'y', 'delta_siki'),
    ('[Exposure] >= 2 measurements on each day', F_MAIN,
     lambda t: (t['n_0_24'] >= 2) & (t['n_24_48'] >= 2), 'y', 'delta_abe'),
    ('[Exposure] First day-1 to last day-2 measurement', f"y ~ abe_ilk_g1 + delta_ilk_son + {COV}",
     None, 'y', 'delta_ilk_son'),
    ('[Cohort] No sodium bicarbonate', f"y ~ abe_ort_0_24 + delta_abe + {drop_term(COV, 'bikarbonat_inf')}",
     lambda t: t['bikarbonat_inf'] == 0, 'y', 'delta_abe'),
    ('[Cohort] No loop diuretics', f"y ~ abe_ort_0_24 + delta_abe + {drop_term(COV, 'loop_diuretik')}",
     lambda t: t['loop_diuretik'] == 0, 'y', 'delta_abe'),
    ('[Cohort] Arterial samples only', F_MAIN, lambda t: t['n_venoz'] == 0, 'y', 'delta_abe'),
    ('[Cohort] Complete non-renal SOFA data (+ SOFA)', F_MAIN + " + sofa_nonrenal",
     lambda t: t['pf_min'].notna() & t['gcs_min'].notna(), 'y', 'delta_abe'),
    ('[Cohort] Albumin available (+ albumin)', F_MAIN + " + log_albumin",
     lambda t: t['log_albumin'].notna(), 'y', 'delta_abe'),
    ('[Baseline] Pre-ICU baseline creatinine only', f"y ~ abe_ort_0_24 + delta_abe + {drop_term(COV, 'bazal_onceki')}",
     lambda t: t['bazal_onceki'] == 1, 'y', 'delta_abe'),
    ('[Baseline] Last pre-ICU value', f"y_bazal_a ~ abe_ort_0_24 + delta_abe + {COV_A}",
     lambda t: ~t['dis_bazal_a'].astype(bool), 'y_bazal_a', 'delta_abe'),
    ('[Baseline] Median value 8-365 days before ICU', f"y_bazal_b ~ abe_ort_0_24 + delta_abe + {COV_B}",
     lambda t: t['kre_bazal_b'].notna() & ~t['dis_bazal_b'].astype(bool), 'y_bazal_b', 'delta_abe'),
    ('[Outcome] Composite: primary outcome or death within 7 days', f"y_bilesik ~ abe_ort_0_24 + delta_abe + {COV}",
     None, 'y_bilesik', 'delta_abe'),
    ('[Outcome] KDIGO stage 3 or KRT', f"y_evre3 ~ abe_ort_0_24 + delta_abe + {COV}", None, 'y_evre3', 'delta_abe'),
    ('[Outcome] Version-2 fixed-baseline definition', f"y_sabit ~ abe_ort_0_24 + delta_abe + {COV}",
     None, 'y_sabit', 'delta_abe'),
    ('[Cohort] No stage 1 AKI by 48 h (fixed-baseline definition)', F_MAIN,
     lambda t: t['evre1_sabit'] == 0, 'y', 'delta_abe'),
]
for lab, f_, sub, outc, term in scen:
    x = d0 if sub is None else d0[sub(d0)]
    try:
        r = pool(f_, [term], sub)
        P(f"SENS | {lab} | n = {len(x)} | events = {int(x[outc].sum())} | {fmt(*r[term])}")
    except Exception as e:
        P(f"SENS | {lab} | FAILED: {str(e)[:150]}")

cc = d0.dropna(subset=['map_min_24s', 'trombosit_min', 'log_bilirubin', 'log_bun', 'hemoglobin_min', 'idrar'])
m_cc = fit(F_MAIN, cc)
P(f"SENS | [Model] Complete-case analysis | n = {len(cc)} | events = {int(cc['y'].sum())} | "
  f"{fmt(m_cc.params['delta_abe'], m_cc.bse['delta_abe'], m_cc.pvalues['delta_abe'])}")

for lab, yc, ex, need in [('last pre-ICU value', 'y_bazal_a', 'dis_bazal_a', 'kre_bazal_a'),
                          ('median 8-365 days before ICU', 'y_bazal_b', 'dis_bazal_b', 'kre_bazal_b')]:
    s_ = d0[d0[need].notna() & ~d0[ex].astype(bool)]
    P(f"Baseline = {lab}: n = {len(s_)}; events {int(s_[yc].sum())}; both definitions "
      f"{int(((s_['y'] == 1) & (s_[yc] == 1)).sum())}, primary only {int(((s_['y'] == 1) & (s_[yc] == 0)).sum())}, "
      f"alternative only {int(((s_['y'] == 0) & (s_[yc] == 1)).sum())}; excluded as early AKI under this baseline "
      f"{int((d0[need].notna() & d0[ex].astype(bool)).sum())}")

H('7b. COMPETING RISK OF DEATH')
try:
    from lifelines import CoxPHFitter, AalenJohansenFitter
    start = d0['intime'] + pd.Timedelta(hours=48)
    end7 = d0['intime'] + pd.Timedelta(days=7)
    end = d0['olum_zaman'].where(d0['olum_zaman'] < end7, end7).fillna(end7)
    ev_t = d0['aki_zaman'].where(d0['y'] == 1)
    aki_obs = ev_t.notna() & (ev_t <= end)          # events at the same time as death count as events
    P(f"Primary events: {int((d0['y'] == 1).sum())}; included as events in the time-to-event analysis: {int(aki_obs.sum())}; "
      f"events after the recorded time of death: {int(((d0['y'] == 1) & ~aki_obs).sum())}; "
      f"events without a time stamp: {int(((d0['y'] == 1) & d0['aki_zaman'].isna()).sum())}")
    stop = ev_t.where(aki_obs, end)
    dur = ((stop - start).dt.total_seconds() / 86400).clip(lower=1e-3).values
    death_obs = (~aki_obs) & d0['olum_zaman'].notna() & (d0['olum_zaman'] <= end7)
    code = np.where(aki_obs, 1, np.where(death_obs, 2, 0))
    P(f"Follow-up from 48 h to day 7: AKI/KRT {int((code == 1).sum())}, death without prior AKI/KRT "
      f"{int((code == 2).sum())}, censored at day 7 {int((code == 0).sum())}")
    bb, ww = [], []
    for t in imps:
        X = patsy.dmatrix(f"abe_ort_0_24 + delta_abe + {COV}", t, return_type='dataframe').drop(columns='Intercept')
        names = list(X.columns); X.columns = [f"v{j}" for j in range(X.shape[1])]
        key = f"v{names.index('delta_abe')}"
        X['T'] = dur; X['E'] = (code == 1).astype(int)
        cph = CoxPHFitter(penalizer=0.0).fit(X, 'T', 'E')
        bb.append(cph.params_[key]); ww.append(cph.standard_errors_[key] ** 2)
    bb, ww = np.array(bb), np.array(ww)
    qb = bb.mean(); se = np.sqrt(ww.mean() + (1 + 1 / len(bb)) * bb.var(ddof=1))
    P('Cause-specific HR for AKI/KRT per 1 mmol/L increase in dABE (death treated as competing event):',
      fmt(qb, se, 2 * norm.sf(abs(qb / se))))
    try:
        from lifelines.statistics import proportional_hazard_test
        ph = proportional_hazard_test(cph, X, time_transform='rank')
        P(f"Proportional hazards (scaled Schoenfeld residuals, rank time, imputation {len(imps)}): dABE P = "
          f"{pfmt(float(ph.summary.loc[key, 'p']))}; global P = {pfmt(float(ph.summary['p'].min()))} (smallest covariate P)")
    except Exception as e:
        P('PH test failed:', str(e)[:120])
    for lab_, sel_, tf_, ef_ in [('48-96 h', np.ones(len(dur), bool), np.minimum(dur, 2.0), (code == 1) & (dur <= 2.0)),
                                 ('96 h-day 7', dur > 2.0, dur - 2.0, (code == 1))]:
        bb2, ww2 = [], []
        for t in imps:
            X = patsy.dmatrix(f"abe_ort_0_24 + delta_abe + {COV}", t, return_type='dataframe').drop(columns='Intercept')
            names = list(X.columns); X.columns = [f"v{j}" for j in range(X.shape[1])]
            key2 = f"v{names.index('delta_abe')}"
            X['T'] = tf_; X['E'] = ef_.astype(int)
            X = X[sel_]
            c2 = CoxPHFitter(penalizer=0.0).fit(X, 'T', 'E')
            bb2.append(c2.params_[key2]); ww2.append(c2.standard_errors_[key2] ** 2)
        bb2, ww2 = np.array(bb2), np.array(ww2)
        q2 = bb2.mean(); se2 = np.sqrt(ww2.mean() + (1 + 1 / len(bb2)) * bb2.var(ddof=1))
        P(f"Cause-specific HR, {lab_} (n at risk {int(sel_.sum())}, events {int(ef_[sel_].sum())}): {fmt(q2, se2, 2 * norm.sf(abs(q2 / se2)))}")
    for k in ['dec_ge2', 'dec_lt2', 'rise_0_2', 'rise_gt2']:
        msk = (d0['delta_kat'] == k).values
        aj1 = AalenJohansenFitter(calculate_variance=False).fit(dur[msk], code[msk], event_of_interest=1)
        aj2 = AalenJohansenFitter(calculate_variance=False).fit(dur[msk], code[msk], event_of_interest=2)
        P(f"Cumulative incidence by day 7, {k}: AKI/KRT {aj1.cumulative_density_.iloc[-1, 0]*100:.1f}%, "
          f"death without AKI/KRT {aj2.cumulative_density_.iloc[-1, 0]*100:.1f}%")
except ImportError:
    P('lifelines not installed: run !pip install lifelines')


# =====================================================================
# 8. COMPONENTS OF ABE
# =====================================================================
H('8. COMPONENTS: dABE vs dBE vs dLactate (same covariates; all imputations)')
for lab, f_, terms in [
        ('dABE (with day-1 ABE): primary model', F_MAIN, ['delta_abe']),
        ('dBE (with day-1 BE)', f"y ~ be_ort_0_24 + delta_be + {COV}", ['delta_be']),
        ('dLactate (day-1 lactate in covariates)', f"y ~ delta_lak + {COV}", ['delta_lak']),
        ('dABE and dLactate together (with day-1 ABE)', f"y ~ abe_ort_0_24 + delta_abe + delta_lak + {COV}",
         ['delta_abe', 'delta_lak']),
        ('dBE and dLactate together (with day-1 BE)', f"y ~ be_ort_0_24 + delta_be + delta_lak + {COV}",
         ['delta_be', 'delta_lak'])]:
    r = pool(f_, terms)
    aucs, briers = [], []
    for t in imps:
        p_ = fit(f_, t).predict(t)
        aucs.append(roc_auc_score(t['y'], p_)); briers.append(np.mean((t['y'] - p_) ** 2))
    P(f"{lab}: " + '; '.join(f"{k} {fmt(*r[k])}" for k in terms) +
      f" | AUROC {np.mean(aucs):.4f} | Brier {np.mean(briers):.5f}")
compare_models(F_MAIN, f"y ~ abe_ort_0_24 + delta_abe + delta_lak + {COV}", 'Primary model -> + dLactate')
compare_models(F_MAIN, f"y ~ be_ort_0_24 + delta_be + {COV}", 'dABE model vs dBE model')


# =====================================================================
# 9. NEW-ONSET AKI vs PROGRESSION
# =====================================================================
H('9. GROUPS BY KDIGO STAGE 1 CRITERIA AT 48 h')
for g, lab in [(0, 'Group A: no stage 1 at 48 h'), (1, 'Group B: stage 1 at 48 h')]:
    s_ = d0[d0['evre1_48s'] == g]
    P(f"{lab}: n = {len(s_)}, events = {int(s_['y'].sum())} ({s_['y'].mean()*100:.1f}%)")
    for k in ['dec_ge2', 'dec_lt2', 'rise_0_2', 'rise_gt2']:
        sk = s_[s_['delta_kat'] == k]
        P(f"   {k}: n = {len(sk)}, events = {int(sk['y'].sum())} ({sk['y'].mean()*100:.1f}%)")
    sub = (lambda t, g=g: t['evre1_48s'] == g)
    for mlab, f_ in [('parsimonious (+ vasopressors, sepsis)',
                      f"y ~ abe_ort_0_24 + delta_abe + {PARS} + vazopressor_48s + sepsis_icd"),
                     ('full model', F_MAIN), ('full model + creatinine change', F_CD)]:
        try:
            P(f"   {mlab}: {fmt(*pool(f_, ['delta_abe'], sub)['delta_abe'])}")
        except Exception as e:
            P(f"   {mlab}: FAILED {str(e)[:100]}")
ri = pool(F_MAIN + " + evre1_48s + delta_abe:evre1_48s", ['delta_abe:evre1_48s'])
P('Interaction dABE x stage 1 at 48 h: P =', pfmt(ri['delta_abe:evre1_48s'][2]))


# =====================================================================
# 10. SUBGROUPS
# =====================================================================
H('10. SUBGROUPS (Bonferroni threshold 0.05/7 = 0.0071)')
subgroups = [
    ('Sepsis', 'sepsis_icd', 'sepsis_icd', ['No', 'Yes']),
    ('Chronic kidney disease', 'kbh', 'kbh', ['No', 'Yes']),
    ('Invasive mechanical ventilation', 'imv_48s', 'imv_48s', ['No', 'Yes']),
    ('Vasopressors', 'vazopressor_48s', 'vazopressor_48s', ['No', 'Yes']),
    ('Loop diuretics', 'loop_diuretik', 'loop_diuretik', ['No', 'Yes']),
    ('ICU type', 'cvicu', 'C(yb_grup)', ['Other ICUs', 'Cardiovascular ICU']),
    ('Cardiac surgery service', 'kalp_cerrahisi', None, ['No', 'Yes']),
]
overall = prim['delta_abe']
fp = []
for name, var, term, labels in subgroups:
    c_ = drop_term(COV, term) if term else COV
    p_int = pool(f"y ~ abe_ort_0_24 + delta_abe * {var} + {c_}", [f"delta_abe:{var}"])[f"delta_abe:{var}"][2]
    fp.append(dict(label=name, header=True))
    for lev in [0, 1]:
        x = d0[d0[var] == lev]
        r = pool(f"y ~ abe_ort_0_24 + delta_abe + {c_}", ['delta_abe'], lambda t, var=var, lev=lev: t[var] == lev)['delta_abe']
        fp.append(dict(label=labels[lev], header=False, n=len(x), ev=int(x['y'].sum()), q=r[0], se=r[1],
                       pint=(p_int if lev == 0 else None)))
        P(f"SUB | {name}: {labels[lev]} | n = {len(x)} | events = {int(x['y'].sum())} | {fmt(*r)}"
          + (f" | P interaction = {pfmt(p_int)}" if lev == 0 else ''))

H('10b. CARDIOVASCULAR ICU vs OTHER ICUs')
for lab, msk in [('Cardiovascular ICU', d0['cvicu'] == 1), ('Other ICUs', d0['cvicu'] == 0)]:
    s_ = d0[msk]
    P(f"{lab}: n = {len(s_)}; cardiac surgery service {s_['kalp_cerrahisi'].mean()*100:.1f}%; "
      f"baseline creatinine {med_iqr(s_['kre_bazal'], 2)}; day-1 lactate {med_iqr(s_['laktat_ort_0_24'])}; "
      f"fluid balance (L) {med_iqr(s_['sivi_dengesi_48s_w'])}; saline (L) {med_iqr(s_['salin_L'])}; "
      f"loop diuretics {s_['loop_diuretik'].mean()*100:.1f}%; sodium bicarbonate {s_['bikarbonat_inf'].mean()*100:.1f}%; "
      f"RBC transfusion {s_['eritrosit_tx'].mean()*100:.1f}%; events {s_['y'].mean()*100:.1f}%")


# =====================================================================
# 11. SELECTION (Supplementary Table S4)
# =====================================================================
H('11. PATIENTS ASSESSED DURING COHORT CONSTRUCTION (Supplementary Table S4)')
sel = q(read_sql('06_selection_comparison.sql'))
sel['male'] = (sel['gender'] == 'M').astype(int)
sel['cvicu'] = sel['first_careunit'].str.contains('CVICU', na=False).astype(int)
sel['micu'] = sel['first_careunit'].str.startswith('Medical Intensive', na=False).astype(int)
sel['tsicu'] = sel['first_careunit'].str.contains('TSICU', na=False).astype(int)


def smd(a, b, binary=False):
    a, b = pd.Series(a).dropna().astype(float), pd.Series(b).dropna().astype(float)
    if binary:
        p1, p2 = a.mean(), b.mean(); sd = np.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
    else:
        sd = np.sqrt((a.var() + b.var()) / 2)
    return (a.mean() - b.mean()) / sd if sd > 0 else np.nan


for g in sorted(sel['grup'].unique()):
    s_ = sel[sel['grup'] == g]
    P(f"{g}: n = {len(s_)}; age {med_iqr(s_['yas'], 0)}; male {s_['male'].mean()*100:.1f}%; "
      f"CVICU {s_['cvicu'].mean()*100:.1f}%; MICU {s_['micu'].mean()*100:.1f}%; TSICU {s_['tsicu'].mean()*100:.1f}%; "
      f"day-1 ABE {med_iqr(s_['abe_d1'])}; day-1 lactate {med_iqr(s_['laktat_d1'])}; ICU stay {med_iqr(s_['los'])}; "
      f"mortality {s_['hospital_expire_flag'].mean()*100:.1f}%")
gs_ = sorted(sel['grup'].unique())
g1, g2 = sel[sel['grup'] == gs_[0]], sel[sel['grup'] == gs_[1]]
P(f'SMD {gs_[0]} vs {gs_[1]}: ' + '; '.join(
    f"{v} {smd(g1[v], g2[v], b):.2f}" for v, b in
    [('yas', False), ('male', True), ('cvicu', True), ('micu', True), ('tsicu', True), ('abe_d1', False),
     ('laktat_d1', False), ('los', False), ('hospital_expire_flag', True)]))


# =====================================================================
# 12. FIGURES
# =====================================================================
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                     'font.size': 10, 'axes.linewidth': 0.8, 'pdf.fonttype': 42})


def save(fig, name):
    fig.savefig(FIG_DIR / f'{name}.png', dpi=600, bbox_inches='tight', facecolor='white')
    fig.savefig(FIG_DIR / f'{name}.pdf', bbox_inches='tight', facecolor='white')
    fig.savefig(FIG_DIR / f'{name}.tiff', dpi=600, bbox_inches='tight', facecolor='white',
                pil_kwargs={'compression': 'tiff_lzw'})


fig = plt.figure(figsize=(6.5, 5.2))
gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1], hspace=0.06)
ax = fig.add_subplot(gs[0]); axh = fig.add_subplot(gs[1], sharex=ax)
ax.fill_between(grid, np.exp(lor - 1.96 * lse), np.exp(lor + 1.96 * lse), color='#9ecae1', alpha=0.55, lw=0)
ax.plot(grid, np.exp(lor), color='#08519c', lw=2)
ax.axhline(1, color='black', lw=0.8, ls='--'); ax.axvline(0, color='gray', lw=0.6, ls=':')
ax.set_yscale('log')
yt = [v for v in [0.3, 0.4, 0.5, 0.7, 1, 1.5, 2, 3] if np.exp(lor - 1.96 * lse).min() * 0.9 <= v <= np.exp(lor + 1.96 * lse).max() * 1.1]
ax.set_yticks(yt); ax.set_yticklabels([f'{v:.1f}' for v in yt]); ax.minorticks_off()
ax.set_ylabel('Adjusted odds ratio (95% CI)'); ax.tick_params(labelbottom=False)
ymin = ax.get_ylim()[0]
ax.plot(KNOTS, [ymin] * len(KNOTS), linestyle='none', marker='^', color='dimgray', ms=5, clip_on=False)
for s in ['top', 'right']:
    ax.spines[s].set_visible(False); axh.spines[s].set_visible(False)
ax.text(0.98, 0.95, 'Reference: ΔABE = 0 mmol/L', transform=ax.transAxes, ha='right', va='top', fontsize=9, color='dimgray')
inside = d0['delta_abe'][(d0['delta_abe'] >= grid_lo) & (d0['delta_abe'] <= grid_hi)]
axh.hist(inside, bins=50, range=(grid_lo, grid_hi), color='#bdbdbd', edgecolor='white', lw=0.3)
axh.set_ylabel('Patients, n'); axh.set_xlabel('ΔABE (mmol/L), mean ABE on day 2 minus mean ABE on day 1')
save(fig, 'Figure2_dABE_OR_curve'); plt.close(fig)

rows = [dict(label='Overall', header=False, overall=True, n=len(d0), ev=int(d0['y'].sum()),
             q=overall[0], se=overall[1], pint=None)] + fp
nR = len(rows)
fig = plt.figure(figsize=(10, 0.36 * nR + 1.2))
gs = gridspec.GridSpec(1, 3, width_ratios=[1.55, 1.15, 1.2], wspace=0.02)
axL, axP, axR = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
ys = np.arange(nR)[::-1]
for a in (axL, axR):
    a.set_xlim(0, 1); a.set_ylim(-0.7, nR + 0.3); a.axis('off')
axP.set_ylim(-0.7, nR + 0.3)
Hh = nR - 0.1
axL.text(0.00, Hh, 'Subgroup', weight='bold', va='bottom')
axL.text(0.70, Hh, 'n', weight='bold', ha='right', va='bottom')
axL.text(0.95, Hh, 'Events', weight='bold', ha='right', va='bottom')
axR.text(0.05, Hh, 'OR (95% CI)', weight='bold', va='bottom')
axR.text(0.98, Hh, 'P for\ninteraction', weight='bold', ha='right', va='bottom', linespacing=1.0)
for yv, r in zip(ys, rows):
    if r['header']:
        axL.text(0.00, yv, r['label'], weight='bold', va='center'); continue
    axL.text(0.0 if r.get('overall') else 0.04, yv, r['label'], va='center', weight='bold' if r.get('overall') else 'normal')
    axL.text(0.70, yv, f"{r['n']:,}", ha='right', va='center')
    axL.text(0.95, yv, f"{r['ev']:,}", ha='right', va='center')
    o, l, u = np.exp(r['q']), np.exp(r['q'] - 1.96 * r['se']), np.exp(r['q'] + 1.96 * r['se'])
    axP.plot([l, u], [yv, yv], color='black', lw=1)
    axP.plot(o, yv, marker='D' if r.get('overall') else 's', color='black', ms=6 if r.get('overall') else 5)
    axR.text(0.05, yv, f"{o:.2f} ({l:.2f}-{u:.2f})", va='center')
    if r['pint'] is not None:
        axR.text(0.98, yv - 0.5, pfmt(r['pint']), ha='right', va='center')
axP.axvline(1, color='black', lw=0.8, ls='--'); axP.axvline(np.exp(overall[0]), color='#08519c', lw=0.8, ls=':')
xl = [np.exp(r['q'] + s * 1.96 * r['se']) for r in rows if not r['header'] for s in (-1, 1)]
axP.set_xlim(min(0.75, min(xl) - 0.02), max(1.15, max(xl) + 0.02)); axP.set_yticks([])
for s in ['top', 'right', 'left']:
    axP.spines[s].set_visible(False)
axP.set_xlabel('Adjusted OR per 1 mmol/L increase in ΔABE')
save(fig, 'Figure3_forest_subgroups'); plt.close(fig)


# =====================================================================
# 13. SOFTWARE VERSIONS
# =====================================================================
H('13. SOFTWARE')
import sys, sklearn, scipy, statsmodels, matplotlib  # noqa: E401,E402
P('Python', sys.version.split()[0], '| pandas', pd.__version__, '| numpy', np.__version__,
  '| statsmodels', statsmodels.__version__, '| scikit-learn', sklearn.__version__,
  '| patsy', patsy.__version__, '| scipy', scipy.__version__, '| matplotlib', matplotlib.__version__)
try:
    import lifelines, tableone  # noqa: E401
    P('lifelines', lifelines.__version__, '| tableone', tableone.__version__)
except Exception:
    pass
P(f'Imputations: {N_IMP}; bootstrap replicates: {N_BOOT}; seed: {SEED}')
OUT.close()
print('\nDone. Results written to results_v2_2.txt; figures saved in', FIG_DIR)
