"""
Additional sensitivity analyses (non-renal SOFA, early creatinine change,
chloride, number of blood gases, secondary outcome).

Run AFTER analysis/abe_aki_analysis.py in the same session (e.g. with
`%run analysis/abe_aki_analysis.py` in Google Colaboratory), and after
sql/05_sofa_day1.sql has created `abe.sofa_d1`.
"""
from scipy.stats import mannwhitneyu

sofa_all = client.query(f"""
SELECT stay_id, pf_min, gcs_min, sofa_resp, sofa_coag, sofa_liver, sofa_cv, sofa_cns,
       sofa_resp + sofa_coag + sofa_liver + sofa_cv + sofa_cns AS sofa_nonrenal
FROM `{PROJECT_ID}.abe.sofa_d1`
""").to_dataframe()
print('=== Non-renal SOFA ===')
print('n:', len(sofa_all), '| PaO2/FiO2 missing:', int(sofa_all['pf_min'].isna().sum()),
      '| GCS not assessable/missing:', int(sofa_all['gcs_min'].isna().sum()))
sofa = sofa_all[['stay_id', 'sofa_nonrenal']]


def add_vars(t):
    t = t.drop(columns=[c for c in ['sofa_nonrenal'] if c in t.columns]).merge(sofa, on='stay_id', how='left')
    t['log_kre_oran_48'] = np.log(t['kre_max_0_48'] / t['kre_bazal'])
    t['evre1_48s'] = ((t['kre_max_0_48'] >= 1.5 * t['kre_bazal']) |
                      (t['kre_max_0_48'] - t['kre_bazal'] >= 0.3)).astype(int)
    t['log_n_gaz'] = np.log(t['n_0_24'] + t['n_24_48'])
    t['klor'] = t['klor_max'].fillna(t['klor_max'].median())
    krt = (t['rrt_ilk'] > t['intime'] + pd.Timedelta(hours=48)) & \
          (t['rrt_ilk'] <= t['intime'] + pd.Timedelta(days=7))
    t['y_evre3'] = ((t['kre_max_48s_7g'] >= 3 * t['kre_bazal']) | krt).astype(int)
    return t


d0 = add_vars(d0)
imps = [add_vars(t) for t in imps]
pool.__globals__['imps'] = imps   # make pool() use the updated imputed datasets

print('\n=== Descriptive ===')
for name, var in [('Non-renal SOFA', 'sofa_nonrenal'), ('Creatinine ratio 0-48 h', 'log_kre_oran_48')]:
    a, b = d0.loc[d0['y'] == 0, var], d0.loc[d0['y'] == 1, var]
    tr = np.exp if var.startswith('log') else (lambda v: v)
    print(f"{name}: overall {tr(d0[var].median()):.2f} | no outcome {tr(a.median()):.2f} "
          f"[{tr(a.quantile(.25)):.2f}-{tr(a.quantile(.75)):.2f}] | outcome {tr(b.median()):.2f} "
          f"[{tr(b.quantile(.25)):.2f}-{tr(b.quantile(.75)):.2f}] | P = {mannwhitneyu(a, b).pvalue:.4f}")
for g in [0, 1]:
    s_ = d0[d0['y'] == g]
    print(f"KDIGO stage 1 by 48 h, outcome={g}: {int(s_['evre1_48s'].sum())} ({s_['evre1_48s'].mean()*100:.1f}%)")
print(f"Secondary outcome (KDIGO stage 3 or KRT): {int(d0['y_evre3'].sum())}")

print('\n=== Additional sensitivity analyses (OR per 1 mmol/L increase in dABE) ===')
scenarios = [
    ('+ non-renal SOFA', f"y ~ abe_ort_0_24 + delta_abe + {COV} + sofa_nonrenal", None, 'y'),
    ('Parsimonious model', "y ~ abe_ort_0_24 + delta_abe + yas + erkek + log_kre_bazal + log_laktat + sofa_nonrenal", None, 'y'),
    ('+ creatinine change 0-48 h', f"y ~ abe_ort_0_24 + delta_abe + {COV} + log_kre_oran_48", None, 'y'),
    ('No KDIGO stage 1 AKI by 48 h', f"y ~ abe_ort_0_24 + delta_abe + {COV}", lambda t: t['evre1_48s'] == 0, 'y'),
    ('+ day-1 chloride', f"y ~ abe_ort_0_24 + delta_abe + {COV} + klor", None, 'y'),
    ('+ number of blood gases', f"y ~ abe_ort_0_24 + delta_abe + {COV} + log_n_gaz", None, 'y'),
    ('+ all four variables', f"y ~ abe_ort_0_24 + delta_abe + {COV} + sofa_nonrenal + log_kre_oran_48 + klor + log_n_gaz", None, 'y'),
    ('Secondary outcome: KDIGO stage 3 or KRT', f"y_evre3 ~ abe_ort_0_24 + delta_abe + {COV}", None, 'y_evre3'),
]
for name, f, sub, outcome in scenarios:
    x = d0 if sub is None else d0[sub(d0)]
    r = pool(f, ['delta_abe'], sub)
    print(f"{name} | n = {len(x)} | events = {int(x[outcome].sum())} | {fmt(*r['delta_abe'])}")
