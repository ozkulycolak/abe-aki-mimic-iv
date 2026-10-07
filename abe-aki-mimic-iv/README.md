# Early changes in alactic base excess and subsequent acute kidney injury (MIMIC-IV)

[![DOI](https://zenodo.org/badge/1407257460.svg)](https://doi.org/10.5281/zenodo.23188220)

Code for the study:

> Kılıç Ö, Yılmaz Çolak Ö. *Early changes in alactic base excess and subsequent acute kidney injury in critically ill adults: a retrospective cohort study using MIMIC-IV.* (manuscript in preparation)

The study examines whether the change in alactic base excess (ABE = base excess + lactate) between day 1 (6 h before to 24 h after ICU admission) and day 2 (24 to 48 h) (ΔABE) is associated with KDIGO stage 2-3 acute kidney injury (creatinine criteria with rolling 48-h and 7-day reference values, as in the MIMIC Code Repository) or kidney replacement therapy between 48 h and day 7.

**This repository contains code only. It does not contain any patient-level data.**

## Data access

The analysis uses [MIMIC-IV v3.1](https://doi.org/10.13026/kpb9-mt58), which is available through PhysioNet to credentialed users who complete the required training and sign the data use agreement. The SQL scripts query MIMIC-IV on Google BigQuery (project `physionet-data`, datasets `mimiciv_3_1_hosp` and `mimiciv_3_1_icu`).

## How to reproduce

1. Obtain credentialed access to MIMIC-IV on PhysioNet and request BigQuery access.
2. Create your own Google Cloud project (the free BigQuery sandbox is sufficient).
3. In the BigQuery console, replace `YOUR_PROJECT_ID` with your project ID and run, in this order:
   - `sql/01_cohort.sql` (creates `abe.kohort_genis` and `abe.kohort`)
   - `sql/02_covariates.sql` (creates `abe.kovaryat`)
   - `sql/05_sofa_day1.sql` (creates `abe.sofa_d1`)
   - `sql/04_flowchart_counts.sql` (counts for Figure 1; the last column equals the number of rows in `abe.kohort`)
   - optional: `sql/08_itemid_check.sql` (labels of the itemids used)
4. In Google Colaboratory, upload or clone this repository and run:
   ```
   !pip install -q tableone lifelines
   import os; os.environ['ABE_PROJECT_ID'] = 'your-project-id'
   %run analysis/abe_aki_analysis.py
   ```
   The script reads `sql/03`, `sql/06`, `sql/07` and `sql/09` directly, prints all results reported in the manuscript, writes them to `results_v2_2.txt` (aggregate results only) and saves Figures 2 and 3 to `figures/`. With 40 imputations and 200 bootstrap replicates it takes about 20-30 minutes.
5. `analysis/figure1_flowchart.py` draws Figure 1 from the aggregate counts.

Multiple imputation and bootstrap use fixed random seeds, so results are reproducible to the reported precision.

## Repository structure

```
sql/
  01_cohort.sql                  cohort, exposure (ABE) and outcome
  02_covariates.sql              covariates
  03_additional_covariates.sql   sodium bicarbonate, blood products, acetazolamide
  04_flowchart_counts.sql        flow chart counts
  05_sofa_day1.sql               non-renal SOFA score, first 24 h
  06_selection_comparison.sql    patients assessed during cohort construction
  07_extra_covariates.sql        cardiac surgery service, crystalloid type
  08_itemid_check.sql            labels of itemids used
  09_kdigo_reconciliation.sql    comparison of KDIGO definitions (rolling vs fixed baseline)
analysis/
  abe_aki_analysis.py            all analyses and figures 2-3
  figure1_flowchart.py           figure 1
figures/                         output folder
```

## Glossary of variable names (Turkish → English)

| Variable | Meaning |
|---|---|
| kohort / kovaryat | cohort / covariates |
| yas, cinsiyet, kabul_tipi, yb_tipi | age, sex, admission type, ICU type |
| abe_ort_0_24, abe_ort_24_48 | mean ABE on day 1 (−6 to 24 h) and day 2 (> 24 to 48 h) |
| abe_ort_0_24_siki, delta_siki | mean ABE in the strict 0-24 h window and corresponding ΔABE |
| abe_ilk_g1, abe_son_g2 | first ABE on day 1, last ABE on day 2 |
| be_ort_24_48, laktat_ort_24_48 | mean base excess and lactate on day 2 |
| klor_g1, klor_g2, sodyum_g1, sodyum_g2 | mean chloride and sodium on day 1 and day 2 |
| abe_min_0_24, abe_min_24_48 | lowest ABE on day 1 and day 2 |
| delta_abe | ΔABE = abe_ort_24_48 − abe_ort_0_24 |
| laktat_ort_0_24, be_ort_0_24 | mean lactate and base excess, 0-24 h |
| n_arteriyel, n_venoz | number of arterial / venous specimens |
| kre_bazal, kre_ilk, bazal_kaynak | baseline creatinine, admission creatinine, source of baseline |
| kre_onceki_min | lowest creatinine in the 7 days before ICU admission (primary baseline) |
| kre_son_onceki, kre_uzak_medyan | last pre-ICU creatinine; median creatinine 8-365 days before ICU (alternative baselines) |
| aki_zaman, krt_sonlanim, kre_sonlanim | time of the primary outcome; KRT component; creatinine component |
| sonlanim_evre3, sonlanim_sabit | KDIGO stage 3 or KRT; version-2 outcome with a fixed baseline (comparison only) |
| evre1_48s_kdigo, dis_kdigo | KDIGO stage 1 by 48 h; KDIGO stage 2-3 within 48 h (exclusion), rolling references |
| kre_oran_48 | highest creatinine 0-48 h / baseline |
| kalp_cerrahisi, salin_ml, dengeli_ml | cardiac surgery service; 0.9% saline and balanced crystalloid volumes 0-48 h |
| kre_max_0_48, kre_max_48s_7g | highest creatinine, −12 to 48 h and 48 h to day 7 |
| rrt_ilk | first kidney replacement therapy start time |
| sonlanim_aki23_rrt | primary outcome (KDIGO stage 2-3 AKI by creatinine or KRT, 48 h to day 7) |
| vazopressor_48s, imv_48s | vasopressors, invasive mechanical ventilation within 48 h |
| sivi_dengesi_48s, idrar_ml_kg_saat | fluid balance (mL), urine output (mL/kg/h), 0-48 h |
| map_min_24s, kilo | lowest mean arterial pressure 0-24 h, body weight |
| trombosit_min, bilirubin_max, bun_max, hemoglobin_min | laboratory values, 0-24 h |
| vankomisin_iv, pip_tazo, loop_diuretik, aminoglikozid | medications |
| bikarbonat_inf, eritrosit_tx, taze_plazma_tx, asetazolamid | sodium bicarbonate, red cells, fresh frozen plasma, acetazolamide |
| kbh, diyabet, kky, siroz, sepsis_icd | chronic kidney disease, diabetes, heart failure, cirrhosis, sepsis |
| olum_7g, y_bilesik | death within 7 days, composite outcome |

## Citation

If you use this code, please cite the article (when published) and this repository (see `CITATION.cff`). Please also cite MIMIC-IV and PhysioNet as requested on the [MIMIC-IV project page](https://physionet.org/content/mimiciv/3.1/).

## Licence

Code is released under the MIT Licence. MIMIC-IV data are subject to the PhysioNet Credentialed Health Data Use Agreement and are not covered by this licence.
