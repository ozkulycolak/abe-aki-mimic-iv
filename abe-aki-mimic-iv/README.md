# Early trajectory of alactic base excess and acute kidney injury (MIMIC-IV)

Code for the study:

> Kılıç Ö, Yılmaz Çolak Ö. *Early trajectory of alactic base excess and acute kidney injury in critically ill adults: a retrospective cohort study using the MIMIC-IV database.* (manuscript submitted)

The study examines whether the change in alactic base excess (ABE = base excess + lactate) between the first and second 24 h after ICU admission (ΔABE) is associated with KDIGO stage 2-3 acute kidney injury or kidney replacement therapy between 48 h and day 7.

**This repository contains code only. It does not contain any patient-level data.**

## Data access

The analysis uses [MIMIC-IV v3.1](https://doi.org/10.13026/kpb9-mt58), which is available through PhysioNet to credentialed users who complete the required training and sign the data use agreement. The SQL scripts query MIMIC-IV on Google BigQuery (project `physionet-data`, datasets `mimiciv_3_1_hosp` and `mimiciv_3_1_icu`).

## How to reproduce

1. Obtain credentialed access to MIMIC-IV on PhysioNet and request BigQuery access.
2. Create your own Google Cloud project (the free BigQuery sandbox is sufficient).
3. In every file under `sql/`, replace `YOUR_PROJECT_ID` with your project ID.
4. In the BigQuery console, run, in this order:
   - `sql/01_cohort.sql` (creates `abe.kohort`; expected n = 6,173, events = 629)
   - `sql/02_covariates.sql` (creates `abe.kovaryat`)
   - `sql/04_flowchart_counts.sql` (counts for Figure 1; last column = 6,173)
5. Open Google Colaboratory, upload or clone this repository, set `PROJECT_ID` in `analysis/abe_aki_analysis.py` and run it. It reads `sql/03_additional_covariates.sql` directly and prints all results reported in the manuscript; Figures 2 and 3 are written to `figures/`.
6. `analysis/figure1_flowchart.py` draws Figure 1 from the aggregate counts.

Multiple imputation uses fixed random seeds (0 to 19), so results should be reproducible to the reported precision.

## Repository structure

```
sql/
  01_cohort.sql                  cohort, exposure (ABE) and outcome
  02_covariates.sql              covariates
  03_additional_covariates.sql   sodium bicarbonate, blood products, acetazolamide
  04_flowchart_counts.sql        flow chart counts
analysis/
  abe_aki_analysis.py            descriptive statistics, models, figures 2-3
  figure1_flowchart.py           figure 1
figures/                         output folder
```

## Glossary of variable names (Turkish → English)

| Variable | Meaning |
|---|---|
| kohort / kovaryat | cohort / covariates |
| yas, cinsiyet, kabul_tipi, yb_tipi | age, sex, admission type, ICU type |
| abe_ort_0_24, abe_ort_24_48 | mean ABE, 0-24 h and 24-48 h |
| abe_min_0_24, abe_min_24_48 | lowest ABE, 0-24 h and 24-48 h |
| delta_abe | ΔABE = abe_ort_24_48 − abe_ort_0_24 |
| laktat_ort_0_24, be_ort_0_24 | mean lactate and base excess, 0-24 h |
| n_arteriyel, n_venoz | number of arterial / venous specimens |
| kre_bazal, kre_ilk, bazal_kaynak | baseline creatinine, admission creatinine, source of baseline |
| kre_onceki_min | lowest creatinine in the 7 days before ICU admission |
| kre_max_0_48, kre_max_48s_7g | highest creatinine, −12 to 48 h and 48 h to day 7 |
| rrt_ilk | first kidney replacement therapy start time |
| sonlanim_aki23_rrt | primary outcome (KDIGO stage 2-3 AKI or KRT, 48 h to day 7) |
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
