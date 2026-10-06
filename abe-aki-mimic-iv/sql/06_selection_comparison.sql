-- =====================================================================
-- 06_selection_comparison.sql
-- Included patients vs eligible patients excluded (Supplementary Table S4).
-- Replace YOUR_PROJECT_ID with your own Google Cloud project ID.
-- =====================================================================
WITH esrd AS (
  SELECT DISTINCT hadm_id
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd`
  WHERE (icd_version = 9  AND icd_code IN ('5856', 'V4511'))
     OR (icd_version = 10 AND (icd_code = 'N186' OR icd_code LIKE 'Z992%'))
),
icu AS (
  SELECT i.subject_id, i.hadm_id, i.stay_id, i.intime, i.los, i.first_careunit,
         p.gender, p.anchor_age + EXTRACT(YEAR FROM i.intime) - p.anchor_year AS yas,
         a.deathtime, a.hospital_expire_flag
  FROM `physionet-data.mimiciv_3_1_icu.icustays` i
  JOIN `physionet-data.mimiciv_3_1_hosp.patients` p ON i.subject_id = p.subject_id
  JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a ON i.hadm_id = a.hadm_id
  WHERE p.anchor_age >= 18
  QUALIFY ROW_NUMBER() OVER (PARTITION BY i.subject_id ORDER BY i.intime) = 1
),
uygun AS (
  SELECT c.* FROM icu c
  LEFT JOIN esrd e ON c.hadm_id = e.hadm_id
  WHERE c.los >= 2
    AND (c.deathtime IS NULL OR c.deathtime > DATETIME_ADD(c.intime, INTERVAL 48 HOUR))
    AND e.hadm_id IS NULL
),
gaz AS (
  SELECT hadm_id, specimen_id, MIN(charttime) AS charttime,
    MAX(IF(itemid = 50802, valuenum, NULL)) AS be,
    MAX(IF(itemid = 50813, valuenum, NULL)) AS laktat
  FROM `physionet-data.mimiciv_3_1_hosp.labevents`
  WHERE itemid IN (50802, 50813) AND valuenum IS NOT NULL
  GROUP BY hadm_id, specimen_id
  HAVING COUNT(DISTINCT itemid) = 2
),
abe AS (
  SELECT u.stay_id,
    COUNTIF(g.charttime <= DATETIME_ADD(u.intime, INTERVAL 24 HOUR)) AS n_0_24,
    COUNTIF(g.charttime >  DATETIME_ADD(u.intime, INTERVAL 24 HOUR)) AS n_24_48,
    AVG(IF(g.charttime <= DATETIME_ADD(u.intime, INTERVAL 24 HOUR), g.be + g.laktat, NULL)) AS abe_d1,
    AVG(IF(g.charttime <= DATETIME_ADD(u.intime, INTERVAL 24 HOUR), g.laktat, NULL)) AS laktat_d1
  FROM uygun u
  JOIN gaz g ON u.hadm_id = g.hadm_id
  WHERE g.charttime BETWEEN DATETIME_SUB(u.intime, INTERVAL 6 HOUR)
                        AND DATETIME_ADD(u.intime, INTERVAL 48 HOUR)
  GROUP BY u.stay_id
),
gruplu AS (
  SELECT u.*, a.abe_d1, a.laktat_d1,
    CASE
      WHEN u.stay_id IN (SELECT stay_id FROM `YOUR_PROJECT_ID.abe.kohort`) THEN '1_included'
      WHEN a.n_24_48 = 0 THEN '2_excluded_no_day2_blood_gas'
      ELSE '3_excluded_other_reasons'
    END AS grup
  FROM uygun u
  JOIN abe a ON u.stay_id = a.stay_id
  WHERE a.n_0_24 >= 1
)
SELECT
  grup,
  COUNT(*) AS n,
  APPROX_QUANTILES(yas, 4) AS age_quartiles,
  ROUND(AVG(IF(gender = 'M', 1, 0)) * 100, 1) AS male_pct,
  ROUND(AVG(IF(first_careunit LIKE '%CVICU%', 1, 0)) * 100, 1) AS cvicu_pct,
  ROUND(AVG(IF(first_careunit LIKE 'Medical Intensive%', 1, 0)) * 100, 1) AS micu_pct,
  ROUND(AVG(IF(first_careunit LIKE '%TSICU%', 1, 0)) * 100, 1) AS tsicu_pct,
  APPROX_QUANTILES(abe_d1, 4) AS abe_day1_quartiles,
  APPROX_QUANTILES(laktat_d1, 4) AS lactate_day1_quartiles,
  APPROX_QUANTILES(los, 4) AS icu_los_days_quartiles,
  ROUND(AVG(hospital_expire_flag) * 100, 1) AS hospital_mortality_pct
FROM gruplu
GROUP BY grup
ORDER BY grup;
