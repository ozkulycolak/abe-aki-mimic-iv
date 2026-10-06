-- =====================================================================
-- 01_cohort.sql
-- Builds the study cohort with exposure (ABE) and outcome variables.
-- Source: MIMIC-IV v3.1 on Google BigQuery (project physionet-data).
-- Replace YOUR_PROJECT_ID with your own Google Cloud project ID.
-- Time zero = ICU admission (icustays.intime).
-- =====================================================================

CREATE SCHEMA IF NOT EXISTS `YOUR_PROJECT_ID.abe` OPTIONS(location = 'US');

CREATE OR REPLACE TABLE `YOUR_PROJECT_ID.abe.kohort` AS
WITH esrd AS (
  -- End-stage kidney disease or chronic dialysis (exclusion)
  SELECT DISTINCT hadm_id
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd`
  WHERE (icd_version = 9  AND icd_code IN ('5856', 'V4511'))
     OR (icd_version = 10 AND (icd_code = 'N186' OR icd_code LIKE 'Z992%'))
),
icu AS (
  -- First ICU stay of adult patients
  SELECT i.subject_id, i.hadm_id, i.stay_id, i.intime, i.outtime, i.los,
         p.gender AS cinsiyet,
         p.anchor_age + EXTRACT(YEAR FROM i.intime) - p.anchor_year AS yas,
         a.admission_type AS kabul_tipi,
         a.deathtime, a.hospital_expire_flag
  FROM `physionet-data.mimiciv_3_1_icu.icustays` i
  JOIN `physionet-data.mimiciv_3_1_hosp.patients` p ON i.subject_id = p.subject_id
  JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a ON i.hadm_id = a.hadm_id
  WHERE p.anchor_age >= 18
  QUALIFY ROW_NUMBER() OVER (PARTITION BY i.subject_id ORDER BY i.intime) = 1
),
icu48 AS (
  -- ICU stay >= 48 h, alive at 48 h, no ESKD
  SELECT * FROM icu
  WHERE los >= 2
    AND (deathtime IS NULL OR deathtime > DATETIME_ADD(intime, INTERVAL 48 HOUR))
    AND hadm_id NOT IN (SELECT hadm_id FROM esrd)
),
gaz AS (
  -- Base excess (50802) and lactate (50813) paired by specimen_id; specimen type (52033)
  SELECT hadm_id, specimen_id,
    MIN(charttime) AS charttime,
    MAX(IF(itemid = 50802, valuenum, NULL)) AS be,
    MAX(IF(itemid = 50813, valuenum, NULL)) AS laktat,
    MAX(IF(itemid = 52033, value, NULL)) AS ornek_tipi
  FROM `physionet-data.mimiciv_3_1_hosp.labevents`
  WHERE itemid IN (52033, 50802, 50813)
  GROUP BY hadm_id, specimen_id
  HAVING MAX(IF(itemid = 50802, valuenum, NULL)) IS NOT NULL
     AND MAX(IF(itemid = 50813, valuenum, NULL)) IS NOT NULL
),
abe_pencere AS (
  -- Day 1: -6 h to 24 h; day 2: > 24 h to 48 h
  SELECT c.stay_id,
    AVG(IF(g.charttime <= DATETIME_ADD(c.intime, INTERVAL 24 HOUR), g.be + g.laktat, NULL)) AS abe_ort_0_24,
    MIN(IF(g.charttime <= DATETIME_ADD(c.intime, INTERVAL 24 HOUR), g.be + g.laktat, NULL)) AS abe_min_0_24,
    AVG(IF(g.charttime >  DATETIME_ADD(c.intime, INTERVAL 24 HOUR), g.be + g.laktat, NULL)) AS abe_ort_24_48,
    MIN(IF(g.charttime >  DATETIME_ADD(c.intime, INTERVAL 24 HOUR), g.be + g.laktat, NULL)) AS abe_min_24_48,
    AVG(IF(g.charttime <= DATETIME_ADD(c.intime, INTERVAL 24 HOUR), g.laktat, NULL)) AS laktat_ort_0_24,
    AVG(IF(g.charttime <= DATETIME_ADD(c.intime, INTERVAL 24 HOUR), g.be, NULL)) AS be_ort_0_24,
    COUNTIF(g.charttime <= DATETIME_ADD(c.intime, INTERVAL 24 HOUR)) AS n_0_24,
    COUNTIF(g.charttime >  DATETIME_ADD(c.intime, INTERVAL 24 HOUR)) AS n_24_48,
    COUNTIF(g.ornek_tipi = 'ART.') AS n_arteriyel,
    COUNTIF(g.ornek_tipi IN ('VEN.', 'MIX.', 'CENTRAL VENOUS.')) AS n_venoz
  FROM icu48 c
  JOIN gaz g ON c.hadm_id = g.hadm_id
  WHERE g.charttime BETWEEN DATETIME_SUB(c.intime, INTERVAL 6 HOUR)
                        AND DATETIME_ADD(c.intime, INTERVAL 48 HOUR)
  GROUP BY c.stay_id
),
kre AS (
  -- Serum creatinine (50912), patient level, -7 days to day 7
  SELECT c.stay_id, c.intime, l.charttime, l.valuenum
  FROM icu48 c
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l ON c.subject_id = l.subject_id
  WHERE l.itemid = 50912
    AND l.valuenum IS NOT NULL AND l.valuenum < 30
    AND l.charttime BETWEEN DATETIME_SUB(c.intime, INTERVAL 7 DAY)
                        AND DATETIME_ADD(c.intime, INTERVAL 7 DAY)
),
kre_ozet AS (
  SELECT stay_id,
    MIN(IF(charttime < intime, valuenum, NULL)) AS kre_onceki_min,
    ARRAY_AGG(IF(charttime BETWEEN DATETIME_SUB(intime, INTERVAL 12 HOUR)
                               AND DATETIME_ADD(intime, INTERVAL 6 HOUR), valuenum, NULL)
              IGNORE NULLS ORDER BY charttime LIMIT 1)[SAFE_OFFSET(0)] AS kre_ilk,
    MAX(IF(charttime BETWEEN DATETIME_SUB(intime, INTERVAL 12 HOUR)
                         AND DATETIME_ADD(intime, INTERVAL 48 HOUR), valuenum, NULL)) AS kre_max_0_48,
    MAX(IF(charttime >  DATETIME_ADD(intime, INTERVAL 48 HOUR)
       AND charttime <= DATETIME_ADD(intime, INTERVAL 7 DAY), valuenum, NULL)) AS kre_max_48s_7g
  FROM kre
  GROUP BY stay_id
),
rrt AS (
  -- Kidney replacement therapy (procedureevents)
  SELECT stay_id, MIN(starttime) AS rrt_ilk
  FROM `physionet-data.mimiciv_3_1_icu.procedureevents`
  WHERE itemid IN (225441, 225802, 225803, 225805, 225809, 225955)
  GROUP BY stay_id
)
SELECT
  c.*,
  a.* EXCEPT(stay_id),
  a.abe_ort_24_48 - a.abe_ort_0_24 AS delta_abe,
  k.* EXCEPT(stay_id),
  COALESCE(k.kre_onceki_min, k.kre_ilk) AS kre_bazal,
  IF(k.kre_onceki_min IS NULL, 'ilk_deger', 'onceki_7gun_min') AS bazal_kaynak,
  r.rrt_ilk,
  -- Primary outcome: creatinine >= 2.0 x baseline or KRT start, > 48 h to day 7
  IF(
    (r.rrt_ilk >  DATETIME_ADD(c.intime, INTERVAL 48 HOUR)
     AND r.rrt_ilk <= DATETIME_ADD(c.intime, INTERVAL 7 DAY))
    OR k.kre_max_48s_7g >= 2 * COALESCE(k.kre_onceki_min, k.kre_ilk),
    1, 0) AS sonlanim_aki23_rrt
FROM icu48 c
JOIN abe_pencere a ON c.stay_id = a.stay_id
JOIN kre_ozet k ON c.stay_id = k.stay_id
LEFT JOIN rrt r ON c.stay_id = r.stay_id
WHERE a.n_0_24 >= 1
  AND a.n_24_48 >= 1
  AND k.kre_ilk IS NOT NULL
  AND (r.rrt_ilk IS NULL OR r.rrt_ilk > DATETIME_ADD(c.intime, INTERVAL 48 HOUR))
  AND k.kre_max_0_48 < 2 * COALESCE(k.kre_onceki_min, k.kre_ilk);
