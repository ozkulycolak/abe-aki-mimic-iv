-- =====================================================================
-- 01_cohort.sql  (version 2.1)
-- Study cohort, exposure (ABE) and outcomes.
-- Source: MIMIC-IV v3.1 on Google BigQuery (project physionet-data).
-- Replace YOUR_PROJECT_ID with your own Google Cloud project ID.
-- Time zero = ICU admission (icustays.intime).
-- Exposure windows: day 1 = from 6 h before to 24 h after ICU admission
-- (inclusive); day 2 = > 24 h to 48 h. A strict day-1 window (0-24 h) is
-- also stored for sensitivity analyses.
-- KDIGO creatinine staging follows the MIMIC Code Repository algorithm:
-- for every creatinine value, reference values are the lowest creatinine in
-- the preceding 48 h and in the preceding 7 days (rolling references).
--   stage >= 1: value >= lowest 48 h + 0.3 mg/dL, or >= 1.5 x lowest 7 days
--   stage >= 2: value >= 2.0 x lowest 7 days
--   stage 3   : value >= 3.0 x lowest 7 days, or value >= 4.0 mg/dL with an
--               acute rise (lowest 48 h <= 3.7 mg/dL or >= 1.5 x lowest 7 days)
-- The previous fixed-baseline definition (version 2) is kept for comparison.
-- Two tables are created: abe.kohort_genis (before the early-AKI exclusion,
-- with both definitions) and abe.kohort (final cohort).
-- =====================================================================

CREATE SCHEMA IF NOT EXISTS `YOUR_PROJECT_ID.abe` OPTIONS(location = 'US');

CREATE OR REPLACE TABLE `YOUR_PROJECT_ID.abe.kohort_genis` AS
WITH esrd AS (
  SELECT DISTINCT hadm_id
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd`
  WHERE (icd_version = 9  AND icd_code IN ('5856', 'V4511'))
     OR (icd_version = 10 AND (icd_code = 'N186' OR icd_code LIKE 'Z992%'))
),
icu AS (
  SELECT i.subject_id, i.hadm_id, i.stay_id, i.intime, i.outtime, i.los,
         p.gender AS cinsiyet,
         p.anchor_age + EXTRACT(YEAR FROM i.intime) - p.anchor_year AS yas,
         p.dod,
         a.admission_type AS kabul_tipi,
         a.dischtime, a.deathtime, a.hospital_expire_flag
  FROM `physionet-data.mimiciv_3_1_icu.icustays` i
  JOIN `physionet-data.mimiciv_3_1_hosp.patients` p ON i.subject_id = p.subject_id
  JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a ON i.hadm_id = a.hadm_id
  WHERE p.anchor_age >= 18
  QUALIFY ROW_NUMBER() OVER (PARTITION BY i.subject_id ORDER BY i.intime) = 1
),
icu48 AS (
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
gaz_p AS (
  SELECT c.stay_id, c.intime, g.charttime, g.be, g.laktat, g.be + g.laktat AS abe, g.ornek_tipi,
    g.charttime <= DATETIME_ADD(c.intime, INTERVAL 24 HOUR) AS g1,
    g.charttime >= c.intime AND g.charttime <= DATETIME_ADD(c.intime, INTERVAL 24 HOUR) AS g1_siki,
    g.charttime >  DATETIME_ADD(c.intime, INTERVAL 24 HOUR) AS g2
  FROM icu48 c
  JOIN gaz g ON c.hadm_id = g.hadm_id
  WHERE g.charttime BETWEEN DATETIME_SUB(c.intime, INTERVAL 6 HOUR)
                        AND DATETIME_ADD(c.intime, INTERVAL 48 HOUR)
),
abe_pencere AS (
  SELECT stay_id,
    AVG(IF(g1, abe, NULL)) AS abe_ort_0_24,
    MIN(IF(g1, abe, NULL)) AS abe_min_0_24,
    AVG(IF(g2, abe, NULL)) AS abe_ort_24_48,
    MIN(IF(g2, abe, NULL)) AS abe_min_24_48,
    AVG(IF(g1, laktat, NULL)) AS laktat_ort_0_24,
    AVG(IF(g2, laktat, NULL)) AS laktat_ort_24_48,
    AVG(IF(g1, be, NULL)) AS be_ort_0_24,
    AVG(IF(g2, be, NULL)) AS be_ort_24_48,
    AVG(IF(g1_siki, abe, NULL)) AS abe_ort_0_24_siki,
    COUNTIF(g1_siki) AS n_0_24_siki,
    ARRAY_AGG(IF(g1, abe, NULL) IGNORE NULLS ORDER BY charttime LIMIT 1)[SAFE_OFFSET(0)] AS abe_ilk_g1,
    ARRAY_AGG(IF(g2, abe, NULL) IGNORE NULLS ORDER BY charttime DESC LIMIT 1)[SAFE_OFFSET(0)] AS abe_son_g2,
    COUNTIF(g1) AS n_0_24,
    COUNTIF(g2) AS n_24_48,
    COUNTIF(ornek_tipi = 'ART.') AS n_arteriyel,
    COUNTIF(ornek_tipi IN ('VEN.', 'MIX.', 'CENTRAL VENOUS.')) AS n_venoz
  FROM gaz_p
  GROUP BY stay_id
),
elek AS (
  -- Serum chloride (50902) and sodium (50983), mean per exposure day
  SELECT c.stay_id,
    AVG(IF(l.itemid = 50902 AND l.charttime <= DATETIME_ADD(c.intime, INTERVAL 24 HOUR), l.valuenum, NULL)) AS klor_g1,
    AVG(IF(l.itemid = 50902 AND l.charttime >  DATETIME_ADD(c.intime, INTERVAL 24 HOUR), l.valuenum, NULL)) AS klor_g2,
    AVG(IF(l.itemid = 50983 AND l.charttime <= DATETIME_ADD(c.intime, INTERVAL 24 HOUR), l.valuenum, NULL)) AS sodyum_g1,
    AVG(IF(l.itemid = 50983 AND l.charttime >  DATETIME_ADD(c.intime, INTERVAL 24 HOUR), l.valuenum, NULL)) AS sodyum_g2
  FROM icu48 c
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l ON c.hadm_id = l.hadm_id
  WHERE l.itemid IN (50902, 50983) AND l.valuenum IS NOT NULL
    AND l.charttime BETWEEN DATETIME_SUB(c.intime, INTERVAL 6 HOUR)
                        AND DATETIME_ADD(c.intime, INTERVAL 48 HOUR)
  GROUP BY c.stay_id
),
kre AS (
  -- Serum creatinine (50912), patient level, -365 days to day 7
  SELECT c.stay_id, c.intime, l.charttime, l.valuenum
  FROM icu48 c
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l ON c.subject_id = l.subject_id
  WHERE l.itemid = 50912
    AND l.valuenum IS NOT NULL AND l.valuenum < 30
    AND l.charttime BETWEEN DATETIME_SUB(c.intime, INTERVAL 365 DAY)
                        AND DATETIME_ADD(c.intime, INTERVAL 7 DAY)
),
kre_ozet AS (
  SELECT stay_id,
    -- primary baseline: lowest value in the 7 days before ICU admission
    MIN(IF(charttime < intime AND charttime >= DATETIME_SUB(intime, INTERVAL 7 DAY), valuenum, NULL)) AS kre_onceki_min,
    ARRAY_AGG(IF(charttime < intime AND charttime >= DATETIME_SUB(intime, INTERVAL 7 DAY),
                 STRUCT(valuenum AS v, charttime AS t), NULL)
              IGNORE NULLS ORDER BY valuenum, charttime DESC LIMIT 1)[SAFE_OFFSET(0)].t AS kre_bazal_zaman,
    -- alternative baselines
    ARRAY_AGG(IF(charttime < intime AND charttime >= DATETIME_SUB(intime, INTERVAL 7 DAY), valuenum, NULL)
              IGNORE NULLS ORDER BY charttime DESC LIMIT 1)[SAFE_OFFSET(0)] AS kre_son_onceki,
    APPROX_QUANTILES(IF(charttime BETWEEN DATETIME_SUB(intime, INTERVAL 365 DAY)
                                      AND DATETIME_SUB(intime, INTERVAL 8 DAY), valuenum, NULL),
                     2 IGNORE NULLS)[SAFE_OFFSET(1)] AS kre_uzak_medyan,
    ARRAY_AGG(IF(charttime BETWEEN DATETIME_SUB(intime, INTERVAL 12 HOUR)
                               AND DATETIME_ADD(intime, INTERVAL 6 HOUR), valuenum, NULL)
              IGNORE NULLS ORDER BY charttime LIMIT 1)[SAFE_OFFSET(0)] AS kre_ilk,
    MAX(IF(charttime BETWEEN DATETIME_SUB(intime, INTERVAL 12 HOUR)
                         AND DATETIME_ADD(intime, INTERVAL 48 HOUR), valuenum, NULL)) AS kre_max_0_48,
    MAX(IF(charttime >  DATETIME_ADD(intime, INTERVAL 48 HOUR)
       AND charttime <= DATETIME_ADD(intime, INTERVAL 7 DAY), valuenum, NULL)) AS kre_max_48s_7g,
    COUNTIF(charttime >  DATETIME_ADD(intime, INTERVAL 48 HOUR)
        AND charttime <= DATETIME_ADD(intime, INTERVAL 7 DAY)) AS n_kre_48s_7g
  FROM kre
  GROUP BY stay_id
),
kre_olay AS (
  -- time of the first creatinine meeting KDIGO stage >= 2 between 48 h and day 7
  SELECT k.stay_id,
    MIN(IF(k.charttime > DATETIME_ADD(k.intime, INTERVAL 48 HOUR)
       AND k.charttime <= DATETIME_ADD(k.intime, INTERVAL 7 DAY)
       AND (k.valuenum >= 2 * COALESCE(o.kre_onceki_min, o.kre_ilk)
            OR (k.valuenum >= 4.0 AND k.valuenum >= COALESCE(o.kre_onceki_min, o.kre_ilk) + 0.3)),
       k.charttime, NULL)) AS aki_kre_zaman
  FROM kre k
  JOIN kre_ozet o ON k.stay_id = o.stay_id
  GROUP BY k.stay_id
),
kre_roll AS (
  -- rolling KDIGO references for every creatinine value (MIMIC Code Repository logic)
  SELECT stay_id, intime, charttime, valuenum,
    MIN(valuenum) OVER (PARTITION BY stay_id ORDER BY UNIX_SECONDS(TIMESTAMP(charttime))
                        RANGE BETWEEN 172800 PRECEDING AND 1 PRECEDING) AS low48,
    MIN(valuenum) OVER (PARTITION BY stay_id ORDER BY UNIX_SECONDS(TIMESTAMP(charttime))
                        RANGE BETWEEN 604800 PRECEDING AND 1 PRECEDING) AS low7
  FROM kre
),
kre_evre AS (
  SELECT stay_id, intime, charttime, valuenum, low48, low7,
    (valuenum >= low48 + 0.3 OR valuenum >= 1.5 * low7) AS st1,
    (valuenum >= 2 * low7) AS st2,
    (valuenum >= 3 * low7 OR (valuenum >= 4.0 AND (low48 <= 3.7 OR valuenum >= 1.5 * low7))) AS st3
  FROM kre_roll
),
kdigo AS (
  SELECT stay_id,
    -- early window (-12 h to 48 h)
    LOGICAL_OR(charttime > DATETIME_SUB(intime, INTERVAL 12 HOUR) AND charttime <= DATETIME_ADD(intime, INTERVAL 48 HOUR)
               AND (IFNULL(st2, FALSE) OR IFNULL(st3, FALSE))) AS erken_aki23,
    LOGICAL_OR(charttime > DATETIME_SUB(intime, INTERVAL 12 HOUR) AND charttime <= DATETIME_ADD(intime, INTERVAL 48 HOUR)
               AND IFNULL(st1, FALSE)) AS evre1_48s_b,
    -- outcome window (> 48 h to day 7)
    MIN(IF(charttime > DATETIME_ADD(intime, INTERVAL 48 HOUR) AND charttime <= DATETIME_ADD(intime, INTERVAL 7 DAY)
           AND (IFNULL(st2, FALSE) OR IFNULL(st3, FALSE)), charttime, NULL)) AS aki_kdigo_zaman,
    LOGICAL_OR(charttime > DATETIME_ADD(intime, INTERVAL 48 HOUR) AND charttime <= DATETIME_ADD(intime, INTERVAL 7 DAY)
               AND IFNULL(st3, FALSE)) AS evre3_kdigo,
    -- acute absolute criterion (>= 4.0 mg/dL with acute rise) in the outcome window
    LOGICAL_OR(charttime > DATETIME_ADD(intime, INTERVAL 48 HOUR) AND charttime <= DATETIME_ADD(intime, INTERVAL 7 DAY)
               AND valuenum >= 4.0 AND (IFNULL(low48 <= 3.7, FALSE) OR IFNULL(valuenum >= 1.5 * low7, FALSE))) AS mutlak_akut
  FROM kre_evre
  GROUP BY stay_id
),
rrt AS (
  -- Kidney replacement therapy (procedureevents)
  SELECT stay_id, MIN(starttime) AS rrt_ilk
  FROM `physionet-data.mimiciv_3_1_icu.procedureevents`
  WHERE itemid IN (225441, 225802, 225803, 225805, 225809, 225955)
  GROUP BY stay_id
),
birlesik AS (
  SELECT
    c.*,
    a.* EXCEPT(stay_id),
    a.abe_ort_24_48 - a.abe_ort_0_24 AS delta_abe,
    e.* EXCEPT(stay_id),
    k.* EXCEPT(stay_id),
    ko.aki_kre_zaman,
    kd.* EXCEPT(stay_id),
    COALESCE(k.kre_onceki_min, k.kre_ilk) AS kre_bazal,
    IF(k.kre_onceki_min IS NULL, 'ilk_deger', 'onceki_7gun_min') AS bazal_kaynak,
    r.rrt_ilk,
    (r.rrt_ilk > DATETIME_ADD(c.intime, INTERVAL 48 HOUR)
     AND r.rrt_ilk <= DATETIME_ADD(c.intime, INTERVAL 7 DAY)) AS krt_pencere
  FROM icu48 c
  JOIN abe_pencere a ON c.stay_id = a.stay_id
  JOIN kre_ozet k ON c.stay_id = k.stay_id
  LEFT JOIN kre_olay ko ON c.stay_id = ko.stay_id
  LEFT JOIN kdigo kd ON c.stay_id = kd.stay_id
  LEFT JOIN elek e ON c.stay_id = e.stay_id
  LEFT JOIN rrt r ON c.stay_id = r.stay_id
)
SELECT
  b.* EXCEPT(krt_pencere),
  IF(IFNULL(b.krt_pencere, FALSE), 1, 0) AS krt_sonlanim,
  IF(b.aki_kdigo_zaman IS NOT NULL, 1, 0) AS kre_sonlanim,
  -- Primary outcome: KDIGO stage 2-3 by creatinine (rolling references) or KRT, > 48 h to day 7
  IF(IFNULL(b.krt_pencere, FALSE) OR b.aki_kdigo_zaman IS NOT NULL, 1, 0) AS sonlanim_aki23_rrt,
  -- Secondary outcome: KDIGO stage 3 by creatinine or KRT
  IF(IFNULL(b.krt_pencere, FALSE) OR IFNULL(b.evre3_kdigo, FALSE), 1, 0) AS sonlanim_evre3,
  IF(IFNULL(b.evre1_48s_b, FALSE), 1, 0) AS evre1_48s_kdigo,
  -- Version-2 definitions (fixed baseline), for comparison
  IF(IFNULL(b.krt_pencere, FALSE) OR b.aki_kre_zaman IS NOT NULL, 1, 0) AS sonlanim_sabit,
  IF(b.kre_max_48s_7g >= 4.0 AND b.kre_max_48s_7g >= b.kre_bazal + 0.3
     AND b.kre_max_48s_7g < 2 * b.kre_bazal, 1, 0) AS mutlak_sadece_sabit,
  IF(IFNULL(b.mutlak_akut, FALSE), 1, 0) AS mutlak_akut_flag,
  (b.kre_max_0_48 >= 2 * b.kre_bazal
   OR (b.kre_max_0_48 >= 4.0 AND b.kre_max_0_48 >= b.kre_bazal + 0.3)) AS dis_sabit,
  IFNULL(b.erken_aki23, FALSE) AS dis_kdigo,
  -- Time of the primary outcome (first creatinine meeting criteria or KRT start)
  CASE
    WHEN IFNULL(b.krt_pencere, FALSE) AND b.aki_kdigo_zaman IS NOT NULL THEN LEAST(b.rrt_ilk, b.aki_kdigo_zaman)
    WHEN IFNULL(b.krt_pencere, FALSE) THEN b.rrt_ilk
    ELSE b.aki_kdigo_zaman
  END AS aki_zaman
FROM birlesik b
WHERE b.n_0_24 >= 1
  AND b.n_24_48 >= 1
  AND b.kre_ilk IS NOT NULL
  AND (b.rrt_ilk IS NULL OR b.rrt_ilk > DATETIME_ADD(b.intime, INTERVAL 48 HOUR));

CREATE OR REPLACE TABLE `YOUR_PROJECT_ID.abe.kohort` AS
SELECT * FROM `YOUR_PROJECT_ID.abe.kohort_genis`
WHERE NOT dis_kdigo;
