-- =====================================================================
-- 04_flowchart_counts.sql
-- Sequential exclusion counts for the study flow chart (Figure 1).
-- The last column must equal the number of rows in abe.kohort (6,173).
-- =====================================================================
WITH esrd AS (
  SELECT DISTINCT hadm_id
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd`
  WHERE (icd_version = 9  AND icd_code IN ('5856', 'V4511'))
     OR (icd_version = 10 AND (icd_code = 'N186' OR icd_code LIKE 'Z992%'))
),
icu AS (
  SELECT i.subject_id, i.hadm_id, i.stay_id, i.intime, i.los, a.deathtime
  FROM `physionet-data.mimiciv_3_1_icu.icustays` i
  JOIN `physionet-data.mimiciv_3_1_hosp.patients` p ON i.subject_id = p.subject_id
  JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a ON i.hadm_id = a.hadm_id
  WHERE p.anchor_age >= 18
  QUALIFY ROW_NUMBER() OVER (PARTITION BY i.subject_id ORDER BY i.intime) = 1
),
gaz AS (
  SELECT hadm_id, specimen_id, MIN(charttime) AS charttime
  FROM `physionet-data.mimiciv_3_1_hosp.labevents`
  WHERE itemid IN (50802, 50813) AND valuenum IS NOT NULL
  GROUP BY hadm_id, specimen_id
  HAVING COUNT(DISTINCT itemid) = 2
),
abe AS (
  SELECT c.stay_id,
    COUNTIF(g.charttime <= DATETIME_ADD(c.intime, INTERVAL 24 HOUR)) AS n_0_24,
    COUNTIF(g.charttime >  DATETIME_ADD(c.intime, INTERVAL 24 HOUR)) AS n_24_48
  FROM icu c
  JOIN gaz g ON c.hadm_id = g.hadm_id
  WHERE g.charttime BETWEEN DATETIME_SUB(c.intime, INTERVAL 6 HOUR)
                        AND DATETIME_ADD(c.intime, INTERVAL 48 HOUR)
  GROUP BY c.stay_id
),
kre AS (
  SELECT c.stay_id, c.intime, l.charttime, l.valuenum
  FROM icu c
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l ON c.subject_id = l.subject_id
  WHERE l.itemid = 50912 AND l.valuenum IS NOT NULL AND l.valuenum < 30
    AND l.charttime BETWEEN DATETIME_SUB(c.intime, INTERVAL 7 DAY)
                        AND DATETIME_ADD(c.intime, INTERVAL 48 HOUR)
),
kre_ozet AS (
  SELECT stay_id,
    MIN(IF(charttime < intime, valuenum, NULL)) AS kre_onceki_min,
    ARRAY_AGG(IF(charttime BETWEEN DATETIME_SUB(intime, INTERVAL 12 HOUR)
                               AND DATETIME_ADD(intime, INTERVAL 6 HOUR), valuenum, NULL)
              IGNORE NULLS ORDER BY charttime LIMIT 1)[SAFE_OFFSET(0)] AS kre_ilk,
    MAX(IF(charttime >= DATETIME_SUB(intime, INTERVAL 12 HOUR), valuenum, NULL)) AS kre_max_0_48
  FROM kre
  GROUP BY stay_id
),
rrt AS (
  SELECT stay_id, MIN(starttime) AS rrt_ilk
  FROM `physionet-data.mimiciv_3_1_icu.procedureevents`
  WHERE itemid IN (225441, 225802, 225803, 225805, 225809, 225955)
  GROUP BY stay_id
),
f AS (
  SELECT
    c.los >= 2 AS s1,
    (c.deathtime IS NULL OR c.deathtime > DATETIME_ADD(c.intime, INTERVAL 48 HOUR)) AS s2,
    e.hadm_id IS NULL AS s3,
    IFNULL(a.n_0_24, 0) >= 1 AS s4,
    IFNULL(a.n_24_48, 0) >= 1 AS s5,
    k.kre_ilk IS NOT NULL AS s6,
    (r.rrt_ilk IS NULL OR r.rrt_ilk > DATETIME_ADD(c.intime, INTERVAL 48 HOUR)) AS s7,
    IFNULL(k.kre_max_0_48 < 2 * COALESCE(k.kre_onceki_min, k.kre_ilk), FALSE) AS s8
  FROM icu c
  LEFT JOIN esrd e ON c.hadm_id = e.hadm_id
  LEFT JOIN abe a ON c.stay_id = a.stay_id
  LEFT JOIN kre_ozet k ON c.stay_id = k.stay_id
  LEFT JOIN rrt r ON c.stay_id = r.stay_id
)
SELECT
  COUNT(*) AS a0_first_adult_icu_stays,
  COUNTIF(s1) AS a1_icu_los_48h,
  COUNTIF(s1 AND s2) AS a2_alive_at_48h,
  COUNTIF(s1 AND s2 AND s3) AS a3_no_eskd,
  COUNTIF(s1 AND s2 AND s3 AND s4) AS a4_paired_be_lactate_0_24h,
  COUNTIF(s1 AND s2 AND s3 AND s4 AND s5) AS a5_paired_be_lactate_24_48h,
  COUNTIF(s1 AND s2 AND s3 AND s4 AND s5 AND s6) AS a6_creatinine_available,
  COUNTIF(s1 AND s2 AND s3 AND s4 AND s5 AND s6 AND s7) AS a7_no_early_krt,
  COUNTIF(s1 AND s2 AND s3 AND s4 AND s5 AND s6 AND s7 AND s8) AS a8_no_early_aki_final
FROM f;
