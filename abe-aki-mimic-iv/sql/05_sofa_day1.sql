-- =====================================================================
-- 05_sofa_day1.sql
-- Non-renal SOFA components during the first 24 h (sensitivity analyses).
-- Components that cannot be assessed are scored as normal (0).
-- GCS values recorded with verbal response 'No Response-ETT' are not assessable.
-- Replace YOUR_PROJECT_ID with your own Google Cloud project ID.
-- =====================================================================
CREATE OR REPLACE TABLE `YOUR_PROJECT_ID.abe.sofa_d1` AS
WITH k AS (
  SELECT stay_id, subject_id, hadm_id, intime FROM `YOUR_PROJECT_ID.abe.kohort`
),
bg AS (
  SELECT k.stay_id, l.specimen_id,
    MIN(l.charttime) AS charttime,
    MAX(IF(l.itemid = 50821, l.valuenum, NULL)) AS po2,
    MAX(IF(l.itemid = 50816, l.valuenum, NULL)) AS fio2_lab,
    MAX(IF(l.itemid = 52033, l.value, NULL)) AS tip
  FROM k
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l ON k.hadm_id = l.hadm_id
  WHERE l.itemid IN (50821, 50816, 52033)
    AND l.charttime BETWEEN DATETIME_SUB(k.intime, INTERVAL 6 HOUR)
                        AND DATETIME_ADD(k.intime, INTERVAL 24 HOUR)
  GROUP BY k.stay_id, l.specimen_id
),
fio2_ch AS (
  SELECT k.stay_id, c.charttime,
    IF(c.valuenum <= 1, c.valuenum * 100, c.valuenum) AS fio2
  FROM k
  JOIN `physionet-data.mimiciv_3_1_icu.chartevents` c ON k.stay_id = c.stay_id
  WHERE c.itemid = 223835
    AND c.valuenum IS NOT NULL
    AND c.charttime BETWEEN DATETIME_SUB(k.intime, INTERVAL 10 HOUR)
                        AND DATETIME_ADD(k.intime, INTERVAL 24 HOUR)
),
pf_join AS (
  SELECT b.stay_id, b.specimen_id, b.charttime, b.po2, b.fio2_lab, f.fio2 AS fio2_ch,
    ROW_NUMBER() OVER (PARTITION BY b.stay_id, b.specimen_id ORDER BY f.charttime DESC) AS rn
  FROM bg b
  LEFT JOIN fio2_ch f
    ON f.stay_id = b.stay_id
   AND f.charttime BETWEEN DATETIME_SUB(b.charttime, INTERVAL 4 HOUR) AND b.charttime
   AND f.fio2 BETWEEN 21 AND 100
  WHERE b.tip = 'ART.' AND b.po2 IS NOT NULL
),
pf AS (
  SELECT stay_id, charttime, po2,
    COALESCE(IF(fio2_lab BETWEEN 21 AND 100, fio2_lab, NULL), fio2_ch) AS fio2
  FROM pf_join
  WHERE rn = 1
),
vent AS (
  SELECT p.stay_id, p.starttime, p.endtime
  FROM `physionet-data.mimiciv_3_1_icu.procedureevents` p
  JOIN k ON p.stay_id = k.stay_id
  WHERE p.itemid = 225792
),
pf_v AS (
  SELECT p.stay_id, p.charttime,
    100 * p.po2 / p.fio2 AS pf_ratio,
    MAX(IF(v.stay_id IS NULL, 0, 1)) AS vent
  FROM pf p
  LEFT JOIN vent v
    ON v.stay_id = p.stay_id AND p.charttime BETWEEN v.starttime AND v.endtime
  WHERE p.fio2 IS NOT NULL
  GROUP BY p.stay_id, p.charttime, p.po2, p.fio2
),
resp AS (
  SELECT stay_id,
    MIN(pf_ratio) AS pf_min,
    MAX(CASE
          WHEN pf_ratio < 100 AND vent = 1 THEN 4
          WHEN pf_ratio < 200 AND vent = 1 THEN 3
          WHEN pf_ratio < 300 THEN 2
          WHEN pf_ratio < 400 THEN 1
          ELSE 0 END) AS sofa_resp
  FROM pf_v
  GROUP BY stay_id
),
vp AS (
  SELECT k.stay_id,
    MAX(IF(ie.itemid = 221906, ie.rate, NULL)) AS ne_max,
    MAX(IF(ie.itemid = 221289, ie.rate, NULL)) AS epi_max,
    MAX(IF(ie.itemid = 221662, ie.rate, NULL)) AS dopa_max,
    MAX(IF(ie.itemid = 221653, ie.rate, NULL)) AS dobu_max
  FROM k
  JOIN `physionet-data.mimiciv_3_1_icu.inputevents` ie ON k.stay_id = ie.stay_id
  WHERE ie.itemid IN (221906, 221289, 221662, 221653)
    AND ie.rateuom = 'mcg/kg/min'
    AND ie.starttime < DATETIME_ADD(k.intime, INTERVAL 24 HOUR)
    AND ie.endtime > k.intime
  GROUP BY k.stay_id
),
gcs AS (
  SELECT k.stay_id, c.charttime,
    MAX(IF(c.itemid = 220739, c.valuenum, NULL)) AS eye,
    MAX(IF(c.itemid = 223900, c.valuenum, NULL)) AS verbal,
    MAX(IF(c.itemid = 223900 AND c.value = 'No Response-ETT', 1, 0)) AS ett,
    MAX(IF(c.itemid = 223901, c.valuenum, NULL)) AS motor
  FROM k
  JOIN `physionet-data.mimiciv_3_1_icu.chartevents` c ON k.stay_id = c.stay_id
  WHERE c.itemid IN (220739, 223900, 223901)
    AND c.charttime BETWEEN DATETIME_SUB(k.intime, INTERVAL 6 HOUR)
                        AND DATETIME_ADD(k.intime, INTERVAL 24 HOUR)
  GROUP BY k.stay_id, c.charttime
),
gcs_min AS (
  SELECT stay_id,
    MIN(IF(ett = 0 AND eye IS NOT NULL AND verbal IS NOT NULL AND motor IS NOT NULL,
           eye + verbal + motor, NULL)) AS gcs_min
  FROM gcs
  GROUP BY stay_id
)
SELECT
  k.stay_id,
  r.pf_min,
  g.gcs_min,
  IFNULL(r.sofa_resp, 0) AS sofa_resp,
  CASE WHEN v.trombosit_min < 20 THEN 4 WHEN v.trombosit_min < 50 THEN 3
       WHEN v.trombosit_min < 100 THEN 2 WHEN v.trombosit_min < 150 THEN 1 ELSE 0 END AS sofa_coag,
  CASE WHEN v.bilirubin_max >= 12 THEN 4 WHEN v.bilirubin_max >= 6 THEN 3
       WHEN v.bilirubin_max >= 2 THEN 2 WHEN v.bilirubin_max >= 1.2 THEN 1 ELSE 0 END AS sofa_liver,
  CASE WHEN vp.dopa_max > 15 OR vp.epi_max > 0.1 OR vp.ne_max > 0.1 THEN 4
       WHEN vp.dopa_max > 5 OR vp.epi_max > 0 OR vp.ne_max > 0 THEN 3
       WHEN vp.dopa_max > 0 OR vp.dobu_max > 0 THEN 2
       WHEN v.map_min_24s < 70 THEN 1 ELSE 0 END AS sofa_cv,
  CASE WHEN g.gcs_min < 6 THEN 4 WHEN g.gcs_min < 10 THEN 3
       WHEN g.gcs_min < 13 THEN 2 WHEN g.gcs_min < 15 THEN 1 ELSE 0 END AS sofa_cns
FROM k
JOIN `YOUR_PROJECT_ID.abe.kovaryat` v ON k.stay_id = v.stay_id
LEFT JOIN resp r ON k.stay_id = r.stay_id
LEFT JOIN vp ON k.stay_id = vp.stay_id
LEFT JOIN gcs_min g ON k.stay_id = g.stay_id;
