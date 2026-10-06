-- =====================================================================
-- 02_covariates.sql
-- Covariates for the cohort created by 01_cohort.sql.
-- Replace YOUR_PROJECT_ID with your own Google Cloud project ID.
-- =====================================================================

CREATE OR REPLACE TABLE `YOUR_PROJECT_ID.abe.kovaryat` AS
WITH k AS (
  SELECT stay_id, subject_id, hadm_id, intime
  FROM `YOUR_PROJECT_ID.abe.kohort`
),
vazo AS (
  -- Vasopressor infusion overlapping 0-48 h
  SELECT k.stay_id,
    1 AS vazopressor_48s,
    COUNT(DISTINCT ie.itemid) AS vazopressor_sayisi,
    MAX(IF(ie.itemid = 221906, 1, 0)) AS norepinefrin
  FROM k
  JOIN `physionet-data.mimiciv_3_1_icu.inputevents` ie ON k.stay_id = ie.stay_id
  WHERE ie.itemid IN (221906, 221289, 222315, 221749, 221662)
    AND ie.starttime < DATETIME_ADD(k.intime, INTERVAL 48 HOUR)
    AND ie.endtime > k.intime
  GROUP BY k.stay_id
),
ventil AS (
  -- Invasive mechanical ventilation overlapping 0-48 h
  SELECT DISTINCT k.stay_id, 1 AS imv_48s
  FROM k
  JOIN `physionet-data.mimiciv_3_1_icu.procedureevents` p ON k.stay_id = p.stay_id
  WHERE p.itemid = 225792
    AND p.starttime < DATETIME_ADD(k.intime, INTERVAL 48 HOUR)
    AND p.endtime > k.intime
),
giris AS (
  -- Fluid intake (mL), starttime 0-48 h
  SELECT k.stay_id, SUM(ie.amount) AS sivi_giris_48s
  FROM k
  JOIN `physionet-data.mimiciv_3_1_icu.inputevents` ie ON k.stay_id = ie.stay_id
  WHERE LOWER(ie.amountuom) = 'ml'
    AND ie.starttime >= k.intime
    AND ie.starttime < DATETIME_ADD(k.intime, INTERVAL 48 HOUR)
  GROUP BY k.stay_id
),
cikis AS (
  -- Total output and urine output, 0-48 h (GU irrigant 227488 subtracted)
  SELECT k.stay_id,
    SUM(IF(o.itemid = 227488, -o.value, o.value)) AS toplam_cikis_48s,
    SUM(CASE
          WHEN o.itemid = 227488 THEN -o.value
          WHEN o.itemid IN (226559, 226560, 226561, 226584, 226563, 226564,
                            226565, 226567, 226557, 226558, 227489) THEN o.value
          ELSE 0 END) AS idrar_48s
  FROM k
  JOIN `physionet-data.mimiciv_3_1_icu.outputevents` o ON k.stay_id = o.stay_id
  WHERE o.charttime >= k.intime
    AND o.charttime < DATETIME_ADD(k.intime, INTERVAL 48 HOUR)
  GROUP BY k.stay_id
),
ce AS (
  -- Lowest MAP (0-24 h) and first admission weight (-6 to 24 h)
  SELECT k.stay_id,
    MIN(IF(c.itemid IN (220052, 220181, 225312)
           AND c.valuenum BETWEEN 20 AND 200
           AND c.charttime >= k.intime, c.valuenum, NULL)) AS map_min_24s,
    ARRAY_AGG(IF(c.itemid = 226512 AND c.valuenum BETWEEN 30 AND 300, c.valuenum, NULL)
              IGNORE NULLS ORDER BY c.charttime LIMIT 1)[SAFE_OFFSET(0)] AS kilo
  FROM k
  JOIN `physionet-data.mimiciv_3_1_icu.chartevents` c ON k.stay_id = c.stay_id
  WHERE c.itemid IN (220052, 220181, 225312, 226512)
    AND c.charttime BETWEEN DATETIME_SUB(k.intime, INTERVAL 6 HOUR)
                        AND DATETIME_ADD(k.intime, INTERVAL 24 HOUR)
  GROUP BY k.stay_id
),
lab AS (
  -- Laboratory values, -6 h to 24 h
  SELECT k.stay_id,
    MIN(IF(l.itemid = 51265, l.valuenum, NULL)) AS trombosit_min,
    MAX(IF(l.itemid = 50885, l.valuenum, NULL)) AS bilirubin_max,
    MAX(IF(l.itemid = 51006, l.valuenum, NULL)) AS bun_max,
    MIN(IF(l.itemid = 50862, l.valuenum, NULL)) AS albumin_min,
    MAX(IF(l.itemid = 50902, l.valuenum, NULL)) AS klor_max,
    MAX(IF(l.itemid = 50983, l.valuenum, NULL)) AS sodyum_max,
    MAX(IF(l.itemid = 50971, l.valuenum, NULL)) AS potasyum_max,
    MIN(IF(l.itemid = 50882, l.valuenum, NULL)) AS bikarbonat_min,
    MIN(IF(l.itemid = 51222, l.valuenum, NULL)) AS hemoglobin_min,
    MAX(IF(l.itemid = 51301, l.valuenum, NULL)) AS lokosit_max
  FROM k
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l ON k.hadm_id = l.hadm_id
  WHERE l.itemid IN (51265, 50885, 51006, 50862, 50902, 50983, 50971, 50882, 51222, 51301)
    AND l.valuenum IS NOT NULL
    AND l.charttime BETWEEN DATETIME_SUB(k.intime, INTERVAL 6 HOUR)
                        AND DATETIME_ADD(k.intime, INTERVAL 24 HOUR)
  GROUP BY k.stay_id
),
ilac AS (
  -- Medications, starttime -24 h to 48 h
  SELECT k.stay_id,
    MAX(IF(LOWER(pr.drug) LIKE '%vancomycin%' AND pr.route LIKE 'IV%', 1, 0)) AS vankomisin_iv,
    MAX(IF(REGEXP_CONTAINS(LOWER(pr.drug), r'gentamicin|tobramycin|amikacin'), 1, 0)) AS aminoglikozid,
    MAX(IF(LOWER(pr.drug) LIKE '%piperacillin%', 1, 0)) AS pip_tazo,
    MAX(IF(REGEXP_CONTAINS(LOWER(pr.drug), r'furosemide|bumetanide|torsemide'), 1, 0)) AS loop_diuretik,
    MAX(IF(REGEXP_CONTAINS(LOWER(pr.drug), r'ibuprofen|ketorolac|naproxen'), 1, 0)) AS nsaid
  FROM k
  JOIN `physionet-data.mimiciv_3_1_hosp.prescriptions` pr ON k.hadm_id = pr.hadm_id
  WHERE pr.starttime BETWEEN DATETIME_SUB(k.intime, INTERVAL 24 HOUR)
                         AND DATETIME_ADD(k.intime, INTERVAL 48 HOUR)
  GROUP BY k.stay_id
),
tani AS (
  -- Comorbidities from ICD discharge codes (index hospitalisation)
  SELECT k.stay_id,
    MAX(IF((d.icd_version = 9 AND d.icd_code LIKE '585%')
        OR (d.icd_version = 10 AND d.icd_code LIKE 'N18%'), 1, 0)) AS kbh,
    MAX(IF((d.icd_version = 9 AND d.icd_code LIKE '250%')
        OR (d.icd_version = 10 AND REGEXP_CONTAINS(d.icd_code, r'^E1[0-4]')), 1, 0)) AS diyabet,
    MAX(IF((d.icd_version = 9 AND d.icd_code LIKE '428%')
        OR (d.icd_version = 10 AND d.icd_code LIKE 'I50%'), 1, 0)) AS kky,
    MAX(IF((d.icd_version = 9 AND d.icd_code IN ('5712', '5715', '5716'))
        OR (d.icd_version = 10 AND REGEXP_CONTAINS(d.icd_code, r'^(K703|K717|K74[3-6])')), 1, 0)) AS siroz,
    MAX(IF((d.icd_version = 9 AND (d.icd_code LIKE '038%' OR d.icd_code IN ('99591', '99592', '78552')))
        OR (d.icd_version = 10 AND REGEXP_CONTAINS(d.icd_code, r'^(A40|A41|R652)')), 1, 0)) AS sepsis_icd
  FROM k
  JOIN `physionet-data.mimiciv_3_1_hosp.diagnoses_icd` d ON k.hadm_id = d.hadm_id
  GROUP BY k.stay_id
)
SELECT
  k.stay_id,
  y.first_careunit AS yb_tipi,
  IFNULL(v.vazopressor_48s, 0) AS vazopressor_48s,
  IFNULL(v.vazopressor_sayisi, 0) AS vazopressor_sayisi,
  IFNULL(v.norepinefrin, 0) AS norepinefrin,
  IFNULL(m.imv_48s, 0) AS imv_48s,
  g.sivi_giris_48s,
  c.toplam_cikis_48s,
  c.idrar_48s,
  IFNULL(g.sivi_giris_48s, 0) - IFNULL(c.toplam_cikis_48s, 0) AS sivi_dengesi_48s,
  SAFE_DIVIDE(c.idrar_48s, ce.kilo * 48) AS idrar_ml_kg_saat,
  ce.map_min_24s,
  ce.kilo,
  lab.* EXCEPT(stay_id),
  IFNULL(ilac.vankomisin_iv, 0) AS vankomisin_iv,
  IFNULL(ilac.aminoglikozid, 0) AS aminoglikozid,
  IFNULL(ilac.pip_tazo, 0) AS pip_tazo,
  IFNULL(ilac.loop_diuretik, 0) AS loop_diuretik,
  IFNULL(ilac.nsaid, 0) AS nsaid,
  IFNULL(t.kbh, 0) AS kbh,
  IFNULL(t.diyabet, 0) AS diyabet,
  IFNULL(t.kky, 0) AS kky,
  IFNULL(t.siroz, 0) AS siroz,
  IFNULL(t.sepsis_icd, 0) AS sepsis_icd
FROM k
JOIN `physionet-data.mimiciv_3_1_icu.icustays` y ON k.stay_id = y.stay_id
LEFT JOIN vazo v ON k.stay_id = v.stay_id
LEFT JOIN ventil m ON k.stay_id = m.stay_id
LEFT JOIN giris g ON k.stay_id = g.stay_id
LEFT JOIN cikis c ON k.stay_id = c.stay_id
LEFT JOIN ce ON k.stay_id = ce.stay_id
LEFT JOIN lab ON k.stay_id = lab.stay_id
LEFT JOIN ilac ON k.stay_id = ilac.stay_id
LEFT JOIN tani t ON k.stay_id = t.stay_id;
