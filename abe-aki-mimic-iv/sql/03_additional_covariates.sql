-- =====================================================================
-- 03_additional_covariates.sql
-- Sodium bicarbonate, blood products (-6 h to 48 h) and acetazolamide
-- (-24 h to 48 h). Returns one row per stay_id; read directly in Python.
-- Replace YOUR_PROJECT_ID with your own Google Cloud project ID.
-- =====================================================================
WITH k AS (
  SELECT stay_id, hadm_id, intime FROM `YOUR_PROJECT_ID.abe.kohort`
),
ie AS (
  SELECT k.stay_id,
    MAX(IF(LOWER(d.label) LIKE '%bicarbonate%', 1, 0)) AS bikarbonat_inf,
    MAX(IF(REGEXP_CONTAINS(LOWER(d.label), r'packed red|prbc'), 1, 0)) AS eritrosit_tx,
    MAX(IF(REGEXP_CONTAINS(LOWER(d.label), r'fresh frozen|ffp'), 1, 0)) AS taze_plazma_tx
  FROM k
  JOIN `physionet-data.mimiciv_3_1_icu.inputevents` x ON k.stay_id = x.stay_id
  JOIN `physionet-data.mimiciv_3_1_icu.d_items` d ON x.itemid = d.itemid
  WHERE x.starttime >= DATETIME_SUB(k.intime, INTERVAL 6 HOUR)
    AND x.starttime <  DATETIME_ADD(k.intime, INTERVAL 48 HOUR)
  GROUP BY k.stay_id
),
asz AS (
  SELECT DISTINCT k.stay_id, 1 AS asetazolamid
  FROM k
  JOIN `physionet-data.mimiciv_3_1_hosp.prescriptions` pr ON k.hadm_id = pr.hadm_id
  WHERE LOWER(pr.drug) LIKE '%acetazolamide%'
    AND pr.starttime BETWEEN DATETIME_SUB(k.intime, INTERVAL 24 HOUR)
                         AND DATETIME_ADD(k.intime, INTERVAL 48 HOUR)
)
SELECT k.stay_id,
  IFNULL(ie.bikarbonat_inf, 0) AS bikarbonat_inf,
  IFNULL(ie.eritrosit_tx, 0) AS eritrosit_tx,
  IFNULL(ie.taze_plazma_tx, 0) AS taze_plazma_tx,
  IFNULL(asz.asetazolamid, 0) AS asetazolamid
FROM k
LEFT JOIN ie ON k.stay_id = ie.stay_id
LEFT JOIN asz ON k.stay_id = asz.stay_id;
