-- =====================================================================
-- 07_extra_covariates.sql
-- Cardiac surgery service, surgical service and crystalloid type
-- (0.9% saline vs balanced crystalloids, 0-48 h). One row per stay_id;
-- read directly in Python. Replace YOUR_PROJECT_ID.
-- =====================================================================
WITH k AS (
  SELECT stay_id, hadm_id, intime FROM `YOUR_PROJECT_ID.abe.kohort`
),
srv AS (
  SELECT k.stay_id,
    MAX(IF(s.curr_service = 'CSURG', 1, 0)) AS kalp_cerrahisi,
    MAX(IF(s.curr_service IN ('SURG', 'CSURG', 'NSURG', 'TSURG', 'VSURG', 'PSURG',
                              'ORTHO', 'TRAUM', 'GU', 'ENT', 'DENT', 'EYE', 'GYN', 'OBS'), 1, 0)) AS cerrahi_servis
  FROM k
  JOIN `physionet-data.mimiciv_3_1_hosp.services` s ON k.hadm_id = s.hadm_id
  WHERE s.transfertime <= DATETIME_ADD(k.intime, INTERVAL 24 HOUR)
  GROUP BY k.stay_id
),
sivi AS (
  SELECT k.stay_id,
    SUM(IF(LOWER(d.label) = 'nacl 0.9%', x.amount, 0)) AS salin_ml,
    SUM(IF(LOWER(d.label) IN ('lr', 'lactated ringers')
           OR REGEXP_CONTAINS(LOWER(d.label), r'plasma ?-?lyte'), x.amount, 0)) AS dengeli_ml
  FROM k
  JOIN `physionet-data.mimiciv_3_1_icu.inputevents` x ON k.stay_id = x.stay_id
  JOIN `physionet-data.mimiciv_3_1_icu.d_items` d ON x.itemid = d.itemid
  WHERE LOWER(x.amountuom) = 'ml'
    AND x.starttime >= k.intime
    AND x.starttime <  DATETIME_ADD(k.intime, INTERVAL 48 HOUR)
  GROUP BY k.stay_id
)
SELECT k.stay_id,
  IFNULL(srv.kalp_cerrahisi, 0) AS kalp_cerrahisi,
  IFNULL(srv.cerrahi_servis, 0) AS cerrahi_servis,
  IFNULL(sivi.salin_ml, 0) AS salin_ml,
  IFNULL(sivi.dengeli_ml, 0) AS dengeli_ml
FROM k
LEFT JOIN srv ON k.stay_id = srv.stay_id
LEFT JOIN sivi ON k.stay_id = sivi.stay_id;
