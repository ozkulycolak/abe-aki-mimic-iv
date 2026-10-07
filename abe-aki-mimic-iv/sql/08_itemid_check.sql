-- =====================================================================
-- 08_itemid_check.sql
-- Labels of the itemids used for KRT, invasive ventilation, vasopressors,
-- GCS and FiO2, and of the crystalloid items, for verification.
-- =====================================================================
SELECT itemid, label, category, linksto
FROM `physionet-data.mimiciv_3_1_icu.d_items`
WHERE itemid IN (225441, 225802, 225803, 225805, 225809, 225955, 225792,
                 221906, 221289, 222315, 221749, 221662, 221653,
                 220739, 223900, 223901, 223835, 226512)
   OR LOWER(label) = 'nacl 0.9%'
   OR LOWER(label) IN ('lr', 'lactated ringers')
   OR REGEXP_CONTAINS(LOWER(label), r'plasma ?-?lyte')
ORDER BY itemid;
