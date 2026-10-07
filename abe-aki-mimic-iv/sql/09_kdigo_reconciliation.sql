-- =====================================================================
-- 09_kdigo_reconciliation.sql
-- Compares the version-2 fixed-baseline definitions with the KDIGO rolling
-- definitions (version 2.1) in abe.kohort_genis. Aggregate counts only.
-- Replace YOUR_PROJECT_ID with your own Google Cloud project ID.
-- =====================================================================
WITH g AS (
  SELECT *,
    (kre_max_0_48 >= 1.5 * kre_bazal OR kre_max_0_48 - kre_bazal >= 0.3) AS evre1_sabit,
    (sonlanim_sabit = 1 AND krt_sonlanim = 0 AND mutlak_sadece_sabit = 1) AS mutlak_olay_sabit
  FROM `YOUR_PROJECT_ID.abe.kohort_genis`
)
SELECT
  COUNTIF(NOT dis_sabit) AS n_kohort_v2,
  COUNTIF(NOT dis_kdigo) AS n_kohort_v21,
  COUNTIF(NOT dis_sabit AND NOT dis_kdigo) AS n_ortak,
  COUNTIF(NOT dis_sabit AND dis_kdigo) AS yalniz_v2,
  COUNTIF(dis_sabit AND NOT dis_kdigo) AS yalniz_v21,
  COUNTIF(NOT dis_sabit AND sonlanim_sabit = 1) AS olay_v2,
  COUNTIF(NOT dis_sabit AND mutlak_olay_sabit) AS v2_yalniz_mutlak_olcut_olay,
  COUNTIF(NOT dis_sabit AND mutlak_olay_sabit AND mutlak_akut_flag = 1) AS bunlardan_akut_mutlak_olcut,
  COUNTIF(NOT dis_sabit AND mutlak_olay_sabit AND NOT dis_kdigo AND sonlanim_aki23_rrt = 1) AS bunlardan_v21_olay,
  COUNTIF(NOT dis_kdigo AND evre1_sabit AND evre1_48s_kdigo = 1) AS evre1_ikisinde,
  COUNTIF(NOT dis_kdigo AND evre1_sabit AND evre1_48s_kdigo = 0) AS evre1_yalniz_sabit,
  COUNTIF(NOT dis_kdigo AND NOT evre1_sabit AND evre1_48s_kdigo = 1) AS evre1_yalniz_kdigo,
  COUNTIF(NOT dis_kdigo AND sonlanim_aki23_rrt = 1 AND sonlanim_sabit = 1) AS olay_ikisinde,
  COUNTIF(NOT dis_kdigo AND sonlanim_aki23_rrt = 1 AND sonlanim_sabit = 0) AS olay_yalniz_kdigo,
  COUNTIF(NOT dis_kdigo AND sonlanim_aki23_rrt = 0 AND sonlanim_sabit = 1) AS olay_yalniz_sabit
FROM g;
