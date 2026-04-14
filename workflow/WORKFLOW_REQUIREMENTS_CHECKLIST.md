# Workflow Requirements Checklist

Date: 2026-04-15

This note is a delivery checklist for the grouped biomarker workflow under `Meta_project/workflow`.

## Requested Requirements

1. Reuse `DB_project/DROMA_R/R` as much as possible.
Status: implemented
Notes: Stage 2 now reuses `DROMA.R::batchFindSignificantFeatures()`. Stage 3 is now a formatter around batch meta output instead of a separate handwritten meta engine.

2. Batch meta-analysis should directly reuse `batchFindSignificantFeatures`.
Status: implemented
Notes: `runWithinStudyScreen()` now delegates screening and meta-analysis to `batchFindSignificantFeatures()`.

3. New code in `Meta_project/DROMA_R2/R` should follow `DROMA_R/R` style and add function comments.
Status: partially implemented
Notes: new helpers were rewritten with roxygen comments and closer naming/style. A full package-wide style sweep is still possible if you want stricter normalization.

4. Biomarker logic should follow `batchFindSignificantFeatures`: expression can come from all matching projects in the same model group, not only the focal drug project.
Status: implemented
Notes: Stage 1 coverage now pools expression/model availability across the full group; Stage 2 uses `batchFindSignificantFeatures()` on the grouped `MultiDromaSet`, so expression traversal follows the same matching logic.

5. Complex functions such as Anderson-Darling filtering should support multithreading similar to `batchFindSignificantFeatures`.
Status: implemented
Notes: TCGA translation filtering now accepts `cores` and uses a parallel row runner when `future` and `furrr` are available.

6. `Meta_project/workflow` scripts should resemble `DB_project/DROMA_Example3/04-Invitro&invivo_flow.R` with more comments.
Status: partially implemented
Notes: bootstrap/common and stages 3-5 were rewritten with section headers and step comments. The remaining stages can be normalized further in a follow-up sweep.

7. `03-Stage2_WithinStudy_Screen.R` should use Spearman instead of Wilcoxon.
Status: implemented
Notes: Stage 2 no longer uses Wilcoxon. It inherits Spearman-based continuous meta-analysis from `batchFindSignificantFeatures()`.

8. Verification should cover `test_top_n`, single-core, and multi-core paths where relevant.
Status: in progress for this delivery
Notes: targeted verification should exercise `batchFindSignificantFeatures(..., test_top_n = ...)` with both `cores = 1` and `cores > 1`.

9. Delivery should include a full rerun as far as the local environment allows.
Status: in progress for this delivery
Notes: the workflow is being rerun stage by stage after logic fixes.
