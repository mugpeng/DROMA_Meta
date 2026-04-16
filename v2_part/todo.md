# Figure 2 TODO

This file tracks the parts of Smirnov-style Figure 2 that are intentionally not included in the current implementation scope.

## In Scope Now

- Figure 2A: biomarker-to-target distance distribution
- Figure 2B: distance versus effect size
- Figure 2C: distance versus validation stage

## Not Yet Implemented

### Figure 2D to Figure 2G

These panels require perturbation-layer data that are not part of the current DROMA meta workflow.

Missing components:

- drug-to-target mappings verified at perturbation level
- DepMap CRISPR gene dependency data
- DepMap RNAi gene dependency data
- logic linking biomarker expression to target knock-out or knock-down response
- comparison framework for:
  - biomarker groups associated versus not associated with CRISPR/RNAi target response
  - distance distributions under perturbation support
  - effect-size differences under perturbation support

## Recommended Next Steps

1. Curate a stable drug-target reference table under `v2_part/data`.
2. Add DepMap CRISPR and RNAi input resources under `v2_part/data`.
3. Define one harmonized target symbol normalization layer shared by:
   - meta biomarker results
   - Reactome graph
   - target annotation
   - DepMap dependency resources
4. Implement a perturbation support table for each biomarker:
   - `supports_crispr_target_response`
   - `supports_rnai_target_response`
   - supporting statistics and effect directions
5. Add workflow scripts for Figure 2D to 2G, likely:
   - `27-Collect_Perturbation_Inputs.R`
   - `28-Compute_Perturbation_Support.R`
   - `29-Plot_Figure2D_2G.R`

## Open Technical Risks

- drug naming may not match cleanly across target and perturbation resources
- target genes may be absent from one or more resources after symbol harmonization
- multi-target drugs need explicit minimum-distance and multi-perturbation rules
- CRISPR and RNAi may disagree in direction or significance
- some final biomarkers may have no usable perturbation-layer counterpart

## Acceptance Condition For Future Work

Only start Figure 2D to Figure 2G after the following are available:

- a reviewed drug-target map
- a reviewed gene-symbol harmonization layer
- local CRISPR and RNAi tables with reproducible load scripts
- at least one successful dry run on a subset of drugs
