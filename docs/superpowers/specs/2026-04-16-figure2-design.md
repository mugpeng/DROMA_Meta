# Figure 2 Design

## Goal

Build a DROMA-specific Figure 2 analysis pipeline under `v2_part/workflow` that measures whether meta-derived biomarkers are functionally close to known drug targets and whether target proximity is associated with biomarker strength or validation depth.

The pipeline will use existing meta outputs from:

- `/Users/peng/Desktop/Project/DROMA/Meta_project/Meta_Example/Output/meta_batch`

and will compare three biomarker layers:

- `selected_genes`
- `selected_genes_ad_filtered`
- `final_biomarkers`

The main narrative should prioritize `final_biomarkers`, while always retaining the other two layers as supporting context.

## Scope

This design covers only Figure 2A to Figure 2C.

Included:

- standardizing figure-2 inputs from meta batch outputs
- collecting or normalizing drug-target annotations under `v2_part`
- building a Reactome-derived gene network
- computing shortest-path target proximity for biomarkers and matched background genes
- generating three figure panels and intermediate tables
- writing workflow-style scripts under `v2_part/workflow`

Excluded:

- Figure 2D to Figure 2G
- DepMap CRISPR/RNAi integration
- manuscript assembly outside figure-2 outputs

## Outputs

All figure-2 assets will live under `v2_part`.

Recommended structure:

- `v2_part/workflow`
- `v2_part/data`
- `v2_part/output/figure2/tables`
- `v2_part/output/figure2/figures`

Primary figure outputs:

- `figure2A_distance_distribution.pdf`
- `figure2B_distance_vs_effect.pdf`
- `figure2C_distance_vs_validation_stage.pdf`

Primary table outputs:

- `figure2_input_features.csv`
- `drug_target_map.csv`
- `reactome_edges.csv`
- `distance_results.csv`
- `background_distance_results.csv`
- `distance_summary_by_stage.csv`

## Data Model

### Biomarker input table

The pipeline will build one unified long-format table with at least:

- `drug`
- `tumor_type`
- `stage` with values `selected_genes`, `selected_genes_ad_filtered`, `final_biomarkers`
- `name`
- `effect_size`
- `effect_source`
- `has_target_mapping`
- `has_network_mapping`

Stage-specific effect-size selection:

- `selected_genes`: prefer `effect_size_pdcpdx`, fallback `effect_size_cell`
- `selected_genes_ad_filtered`: prefer `effect_size_pdcpdx`, fallback `effect_size_cell`
- `final_biomarkers`: prefer `effect_size_ctrdb`, fallback `effect_size_pdcpdx_invitro`, then `effect_size_cell_invitro`

This preserves comparability while keeping `final_biomarkers` tied to the clinical layer whenever available.

### Drug-target table

The pipeline will standardize a single target table under `v2_part/data` with at least:

- `drug`
- `target_gene`
- `target_source`
- `evidence_level`

Drugs may map to multiple targets. The distance metric will use the minimum shortest path from a biomarker gene to any target gene for that drug.

### Reactome network

The network will be built as an undirected gene-gene graph where two genes are connected if they co-occur in at least one Reactome pathway.

Stored outputs:

- pathway membership cache
- edge list
- optional node list

The graph must be deterministic and reusable across runs.

## Figure Definitions

### Figure 2A

Distance distribution panel comparing:

- biomarker to target distances for each stage
- random gene to target distances
- random gene to random gene distances

Purpose:

- show whether DROMA-prioritized biomarkers are closer to known targets than background genes

Suggested display:

- violin plus boxplot or density plus boxplot
- stage-aware facets or grouped comparisons

### Figure 2B

Distance versus effect-size panel.

Primary emphasis:

- `final_biomarkers`

Supporting comparisons:

- `selected_genes`
- `selected_genes_ad_filtered`

Suggested display:

- boxplot or jittered boxplot of effect size by distance bin
- optional facet by stage

Purpose:

- test whether functionally proximal biomarkers tend to show stronger signal

### Figure 2C

Distance versus validation-stage panel.

Compare distance distributions across:

- `selected_genes`
- `selected_genes_ad_filtered`
- `final_biomarkers`

Purpose:

- test whether target proximity is enriched in biomarkers that survive deeper validation

Suggested display:

- violin or boxplot by stage
- annotate counts and unmapped fractions

## Workflow Scripts

The implementation should follow the existing repository pattern of numbered scripts.

Planned scripts:

1. `21-Collect_Figure2_Inputs.R`
2. `22-Build_Target_Network.R`
3. `23-Compute_Target_Distance.R`
4. `24-Plot_Figure2A_Distance_Distribution.R`
5. `25-Plot_Figure2B_Distance_vs_Effect.R`
6. `26-Plot_Figure2C_Distance_vs_ValidationStage.R`

Optional helper files if the logic grows:

- `v2_part/R/FuncFigure2Helper.R`
- `v2_part/R/FuncFigure2Plot.R`

## Processing Logic

### Step 1: collect inputs

Scan `meta_batch` recursively and load:

- `selected_genes.csv`
- `selected_genes_ad_filtered.csv`
- `final_biomarkers.csv`

Normalize column names into one long-format figure-2 input table.

### Step 2: normalize drug targets

Load user-provided or locally curated target resources from `v2_part`.

Normalize drug names conservatively:

- exact match first
- controlled manual alias map second

No fuzzy auto-match should be used without an explicit alias table.

### Step 3: build Reactome graph

Construct a reusable undirected edge list from Reactome pathway membership.

Only genes present in the chosen Reactome resource should appear in the graph.

### Step 4: compute distances

For each `(drug, tumor_type, stage, gene)` row:

- obtain all mapped targets for the drug
- compute shortest path from the gene to all mapped targets
- keep the minimum finite distance
- store mapping status separately from distance value

Also compute background distributions with a fixed random seed.

### Step 5: summarize and plot

Generate per-stage and pooled summaries, then render Figure 2A to 2C.

## Error Handling

The pipeline should not silently drop hard cases.

Track and report:

- drugs with no target annotation
- genes not found in the graph
- targets not found in the graph
- pairs with empty stage files
- figure panels skipped due to no valid data

Missing mappings must be represented explicitly in summary tables.

## Reproducibility

- fixed random seed for all background sampling
- deterministic graph build
- workflow scripts rerunnable independently
- intermediate CSV caches written before plotting

## Testing and Verification

Verification should be lightweight but explicit:

- check that figure-2 input collection reads all three stage file types
- check that multi-target drugs produce minimum distance, not first-match distance
- check that unmapped genes and drugs are counted correctly
- check that each plotting script fails clearly when required input tables are missing

At minimum, run the workflow on a small subset of drug-tumor pairs before full execution.

## Recommendation

Proceed with the six-script workflow under `v2_part/workflow`, with `final_biomarkers` as the main figure story and the other two layers as context. This gives the strongest result now and leaves a clean path to add Figure 2D to 2G later without rewriting the early pipeline.
