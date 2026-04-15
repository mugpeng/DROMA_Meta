# DROMA.Meta: Application Package Built on DROMA.R

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: MPL-2.0](https://img.shields.io/badge/License-MPL--2.0-yellow.svg)](https://opensource.org/licenses/MPL-2.0)

## Overview

**DROMA.Meta** is an R application package built on top of **DROMA.R** and **DROMA.Set**.  
Its role is not to replace the core analysis packages, but to package a concrete
multi-step biomarker discovery workflow that combines:

1. preclinical batch screening,
2. preclinical biomarker selection,
3. TCGA/TARGET Anderson-Darling concordance filtering,
4. clinical validation on CTRDB.

In the DROMA ecosystem:

- **DROMA.Set** provides the data structures and database access layer.
- **DROMA.R** provides the statistical analysis functions.
- **DROMA.Meta** is the application-layer workflow package that orchestrates
  those lower-level capabilities into an end-to-end biomarker discovery pipeline.

This makes **DROMA.Meta** an application package for **DROMA.R**, focused on
running a standardized translational meta-workflow rather than exposing a broad
new analysis API.

## Position in the DROMA Stack

```text
DROMA.Set  ->  DROMA.R  ->  DROMA.Meta
data layer     analysis     workflow application layer
```

- **DROMA.Set**: builds and manages `DromaSet` / `MultiDromaSet` objects.
- **DROMA.R**: performs drug-omics association analysis and batch screening.
- **DROMA.Meta**: turns those capabilities into a reusable workflow package for
  preclinical-to-clinical biomarker discovery.

## What This Package Does

- Wraps the four-step DROMA meta workflow into reusable R functions.
- Keeps workflow logic in `R/` as package-style function files.
- Lets users pass explicit runtime parameters such as:
  - `drug`
  - `tumor_type`
  - `data_type`
  - `cores`
  - `db_path`
  - `tcga_rna_counts_dir`
  - `gene_probe_map_path`
  - threshold settings for each workflow step
- Supports script-level looping over drug and tumor-type combinations through
  [one.R](one.R).

## Core Entry Points

### Package Function

The main workflow entry point is:

```r
runMetaWorkflow(
  drug,
  tumor_type,
  feature2_type = "mRNA",
  data_type = "all",
  cores = 3,
  cell_min_intersected_cells = 20,
  pdcpdx_min_intersected_cells = 8,
  db_path,
  ctrdb_path,
  tcga_rna_counts_dir,
  gene_probe_map_path,
  output_base,
  override = FALSE
)
```

This function runs one full workflow for one `drug` and one `tumor_type`, and
returns a one-row summary table while writing intermediate outputs to the
workflow output directory. When `override = FALSE`, existing stage outputs under
`output_base/<drug>/<tumor_type>/` are reused and those stages are skipped.

### Driver Script

[one.R](one.R) is the application driver script.

It is responsible for:

- loading package dependencies,
- loading local function files from `R/`,
- defining runtime parameters,
- reading `valid_drugs.csv` and `valid_tumor_types.csv`,
- looping over all selected combinations,
- writing a batch summary CSV.

This separation keeps the files under `R/` closer to package-style reusable
functions, while leaving actual experiment execution in a script.

## Installation

### Prerequisites

Install the lower-level DROMA packages first:

```r
# devtools::install_github("mugpeng/DROMA_Set")
# devtools::install_github("mugpeng/DROMA_R")
```

### Install DROMA.Meta

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("mugpeng/DROMA_Meta")
```

## Quick Start

### 1. Load the DROMA Stack

```r
library(DROMA.Set)
library(DROMA.R)
```

### 2. Load DROMA.Meta Workflow Functions

```r
source("R/FuncHelper.R")
source("R/FuncValidCheck.R")
source("R/FuncTcgaAD.R")
source("R/FuncMetaWorkflow.R")
```

### 3. Run One Drug and Tumor-Type Workflow

```r
result <- runMetaWorkflow(
  drug = "Paclitaxel",
  tumor_type = "breast cancer",
  feature2_type = "mRNA",
  data_type = "all",
  cores = 3,
  cell_min_intersected_cells = 20,
  pdcpdx_min_intersected_cells = 8,
  db_path = "/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite",
  ctrdb_path = "/Users/peng/Desktop/Project/DROMA/Data/ctrdb.sqlite",
  tcga_rna_counts_dir = "/path/to/tcga/rna_counts",
  gene_probe_map_path = "/path/to/gencode.human.v49.annotation.gene.probeMap",
  output_base = "workflow/Output",
  override = FALSE
)

print(result)
```

### 4. Batch Run by Script-Level Looping

Use [one.R](one.R) as the application script.  
Edit the runtime parameters in that file, then run:

```r
Rscript one.R
```

## Repository Structure

```text
R/
  FuncHelper.R
  FuncValidCheck.R
  FuncTcgaAD.R
  FuncMetaWorkflow.R
workflow/
  00-Eligible_Drug_Tumor.R
  01-Batch_Preclinical.R
  02-Select_Preclinical.R
  03-TCGA_AD_Filter.R
  04-Clinical_Validation.R
one.R
DESCRIPTION
NAMESPACE
```

## Intended Use

Use **DROMA.Meta** when you already rely on **DROMA.R** for analysis and want a
reproducible application package for a fixed biomarker discovery workflow.

Use **DROMA.R** directly when you need lower-level analysis functions, custom
plots, or exploratory workflows.

## Related Packages

- [DROMA.Set](https://github.com/mugpeng/DROMA_Set): data structures and database layer
- [DROMA.R](https://github.com/mugpeng/DROMA_R): statistical analysis layer
- [DROMA.Meta](https://github.com/mugpeng/DROMA_Meta): workflow application layer
