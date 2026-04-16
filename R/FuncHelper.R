# Helper functions for DROMA.Meta ----

fallbackIfMissing <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "") || (is.character(x) && !nzchar(x[1]))) {
    y
  } else {
    x
  }
}

#' Sanitize a Label for File-System Use
#'
#' @description Converts arbitrary text into a conservative slug composed of
#' letters, numbers, and underscores. Used for stable workflow output folders.
#' @param x Character vector of labels to sanitize.
#' @return A length-one character string safe to use in file paths.
#' @export
sanitizeName <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  if (!nzchar(out)) {
    "unknown"
  } else {
    out
  }
}

#' Create an Empty Meta-analysis Result Frame
#'
#' @description Builds an empty result data frame with the expected columns used
#' by downstream workflow helpers when a step yields no usable result.
#' @return A zero-row data.frame with standard meta-analysis columns.
#' @keywords internal
createEmptyMetaDf <- function() {
  data.frame(
    name = character(0),
    effect_size = numeric(0),
    p_value = numeric(0),
    q_value = numeric(0),
    n_datasets = integer(0),
    stringsAsFactors = FALSE
  )
}

#' Create an Empty Final-biomarker Result Frame
#'
#' @description Builds an empty result data frame for the final biomarker export
#' so downstream scripts can still read headers when no biomarker survives.
#' @return A zero-row data.frame with standard final-biomarker columns.
#' @keywords internal
createEmptyFinalBiomarkersDf <- function() {
  data.frame(
    name = character(0),
    drug = character(0),
    tumor_type = character(0),
    direction = character(0),
    ctrdb_fallback = logical(0),
    ctrdb_status = character(0),
    stringsAsFactors = FALSE
  )
}
