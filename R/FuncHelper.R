# Helper functions for DROMA.Meta ----

#' Null-Coalescing Helper
#'
#' @description Returns `y` when `x` is NULL, empty, or a blank string.
#' This is useful when package helpers need a compact fallback operator.
#' @param x Primary value.
#' @param y Fallback value.
#' @return Either `x` or `y`.
#' @export
`%||%` <- function(x, y) {
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
#' @export
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
