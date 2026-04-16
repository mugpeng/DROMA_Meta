#' @keywords internal
#' @importFrom grDevices dev.off pdf
#' @importFrom rlang .data
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".", ":=", ".N", ".SD", "N", "Stage", "concordant", "count",
    "direction", "drug", "effect_size", "end", "es", "group", "name", "pair", "pct",
    "pct_retained", "stage", "start", "total", "tumor_type", "x"
  ))
}
