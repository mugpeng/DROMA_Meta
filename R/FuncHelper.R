`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "") || (is.character(x) && !nzchar(x[1]))) {
    y
  } else {
    x
  }
}

sanitizeName <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  if (!nzchar(out)) {
    "unknown"
  } else {
    out
  }
}

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
