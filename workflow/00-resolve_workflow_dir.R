#' Resolve the directory that contains 00-Workflow_Common.R (the workflow folder).
resolve_workflow_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }
  for (f in rev(sys.frames())) {
    if (!exists("ofile", envir = f, inherits = FALSE)) {
      next
    }
    ofile <- tryCatch(get("ofile", envir = f, inherits = FALSE), error = function(e) NULL)
    if (is.null(ofile) || !length(ofile) || !nzchar(ofile[1])) {
      next
    }
    if (file.exists(ofile[1])) {
      return(dirname(normalizePath(ofile[1], mustWork = TRUE)))
    }
  }
  wd <- getwd()
  for (d in c(file.path(wd, "workflow"), wd)) {
    common <- file.path(d, "00-Workflow_Common.R")
    if (file.exists(common)) {
      return(normalizePath(d, mustWork = TRUE))
    }
  }
  stop(
    "Could not find 00-Workflow_Common.R. setwd() to .../Meta_project or .../Meta_project/workflow, ",
    "or run with Rscript --file=path/to/script.R",
    call. = FALSE
  )
}
