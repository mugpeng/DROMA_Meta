library(testthat)

test_root <- normalizePath("DROMA_R2/tests", mustWork = TRUE)
r_dir <- normalizePath(file.path(test_root, "..", "R"), mustWork = TRUE)

for (path in sort(list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE))) {
  source(path, local = FALSE)
}

test_dir(file.path(test_root, "testthat"))
