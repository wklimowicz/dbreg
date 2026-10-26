# Installed tests also check the native validation inventory without credentials.
local({
  directory = system.file("tinytest", package = "dbreg")
  manifest = file.path(directory, "native-suites.txt")
  expect_true(file.exists(manifest))
  if (file.exists(manifest)) {
    suites = readLines(manifest)
    expect_true(length(suites) > 0L)
    expect_identical(anyDuplicated(suites), 0L)
    expect_true(all(grepl("^test_[a-z_]+[.]R$", suites)))
    expect_true(all(file.exists(file.path(directory, suites))))
    live = list.files(directory, pattern = "^test_.*_live[.]R$")
    expect_equal(setdiff(live, suites), character())
    root = normalizePath(file.path(directory, "..", ".."))
    script = file.path(root, ".github", "scripts", "run-databricks-tests.R")
    if (file.exists(script)) {
      previous = setwd(root)
      on.exit(setwd(previous), add = TRUE)
      listed = system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), "--list-suites"), stdout = TRUE, stderr = TRUE)
      expect_equal(listed, suites)
      expect_identical(attr(listed, "status"), NULL)
    }
  }
})
