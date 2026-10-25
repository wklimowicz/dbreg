root = "_benchmarks/formula_i_20261024"
full = read.csv(file.path(root,"local-final-preflight/full-summary.csv"))
native = read.csv(file.path(root,"native-acceptance/summary.csv"))
stopifnot(nrow(full)==1L,full[["assertions"]]==5418L,full[["failures"]]==0L,
 nrow(native)==2L,all(native[["assertions"]]>0L),all(native[["failures"]]==0L))
check = readLines(file.path(root,"package-check-final.txt"))
stopifnot(any(grepl("Status: 1 NOTE",check,fixed=TRUE)),!any(grepl("ERROR|WARNING",check)))
performance = read.csv(file.path(root,"performance-summary.csv"))
stopifnot(nrow(performance)==3L)
lines = c(
 "Supported i() formula coding and final compatibility acceptance",
 "=============================================================",
 paste("UTC:",format(Sys.time(),tz="UTC")),
 paste("Parent revision:",system2("git",c("rev-parse","HEAD"),stdout=TRUE)),
 "Final source and the native validation delta are identified by local-final-preflight/native-acceptance/preflight manifests and audit-source.R.",
 "Reference: installed fixest 0.14.0; adjacent 0.14.2 source 789f888416e97673d154c840eec2d99101fc07b1.",
 "Development uses pkgload::load_all(). No additional essential Imports.",
 "",
 paste0("Full local: ",full[["assertions"]]," assertions, zero failures; ",round(full[["elapsed_seconds"]],2)," seconds."),
 paste0("Native i(): ",native[["assertions"]][1]," assertions, zero failures; ",round(native[["elapsed_seconds"]][1],2)," seconds."),
 paste0("Native adapters: ",native[["assertions"]][2]," assertions, zero failures; ",round(native[["elapsed_seconds"]][2],2)," seconds."),
 "Package check: zero errors/warnings; existing optional binsreg cross-reference NOTE.",
 "Package tests were skipped in that check; the complete source tinytest suite ran separately.",
 "Tinytest's progress label incorrectly displays minutes as seconds for these native runs.",
 "The elapsed_seconds fields above are the runner's wall-clock seconds.",
 "",
 "The shared fixture compares named coefficients, full covariance, likelihood,
 R/SQL/fitted-sample predictions, post-fit inference and cleanup with fixest.",
 "It includes indicators, numeric/logical group slopes, literal controls, two FEs,
 varying slopes, exact BIGINT values above 2^53 and adjacent DOUBLE categories.",
 "R/SQL prediction agrees when an unknown category accompanies an unknown FE or
 missing regressor. Otherwise valid unseen/malformed categories error; lossy
 integer inputs are rejected. Typed Date/enum/integer64 and at-budget probes pass.",
 "Local-only tests add Gaussian/Poisson/logit/probit, no-intercept/multiple terms,
 serialization/lean release, unsupported forms, resource limits, ambiguous public
 names and a nonsyntactic first-column name. Native shared data includes quoted
 text containing an internal-looking alias, which must remain unchanged in labels.",
 "",
 "Initial full validation caught internal FE aliases being passed as formula terms.",
 "The fix separates the i() absorbed-intercept context from ordinary factor coding.",
 "Later review caught SQL data-frame prediction rejecting category keys before
 applying native R missing/unknown-FE eligibility. The saved marker now lets the
 existing SQL validity checks decide in the same order; unsupported precision
 cannot be rounded into an apparently valid fitted category.",
 "All earlier sources are archived exactly where they differ. audit-source.R
 verifies historical manifests against those files. Final local and preflight
 manifests match the current tree. Native source differs only by a later guard
 rejecting unsupported i interactions before setup, plus two local-only tests;
 the audit proves that exact delta and current preflight probes cover both APIs. Earlier native 92+118 and 116+118 passes are historical gates,
 not a claim that those runs exercised the later eligibility correction.",
 "",
 "Million-row resource fixture: one warmup and three measured fits per mode;
 fresh process per mode, DuckDB threads=2/memory_limit=2GB, explicit id in all modes.",
 "Native baseline uses y~x+z+f moments; adapters use equivalent y~x+z+i(f).
 The source has three f levels derived from id and the existing paired-noise panel.",
 capture.output(print(performance,row.names=FALSE)),
 "",
 "Timing includes fit-time lean release; explicit full close is outside fit timing.",
 "API read/execute counts measure DBI calls during fitting, not every driver's
 internal network operation. R query object bytes are not transport byte counts.",
 "Every completed fit matches native coefficients, full V, N and likelihood.
 At most three rows/query cross into R. Closing each fit restores the table
 inventory and preserves the caller connection. Post-close RSS varies with
 embedded DuckDB/allocator caching; three repetitions do not prove zero long-run
 memory growth. Peak RSS includes DuckDB and temporary row-sized state.",
 "i()/factor-coercion/varying-slope OLS uses the shared Gaussian projector.
 It can cost substantially more than ordinary compact moments. Full compatibility
 retains rows/scores in SQL; lean transiently retains inference rows and avoids
 learned-FE recovery. Existing native compact fast paths are preserved.",
 "These warmed local timings do not establish warehouse or disk-bound costs.",
 "",
 "Reproduction from repository root (use a new output database path):",
 "Rscript _benchmarks/two_fe_prediction_20261024/large-state.R prepare /tmp/dbreg-i-repro.duckdb",
 "Rscript _benchmarks/formula_i_20261024/large-fit.R /tmp/dbreg-i-repro.duckdb baseline",
 "Rscript _benchmarks/formula_i_20261024/large-fit.R /tmp/dbreg-i-repro.duckdb lean",
 "Rscript _benchmarks/formula_i_20261024/large-fit.R /tmp/dbreg-i-repro.duckdb full",
 "Rscript _benchmarks/check_glm_fractional.R full OUTPUT_DIRECTORY",
 "Rscript .github/scripts/run-databricks-tests.R OUTPUT_DIRECTORY test_formula_i_live.R test_compat_api_live.R",
 "Rscript _benchmarks/formula_i_20261024/check-typed.R",
 "Rscript _benchmarks/formula_i_20261024/audit-source.R",
 "",
 "Support table and executable lifecycle examples: vignettes/intro.qmd and
 inst/examples/fixest-compatibility.R. Installed statistical/lifetime contract:
 inst/FIXEST-COMPATIBILITY.txt. The API is the documented regression subset;
 external cluster mappings, slope-only FEs, unsupported i() forms, unrestricted
 R expressions and oversized 3+ FE/varying-slope prediction recovery remain
 explicit boundaries. The separate wide-AP-transfer TODO is unchanged."
)
writeLines(lines,file.path(root,"results.txt"))
