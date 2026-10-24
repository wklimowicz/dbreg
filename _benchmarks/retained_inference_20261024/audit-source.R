root = "_benchmarks/retained_inference_20261024"
files = c("local-final/full-source-md5.csv", "native-numerics/source-md5.csv",
  "prediction-source-md5.csv", "all-source-md5.csv")
audit = do.call(rbind, lapply(files, function(file) {
  manifest = read.csv(file.path(root, file))
  differences = unname(tools::md5sum(manifest[["file"]])) != manifest[["md5"]]
  if (anyNA(differences) || any(differences)) stop("Source differs from accepted gate: ", file)
  data.frame(manifest = file, files = nrow(manifest), differences = sum(differences), scope = "final implementation")
}))
manifest = read.csv(file.path(root, "native/source-md5.csv"))
differences = unname(tools::md5sum(manifest[["file"]])) != manifest[["md5"]]
allowed = c(".github/scripts/run-databricks-tests.R", "R/inference-vcov.R", "R/sample-identity.R",
  "R/score-moments.R", "inst/tinytest/helper_retained_inference.R", "inst/tinytest/test_retained_inference_local.R")
changed = manifest[["file"]][differences]
if (anyNA(differences) || !setequal(changed, allowed)) stop("Unexpected broad-native source differences.")
archived = file.path(root, "native/pre-numerics-source", changed)
if (any(unname(tools::md5sum(archived)) != manifest[["md5"]][differences])) stop("Broad-native source archive differs.")
audit = rbind(audit, data.frame(manifest = "native/source-md5.csv", files = nrow(manifest),
  differences = sum(differences), scope = "before numerical stabilization, NaN key rejection and test namespace qualification"))
write.csv(audit, file.path(root, "source-audit.csv"), row.names = FALSE)
print(audit, row.names = FALSE)
