root = "_benchmarks/sql_prediction_20261024"
files = c("local-final-tail/full-source-md5.csv", "native-tails/source-md5.csv", "large-output-source-md5.csv")
audit = do.call(rbind, lapply(files, function(file) {
  manifest = read.csv(file.path(root, file))
  differences = unname(tools::md5sum(manifest[["file"]])) != manifest[["md5"]]
  if (anyNA(differences) || any(differences)) stop("Source differs from accepted gate: ", file)
  data.frame(manifest = file, files = nrow(manifest), differences = sum(differences), scope = "final implementation")
}))
# The broad native gate preceded only the response-helper/tail correction and
# registration of its focused fixture. Preserve and verify that exact source.
manifest = read.csv(file.path(root, "native-accepted/source-md5.csv"))
differences = unname(tools::md5sum(manifest[["file"]])) != manifest[["md5"]]
changed = manifest[["file"]][differences]
allowed = c("R/sql-prediction.R", "R/sql-prediction-equation.R", ".github/scripts/run-databricks-tests.R")
if (!setequal(changed, allowed)) stop("Unexpected broad-native source differences.")
archived = file.path(root, "native-accepted/pre-tail-source", changed)
if (any(unname(tools::md5sum(archived)) != manifest[["md5"]][differences])) stop("Broad-native source archive differs.")
audit = rbind(audit, data.frame(manifest = "native-accepted/source-md5.csv", files = nrow(manifest),
  differences = sum(differences), scope = "pre-tail response helper; final tails covered separately"))
write.csv(audit, file.path(root, "source-audit.csv"), row.names = FALSE)
print(audit, row.names = FALSE)
