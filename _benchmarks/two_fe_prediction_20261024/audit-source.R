root = "_benchmarks/two_fe_prediction_20261024"
manifests = c("local/full-source-md5.csv", "native-two-fe/source-md5.csv", "native-single-fe/source-md5.csv",
  "performance/none-source-md5.csv", "performance/prediction-source-md5.csv")
results = lapply(manifests, function(file) {
  manifest = read.csv(file.path(root, file))
  actual = unname(tools::md5sum(manifest[["file"]]))
  changed = is.na(actual) | actual != manifest[["md5"]]
  if (any(changed)) stop("Source differs for ", file, ": ", paste(manifest[["file"]][changed], collapse = ", "))
  data.frame(manifest = file, files = nrow(manifest), differences = sum(changed))
})
audit = do.call(rbind, results)
write.csv(audit, file.path(root, "source-audit.csv"), row.names = FALSE)
print(audit, row.names = FALSE)
