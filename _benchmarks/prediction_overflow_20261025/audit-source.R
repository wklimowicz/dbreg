root = "_benchmarks/prediction_overflow_20261025"
manifests = list.files(root,pattern="source-md5[.]csv$",recursive=TRUE,full.names=TRUE)
for (manifest in manifests) {
 rows = read.csv(manifest)
 stopifnot(all(unname(tools::md5sum(rows[["file"]]))==rows[["md5"]]))
}
# The preceding broad native acceptance belongs to the prior committed code.
# Verify it against that commit, including its exact archived preflight delta.
previous = "447f92c8cd802cbf444f7206c1b0f41c90ed7192"
old_root = "_benchmarks/formula_i_20261024"
records = list()
for (manifest in c(file.path(old_root,"local-final-preflight/full-source-md5.csv"),
 file.path(old_root,"native-acceptance/source-md5.csv"))) {
 rows = read.csv(manifest)
 for (j in seq_len(nrow(rows))) {
  file = rows[["file"]][j]
  candidates = c(file,file.path(old_root,"pre-additive-guard-source",file))
  matched = FALSE
  for (candidate in candidates) {
   output = tempfile()
   errors = tempfile()
   status = system2("git",c("show",shQuote(paste0(previous,":",candidate))),stdout=output,stderr=errors)
   if(status==0L && identical(unname(tools::md5sum(output)),rows[["md5"]][j])) matched=TRUE
   unlink(c(output,errors))
   if(matched) break
  }
  stopifnot(matched)
  records[[length(records)+1L]] = data.frame(manifest=manifest,file=file,verified_commit=previous,source=candidate)
 }
}
write.csv(do.call(rbind,records),file.path(root,"previous-source-audit.csv"),row.names=FALSE)
message("Current follow-up manifests match; preceding acceptance is verified against commit 447f92c and its exact archives.")
