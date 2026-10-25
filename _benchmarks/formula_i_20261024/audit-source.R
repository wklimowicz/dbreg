# Verify each recorded source against the final tree or an exact archived file.
root = "_benchmarks/formula_i_20261024"
archives = list.files(root,recursive=TRUE,full.names=TRUE)
archives = archives[grepl("/pre-[^/]+/",archives) & !dir.exists(archives)]
archived_hashes = tools::md5sum(archives)
manifests = list.files(root,pattern="source-md5[.]csv$",recursive=TRUE,full.names=TRUE)
records = list()
for (manifest in manifests) {
 rows = read.csv(manifest,stringsAsFactors=FALSE)
 for (j in seq_len(nrow(rows))) {
  file = rows[["file"]][j]
  tested = rows[["md5"]][j]
  current = if(file.exists(file)) unname(tools::md5sum(file)) else NA_character_
  same = !is.na(current) && identical(tested,current)
  matched = names(archived_hashes)[which(archived_hashes==tested)]
  if(!same && !length(matched)) stop("Missing exact source archive: ",manifest," / ",file)
  records[[length(records)+1L]] = data.frame(manifest=manifest,file=file,tested_md5=tested,current_md5=current,
    verified=if(same) "current" else "archived",archive=if(same) "" else matched[[1L]])
 }
}
result = do.call(rbind,records)
write.csv(result,file.path(root,"source-audit.csv"),row.names=FALSE)
print(with(result,table(manifest,verified)))
final = grepl("/(local-final-preflight|preflight)/",result[["manifest"]])
stopifnot(any(final),all(result[["verified"]][final]=="current"))
# The last change only rejects unsupported i interactions before setup. It is
# inactive for every valid native model. Verify that exact, narrowly scoped delta.
native = grepl("/native-acceptance/",result[["manifest"]])
changed = result[["file"]][native & result[["verified"]]=="archived"]
stopifnot(setequal(changed,c("R/glm-sample.R","inst/tinytest/test_formula_i.R")))
current = readLines("R/glm-sample.R")
line = which(grepl("contains_i_term(e))",current,fixed=TRUE))
stopifnot(length(line)==1L,identical(current[-(line:(line+2L))],
 readLines(file.path(root,"pre-additive-guard-source/R/glm-sample.R"))))
current = readLines("inst/tinytest/test_formula_i.R")
added = grepl("expect_error(dbreg::feols(y~i(f)",current,fixed=TRUE) & grepl('table="absent"',current,fixed=TRUE)
stopifnot(sum(added)==2L,identical(current[!added],
 readLines(file.path(root,"pre-additive-guard-source/inst/tinytest/test_formula_i.R"))))
message("Native delta is exactly the pre-query interaction guard plus two local-only regressions; current preflight tests cover both APIs.")
