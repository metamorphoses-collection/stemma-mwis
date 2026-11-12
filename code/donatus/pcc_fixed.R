# Camps & Cafiero PCC, automatic mode, with the crash path repaired (local wrapper; package untouched).
# stemmatology 0.3.2 crashes when PCC.equipollent returns no database (1:length() on an empty list).
# Instead of the equipollent branch, this applies the method's own steps until the base is conflict-free:
#   1. PCC.overconflicting + PCC.elimination at the same threshold, repeated while it removes something;
#   2. if conflicts remain and nothing more is removed: drop every variant location still in conflict
#      (Camps & Cafiero 2015, step 3a: conflicting 'potentially reliable' locations are both deleted).
# Then PCC.Stemma on the conflict-free base. Usage:
#   Rscript pcc_fixed.R <matrix.num> <quantile> <out.edges>      (Donatus .num, threshold = quantile of index)
#   Rscript pcc_fixed.R heinrichi <quantile> <out.edges>          (package dataset)
suppressMessages(library(stemmatology))
a <- commandArgs(TRUE)
if (a[1] == "heinrichi") {
  data(heinrichi); X <- heinrichi; rownames(X) <- as.character(seq_len(nrow(X)))
} else {
  L <- readLines(a[1]); wit <- strsplit(L[1], " +")[[1]]
  M <- do.call(rbind, lapply(L[-1][nzchar(L[-1])], function(l) strsplit(l, " +")[[1]]))
  M[M == "?"] <- NA; X <- matrix(as.integer(M) + 1L, nrow = nrow(M))
  colnames(X) <- wit; rownames(X) <- as.character(seq_len(nrow(X)))
}
cf <- PCC.conflicts(X)
ci <- cf$conflictsTotal[, 2]
thr <- as.numeric(quantile(ci[ci > 0], as.numeric(a[2])))
cat("locations", nrow(X), "conflicts", nrow(cf$edgelist), "threshold", thr, "\n")
db <- X; it <- 0; step3a <- 0
repeat {
  cf <- PCC.conflicts(db)
  if (sum(cf$conflictsTotal[, 1]) == 0) break
  it <- it + 1
  oc <- PCC.overconflicting(cf, ask = FALSE, threshold = thr)
  db2 <- PCC.elimination(oc)
  if (nrow(db2) < nrow(db)) { db <- db2; next }
  bad <- unique(as.vector(cf$edgelist))               # step 3a: drop all locations still in conflict
  db <- db[!(rownames(db) %in% bad), , drop = FALSE]; step3a <- length(bad); break
}
cat("iterations", it, "step3a_dropped", step3a, "kept_locations", nrow(db), "\n")
pdf(NULL)
st <- PCC.Stemma(db, ask = FALSE, verbose = FALSE)
el <- if (!is.null(st$edgelist)) st$edgelist else st[[1]]$edgelist
write.table(el, a[3], quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
cat("edges", nrow(el), "\n")
