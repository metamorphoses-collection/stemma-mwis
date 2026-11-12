suppressMessages(library(phangorn))
args <- commandArgs(TRUE)
L <- readLines(args[1]); wit <- strsplit(L[1], " +")[[1]]
M <- do.call(rbind, lapply(L[-1][nzchar(L[-1])], function(l) strsplit(l, " +")[[1]]))
M[M == "?"] <- NA
X <- t(M); rownames(X) <- wit
lev <- sort(unique(as.vector(X[!is.na(X)])))
X[is.na(X)] <- "?"
pd <- phyDat(X, type = "USER", levels = lev, ambiguity = "?")
set.seed(1)
tr <- pratchet(pd, trace = 0, minit = 100, maxit = 1000, all = FALSE)
cat("parsimony score", parsimony(tr, pd), "\n")
write.tree(unroot(tr), file = args[2])
