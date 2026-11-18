source("analyze.R")

observations <- read_measure(c(0.05, NA, 2, 0.8, NA, 2, 0), c(0.1, 0.12, 0.1, 0.1, NA, 0.1, 0.1),
                             c(0.3, 0.3, 0.3, 0.3, 0.25, 0.3, 0.3),
                             c("", "ND", "L", "H", "S;ND", "Q", ""), factor = 0.001)
stopifnot(isTRUE(all.equal(observations$half, c(0.00005, 0.00006, 0.002, NA, 0.000125, NA, 0.00005))))
stopifnot(identical(observations$censored, c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE)))
stopifnot(isTRUE(all.equal(weighted_quantile(c(8, 1, 4), c(1, 9, 1), 0.50), 1)))
assertion <- try(assert(FALSE, "fixture"), silent = TRUE)
stopifnot(inherits(assertion, "try-error"))
pca <- weighted_pca(data.frame(a = c(1, 2, 3), b = c(2, 3, 4)), c(1, 3, 2))
stopifnot(abs(sum(pca$explained) - 1) < 1e-12, pca$explained[1] > 0.99)
cat("Censoring, quality flags, unit conversion, weighting, audits and PCA verified.\n")
