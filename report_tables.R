# Write Quarto include fragments from the audited R tables.
table_fragment <- function(df, path, caption) {
    cleaned <- lapply(df, function(column) {
        if (is.numeric(column)) {
            return(ifelse(is.na(column), "Unavailable", formatC(column, digits = 4, format = "fg")))
        }
        gsub("|", "\\|", column, fixed = TRUE)
    })
    row <- function(cells) paste0("| ", paste(trimws(cells), collapse = " | "), " |")
    content <- c(row(names(df)), row(rep("---", length(df))),
                 apply(as.data.frame(cleaned), 1, row), "", paste(":", caption), "")
    writeLines(content, path, useBytes = TRUE)
}

read_run <- function(name) read.csv(file.path("runs", paste0(name, ".csv")))
estimates <- read_run("nutrient-estimates")
table_fragment(subset(estimates, zone == "National")[c("analyte", "stations", "censored", "weight_coverage",
               "unweighted", "weighted", "ci_low", "ci_high", "median", "p90", "cv")],
               "runs/national.qmd",
               "National mean and 95% interval using half-limit substitution. Phosphorus is in mg/L. Nitrogen and nitrate are in mg N/L.")
table_fragment(subset(estimates, zone != "National")[c("zone", "analyte", "stations", "lower_sub_mean",
               "weighted", "upper_sub_mean", "weight_coverage", "kish_n")],
               "runs/ecoregions.qmd", "Detection-limit sensitivity and retained weight by reporting domain.")
table_fragment(read_run("pca-explained"), "runs/pca.qmd", "Weighted PCA variance explained.")
table_fragment(read_run("glm-estimates"), "runs/glm.qmd", "Survey-weighted phosphorus model.")
