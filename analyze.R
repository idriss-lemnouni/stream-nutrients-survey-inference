# Estimate nutrient concentrations using the released NRSA population weights.

assert <- function(ok, message) {
    if (!isTRUE(ok)) stop(message, call. = FALSE)
}

load_nrsa <- function(input) {
    read <- function(name) read.csv(file.path(input, paste0(name, ".csv")), fileEncoding = "UTF-8-BOM",
                                    colClasses = "character", na.strings = c("", "NA"), check.names = FALSE)
    population <- read("population")
    site <- read("site")
    chemistry <- read("chemistry")
    for (table in list(population, site, chemistry)) {
        assert(!anyNA(table$UID) && !anyDuplicated(table$UID), "Missing or duplicate UID")
    }
    selected <- population$VISIT_NO == "1" & population$EVAL_CAT == "Target_Sampled" &
        population$SITESAMP == "Y" & as.numeric(population$WGT_TP) > 0
    population <- population[which(selected), , drop = FALSE]
    site_order <- match(population$UID, site$UID)
    chemistry_order <- match(population$UID, chemistry$UID)
    assert(!anyNA(site_order) && !anyNA(chemistry_order), "Population UID failed to join")
    site <- site[site_order, , drop = FALSE]
    chemistry <- chemistry[chemistry_order, , drop = FALSE]
    assert(all(population$SITE_ID == site$SITE_ID & population$SITE_ID == chemistry$SITE_ID),
           "Station identifier disagrees across UID join")
    assert(!anyDuplicated(population$SITE_ID) && all(site$SITETYPE == "PROB"),
           "Population includes a duplicated or non-probability station")
    result <- data.frame(site = population$SITE_ID, weight = as.numeric(population$WGT_TP),
                         region = population$AG_ECO9, stratum = site$STRATUM,
                         lon = as.numeric(population$LON_DD83), lat = as.numeric(population$LAT_DD83))
    assert(all(is.finite(result$weight) & result$weight > 0) && !anyNA(result$region) && !anyNA(result$stratum),
           "Invalid design fields")
    assert(all(is.finite(result$lon) & result$lon >= -125 & result$lon <= -66 &
               is.finite(result$lat) & result$lat >= 24 & result$lat <= 50), "Unexpected site coordinates")
    list(design = result, chemistry = chemistry)
}

read_measure <- function(result, mdl, rl, flag, factor = 1) {
    below <- function(code) grepl(paste0("(^|;)", code, "($|;)"), flag)
    failed <- below("Q") | below("H")
    limit <- ifelse(is.finite(mdl) & mdl > 0, mdl, ifelse(is.finite(rl) & rl > 0, rl, NA_real_))
    censored <- !failed & (below("ND") |
        (is.finite(result) & is.finite(limit) & result >= 0 & result < limit))
    measured <- !failed & !censored & is.finite(result) & result >= 0
    substitute <- function(level) factor * ifelse(censored & is.finite(limit), level * limit,
                                                   ifelse(measured, result, NA_real_))
    list(zero = substitute(0), half = substitute(0.5), upper = substitute(1),
         censored = censored & is.finite(limit), failed = failed)
}

extract_measure <- function(chemistry, code, units, factor = 1) {
    col <- function(suffix) chemistry[[paste0(code, "_", suffix)]]
    number <- function(suffix) as.numeric(col(suffix))
    flag <- col("NARS_FLAG")
    flag[is.na(flag)] <- ""
    flag <- gsub(",", ";", flag, fixed = TRUE)
    qualifiers <- unique(unlist(strsplit(flag, ";", fixed = TRUE)))
    assert(all(qualifiers %in% c("ND", "L", "S", "Q", "H")), paste("Unmapped EPA flags for", code))
    assert(all(is.na(col("UNITS")) | col("UNITS") %in% units), paste("Unexpected units for", code))
    read_measure(number("RESULT"), number("MDL"), number("RL"), flag, factor)
}

weighted_quantile <- function(values, weights, probability) {
    sorted <- order(values)
    values[sorted][which(cumsum(weights[sorted]) >= probability * sum(weights))[1]]
}

nutrient_summary <- function(replicates, analyte, zone) {
    d <- replicates$variables
    valid <- is.finite(d[[paste0(analyte, "_half")]]) & (d$region == zone | zone == "National")
    in_zone <- d$region == zone | zone == "National"
    assert(sum(valid) > 2, paste("Too few valid observations in", zone))
    domain <- subset(replicates, valid)
    group <- d[valid, , drop = FALSE]
    formula <- reformulate(paste0(analyte, c("_zero", "_half", "_upper")))
    estimates <- survey::svymean(formula, domain)
    means <- as.numeric(coef(estimates))
    center <- group[[paste0(analyte, "_half")]]
    w <- group$weight
    se <- as.numeric(survey::SE(estimates))[2]
    margin <- qt(0.975, survey::degf(domain)) * se
    data.frame(analyte = analyte, zone = zone, stations = sum(valid),
               censored = sum(group[[paste0(analyte, "_censored")]]),
               kish_n = sum(w)^2 / sum(w^2), weight_coverage = sum(w) / sum(d$weight[in_zone]),
               unweighted = mean(center), lower_sub_mean = means[1], weighted = means[2], upper_sub_mean = means[3],
               median = weighted_quantile(center, w, 0.50), p90 = weighted_quantile(center, w, 0.90),
               se = se, cv = se / means[2], ci_low = means[2] - margin, ci_high = means[2] + margin)
}

weighted_pca <- function(values, weights) {
    result <- cov.wt(log1p(values), wt = weights, cor = TRUE)
    axes <- eigen(result$cor, symmetric = TRUE)
    rownames(axes$vectors) <- names(values)
    list(explained = axes$values / sum(axes$values), loadings = axes$vectors)
}

run <- function(input = "data/raw", out = "runs") {
    study <- load_nrsa(input)
    d <- study$design
    c <- study$chemistry
    measures <- list(nitrate = extract_measure(c, "NITRATE_N", "mg N/L"),
                     phosphorus = extract_measure(c, "PTL", "ug/L", 0.001),
                     nitrogen = extract_measure(c, "NTL", "mg N/L"),
                     carbon = extract_measure(c, "DOC", "mg/L"))
    for (analyte in names(measures)) {
        for (stat in names(measures[[analyte]])) {
            d[[paste(analyte, stat, sep = "_")]] <- measures[[analyte]][[stat]]
        }
    }
    design <- survey::svydesign(ids = ~1, weights = ~weight, strata = ~stratum, data = d)
    set.seed(1519)
    replicates <- survey::as.svrepdesign(design, type = "subbootstrap", replicates = 500, mse = TRUE)
    series <- lapply(c("nitrate", "phosphorus", "nitrogen"), function(name) {
        do.call(rbind, lapply(c("National", sort(unique(d$region))), function(zone) nutrient_summary(replicates, name, zone)))
    })
    estimates <- do.call(rbind, series)
    numeric <- d[paste0(names(measures), "_half")]
    pca_ok <- complete.cases(numeric)
    pca <- weighted_pca(numeric[pca_ok, ], d$weight[pca_ok])
    log_n <- update(design, log_phosphorus = log1p(phosphorus_half * 1000))
    regression <- survey::svyglm(log_phosphorus ~ log1p(nitrate_half) + region, design = log_n, family = gaussian())
    diagnostic <- summary(regression)$coefficients
    model <- data.frame(term = rownames(diagnostic), estimate = diagnostic[, 1], se = diagnostic[, 2],
                        statistic = diagnostic[, 3], p = diagnostic[, 4])
    audit <- data.frame(stations = nrow(d), design_strata = length(unique(d$stratum)), domain_count = length(unique(d$region)),
                        total_weight_miles = sum(d$weight), pca_sites = sum(pca_ok),
                        glm_sites = nobs(regression), bootstrap_replicates = 500)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    write.csv(estimates, file.path(out, "nutrient-estimates.csv"), row.names = FALSE)
    write.csv(data.frame(variable = rownames(pca$loadings), pca$loadings), file.path(out, "pca-loadings.csv"), row.names = FALSE)
    write.csv(data.frame(component = seq_along(pca$explained), explained = pca$explained),
              file.path(out, "pca-explained.csv"), row.names = FALSE)
    write.csv(model, file.path(out, "glm-estimates.csv"), row.names = FALSE)
    write.csv(audit, file.path(out, "design-audit.csv"), row.names = FALSE)
    print(audit)
    print(estimates[estimates$zone == "National", ])
}

if (sys.nframe() == 0) run()
