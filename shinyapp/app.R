library(shiny)

cramers_v <- function(tab) {
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  n <- sum(tab)
  k <- min(dim(tab))
  v <- if (k > 1) sqrt(as.numeric(chi$statistic) / (n * (k - 1))) else NA_real_
  list(v = v, chi2 = as.numeric(chi$statistic), p_value = chi$p.value, n = n)
}

eta_squared <- function(df, group_col, value_col) {
  f <- as.formula(paste0("`", value_col, "` ~ `", group_col, "`"))
  m <- tryCatch(aov(f, data = df), error = function(e) NULL)
  if (is.null(m)) return(NA_real_)
  ss <- summary(m)[[1]][["Sum Sq"]]
  ss[1] / sum(ss)
}

# ---- Shared severity-tier vocabulary (matches the explainer's "how to read
# the output" section, applied consistently everywhere an association or
# design property gets flagged) ----
tier_label <- function(tier) {
  switch(tier,
    unfixable = "Unfixable (structural)",
    fixable = "Fixable at a cost",
    clean = "Clean",
    "Not checked")
}

tier_cls <- function(tier) {
  switch(tier, unfixable = "poor", fixable = "caution", clean = "good", "neutral")
}

# Red/yellow/green-style risk labels for the top-level Criteria assessment
# cards (a plain-language summary view) -- the detailed tables below keep
# the more precise Unfixable/Fixable/Clean vocabulary via tier_label().
tier_label_risk <- function(tier) {
  switch(tier,
    unfixable = "High risk",
    fixable = "Moderate risk",
    clean = "Low risk",
    "Not checked")
}

# Association-strength tier shared by every Cramer's V / eta-squared based
# check (batch confounding, covariate screening, interpretability,
# collinearity, missingness) so the whole app applies one consistent rule.
severity_from_v <- function(v) {
  if (is.na(v)) return("not_checked")
  if (v >= 0.6) return("unfixable")
  if (v >= 0.3) return("fixable")
  return("clean")
}

diagnose_batch <- function(df, condition_col, batch_col) {
  tab <- table(df[[condition_col]], df[[batch_col]])
  cv <- cramers_v(tab)
  list(batch_col = batch_col, v = cv$v, chi2 = cv$chi2, p_value = cv$p_value,
       severity = severity_from_v(cv$v))
}

diagnose_covariate <- function(df, condition_col, covariate_col) {
  is_num <- is.numeric(df[[covariate_col]])
  if (is_num) {
    assoc <- eta_squared(df, condition_col, covariate_col)
    method <- "eta_squared"
  } else {
    tab <- table(df[[condition_col]], df[[covariate_col]])
    assoc <- cramers_v(tab)$v
    method <- "cramers_v"
  }
  vals <- stats::na.omit(df[[covariate_col]])
  list(covariate = covariate_col, is_numeric = is_num, method = method,
       association = assoc, severity = severity_from_v(assoc),
       n_levels = length(unique(vals)))
}

# Full covariate screen: test EVERY column in the sheet (not just ones the
# user declared as batch/covariate) against condition, so unlabelled fields
# (extraction date, operator ID, etc.) that turn out to be the real batch
# effect don't get missed just because nobody thought to name them.
full_covariate_screen <- function(df, condition_col, sample_id_cols, roles) {
  other_cols <- setdiff(names(df), c(condition_col, sample_id_cols))
  lapply(other_cols, function(col) {
    d <- diagnose_covariate(df, condition_col, col)
    d$declared_role <- if (!is.null(roles[[col]])) roles[[col]] else "ignore"
    d$declared <- d$declared_role %in% c("batch", "covariate")
    d
  })
}

# Undeclared-but-associated columns from the full screen -- used both to
# register their dynamic plots (server) and to lay them out (renderUI), via
# this one shared definition so the two can't drift apart.
screen_undeclared_flagged <- function(full_screen) {
  Filter(function(s) !s$declared && !is.na(s$association) && s$association >= 0.3, full_screen)
}

# Collinearity among covariates: do any two covariates carry (nearly) the
# same information, so including both over-parameterizes the model even
# though neither is aliased with condition itself? |r| > 0.8 is the
# conventional threshold for continuous-continuous redundancy; categorical/
# mixed pairs fall back to the app's general association-tier thresholds,
# since there's no single agreed-upon cutoff for those statistics here.
severity_from_r_collinearity <- function(r) {
  if (is.na(r)) return("not_checked")
  if (r >= 0.8) return("fixable")
  return("clean")
}

covariate_collinearity <- function(df, covariate_cols) {
  if (length(covariate_cols) < 2) return(list())
  pairs <- utils::combn(covariate_cols, 2, simplify = FALSE)
  lapply(pairs, function(p) {
    a <- p[1]; b <- p[2]
    a_num <- is.numeric(df[[a]]); b_num <- is.numeric(df[[b]])
    if (a_num && b_num) {
      val <- suppressWarnings(as.numeric(stats::cor(df[[a]], df[[b]], use = "pairwise.complete.obs")))
      val <- abs(val); method <- "pearson_r"
      severity <- severity_from_r_collinearity(val)
    } else if (!a_num && !b_num) {
      tab <- table(df[[a]], df[[b]])
      val <- cramers_v(tab)$v; method <- "cramers_v"
      severity <- severity_from_v(val)
    } else {
      grp <- if (a_num) b else a
      value_col <- if (a_num) a else b
      val <- eta_squared(df, grp, value_col); method <- "eta_squared"
      severity <- severity_from_v(val)
    }
    list(a = a, b = b, method = method, value = val, severity = severity)
  })
}

# Missingness vs. group: is data missing preferentially within one condition
# group? If so, imputation can't fix the resulting bias -- the missingness
# pattern itself is informative. Only columns that actually have some
# missing data are reported.
missingness_screen <- function(df, condition_col) {
  cols <- setdiff(names(df), condition_col)
  results <- lapply(cols, function(col) {
    v <- df[[col]]
    is_missing <- is.na(v) | (is.character(v) & trimws(v) == "")
    n_missing <- sum(is_missing)
    n <- length(v)
    if (n_missing == 0 || n_missing == n) {
      return(list(column = col, n_missing = n_missing, pct_missing = round(100 * n_missing / n, 1),
                   v = NA_real_, p_value = NA_real_, severity = "clean"))
    }
    tab <- table(is_missing, df[[condition_col]])
    if (nrow(tab) < 2) {
      return(list(column = col, n_missing = n_missing, pct_missing = round(100 * n_missing / n, 1),
                   v = NA_real_, p_value = NA_real_, severity = "clean"))
    }
    cv <- cramers_v(tab)
    list(column = col, n_missing = n_missing, pct_missing = round(100 * n_missing / n, 1),
         v = cv$v, p_value = cv$p_value, severity = severity_from_v(cv$v))
  })
  results[vapply(results, function(r) r$n_missing > 0, logical(1))]
}

# ============================================================================
# Plotting helpers -- base R graphics only (barplot/boxplot/plot), so these
# run under webR/shinylive with no extra package dependencies. Each function
# renders one plot for one check result; wired up as dynamic renderPlot()
# outputs in the server, one per batch/covariate/pair/column found.
# ============================================================================

PALETTE <- c("#1A345E", "#E64626", "#BDDC96", "#FFB800", "#91BDE5", "#ded4b9")

# Shared plot-output-id builder -- used identically when registering a
# renderPlot() in the server and when creating the matching plotOutput() in
# the UI, so the two never drift apart.
plot_id <- function(prefix, ...) {
  parts <- vapply(list(...), make.names, character(1))
  paste0(prefix, "_", paste(parts, collapse = "_"))
}

plot_condition_by_batch <- function(df, condition_col, batch_col) {
  tab <- table(df[[batch_col]], df[[condition_col]])
  m <- t(tab)
  cols <- PALETTE[((seq_len(nrow(m)) - 1) %% length(PALETTE)) + 1]
  op <- graphics::par(mar = c(4, 4, 3, 1)); on.exit(graphics::par(op))
  graphics::barplot(m, beside = TRUE, col = cols, ylab = "samples", las = 1,
                     main = paste0(condition_col, " within each ", batch_col),
                     legend.text = rownames(m),
                     args.legend = list(x = "topright", bty = "n", cex = 0.8, title = condition_col))
}

plot_covariate_by_condition <- function(df, condition_col, covariate_col) {
  if (is.numeric(df[[covariate_col]])) {
    op <- graphics::par(mfrow = c(2, 1), mar = c(2, 4, 3, 1)); on.exit(graphics::par(op))
    f <- stats::as.formula(paste0("`", covariate_col, "` ~ `", condition_col, "`"))
    graphics::boxplot(f, data = df, col = "#91BDE5", ylab = covariate_col, xlab = "",
                       main = paste0(covariate_col, " by ", condition_col), las = 1)
    graphics::par(mar = c(4, 4, 1, 1))
    graphics::stripchart(f, data = df, vertical = TRUE, method = "jitter", jitter = 0.15,
                          pch = 19, col = "#1A345E", ylab = covariate_col, xlab = condition_col, las = 1)
  } else {
    op <- graphics::par(mar = c(4, 4, 3, 1)); on.exit(graphics::par(op))
    tab <- table(df[[condition_col]], df[[covariate_col]])
    m <- t(tab)
    cols <- PALETTE[((seq_len(nrow(m)) - 1) %% length(PALETTE)) + 1]
    graphics::barplot(m, beside = TRUE, col = cols, ylab = "samples", las = 1,
                       main = paste0(covariate_col, " by ", condition_col),
                       legend.text = rownames(m), args.legend = list(x = "topright", bty = "n", cex = 0.8))
  }
}

plot_vif_bar <- function(vif_result) {
  if (length(vif_result) == 0) return(invisible(NULL))
  terms_v <- vapply(vif_result, function(v) v$term, character(1))
  vals <- vapply(vif_result, function(v) if (is.finite(v$gvif_std)) v$gvif_std else NA_real_, numeric(1))
  op <- graphics::par(mar = c(6, 4, 3, 1)); on.exit(graphics::par(op))
  finite_vals <- vals[is.finite(vals)]
  ymax <- if (length(finite_vals)) max(c(finite_vals, sqrt(10))) * 1.15 else sqrt(10) * 1.15
  bar_vals <- ifelse(is.finite(vals), vals, 0)
  bp <- graphics::barplot(bar_vals, names.arg = terms_v, col = "#1A345E", ylim = c(0, ymax),
                           ylab = "GVIF (df-adjusted)", las = 2, main = "VIF per term")
  graphics::abline(h = sqrt(5), col = "#FFB800", lty = 2)
  graphics::abline(h = sqrt(10), col = "#E64626", lty = 2)
  na_idx <- which(!is.finite(vals))
  if (length(na_idx)) graphics::text(bp[na_idx], ymax * 0.5, "Inf", col = "#E64626", font = 2)
}

plot_collinearity_pair <- function(df, a, b) {
  op <- graphics::par(mar = c(4, 4, 3, 1)); on.exit(graphics::par(op))
  a_num <- is.numeric(df[[a]]); b_num <- is.numeric(df[[b]])
  if (a_num && b_num) {
    graphics::plot(df[[a]], df[[b]], pch = 19, col = "#1A345E",
                   xlab = a, ylab = b, main = paste(a, "vs", b), las = 1)
    ok <- stats::complete.cases(df[[a]], df[[b]])
    if (sum(ok) >= 2) {
      fit <- tryCatch(stats::lm(df[[b]][ok] ~ df[[a]][ok]), error = function(e) NULL)
      if (!is.null(fit)) graphics::abline(fit, col = "#E64626", lwd = 2)
    }
  } else if (a_num != b_num) {
    num_col <- if (a_num) a else b
    cat_col <- if (a_num) b else a
    f <- stats::as.formula(paste0("`", num_col, "` ~ `", cat_col, "`"))
    graphics::boxplot(f, data = df, col = "#91BDE5", xlab = cat_col, ylab = num_col,
                       main = paste(a, "vs", b), las = 1)
  } else {
    tab <- table(df[[a]], df[[b]])
    cols <- PALETTE[((seq_len(nrow(tab)) - 1) %% length(PALETTE)) + 1]
    graphics::barplot(tab, beside = TRUE, col = cols, ylab = "samples", las = 1,
                       main = paste(a, "vs", b),
                       legend.text = rownames(tab), args.legend = list(x = "topright", bty = "n", cex = 0.8))
  }
}

plot_missingness_by_group <- function(df, condition_col, column) {
  op <- graphics::par(mar = c(4, 4, 3, 1)); on.exit(graphics::par(op))
  v <- df[[column]]
  is_missing <- is.na(v) | (is.character(v) & trimws(v) == "")
  pct <- tapply(is_missing, df[[condition_col]], function(x) 100 * mean(x))
  graphics::barplot(pct, col = "#E64626", ylab = "% missing", las = 1, ylim = c(0, 100),
                     main = paste0("Missingness of ", column, " by ", condition_col))
}

# PCA on the batch + covariate design (one-hot encoded, zero-variance columns
# dropped), colored by condition -- an exploratory visual complement to the
# rank/VIF/association checks, since the app has no molecular data to run a
# "real" PCA on. Base R only (prcomp).
compute_pca <- function(df, batch_cols, covariate_cols) {
  design_cols <- unique(c(batch_cols, covariate_cols))
  if (length(design_cols) == 0) return(NULL)
  f <- stats::as.formula(paste("~ 0 +", paste(sprintf("`%s`", design_cols), collapse = " + ")))
  mm <- tryCatch(stats::model.matrix(f, data = df), error = function(e) NULL)
  if (is.null(mm) || nrow(mm) < 3) return(NULL)
  vars <- apply(mm, 2, stats::var)
  mm <- mm[, vars > 1e-8, drop = FALSE]
  if (ncol(mm) < 2) return(NULL)
  pca <- tryCatch(stats::prcomp(mm, center = TRUE, scale. = TRUE), error = function(e) NULL)
  if (is.null(pca)) return(NULL)
  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
  list(scores = pca$x, var_explained = var_explained, n = nrow(mm))
}

plot_pca <- function(pca_result, condition_vec) {
  if (is.null(pca_result) || ncol(pca_result$scores) < 2) return(invisible(NULL))
  op <- graphics::par(mar = c(4, 4, 3, 8), xpd = TRUE); on.exit(graphics::par(op))
  scores <- pca_result$scores
  cond_f <- factor(condition_vec)
  cols <- PALETTE[as.integer(cond_f)]
  ve <- pca_result$var_explained
  graphics::plot(scores[, 1], scores[, 2], col = cols, pch = 19, cex = 1.3,
                 xlab = sprintf("PC1 (%.0f%% var)", ve[1] * 100),
                 ylab = sprintf("PC2 (%.0f%% var)", ve[2] * 100),
                 main = "PCA of batch/covariate design", las = 1)
  graphics::legend("topright", inset = c(-0.3, 0), legend = levels(cond_f),
                    col = PALETTE[seq_along(levels(cond_f))], pch = 19, bty = "n",
                    title = "condition", cex = 0.85, xpd = TRUE)
}

replicate_adequacy <- function(df, condition_col) {
  counts <- table(df[[condition_col]])
  data.frame(condition = names(counts), n = as.integer(counts), stringsAsFactors = FALSE)
}

batch_balance <- function(df, batch_col) {
  counts <- table(df[[batch_col]])
  cv <- if (length(counts) > 1) stats::sd(counts) / mean(counts) else NA_real_
  list(batch_col = batch_col, cv = cv)
}

# Rank / aliasing check: is condition separable from batch + covariates, and
# are batch/covariates separable from each other? Uses rank-contribution
# (does adding a term's columns to a model raise the rank by its full column
# count, or does some of that information turn out to be redundant?) rather
# than raw QR pivot order, which is arbitrary w.r.t. formula term order and
# would otherwise "blame" whichever term happens to come later. Base R only
# (model.matrix + qr), so this runs fine under webR too.
term_rank_contrib <- function(df, focus_cols, other_cols) {
  f_full <- as.formula(paste("~", paste(sprintf("`%s`", c(focus_cols, other_cols)), collapse = " + ")))
  mm_full <- tryCatch(stats::model.matrix(f_full, data = df), error = function(e) NULL)
  if (is.null(mm_full)) return(NULL)
  rank_full <- qr(mm_full)$rank

  if (length(other_cols) == 0) {
    rank_without <- 1L; ncol_without <- 1L
  } else {
    f_without <- as.formula(paste("~", paste(sprintf("`%s`", other_cols), collapse = " + ")))
    mm_without <- tryCatch(stats::model.matrix(f_without, data = df), error = function(e) NULL)
    if (is.null(mm_without)) return(NULL)
    rank_without <- qr(mm_without)$rank
    ncol_without <- ncol(mm_without)
  }
  focus_ncol <- ncol(mm_full) - ncol_without
  gained <- rank_full - rank_without
  list(deficient = gained < focus_ncol, focus_ncol = focus_ncol, gained = gained,
       rank_full = rank_full, ncol_full = ncol(mm_full))
}

check_design_rank <- function(df, condition_col, batch_cols, covariate_cols) {
  other_terms <- unique(c(batch_cols, covariate_cols))
  if (length(other_terms) == 0) {
    return(list(status = "single_term"))
  }
  overall <- tryCatch(term_rank_contrib(df, condition_col, other_terms), error = function(e) NULL)
  if (is.null(overall)) {
    return(list(status = "error"))
  }

  if (isTRUE(overall$deficient)) {
    culprits <- character(0)
    for (t in other_terms) {
      res_t <- tryCatch(term_rank_contrib(df, condition_col, t), error = function(e) NULL)
      if (!is.null(res_t) && isTRUE(res_t$deficient)) culprits <- c(culprits, t)
    }
    if (length(culprits) == 0) culprits <- other_terms
    return(list(status = "condition_aliased", aliased_terms = culprits,
                rank = overall$rank_full, ncol = overall$ncol_full))
  }

  if (length(other_terms) < 2) {
    return(list(status = "full_rank", rank = overall$rank_full, ncol = overall$ncol_full))
  }
  f_other <- as.formula(paste("~", paste(sprintf("`%s`", other_terms), collapse = " + ")))
  mm_other <- tryCatch(stats::model.matrix(f_other, data = df), error = function(e) NULL)
  if (is.null(mm_other)) return(list(status = "full_rank", rank = overall$rank_full, ncol = overall$ncol_full))
  qr_obj <- qr(mm_other)
  rank <- qr_obj$rank
  p <- ncol(mm_other)
  if (rank == p) {
    return(list(status = "full_rank", rank = overall$rank_full, ncol = overall$ncol_full))
  }
  dependent_pos <- qr_obj$pivot[(rank + 1):p]
  term_labels <- attr(stats::terms(f_other), "term.labels")
  assign_idx <- attr(mm_other, "assign")
  dependent_term_idx <- sort(unique(assign_idx[dependent_pos]))
  dependent_term_idx <- dependent_term_idx[dependent_term_idx > 0]
  aliased_terms <- gsub("`", "", term_labels[dependent_term_idx])
  list(status = "rank_deficient", rank = rank, ncol = p, aliased_terms = aliased_terms)
}

assess_rank <- function(rank_check, condition_col) {
  if (is.null(rank_check) || rank_check$status %in% c("single_term", "error")) {
    return(list(cls = "neutral", tier = "not_checked",
      reasons = "Mark at least a condition column plus one batch or covariate column to run this check."))
  }
  if (rank_check$status == "condition_aliased") {
    other <- rank_check$aliased_terms
    return(list(cls = "poor", tier = "unfixable",
      reasons = sprintf(
        "%s is perfectly (or near-perfectly) aliased with %s — the design matrix is rank-deficient (rank %d of %d columns). No statistical adjustment can separate a %s effect from %s in this data; the only fix is redesigning sample allocation so these terms vary independently.",
        condition_col, paste(other, collapse = ", "), rank_check$rank, rank_check$ncol, condition_col, paste(other, collapse = " / "))))
  }
  if (rank_check$status == "full_rank") {
    return(list(cls = "good", tier = "clean",
      reasons = sprintf(
        "No term in this design (condition, batch, covariates) is a perfect linear combination of the others — the design matrix has full rank (%d of %d columns). Each term's effect is in principle separable from the rest.",
        rank_check$rank, rank_check$ncol)))
  }
  list(cls = "caution", tier = "fixable",
    reasons = sprintf(
      "%s are perfectly redundant given the rest of the design (rank %d of %d columns among batch/covariate terms). %s itself is not affected, but drop one of these terms or the model will be unidentifiable for them.",
      paste(rank_check$aliased_terms, collapse = " & "), rank_check$rank, rank_check$ncol, condition_col))
}

# Variance Inflation Factor per term (generalized VIF, Fox & Monette 1992),
# computed from the correlation matrix of the design matrix -- base R only
# (cor + det), so this is safe under webR. GVIF^(1/(2*df)) is the
# df-standardized version, comparable across terms with different numbers
# of dummy columns. VIF is undefined (Inf) exactly when a term is aliased,
# which the rank/aliasing check above already reports -- so an Inf/NA VIF
# here just cross-references that check rather than showing a bogus number.
compute_vif <- function(df, condition_col, batch_cols, covariate_cols) {
  terms_in <- unique(c(condition_col, batch_cols, covariate_cols))
  if (length(terms_in) < 2) return(list())
  f <- as.formula(paste("~", paste(sprintf("`%s`", terms_in), collapse = " + ")))
  mm_full <- tryCatch(stats::model.matrix(f, data = df), error = function(e) NULL)
  if (is.null(mm_full)) return(list())
  assign_full <- attr(mm_full, "assign")
  keep <- assign_full != 0
  mm <- mm_full[, keep, drop = FALSE]
  assign_idx <- assign_full[keep]
  term_labels <- gsub("`", "", attr(stats::terms(f), "term.labels"))
  p <- ncol(mm)
  if (p < 2) return(list())
  R <- tryCatch(stats::cor(mm), error = function(e) NULL)
  if (is.null(R) || any(!is.finite(R))) return(list())
  detR <- det(R)

  results <- list()
  for (i in seq_along(term_labels)) {
    cols_j <- which(assign_idx == i)
    cols_rest <- which(assign_idx != i)
    dfj <- length(cols_j)
    if (dfj == 0) next
    if (length(cols_rest) == 0) {
      gvif <- 1
    } else if (!is.finite(detR) || abs(detR) < 1e-10) {
      gvif <- Inf
    } else {
      gvif <- (det(R[cols_j, cols_j, drop = FALSE]) * det(R[cols_rest, cols_rest, drop = FALSE])) / detR
    }
    gvif_std <- if (is.finite(gvif) && gvif > 0) gvif^(1 / (2 * dfj)) else Inf
    n_effective <- if (is.finite(gvif) && gvif >= 1) nrow(df) / gvif else NA_real_
    results[[term_labels[i]]] <- list(term = term_labels[i], df = dfj,
                                       gvif = gvif, gvif_std = gvif_std, n_effective = n_effective)
  }
  results
}

assess_vif <- function(vif_entry) {
  if (is.null(vif_entry) || !is.finite(vif_entry$gvif_std)) {
    return(list(tier = "unfixable", cls = "poor",
      reason = "VIF is undefined for this term because it is (near-)perfectly aliased with another term — see the Rank / aliasing check above."))
  }
  gs <- vif_entry$gvif_std
  if (gs >= sqrt(10)) {
    return(list(tier = "fixable", cls = "poor",
      reason = sprintf("GVIF %.2f (df-standardized %.2f) — severe: equivalent to running this comparison with an effective sample size of about %.0f instead of the %s actually collected.",
                        vif_entry$gvif, gs, vif_entry$n_effective, "full sample")))
  }
  if (gs >= sqrt(5)) {
    return(list(tier = "fixable", cls = "caution",
      reason = sprintf("GVIF %.2f (df-standardized %.2f) — moderate: this term shares meaningful variance with the rest of the design, costing some precision if included as a covariate.",
                        vif_entry$gvif, gs)))
  }
  list(tier = "clean", cls = "good", reason = "Not meaningfully predictable from the rest of the design.")
}

guess_role <- function(col) {
  c <- tolower(col)
  if (grepl("id", c)) return("sample_id")
  if (grepl("condition|group|status|disease|treatment", c)) return("condition")
  if (grepl("batch|site|day|run|lane|operator|kit|plate", c)) return("batch")
  if (grepl("sex|age|cohort|centre|center", c)) return("covariate")
  return("ignore")
}

# ---- Power calculator functions (verbatim from the tested webR version) ----

count_effect_sd <- function(depth, cv) sqrt(2 * (1 / depth + cv^2) / (log(2)^2))

count_power <- function(n, depth, cv, log2fc, alpha, target_power) {
  sd <- count_effect_sd(depth, cv)
  pt <- power.t.test(n = n, delta = abs(log2fc), sd = sd, sig.level = alpha)
  n_needed <- tryCatch(
    power.t.test(delta = abs(log2fc), sd = sd, sig.level = alpha, power = target_power)$n,
    error = function(e) NA_real_)
  list(achieved_power = as.numeric(pt$power), sd_used = sd, n_needed = n_needed)
}

cohend_power <- function(n, cohend, alpha, target_power) {
  pt <- power.t.test(n = n, delta = abs(cohend), sd = 1, sig.level = alpha)
  n_needed <- tryCatch(
    power.t.test(delta = abs(cohend), sd = 1, sig.level = alpha, power = target_power)$n,
    error = function(e) NA_real_)
  list(achieved_power = as.numeric(pt$power), n_needed = n_needed)
}

capture_probability <- function(n_units, freq, k_min) 1 - pbinom(k_min - 1, size = n_units, prob = freq)

capture_n_needed <- function(freq, k_min, target_prob, n_max = 5e6) {
  n <- k_min
  while (n <= n_max) {
    if (capture_probability(n, freq, k_min) >= target_prob) return(n)
    n <- ceiling(n * 1.08) + 1
  }
  NA_real_
}

depth_n_needed <- function(af, k_min, target_prob, n_max = 2000) capture_n_needed(af, k_min, target_prob, n_max)

assoc_power <- function(n_cases, n_controls, p_cases, p_controls, alpha, target_power) {
  n_use <- min(n_cases, n_controls)
  pt <- power.prop.test(n = n_use, p1 = p_cases, p2 = p_controls, sig.level = alpha)
  n_needed <- tryCatch(
    power.prop.test(p1 = p_cases, p2 = p_controls, sig.level = alpha, power = target_power)$n,
    error = function(e) NA_real_)
  list(achieved_power = as.numeric(pt$power), n_per_group_used = n_use, n_needed = n_needed)
}

# ---- Four-criteria assessment (ported from the JS version) ----

assess_resource <- function(worst_cv, rep_df) {
  sizes <- rep_df$n
  max_ratio <- if (min(sizes) > 0) max(sizes) / min(sizes) else Inf
  reasons <- c()
  if (!is.na(worst_cv) && (worst_cv > 0.4 || max_ratio > 3)) {
    cls <- "poor"; tier <- "fixable"
    if (!is.na(worst_cv) && worst_cv > 0.4) reasons <- c(reasons, "Samples are very unevenly spread across batches — some batches carry little unique information.")
    if (max_ratio > 3) reasons <- c(reasons, "Condition group sizes are highly imbalanced, which wastes power relative to the total sample collected.")
  } else if (!is.na(worst_cv) && (worst_cv > 0.15 || max_ratio > 1.5)) {
    cls <- "caution"; tier <- "fixable"
    reasons <- "Some imbalance across batches or condition groups — not fatal, but resources aren't optimally spent."
  } else {
    cls <- "good"; tier <- "clean"
    reasons <- "Samples are evenly spread across batches and condition groups, so every sample is pulling its weight."
  }
  list(cls = cls, tier = tier, reasons = reasons)
}

assess_interpretability <- function(worst_batch_v, worst_covariate_v) {
  worst <- suppressWarnings(max(c(worst_batch_v, worst_covariate_v), na.rm = TRUE))
  if (!is.finite(worst)) worst <- 0
  if (worst >= 0.5) {
    cls <- "poor"; tier <- "fixable"
    reasons <- "Batch and/or a covariate is strongly entangled with condition — an observed effect could belong to any of them. See the Rank / aliasing check above for whether this is exact (unfixable) or partial (adjustable at a cost)."
  } else if (worst >= 0.25) {
    cls <- "caution"; tier <- "fixable"
    reasons <- "Some entanglement present — findings should be reported with an explicit caveat about what can and can't be disentangled."
  } else {
    cls <- "good"; tier <- "clean"
    reasons <- "Condition is reasonably independent of batch and known covariates, so an effect can be attributed with more confidence."
  }
  list(cls = cls, tier = tier, reasons = reasons)
}

assess_generalisability <- function(covariates) {
  if (length(covariates) == 0) {
    return(list(cls = "neutral", tier = "not_checked",
      reasons = "No covariates marked — add columns like sex, age, or site to assess how broadly this cohort is represented."))
  }
  narrow_count <- 0
  reasons <- c()
  for (cov in covariates) {
    if (cov$n_levels <= 1) {
      narrow_count <- narrow_count + 1
      reasons <- c(reasons, sprintf('"%s" has only one value represented — findings may not extend beyond this %s.', cov$covariate, cov$covariate))
    }
  }
  if (narrow_count == 0) {
    cls <- "good"; tier <- "clean"
    reasons <- c(reasons, "The covariates provided show some spread.")
  } else if (narrow_count < length(covariates)) { cls <- "caution"; tier <- "fixable" }
  else { cls <- "poor"; tier <- "fixable" }
  reasons <- c(reasons, "Remember: this can only speak to the covariates you've included, not to populations absent from the metadata entirely.")
  list(cls = cls, tier = tier, reasons = reasons)
}

# ============================================================================
# Example datasets 
# ============================================================================

examples <- list(
  siteConfound = "sample_id,condition,site,sex,age
S01,disease,SiteX,F,58
S02,disease,SiteX,M,62
S03,disease,SiteX,F,49
S04,disease,SiteX,M,55
S05,disease,SiteX,F,67
S06,disease,SiteX,M,71
S07,control,SiteY,F,52
S08,control,SiteY,M,60
S09,control,SiteY,F,45
S10,control,SiteY,M,58
S11,control,SiteY,F,63
S12,control,SiteY,M,50",
  batchConfound = "sample_id,condition,batch,sex,age
S01,early,Batch1,F,34
S02,early,Batch1,M,29
S03,early,Batch1,F,41
S04,early,Batch1,M,37
S05,early,Batch1,F,45
S06,early,Batch1,M,31
S07,late,Batch2,F,38
S08,late,Batch2,M,44
S09,late,Batch2,F,29
S10,late,Batch2,M,50
S11,late,Batch2,F,33
S12,late,Batch2,M,40",
  balanced = "sample_id,condition,batch,sex,age
S01,disease,Batch1,F,55
S02,disease,Batch2,M,61
S03,disease,Batch3,F,49
S04,disease,Batch1,M,58
S05,disease,Batch2,F,63
S06,disease,Batch3,M,52
S07,control,Batch1,M,57
S08,control,Batch2,F,60
S09,control,Batch3,M,50
S10,control,Batch1,F,64
S11,control,Batch2,M,53
S12,control,Batch3,F,59"
)

ASSAY_LABELS <- c(
  bulk_rna = "Bulk RNA-seq", atac = "ATAC-seq", scrna = "scRNA-seq",
  prot_dia = "Proteomics (DIA)", prot_tmt = "Proteomics (TMT)",
  metabolomics = "Metabolomics", wgs = "WGS / WES", other = "Other"
)

# ============================================================================
# Explainer content (Design diagnostics tab)
# ============================================================================

explainer_intro <- "Input: a sample metadata sheet with group assignment, batch/technical variables, and covariates. Every check below runs on metadata alone, before any molecular data exists — it tells you whether your design will let you draw valid conclusions, while the design is still fixable. Hover the ⓘ next to any check for its test, statistic, and how to read it."

# Per-check explanations, surfaced as hover tooltips directly on the result
# that they explain (rather than as a separate wall of text up top).
info_tip <- function(text) {
  tags$span(class = "info-tip", HTML("&#9432;"),
    tags$span(class = "tip-bubble", text))
}

tip_rank <- "Test: builds the design matrix X implied by condition + batch + covariates and computes rank(X) via QR/SVD. Statistic: rank deficiency d = ncol(X) − rank(X). Interpretation: d = 0 means every term is estimable given the others; d > 0 means at least one term is an exact linear combination of another (e.g. batch = A for every case sample). Necessary but not sufficient — full rank means the model is estimable, not that estimates are precise (that's what VIF is for)."
tip_vif <- "Test: regresses each term on all other terms in the design and takes VIF = 1/(1 − R²) (generalized VIF for multi-level factors, Fox & Monette 1992, so it's comparable across terms with different degrees of freedom). Interpretation: VIF = 1 means orthogonal to the rest of the design; VIF > 5 warrants attention, VIF > 10 is severe. Equivalent to estimating that term's effect with an effective sample size of n / VIF."
tip_batch_assoc <- "Test: cross-tabulates condition against each batch/technical variable (chi-square test of independence, or Cramer's V as effect size). Interpretation: a zero cell (a condition/batch combination with no samples) is a structural, unfixable confound — there's no data to estimate the batch effect independently of condition in that stratum. A singleton cell (n = 1) has no within-cell variance to separate batch noise from a genuine data point. A significant result with all cells reasonably populated is the fixable, partial-confounding case."
tip_balance <- "Test: how evenly is condition distributed within each batch, versus what proportional allocation would predict (entropy or chi-square goodness-of-fit). Interpretation: even spread means batch behaves like a random block; condition-homogeneous batches are the imbalance pattern that produces the confounding/VIF problems above — this is a fast per-batch scan before drilling into the full contingency table."
tip_full_screen <- "Test: every column in the sheet — not just ones marked Batch or Covariate — tested against condition (chi-square/Fisher for categorical columns, ANOVA/eta-squared for continuous ones). Interpretation: a flagged column is a candidate confounder whether or not anyone intended to model it — unlabelled fields like processing date or operator ID surface here even when nobody thought to flag them as batch."
tip_covariate_assoc <- "Test: association between each declared covariate and condition (Cramer's V for categorical, eta-squared from ANOVA for continuous). Interpretation: a high association means that covariate is a candidate confounder and should be included as an adjustment in any model, at some cost to power."
tip_collinearity <- "Test: pairwise correlation among covariates (Pearson r for continuous pairs, Cramer's V / eta-squared for categorical or mixed pairs). Interpretation: |r| > 0.8 between two continuous covariates means they carry substantially redundant information — including both inflates variance for both without adding independent explanatory power. This is a special case of the VIF problem, restricted to covariate-covariate pairs rather than covariate-vs-condition."
tip_missingness <- "Test: chi-square test on a 2×k table of (missing/present) × condition, for each column with missing data. Interpretation: a significant association means data isn't missing completely at random with respect to condition — standard imputation does not correct for group-dependent missingness, it will reproduce or amplify the bias rather than remove it."

# ============================================================================
# UI
# ============================================================================

ui <- fluidPage(
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),
  div(class = "top-navbar",
      span(class = "navbar-brand", "Omics experiment calculator")
  ),
  div(class = "app-container",
    tabsetPanel(
      id = "mainTabs",

      tabPanel("Design diagnostics",
        br(),
        div(class = "card-panel",
          div(class = "hint", explainer_intro)
        ),
        div(class = "card-row",
          div(class = "card-panel",
            h4("1. Load metadata"),
            div(
              actionButton("loadSite", "Load: site confound", class = "btn-default example-btn"),
              actionButton("loadBatch", "Load: batch confound", class = "btn-default example-btn"),
              actionButton("loadBalanced", "Load: balanced design", class = "btn-default example-btn")
            ),
            textAreaInput("csvText", NULL, rows = 10, width = "100%",
                          placeholder = "sample_id,condition,batch,sex,age\nS01,disease,batch1,F,54\n..."),
            fileInput("csvFile", "...or upload a CSV file", accept = ".csv"),
            actionButton("parseBtn", "Parse data →", class = "btn-primary")
          ),
          div(class = "card-panel",
            h4("Preview"),
            tableOutput("previewTable")
          )
        ),
        div(class = "card-panel",
          uiOutput("mappingUI"),
          uiOutput("diagnoseButtonUI")
        ),
        uiOutput("diagnosisResults")
      ),

      tabPanel("Power calculator",
        br(),
        div(class = "card-panel",
          div(class = "hint", "No fixed replicate-number thresholds — this runs a real closed-form power calculation (power.t.test, power.prop.test, or exact binomial coverage) against a target power ", tags$b("you"), " set."),
          selectInput("pcAssay", "Sequencing / assay type",
                      choices = c("Bulk RNA-seq" = "bulk_rna", "ATAC-seq" = "atac", "scRNA-seq" = "scrna",
                                  "Proteomics — DIA" = "prot_dia", "Proteomics — TMT" = "prot_tmt",
                                  "Metabolomics" = "metabolomics", "WGS / WES cohort" = "wgs", "Other / not sure" = "other")),
          uiOutput("pcModeUI"),
          uiOutput("pcFormUI"),
          br(),
          actionButton("pcRunBtn", "Check power →", class = "btn-primary"),
          uiOutput("pcResults")
        )
      )
    )
  ),
  tags$footer(class = "app-footer",
    span("© 2026 Sydney Informatics Hub, University of Sydney"),
    a(href = "https://github.com/Sydney-Informatics-Hub/omics-calculator", target = "_blank", "GitHub")
  )
)

# ============================================================================
# Server
# ============================================================================

server <- function(input, output, session) {

  rv <- reactiveValues(df = NULL, col_ids = NULL, diag_result = NULL, pc_result = NULL, pc_mode = list())

  # ---- Load examples ----
  observeEvent(input$loadSite, updateTextAreaInput(session, "csvText", value = examples$siteConfound))
  observeEvent(input$loadBatch, updateTextAreaInput(session, "csvText", value = examples$batchConfound))
  observeEvent(input$loadBalanced, updateTextAreaInput(session, "csvText", value = examples$balanced))

  observeEvent(input$csvFile, {
    txt <- paste(readLines(input$csvFile$datapath, warn = FALSE), collapse = "\n")
    updateTextAreaInput(session, "csvText", value = txt)
  })

  # ---- Parse ----
  observeEvent(input$parseBtn, {
    txt <- input$csvText
    validate(need(nchar(trimws(txt)) > 0, "Paste or load some CSV data first."))
    df <- tryCatch(read.csv(text = txt, stringsAsFactors = FALSE, check.names = FALSE),
                    error = function(e) NULL)
    validate(need(!is.null(df), "Could not parse this as CSV."))
    rv$df <- df
    rv$col_ids <- setNames(make.names(names(df), unique = TRUE), names(df))
    rv$diag_result <- NULL
    rv$pc_result <- NULL
  })

  output$previewTable <- renderTable({
    req(rv$df)
    head(rv$df, 8)
  })

  # ---- Column role mapping UI ----
  output$mappingUI <- renderUI({
    req(rv$df)
    cols <- names(rv$df)
    tagList(
      h4("2. Map columns"),
      div(class = "hint", "Mark exactly one column as Condition (the biological variable of interest); mark technical variables as Batch and demographic/clinical variables as Covariate."),
      lapply(cols, function(col) {
        id <- paste0("role_", rv$col_ids[[col]])
        selectInput(id, label = col, width = "100%",
                    choices = c("Ignore" = "ignore", "Sample ID" = "sample_id",
                                "Condition (biological variable)" = "condition",
                                "Batch / technical variable" = "batch",
                                "Covariate" = "covariate"),
                    selected = guess_role(col))
      })
    )
  })

  output$diagnoseButtonUI <- renderUI({
    req(rv$df)
    tagList(h4("3. Diagnose"), actionButton("diagnoseBtn", "Diagnose this design →", class = "btn-primary"))
  })

  get_roles <- function() {
    req(rv$df)
    cols <- names(rv$df)
    roles <- vapply(cols, function(col) {
      id <- paste0("role_", rv$col_ids[[col]])
      val <- input[[id]]
      if (is.null(val)) "ignore" else val
    }, character(1))
    names(roles) <- cols
    roles
  }

  # ---- Run diagnosis ----
  observeEvent(input$diagnoseBtn, {
    req(rv$df)
    roles <- get_roles()
    condition_col <- names(roles)[roles == "condition"]
    validate(need(length(condition_col) >= 1, 'Please mark one column as "Condition (biological variable)" to run diagnostics.'))
    condition_col <- condition_col[1]
    batch_cols <- names(roles)[roles == "batch"]
    covariate_cols <- names(roles)[roles == "covariate"]
    sample_id_cols <- names(roles)[roles == "sample_id"]

    df <- rv$df
    batches <- lapply(batch_cols, function(b) diagnose_batch(df, condition_col, b))
    covariates <- lapply(covariate_cols, function(cv) diagnose_covariate(df, condition_col, cv))
    replicates <- replicate_adequacy(df, condition_col)
    balance <- lapply(batch_cols, function(b) batch_balance(df, b))
    rank_check <- check_design_rank(df, condition_col, batch_cols, covariate_cols)
    vif_result <- tryCatch(compute_vif(df, condition_col, batch_cols, covariate_cols), error = function(e) list())
    full_screen <- tryCatch(full_covariate_screen(df, condition_col, sample_id_cols, roles), error = function(e) list())
    collinearity <- tryCatch(covariate_collinearity(df, covariate_cols), error = function(e) list())
    missingness <- tryCatch(missingness_screen(df, condition_col), error = function(e) list())
    pca_result <- tryCatch(compute_pca(df, batch_cols, covariate_cols), error = function(e) NULL)
    undeclared_flagged <- screen_undeclared_flagged(full_screen)

    rv$diag_result <- list(condition_col = condition_col, batches = batches,
                            covariates = covariates, replicates = replicates, balance = balance,
                            rank_check = rank_check, vif_result = vif_result,
                            full_screen = full_screen, collinearity = collinearity,
                            missingness = missingness, has_pca = !is.null(pca_result))

    # ---- Dynamic plots: one output per batch/covariate/pair/flagged column,
    # registered here and matched to plotOutput() ids built the same way
    # (via plot_id()) in the renderUI below. ----
    failed_plot <- function() { graphics::plot.new(); graphics::text(0.5, 0.5, "Could not render this plot") }

    lapply(batch_cols, function(b) {
      local({
        bb <- b
        output[[plot_id("plot_batch", bb)]] <- renderPlot({
          tryCatch(plot_condition_by_batch(df, condition_col, bb), error = function(e) failed_plot())
        })
      })
    })

    lapply(covariate_cols, function(cv) {
      local({
        cvcv <- cv
        output[[plot_id("plot_cov", cvcv)]] <- renderPlot({
          tryCatch(plot_covariate_by_condition(df, condition_col, cvcv), error = function(e) failed_plot())
        })
      })
    })

    lapply(undeclared_flagged, function(s) {
      local({
        col <- s$covariate
        output[[plot_id("plot_screen", col)]] <- renderPlot({
          tryCatch(plot_covariate_by_condition(df, condition_col, col), error = function(e) failed_plot())
        })
      })
    })

    output$plot_vif <- renderPlot({
      tryCatch(plot_vif_bar(vif_result), error = function(e) failed_plot())
    })

    lapply(collinearity, function(cpair) {
      local({
        aa <- cpair$a; bb <- cpair$b
        output[[plot_id("plot_collin", aa, bb)]] <- renderPlot({
          tryCatch(plot_collinearity_pair(df, aa, bb), error = function(e) failed_plot())
        })
      })
    })

    lapply(missingness, function(m) {
      local({
        col <- m$column
        output[[plot_id("plot_missing", col)]] <- renderPlot({
          tryCatch(plot_missingness_by_group(df, condition_col, col), error = function(e) failed_plot())
        })
      })
    })

    output$plot_pca <- renderPlot({
      req(pca_result)
      tryCatch(plot_pca(pca_result, df[[condition_col]]), error = function(e) failed_plot())
    })
  })

  # ---- Render diagnosis results ----
  output$diagnosisResults <- renderUI({
    req(rv$diag_result)
    res <- rv$diag_result

    worst_batch_v <- if (length(res$batches) == 0) NA_real_ else max(sapply(res$batches, function(b) b$v), na.rm = TRUE)
    worst_cov_v <- if (length(res$covariates) == 0) NA_real_ else max(sapply(res$covariates, function(c) c$association), na.rm = TRUE)
    worst_cv <- if (length(res$balance) == 0) NA_real_ else max(sapply(res$balance, function(b) b$cv), na.rm = TRUE)

    resource <- assess_resource(worst_cv, res$replicates)
    interpret <- assess_interpretability(worst_batch_v, worst_cov_v)
    general <- assess_generalisability(res$covariates)
    rank_assess <- assess_rank(res$rank_check, res$condition_col)

    card <- function(title, a, tip = NULL) {
      div(class = "result-card",
          h5(title, if (!is.null(tip)) info_tip(tip)),
          span(class = paste0("flag flag-", a$cls), tier_label_risk(a$tier)),
          tags$ul(lapply(a$reasons, tags$li)))
    }

    batch_rows <- if (length(res$batches) > 0) {
      do.call(rbind, lapply(res$batches, function(b) data.frame(
        Batch = b$batch_col, `Cramers V` = round(b$v, 3), Severity = tier_label(b$severity),
        `p value` = signif(b$p_value, 3), check.names = FALSE)))
    } else NULL

    cov_rows <- if (length(res$covariates) > 0) {
      do.call(rbind, lapply(res$covariates, function(c) data.frame(
        Covariate = c$covariate, Method = c$method, Association = round(c$association, 3),
        Severity = tier_label(c$severity), `N levels` = c$n_levels, check.names = FALSE)))
    } else NULL

    vif_rows <- if (length(res$vif_result) > 0) {
      do.call(rbind, lapply(res$vif_result, function(v) {
        a <- assess_vif(v)
        data.frame(Term = v$term, GVIF = round(v$gvif, 2),
                   `GVIF (df-adj.)` = round(v$gvif_std, 2),
                   `Effective n` = if (is.na(v$n_effective)) "—" else round(v$n_effective, 1),
                   Severity = tier_label(a$tier), check.names = FALSE)
      }))
    } else NULL

    screen_rows <- if (length(res$full_screen) > 0) {
      do.call(rbind, lapply(res$full_screen, function(s) data.frame(
        Column = s$covariate, `Declared as` = s$declared_role, Method = s$method,
        Association = round(s$association, 3), Severity = tier_label(s$severity),
        Flag = if (!s$declared && !is.na(s$association) && s$association >= 0.3) "not declared as batch/covariate" else "",
        check.names = FALSE)))
    } else NULL

    collin_rows <- if (length(res$collinearity) > 0) {
      do.call(rbind, lapply(res$collinearity, function(cpair) data.frame(
        `Covariate A` = cpair$a, `Covariate B` = cpair$b, Method = cpair$method,
        Association = round(cpair$value, 3), Severity = tier_label(cpair$severity), check.names = FALSE)))
    } else NULL

    missing_rows <- if (length(res$missingness) > 0) {
      do.call(rbind, lapply(res$missingness, function(m) data.frame(
        Column = m$column, `N missing` = m$n_missing, `% missing` = m$pct_missing,
        `Assoc. w/ condition` = if (is.na(m$v)) "—" else round(m$v, 3),
        Severity = tier_label(m$severity), check.names = FALSE)))
    } else NULL

    undeclared_flagged <- screen_undeclared_flagged(res$full_screen)

    plot_grid <- function(ids, height = "260px") {
      div(class = "plot-grid", lapply(ids, function(id) plotOutput(id, height = height)))
    }

    tagList(
      div(class = "card-panel",
        h4("Criteria assessment"),
        div(class = "assessment-row",
          card("Rank / aliasing check", rank_assess, tip_rank),
          card("Resource usage", resource, tip_balance),
          card("Interpretability", interpret),
          card("Generalisability", general)
        )
      ),

      if (isTRUE(res$has_pca)) div(class = "card-panel",
        h4("PCA of batch/covariate design", info_tip("Principal components of the one-hot-encoded batch + covariate design (not molecular data — the app has none), colored by condition. If condition separates cleanly along PC1/PC2, that's a visual echo of the same confounding the checks above quantify numerically.")),
        plotOutput("plot_pca", height = "320px")
      ),

      div(class = "card-panel",
        h4("Variance Inflation Factor (VIF) per term", info_tip(tip_vif)),
        div(class = "box-with-plot",
          div(class = "box-table",
            if (!is.null(vif_rows)) renderTable(vif_rows, rownames = FALSE) else p("Need at least two design terms (condition + batch/covariate) to compute VIF.")
          ),
          if (!is.null(vif_rows)) plotOutput("plot_vif", height = "280px")
        )
      ),

      div(class = "card-panel",
        h4("Batch / condition association", info_tip(tip_batch_assoc)),
        div(class = "box-with-plot",
          div(class = "box-table",
            if (!is.null(batch_rows)) renderTable(batch_rows, rownames = FALSE) else p("No batch columns marked.")
          ),
          plot_grid(vapply(res$batches, function(b) plot_id("plot_batch", b$batch_col), character(1)))
        )
      ),

      div(class = "card-panel",
        h4("Covariate association", info_tip(tip_covariate_assoc)),
        div(class = "box-with-plot",
          div(class = "box-table",
            if (!is.null(cov_rows)) renderTable(cov_rows, rownames = FALSE) else p("No covariate columns marked.")
          ),
          plot_grid(vapply(res$covariates, function(c) plot_id("plot_cov", c$covariate), character(1)), height = "440px")
        )
      ),

      div(class = "card-panel",
        h4("Full covariate screen (every column, not just declared ones)", info_tip(tip_full_screen)),
        div(class = "box-with-plot",
          div(class = "box-table",
            if (!is.null(screen_rows)) renderTable(screen_rows, rownames = FALSE) else p("No other columns to screen."),
            if (length(undeclared_flagged) > 0) div(class = "hint", "Plots below are only for undeclared columns flagged here — declared batch/covariate columns are already plotted in their own box above.")
          ),
          plot_grid(vapply(undeclared_flagged, function(s) plot_id("plot_screen", s$covariate), character(1)), height = "440px")
        )
      ),

      div(class = "card-panel",
        h4("Collinearity among covariates", info_tip(tip_collinearity)),
        div(class = "box-with-plot",
          div(class = "box-table",
            if (!is.null(collin_rows)) renderTable(collin_rows, rownames = FALSE) else p("Need at least two covariates marked to check collinearity.")
          ),
          plot_grid(vapply(res$collinearity, function(cpair) plot_id("plot_collin", cpair$a, cpair$b), character(1)))
        )
      ),

      div(class = "card-panel",
        h4("Missingness vs. group", info_tip(tip_missingness)),
        div(class = "box-with-plot",
          div(class = "box-table",
            if (!is.null(missing_rows)) renderTable(missing_rows, rownames = FALSE) else p("No missing data detected in this sheet.")
          ),
          plot_grid(vapply(res$missingness, function(m) plot_id("plot_missing", m$column), character(1)))
        )
      ),

      div(class = "card-panel",
        h4("Samples per condition group"),
        renderTable(res$replicates, rownames = FALSE)
      )
    )
  })

  # ==========================================================================
  # Power calculator
  # ==========================================================================

  output$pcModeUI <- renderUI({
    a <- input$pcAssay
    if (a == "scrna") {
      radioButtons("pcScrnaMode", NULL,
                   choices = c("Rare population capture" = "capture", "Differential expression (pseudobulk)" = "de"),
                   selected = "capture", inline = TRUE)
    } else if (a == "wgs") {
      radioButtons("pcWgsMode", NULL,
                   choices = c("Variant-calling depth adequacy" = "depth", "Case-control association power" = "assoc"),
                   selected = "depth", inline = TRUE)
    } else NULL
  })

  output$pcFormUI <- renderUI({
    a <- input$pcAssay
    count_fields <- function(n_label = "Samples per group (n)") {
      tagList(
        fluidRow(
          column(4, numericInput("pc_n", n_label, value = 3)),
          column(4, numericInput("pc_depth", "Per-feature read depth (normalised mean count)", value = 10)),
          column(4, numericInput("pc_cv", "Biological coefficient of variation (CV)", value = 0.4, step = 0.05))
        ),
        fluidRow(
          column(4, numericInput("pc_log2fc", "Minimum log2 fold-change of interest", value = 1, step = 0.1)),
          column(4, numericInput("pc_alpha", "Significance threshold (alpha)", value = 0.05, step = 0.001)),
          column(4, numericInput("pc_target", "Your target power", value = 0.8, step = 0.05))
        ),
        div(class = "hint", "Depth guide: a typical moderately-expressed feature at standard bulk-seq depth sits around 10-20. CV guide: ~0.1-0.2 for cell lines/inbred model organisms, ~0.4 for outbred human cohorts (Hart et al. 2013).")
      )
    }
    cohend_fields <- function() {
      tagList(
        fluidRow(
          column(4, numericInput("pc_n", "Samples per group (n)", value = 4)),
          column(4, numericInput("pc_cohend", "Effect size (Cohen's d)", value = 0.8, step = 0.1)),
          column(4, numericInput("pc_alpha", "Significance threshold (alpha)", value = 0.05, step = 0.001))
        ),
        numericInput("pc_target", "Your target power", value = 0.8, step = 0.05),
        div(class = "hint", "Convention: small d ~0.2, medium ~0.5, large ~0.8.")
      )
    }

    if (a %in% c("bulk_rna", "atac")) {
      count_fields()
    } else if (a == "scrna") {
      mode <- input$pcScrnaMode %||% "capture"
      if (identical(mode, "capture")) {
        tagList(
          fluidRow(
            column(4, numericInput("pc_ncells", "Total cells planned (across the group)", value = 5000)),
            column(4, numericInput("pc_freq", "Expected population frequency (0-1)", value = 0.02, step = 0.005)),
            column(4, numericInput("pc_kmin", "Minimum cells needed to call it detected", value = 10))
          ),
          numericInput("pc_target", "Your target capture probability", value = 0.95, step = 0.01)
        )
      } else {
        tagList(count_fields(n_label = "Independent samples/individuals per group"),
                div(class = "hint", "Aggregate cells to pseudobulk per sample first — cells from the same individual are not independent replicates."))
      }
    } else if (a %in% c("prot_dia", "metabolomics", "other")) {
      cohend_fields()
    } else if (a == "prot_tmt") {
      tagList(cohend_fields(), checkboxInput("pc_multiplex", "Samples split across multiple TMT plexes", value = FALSE))
    } else if (a == "wgs") {
      mode <- input$pcWgsMode %||% "depth"
      if (identical(mode, "depth")) {
        tagList(
          fluidRow(
            column(4, numericInput("pc_depth2", "Mean sequencing depth (x)", value = 30)),
            column(4, numericInput("pc_af", "Expected allele fraction", value = 0.5, step = 0.05)),
            column(4, numericInput("pc_kmin2", "Minimum supporting reads to call confidently", value = 3))
          ),
          numericInput("pc_target", "Your target confidence", value = 0.99, step = 0.005)
        )
      } else {
        tagList(
          fluidRow(
            column(3, numericInput("pc_ncases", "Number of cases", value = 500)),
            column(3, numericInput("pc_ncontrols", "Number of controls", value = 500)),
            column(3, numericInput("pc_pcontrols", "Risk allele frequency in controls", value = 0.2, step = 0.01)),
            column(3, numericInput("pc_pcases", "Risk allele frequency in cases", value = 0.26, step = 0.01))
          ),
          fluidRow(
            column(6, numericInput("pc_alpha", "Significance threshold (alpha)", value = 5e-8)),
            column(6, numericInput("pc_target", "Your target power", value = 0.8, step = 0.05))
          ),
          div(class = "hint", "Genome-wide significance is conventionally 5e-8 — adjust if this is a targeted/candidate-gene test instead.")
        )
      }
    }
  })

  `%||%` <- function(a, b) if (is.null(a)) b else a

  verdict <- function(achieved, target) {
    if (achieved >= target) list(cls = "good", label = "meets your target")
    else if (achieved >= target - 0.15) list(cls = "caution", label = "below target")
    else list(cls = "poor", label = "well below target")
  }

  observeEvent(input$pcRunBtn, {
    a <- input$pcAssay
    result <- NULL

    if (a %in% c("bulk_rna", "atac") || (a == "scrna" && identical(input$pcScrnaMode, "de"))) {
      r <- count_power(input$pc_n, input$pc_depth, input$pc_cv, input$pc_log2fc, input$pc_alpha, input$pc_target)
      v <- verdict(r$achieved_power, input$pc_target)
      result <- list(
        assay_label = ASSAY_LABELS[[a]], achieved_power = r$achieved_power, target_power = input$pc_target,
        cls = v$cls, label = v$label,
        big = sprintf("%.2f", r$achieved_power), big_label = sprintf("achieved power at n=%s/group", input$pc_n),
        lines = c(sprintf("Effective SD of log2 fold-change at this depth/CV: %.3f.", r$sd_used),
                  if (is.finite(r$n_needed)) sprintf("Samples per group needed for power >= %.2f: %d.", input$pc_target, ceiling(r$n_needed))
                  else "Could not solve for a required n at this target — try a larger effect size or lower target."),
        caution = "Approximation in the spirit of Hart et al. (2013, J Comput Biol): Poisson shot noise + biological CV, tested via a t-test on the log2 scale. Not the literal RNASeqPower package."
      )
    } else if (a %in% c("prot_dia", "metabolomics", "other", "prot_tmt")) {
      r <- cohend_power(input$pc_n, input$pc_cohend, input$pc_alpha, input$pc_target)
      v <- verdict(r$achieved_power, input$pc_target)
      multiplex_note <- if (identical(a, "prot_tmt") && isTRUE(input$pc_multiplex))
        " Samples span multiple plexes: this treats variance as one pooled number. If plex assignment lines up with your condition groups, check that in the Design diagnostics tab first." else ""
      result <- list(
        assay_label = ASSAY_LABELS[[a]], achieved_power = r$achieved_power, target_power = input$pc_target,
        cls = v$cls, label = v$label,
        big = sprintf("%.2f", r$achieved_power), big_label = sprintf("achieved power at n=%s/group", input$pc_n),
        lines = if (is.finite(r$n_needed)) sprintf("Samples per group needed for power >= %.2f: %d.", input$pc_target, ceiling(r$n_needed))
                else "Could not solve for a required n at this target.",
        caution = paste0("Standard two-sample t-test power on Cohen's d (Oberg & Vitek 2009; Blaise et al. 2016 for why metabolomics effect sizes are hard to state a priori).", multiplex_note)
      )
    } else if (a == "scrna" && identical(input$pcScrnaMode, "capture")) {
      p <- capture_probability(input$pc_ncells, input$pc_freq, input$pc_kmin)
      nn <- capture_n_needed(input$pc_freq, input$pc_kmin, input$pc_target)
      v <- verdict(p, input$pc_target)
      result <- list(
        assay_label = ASSAY_LABELS[["scrna"]], achieved_power = p, target_power = input$pc_target,
        cls = v$cls, label = v$label,
        big = sprintf("%.3f", p), big_label = sprintf("probability of capturing >= %s cells at n=%s", input$pc_kmin, input$pc_ncells),
        lines = if (is.finite(nn)) sprintf("Total cells needed for >= %.0f%% confidence: %d.", input$pc_target * 100, ceiling(nn))
                else "Could not solve for a required cell count — try a lower confidence target.",
        caution = "Exact binomial coverage calculation, not a hypothesis test — this only asks 'will I see enough of it,' not 'is the difference real.'"
      )
    } else if (a == "wgs" && identical(input$pcWgsMode, "depth")) {
      p <- capture_probability(input$pc_depth2, input$pc_af, input$pc_kmin2)
      nn <- depth_n_needed(input$pc_af, input$pc_kmin2, input$pc_target)
      v <- verdict(p, input$pc_target)
      result <- list(
        assay_label = ASSAY_LABELS[["wgs"]], achieved_power = p, target_power = input$pc_target,
        cls = v$cls, label = v$label,
        big = sprintf("%.3f", p), big_label = sprintf("confidence of >= %s supporting reads at %sx", input$pc_kmin2, input$pc_depth2),
        lines = if (is.finite(nn)) sprintf("Depth needed for >= %.0f%% confidence: %dx.", input$pc_target * 100, ceiling(nn))
                else "Could not solve for a required depth.",
        caution = "Binomial read-sampling model (Sims et al. 2014, Nat Rev Genet) — assumes reads sample alleles independently at the stated fraction; real mapping bias will shift this."
      )
    } else if (a == "wgs" && identical(input$pcWgsMode, "assoc")) {
      r <- assoc_power(input$pc_ncases, input$pc_ncontrols, input$pc_pcases, input$pc_pcontrols, input$pc_alpha, input$pc_target)
      v <- verdict(r$achieved_power, input$pc_target)
      result <- list(
        assay_label = ASSAY_LABELS[["wgs"]], achieved_power = r$achieved_power, target_power = input$pc_target,
        cls = v$cls, label = v$label,
        big = sprintf("%.2f", r$achieved_power), big_label = sprintf("achieved power (n/group used: %d)", r$n_per_group_used),
        lines = if (is.finite(r$n_needed)) sprintf("Per-group n needed for power >= %.2f: %d.", input$pc_target, ceiling(r$n_needed))
                else "Could not solve for a required n at this alpha/effect size.",
        caution = "Two-proportion test on allele frequency — same statistical shape as the confounding chi-square check in the Design diagnostics tab. Uses the smaller of cases/controls as a conservative per-group n."
      )
    }

    rv$pc_result <- result
  })

  output$pcResults <- renderUI({
    req(rv$pc_result)
    r <- rv$pc_result
    div(class = "result-card",
        span(class = paste0("flag flag-", r$cls), r$label),
        span(style = "font-size:28px; margin-left:12px;", r$big),
        span(style = "color:#777; margin-left:8px;", r$big_label),
        tags$ul(lapply(r$lines, tags$li)),
        div(class = "caution-note", r$caution))
  })
}

shinyApp(ui, server)
