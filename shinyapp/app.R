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

# ---- Shared severity-tier vocabulary (a neutral low/medium/high severity
# scale, applied consistently everywhere an association or design property
# gets flagged, rather than judgment-laden words like "unfixable") ----
tier_label <- function(tier) {
  switch(tier,
    unfixable = "High",
    fixable = "Medium",
    clean = "Low",
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
                     ylim = c(0, max(m, 1) * 1.15),
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
                       ylim = c(0, max(m, 1) * 1.15),
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
                       ylim = c(0, max(tab, 1) * 1.15),
                       main = paste(a, "vs", b),
                       legend.text = rownames(tab), args.legend = list(x = "topright", bty = "n", cex = 0.8))
  }
}

plot_missingness_by_group <- function(df, condition_col, column) {
  op <- graphics::par(mar = c(4, 4, 3, 1)); on.exit(graphics::par(op))
  v <- df[[column]]
  is_missing <- is.na(v) | (is.character(v) & trimws(v) == "")
  pct <- tapply(is_missing, df[[condition_col]], function(x) 100 * mean(x))
  graphics::barplot(pct, col = "#E64626", ylab = "% missing", las = 1, ylim = c(0, max(pct, 1) * 1.15),
                     main = paste0("Missingness of ", column, " by ", condition_col))
}

plot_replicate_counts <- function(replicates) {
  op <- graphics::par(mar = c(4, 4, 3, 1)); on.exit(graphics::par(op))
  counts <- stats::setNames(replicates$n, replicates$condition)
  graphics::barplot(counts, col = PALETTE[seq_along(counts)], ylab = "samples", las = 1,
                     ylim = c(0, max(counts, 1) * 1.15), main = "Samples per condition group")
}

# MDS (multidimensional scaling) on a Gower distance over the batch +
# covariate design, colored by condition -- an exploratory visual
# complement to the rank/VIF/association checks, since the app has no
# molecular data to run a "real" ordination on. Gower distance (not plain
# PCA on one-hot dummies) is the statistically appropriate way to combine
# categorical and continuous columns into a single distance: each numeric
# column contributes a scaled absolute difference, each categorical column
# contributes a simple mismatch indicator, and the per-column distances are
# averaged. Base R only (outer/range for Gower, stats::cmdscale for the
# embedding) -- cluster::daisy() would be the usual shortcut but isn't part
# of webR's bundled package set.
gower_distance <- function(df, cols) {
  n <- nrow(df)
  total <- matrix(0, n, n)
  used <- 0
  for (col in cols) {
    v <- df[[col]]
    if (is.numeric(v)) {
      rng <- diff(range(v, na.rm = TRUE))
      if (!is.finite(rng) || rng == 0) next
      d <- abs(outer(v, v, "-")) / rng
    } else {
      v <- as.character(v)
      d <- outer(v, v, function(a, b) as.numeric(a != b))
    }
    total <- total + d
    used <- used + 1
  }
  if (used == 0) return(NULL)
  total / used
}

compute_mds <- function(df, batch_cols, covariate_cols) {
  design_cols <- unique(c(batch_cols, covariate_cols))
  if (length(design_cols) == 0 || nrow(df) < 4) return(NULL)
  d <- tryCatch(gower_distance(df, design_cols), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  fit <- tryCatch(stats::cmdscale(stats::as.dist(d), k = 2, eig = TRUE), error = function(e) NULL)
  if (is.null(fit) || is.null(fit$points) || ncol(fit$points) < 2) return(NULL)
  eig_pos <- pmax(fit$eig, 0)
  var_explained <- if (sum(eig_pos) > 0) eig_pos / sum(eig_pos) else rep(NA_real_, length(fit$eig))
  list(points = fit$points, var_explained = var_explained, n = nrow(df))
}

plot_mds <- function(mds_result, condition_vec) {
  if (is.null(mds_result)) return(invisible(NULL))
  op <- graphics::par(mar = c(4, 4, 3, 8), xpd = TRUE); on.exit(graphics::par(op))
  pts <- mds_result$points
  cond_f <- factor(condition_vec)
  cols <- PALETTE[as.integer(cond_f)]
  ve <- mds_result$var_explained
  xlab <- if (is.finite(ve[1])) sprintf("MDS1 (%.0f%% of variance)", ve[1] * 100) else "MDS1"
  ylab <- if (is.finite(ve[2])) sprintf("MDS2 (%.0f%% of variance)", ve[2] * 100) else "MDS2"
  graphics::plot(pts[, 1], pts[, 2], col = cols, pch = 19, cex = 1.3,
                 xlab = xlab, ylab = ylab,
                 main = "MDS of batch/covariate design (Gower distance)", las = 1)
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

# Two ways a batch/covariate column can break the rank/VIF/MDS design (all
# of which build a model.matrix()):
#  1. A categorical column with fewer than 2 levels (e.g. every sample is
#     the same batch) can't be contrast-encoded -- "contrasts can be
#     applied only to factors with 2 or more levels".
#  2. ANY missing values, even in a different column, can trigger the same
#     error indirectly: model.matrix()'s default na.omit drops every row
#     with an NA in *any* referenced column, which can silently make some
#     OTHER column single-level among what's left (e.g. a covariate that's
#     only ever recorded for one condition group drags that whole group
#     out of the design). This is a real, common case (single-batch
#     studies; a covariate only collected for one arm), not a user error,
#     so such columns are excluded here before those checks run, rather
#     than letting the whole check silently fail with a misleading
#     "mark more columns" message. Missingness itself is still fully
#     covered by the "Missingness vs. group" check.
has_enough_levels <- function(v) {
  if (is.numeric(v)) return(TRUE)
  length(unique(stats::na.omit(v))) >= 2
}

design_col_exclusion_reason <- function(v) {
  if (anyNA(v)) return("has missing values")
  if (!has_enough_levels(v)) return("has only 1 level")
  NA_character_
}

usable_design_cols <- function(df, cols) {
  Filter(function(col) is.na(design_col_exclusion_reason(df[[col]])), cols)
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

excluded_cols_text <- function(excluded_cols) {
  if (length(excluded_cols) == 0) return(NULL)
  paste(vapply(excluded_cols, function(e) sprintf("%s (%s)", e$column, e$reason), character(1)), collapse = ", ")
}

assess_rank <- function(rank_check, condition_col, excluded_cols = list()) {
  if (is.null(rank_check) || rank_check$status %in% c("single_term", "error")) {
    excl_text <- excluded_cols_text(excluded_cols)
    reason <- if (!is.null(excl_text)) {
      sprintf("Excluded from this check: %s. Mark a batch or covariate column with at least 2 levels and no missing values to check against %s.",
              excl_text, condition_col)
    } else {
      "Mark at least a condition column plus one batch or covariate column to run this check."
    }
    return(list(cls = "neutral", tier = "not_checked", reasons = reason))
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
      reason = "Inf: VIF is mathematically undefined because the overall design is exactly rank-deficient (some term is a perfect linear combination of others) — this makes every term's VIF undefined at once, not just the aliased ones. See the Rank / aliasing check above for which specific terms are responsible."))
  }
  gs <- vif_entry$gvif_std
  if (gs >= sqrt(10)) {
    return(list(tier = "fixable", cls = "caution",
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

# ---- Four-criteria assessment (ported from the JS version) ----

assess_resource <- function(worst_cv, rep_df) {
  sizes <- rep_df$n
  max_ratio <- if (min(sizes) > 0) max(sizes) / min(sizes) else Inf
  reasons <- c()
  if (!is.na(worst_cv) && (worst_cv > 0.4 || max_ratio > 3)) {
    cls <- "caution"; tier <- "fixable"
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
    cls <- "caution"; tier <- "fixable"
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
  } else {
    cls <- "caution"; tier <- "fixable"
  }
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

# ============================================================================
# Explainer content (Design diagnostics tab)
# ============================================================================

explainer_intro <- c(
  "This tool checks your experimental design for potential issues that may compromise the validity of your analyses. The checks are based solely on the metadata you provide, no molecular data is needed. Each check evaluates a different aspect of your design, like whether your condition groups are confounded with batch effects or covariates, whether there are missing data patterns that could bias your results, and whether all terms in your model are estimable.",
  "Input: a sample metadata sheet with group assignment, batch/technical variables, and covariates. Every check below runs on metadata alone, before any molecular data exists. Hover the ⓘ next to any check for its test, statistic, and how to read it."
)

# Per-check explanations, surfaced as hover tooltips directly on the result
# that they explain (rather than as a separate wall of text up top).
info_tip <- function(text) {
  tags$span(class = "info-tip", HTML("&#9432;"),
    tags$span(class = "tip-bubble", text))
}

tip_rank <- "Test: builds the design matrix X implied by condition + batch + covariates and computes rank(X) via QR/SVD. Statistic: rank deficiency d = ncol(X) − rank(X). Interpretation: d = 0 means every term is estimable given the others; d > 0 means at least one term is an exact linear combination of another (e.g. batch = A for every case sample). Necessary but not sufficient — full rank means the model is estimable, not that estimates are precise (that's what VIF is for)."
tip_vif <- "Test: regresses each term on all other terms in the design and takes VIF = 1/(1 − R²) (generalized VIF for multi-level factors, Fox & Monette 1992, so it's comparable across terms with different degrees of freedom). Interpretation: VIF = 1 means orthogonal to the rest of the design; VIF > 5 warrants attention, VIF > 10 is severe. Equivalent to estimating that term's effect with an effective sample size of n / VIF. Inf means VIF is mathematically undefined — the overall design is exactly rank-deficient (see the Rank / aliasing check above for which terms are responsible); this shows Inf for every term at once, not just the aliased ones, since they all share the same singular correlation matrix."
tip_batch_assoc <- "Test: cross-tabulates condition against each batch/technical variable (chi-square test of independence, or Cramer's V as effect size). Interpretation: a zero cell (a condition/batch combination with no samples) is a structural, unfixable confound — there's no data to estimate the batch effect independently of condition in that stratum. A singleton cell (n = 1) has no within-cell variance to separate batch noise from a genuine data point. A significant result with all cells reasonably populated is the fixable, partial-confounding case."
tip_balance <- "Test: how evenly is condition distributed within each batch, versus what proportional allocation would predict (entropy or chi-square goodness-of-fit). Interpretation: even spread means batch behaves like a random block; condition-homogeneous batches are the imbalance pattern that produces the confounding/VIF problems above — this is a fast per-batch scan before drilling into the full contingency table."
tip_full_screen <- "Test: every column in the sheet — not just ones marked Batch or Covariate — tested against condition (chi-square/Fisher for categorical columns, ANOVA/eta-squared for continuous ones). Interpretation: a flagged column is a candidate confounder whether or not anyone intended to model it — unlabelled fields like processing date or operator ID surface here even when nobody thought to flag them as batch."
tip_covariate_assoc <- "Test: association between each declared covariate and condition (Cramer's V for categorical, eta-squared from ANOVA for continuous). Interpretation: a high association means that covariate is a candidate confounder and should be included as an adjustment in any model, at some cost to power."
tip_collinearity <- "Test: pairwise correlation among covariates (Pearson r for continuous pairs, Cramer's V / eta-squared for categorical or mixed pairs). Interpretation: |r| > 0.8 between two continuous covariates means they carry substantially redundant information — including both inflates variance for both without adding independent explanatory power. This is a special case of the VIF problem, restricted to covariate-covariate pairs rather than covariate-vs-condition."
tip_missingness <- "Test: chi-square test on a 2×k table of (missing/present) × condition, for each column with missing data. Interpretation: a significant association means data isn't missing completely at random with respect to condition — standard imputation does not correct for group-dependent missingness, it will reproduce or amplify the bias rather than remove it."

# ============================================================================
# UI
# ============================================================================

theme_toggle_script <- r"-----(
(function() {
  var STORAGE_KEY = 'omics-calculator-theme';
  function systemPrefersDark() {
    return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
  }
  function applyTheme(theme) {
    document.documentElement.setAttribute('data-theme', theme);
  }
  function initTheme() {
    var saved = null;
    try { saved = localStorage.getItem(STORAGE_KEY); } catch (e) {}
    applyTheme(saved || (systemPrefersDark() ? 'dark' : 'light'));
  }
  initTheme();
  document.addEventListener('click', function(e) {
    var btn = e.target.closest ? e.target.closest('#themeToggleBtn') : null;
    if (!btn) return;
    var current = document.documentElement.getAttribute('data-theme') || 'light';
    var next = current === 'dark' ? 'light' : 'dark';
    applyTheme(next);
    try { localStorage.setItem(STORAGE_KEY, next); } catch (e) {}
  });
})();
)-----"

ui <- fluidPage(
  title = "Omics Experiment Calculator",
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$script(HTML(theme_toggle_script))
  ),
  div(class = "top-navbar",
      span(class = "navbar-brand", "Sydney Informatics Hub - Omics Experiment Calculator"),
      tags$button(id = "themeToggleBtn", class = "theme-toggle-btn", type = "button",
        tags$span(class = "icon-light", "\U0001F319 Dark"),
        tags$span(class = "icon-dark", "\U00002600 Light")
      )
  ),
  div(class = "app-container",
    tabsetPanel(
      id = "mainTabs",

      tabPanel("Design diagnostics",
        br(),
        div(class = "card-panel",
          div(class = "intro-text", lapply(explainer_intro, p))
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

  rv <- reactiveValues(df = NULL, col_ids = NULL, diag_result = NULL)

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

    # model.matrix() needs >= 2 levels per categorical term and no missing
    # values (see design_col_exclusion_reason() above for why), so such
    # columns are excluded from the rank/VIF/MDS design specifically --
    # they're still fully covered by the association/missingness tables
    # above via diagnose_batch/diagnose_covariate/missingness_screen.
    rank_batch_cols <- usable_design_cols(df, batch_cols)
    rank_covariate_cols <- usable_design_cols(df, covariate_cols)
    excluded_cols_raw <- setdiff(c(batch_cols, covariate_cols), c(rank_batch_cols, rank_covariate_cols))
    excluded_cols <- lapply(excluded_cols_raw, function(col) {
      list(column = col, reason = design_col_exclusion_reason(df[[col]]))
    })

    rank_check <- check_design_rank(df, condition_col, rank_batch_cols, rank_covariate_cols)
    vif_result <- tryCatch(compute_vif(df, condition_col, rank_batch_cols, rank_covariate_cols), error = function(e) list())
    full_screen <- tryCatch(full_covariate_screen(df, condition_col, sample_id_cols, roles), error = function(e) list())
    collinearity <- tryCatch(covariate_collinearity(df, covariate_cols), error = function(e) list())
    missingness <- tryCatch(missingness_screen(df, condition_col), error = function(e) list())
    mds_result <- tryCatch(compute_mds(df, rank_batch_cols, rank_covariate_cols), error = function(e) NULL)
    undeclared_flagged <- screen_undeclared_flagged(full_screen)

    rv$diag_result <- list(condition_col = condition_col, batches = batches,
                            covariates = covariates, replicates = replicates, balance = balance,
                            rank_check = rank_check, vif_result = vif_result,
                            full_screen = full_screen, collinearity = collinearity,
                            missingness = missingness, has_mds = !is.null(mds_result),
                            excluded_cols = excluded_cols)

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

    output$plot_mds <- renderPlot({
      req(mds_result)
      tryCatch(plot_mds(mds_result, df[[condition_col]]), error = function(e) failed_plot())
    })

    output$plot_replicates <- renderPlot({
      tryCatch(plot_replicate_counts(replicates), error = function(e) failed_plot())
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
    rank_assess <- assess_rank(res$rank_check, res$condition_col, res$excluded_cols)

    card <- function(title, a, tip = NULL) {
      div(class = "result-card",
          h5(title, if (!is.null(tip)) info_tip(tip)),
          span(class = paste0("flag flag-", a$cls), tier_label_risk(a$tier)),
          lapply(a$reasons, p))
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
                   Severity = tier_label(a$tier), Note = a$reason, check.names = FALSE)
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

      div(class = "card-panel",
        h4("MDS of batch/covariate design", info_tip("Multidimensional scaling on a Gower distance over the batch + covariate design (not molecular data — the app has none), colored by condition. Gower distance properly combines categorical and continuous columns (unlike plain PCA on one-hot dummies): each numeric column contributes a scaled difference, each categorical column a match/mismatch indicator. If condition separates cleanly, that's a visual echo of the same confounding the checks above quantify numerically.")),
        if (isTRUE(res$has_mds)) plotOutput("plot_mds", height = "320px")
        else p(if (!is.null(excluded_cols_text(res$excluded_cols)))
          sprintf("Not enough usable design columns to run MDS — excluded: %s.", excluded_cols_text(res$excluded_cols))
          else "Mark at least one usable batch/covariate column (and at least 4 samples) to run MDS.")
      ),

      div(class = "card-panel",
        h4("Variance Inflation Factor (VIF) per term", info_tip(tip_vif)),
        div(class = "box-with-plot",
          div(class = "box-table",
            if (!is.null(vif_rows)) renderTable(vif_rows, rownames = FALSE) else p("Need at least two design terms (condition + batch/covariate) to compute VIF."),
            if (!is.null(excluded_cols_text(res$excluded_cols))) div(class = "hint",
              sprintf("Excluded from VIF/rank/MDS: %s — still covered by the association and missingness tables elsewhere on this page.",
                      excluded_cols_text(res$excluded_cols)))
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
        div(class = "box-with-plot",
          div(class = "box-table", renderTable(res$replicates, rownames = FALSE)),
          plot_grid("plot_replicates")
        )
      )
    )
  })

}

shinyApp(ui, server)
