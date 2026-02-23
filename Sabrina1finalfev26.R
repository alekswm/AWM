# Script: GLM analysis (Gaussian/Gamma/Inverse-Gaussian) + diagnostics + post-tests
# Author: Aleksander Westphal Muniz
# Date: 2026-02-23
# AI assistance: initial code structure generated with the help of an AI assistant (Adapta ONE),
# then reviewed/edited by the author. All analytical decisions and interpretations are the author's.

# ==========================================================
# COMPLETE SCRIPT - GLM (Gaussian/Gamma/Inverse-Gaussian)
# Read data -> check balance -> boxplots -> fit models -> AIC
# standardized Pearson residuals -> outliers (Cook + Leverage)
# clean dataset -> GLM ANOVA (clean data) -> post-tests (Tukey + LSD)
# for the lowest-AIC model (ΔAIC < 2 => prefer Gaussian)
# + residual boxplots (all variables and models)
# + BAR PLOTS:
#   (i) means (original scale)
#   (ii) "mean residual" bars
#   (iii) mean bars (original scale) WITH LSD letters
# ==========================================================

# =========================
# 0) Packages
# =========================
pkgs <- c(
  "tidyverse", "readr", "janitor",
  "performance", "DHARMa", "emmeans",
  "broom", "statmod", "multcomp"
)

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# =========================
# 1) Read data and define variables
# =========================
dat <- read_table("sabrina1fev26B.txt", col_types = cols()) %>%
  clean_names() %>%
  mutate(trat = factor(trat, levels = c("NC", "PC", "2052", "2053")))

vars <- c("alt", "raiz", "mpas", "msr")

str(dat)
summary(dat)

# =========================
# 2) Check balance (n per treatment)
# =========================
balanceamento <- dat %>%
  count(trat, name = "n") %>%
  arrange(trat)
print(balanceamento)

na_check <- dat %>%
  group_by(trat) %>%
  summarise(across(all_of(vars), ~sum(is.na(.x)), .names = "na_{.col}"),
            .groups = "drop")
print(na_check)

# =========================
# 3) Boxplots of raw variables
# =========================
boxplot_gg <- function(df, var){
  ggplot(df, aes(x = trat, y = .data[[var]], fill = trat)) +
    geom_boxplot(width = 0.65, outlier.alpha = 0.7) +
    geom_jitter(width = 0.12, alpha = 0.30, size = 1.8) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none") +
    labs(x = "Treatment", y = var)
}

plots_raw <- purrr::map(vars, ~boxplot_gg(dat, .x))
names(plots_raw) <- vars

plots_raw$alt
plots_raw$raiz
plots_raw$mpas
plots_raw$msr

# =========================
# 4) Helper functions: GLM fit, AIC, Pearson residuals, outliers, post-tests, plots
# =========================

# 4.1) Fit the 3 GLM families
fit_glm_set <- function(df, y){
  f <- as.formula(paste0(y, " ~ trat"))
  
  m_gauss <- glm(f, data = df, family = gaussian(link = "identity"))
  m_gamma <- glm(f, data = df, family = Gamma(link = "log"))
  m_ig    <- glm(f, data = df, family = inverse.gaussian(link = "log"))
  
  list(gauss = m_gauss, gamma = m_gamma, invgauss = m_ig)
}

# 4.2) AIC comparison table
compare_aic <- function(models, var_name){
  tibble(
    variable = var_name,
    model = c("gauss", "gamma", "invgauss"),
    AIC = c(AIC(models$gauss), AIC(models$gamma), AIC(models$invgauss))
  ) %>%
    arrange(AIC) %>%
    mutate(delta_AIC = AIC - min(AIC))
}

# 4.3) Pick best model by AIC with parsimony:
# if Gaussian is within ΔAIC < 2 from the best, choose Gaussian.
pick_best_model <- function(models, prefer_gauss_if_delta_lt = 2){
  
  aics <- c(
    gauss    = AIC(models$gauss),
    gamma    = AIC(models$gamma),
    invgauss = AIC(models$invgauss)
  )
  
  ord <- sort(aics)
  best_name <- names(ord)[1]
  delta_best <- ord - ord[1]
  
  if(!isFALSE(prefer_gauss_if_delta_lt)){
    if("gauss" %in% names(delta_best) && (delta_best["gauss"] < prefer_gauss_if_delta_lt)){
      best_name <- "gauss"
    }
  }
  
  list(
    best_name = best_name,
    best_model = models[[best_name]],
    aic_table = tibble(
      model = names(aics),
      AIC = unname(aics),
      delta_AIC = unname(aics - min(aics))
    ) %>% arrange(AIC)
  )
}

# 4.4) Standardized Pearson residuals + influence measures
diag_influence <- function(model){
  n <- nobs(model)
  p <- length(coef(model))
  
  tibble(
    .id = seq_len(n),
    fitted = fitted(model),
    resid_pearson_std = rstandard(model, type = "pearson"),
    cook = cooks.distance(model),
    leverage = hatvalues(model),
    cook_thr = 4 / n,
    lev_thr = 2 * p / n
  )
}

flag_outliers <- function(diag_tbl){
  diag_tbl %>%
    mutate(
      out_cook = cook > cook_thr,
      out_lev  = leverage > lev_thr,
      out_any  = out_cook | out_lev
    )
}

# 4.5) Pearson residual boxplot by treatment
boxplot_resid <- function(df_joined, y, model_name){
  ggplot(df_joined, aes(x = trat, y = resid_pearson_std, fill = trat)) +
    geom_boxplot(width = 0.65, outlier.alpha = 0.7) +
    geom_jitter(width = 0.12, alpha = 0.30, size = 1.8) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none") +
    labs(
      x = "Treatment",
      y = "Standardized Pearson residual",
      title = paste0("Pearson residuals - ", y, " (", model_name, ")")
    ) +
    geom_hline(yintercept = 0, col = "red", linetype = 2)
}

# 4.6) Post-tests: Tukey and LSD via emmeans (no pairs() usage)
# Tukey: adjust="tukey"
# LSD: adjust="none"
post_tests_tukey_lsd <- function(model, type_resp = FALSE){
  
  emm <- if(type_resp) emmeans::emmeans(model, ~ trat, type = "response") else emmeans::emmeans(model, ~ trat)
  
  tukey <- emmeans::contrast(emm, method = "pairwise", adjust = "tukey")
  lsd   <- emmeans::contrast(emm, method = "pairwise", adjust = "none")
  
  cld_tukey <- multcomp::cld(emm, Letters = letters, adjust = "tukey")
  cld_lsd   <- multcomp::cld(emm, Letters = letters, adjust = "none")
  
  list(
    emmeans = emm,
    tukey = tukey,
    lsd = lsd,
    cld_tukey = cld_tukey,
    cld_lsd = cld_lsd
  )
}

# 4.7) Bar-plot data: means and "mean residual" (value - treatment mean)
mean_residual_bar_data <- function(df, y){
  df %>%
    select(trat, y = all_of(y)) %>%
    group_by(trat) %>%
    mutate(mean_trat = mean(y, na.rm = TRUE),
           resid_mean = y - mean_trat) %>%
    summarise(
      mean = mean(y, na.rm = TRUE),
      mean_resid = mean(resid_mean, na.rm = TRUE),
      sd_resid = sd(resid_mean, na.rm = TRUE),
      se_resid = sd(resid_mean, na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    )
}

barplot_means <- function(sum_tbl, y){
  ggplot(sum_tbl, aes(x = trat, y = mean, fill = trat)) +
    geom_col(width = 0.7) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none") +
    labs(
      x = "Treatment",
      y = paste0("Mean of ", y, " (original scale)"),
      title = paste0("Bar plot - treatment means (", y, ")")
    )
}

barplot_mean_residuals <- function(sum_tbl, y, use = c("se", "sd")){
  use <- match.arg(use)
  err <- if(use == "se") "se_resid" else "sd_resid"
  
  ggplot(sum_tbl, aes(x = trat, y = mean_resid, fill = trat)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = mean_resid - .data[[err]],
                      ymax = mean_resid + .data[[err]]),
                  width = 0.2, linewidth = 0.6) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none") +
    labs(
      x = "Treatment",
      y = paste0("Mean residual (", if(use == "se") "±SE" else "±SD", ")"),
      title = paste0("Bar plot - residuals around treatment mean (", y, ")")
    ) +
    geom_hline(yintercept = 0, col = "red", linetype = 2)
}

# 4.8) NEW: Mean bars (original scale) with LSD letters (best AIC model)
barplot_means_with_lsd_letters <- function(model, best_name, y, adjust_y = 0.05){
  
  emm_resp <- if(best_name %in% c("gamma", "invgauss")) {
    emmeans::emmeans(model, ~ trat, type = "response")
  } else {
    emmeans::emmeans(model, ~ trat)
  }
  
  cld_lsd <- multcomp::cld(emm_resp, Letters = letters, adjust = "none")
  
  plot_df <- as.data.frame(cld_lsd) %>%
    as_tibble() %>%
    mutate(
      .group = gsub(" ", "", .group),
      trat = as.factor(trat)
    )
  
  mean_col <- if("response" %in% names(plot_df)) "response" else "emmean"
  se_col <- "SE"
  
  ymax <- max(plot_df[[mean_col]], na.rm = TRUE)
  ypad <- ymax * adjust_y
  
  ggplot(plot_df, aes(x = trat, y = .data[[mean_col]], fill = trat)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = .data[[mean_col]] - .data[[se_col]],
                      ymax = .data[[mean_col]] + .data[[se_col]]),
                  width = 0.2, linewidth = 0.6) +
    geom_text(aes(y = .data[[mean_col]] + .data[[se_col]] + ypad, label = .group),
              size = 5) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none") +
    labs(
      x = "Treatment",
      y = paste0("Estimated mean of ", y, " (original scale)"),
      title = paste0("Means (original scale) + LSD letters | ", y, " | model: ", best_name)
    )
}

# =========================
# 5) Fit (full dataset): ANOVA + AIC + diagnostics (Pearson/Cook/Leverage)
# =========================
results <- list()
aic_all <- list()

for(y in vars){
  
  models <- fit_glm_set(dat, y)
  
  aic_tbl <- compare_aic(models, y)
  aic_all[[y]] <- aic_tbl
  
  cat("\n=========================\n")
  cat("FULL DATA | Variable:", y, "\n")
  cat("=========================\n")
  cat("\n--- AIC comparison (lower is better) ---\n")
  print(aic_tbl)
  
  # ANOVA: Gaussian -> F; Gamma/InvGaussian -> LR/Chisq
  cat("\n--- ANOVA Gaussian (F) ---\n")
  print(anova(models$gauss, test = "F"))
  
  cat("\n--- ANOVA Gamma (LR / Chisq) ---\n")
  print(anova(models$gamma, test = "Chisq"))
  
  cat("\n--- ANOVA Inverse-Gaussian (LR / Chisq) ---\n")
  print(anova(models$invgauss, test = "Chisq"))
  
  # Diagnostics
  diag_gauss <- diag_influence(models$gauss) %>% flag_outliers() %>% mutate(var = y, model = "gauss")
  diag_gamma <- diag_influence(models$gamma) %>% flag_outliers() %>% mutate(var = y, model = "gamma")
  diag_ig    <- diag_influence(models$invgauss) %>% flag_outliers() %>% mutate(var = y, model = "invgauss")
  
  results[[y]] <- list(
    models = models,
    diag = bind_rows(diag_gauss, diag_gamma, diag_ig)
  )
}

aic_summary <- bind_rows(aic_all) %>%
  group_by(variable) %>%
  arrange(AIC, .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

cat("\n=========================\n")
cat("AIC summary (full data)\n")
cat("=========================\n")
print(aic_summary)

# =========================
# 6) Residual boxplots (full data) - all variables and models
# =========================
for(y in vars){
  diag_tbl <- results[[y]]$diag
  
  base_join <- dat %>%
    mutate(.id = row_number()) %>%
    select(.id, trat)
  
  diag_tbl2 <- diag_tbl %>%
    left_join(base_join, by = ".id")
  
  print(boxplot_resid(filter(diag_tbl2, model == "gauss"), y, "gauss"))
  print(boxplot_resid(filter(diag_tbl2, model == "gamma"), y, "gamma"))
  print(boxplot_resid(filter(diag_tbl2, model == "invgauss"), y, "inv-gauss"))
}

# =========================
# 6.1) Bar plots (full data): means and mean residual bars
# =========================
for(y in vars){
  sum_tbl <- mean_residual_bar_data(dat, y)
  print(barplot_means(sum_tbl, y))
  print(barplot_mean_residuals(sum_tbl, y, use = "se"))
}

# =========================
# 7) Outliers (Cook + Leverage) and clean dataset
# Rule: remove observation if flagged for ANY variable and ANY model
# =========================
all_diag <- bind_rows(lapply(results, function(x) x$diag))

outlier_ids <- all_diag %>%
  group_by(.id) %>%
  summarise(
    out_any = any(out_any),
    out_cook = any(out_cook),
    out_lev = any(out_lev),
    n_flags = sum(out_any),
    .groups = "drop"
  ) %>%
  filter(out_any) %>%
  arrange(desc(n_flags))

cat("\n=========================\n")
cat("Observations flagged as outliers (Cook OR Leverage)\n")
cat("=========================\n")
print(outlier_ids)

dat_clean <- dat %>%
  mutate(.id = row_number()) %>%
  filter(!.id %in% outlier_ids$.id) %>%
  select(-.id)

cat("\nN original:", nrow(dat), "\n")
cat("N after removal:", nrow(dat_clean), "\n")

balanceamento_clean <- dat_clean %>%
  count(trat, name = "n") %>%
  arrange(trat)

cat("\n=========================\n")
cat("Balance check on clean data\n")
cat("=========================\n")
print(balanceamento_clean)

# =========================
# 8) CLEAN DATA: GLM ANOVA + AIC + best model selection
# =========================
results_clean <- list()
aic_all_clean <- list()

for(y in vars){
  
  models <- fit_glm_set(dat_clean, y)
  
  aic_tbl_clean <- compare_aic(models, y)
  aic_all_clean[[y]] <- aic_tbl_clean
  
  sel <- pick_best_model(models, prefer_gauss_if_delta_lt = 2)
  
  cat("\n====================================================\n")
  cat("CLEAN DATA | Variable:", y, "\n")
  cat("AIC table (clean data):\n")
  print(sel$aic_table)
  cat("Selected model (ΔAIC<2 => Gaussian):", sel$best_name, "\n")
  cat("====================================================\n")
  
  cat("\n--- ANOVA Gaussian (F) ---\n")
  an_gauss <- anova(models$gauss, test = "F")
  print(an_gauss)
  
  cat("\n--- ANOVA Gamma (LR / Chisq) ---\n")
  an_gamma <- anova(models$gamma, test = "Chisq")
  print(an_gamma)
  
  cat("\n--- ANOVA Inverse-Gaussian (LR / Chisq) ---\n")
  an_ig <- anova(models$invgauss, test = "Chisq")
  print(an_ig)
  
  results_clean[[y]] <- list(
    models = models,
    selection = sel,
    anova = list(
      gauss = an_gauss,
      gamma = an_gamma,
      invgauss = an_ig
    )
  )
}

aic_summary_clean <- bind_rows(aic_all_clean) %>%
  group_by(variable) %>%
  arrange(AIC, .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

cat("\n=========================\n")
cat("AIC summary (clean data)\n")
cat("=========================\n")
print(aic_summary_clean)

# =========================
# 9) CLEAN DATA: post-tests (Tukey + LSD) for the selected (best AIC) model
# + NEW: mean bars (original scale) with LSD letters
# =========================
post_best <- list()

for(y in vars){
  
  sel <- results_clean[[y]]$selection
  best_name <- sel$best_name
  best_model <- sel$best_model
  
  cat("\n====================================================\n")
  cat("POST-TESTS (clean data) | Variable:", y, "\n")
  cat("Selected model:", best_name, "\n")
  cat("====================================================\n")
  
  if(best_name %in% c("gamma","invgauss")){
    
    cat("\n--- emmeans on the LINK scale ---\n")
    pt_link <- post_tests_tukey_lsd(best_model, type_resp = FALSE)
    print(pt_link$emmeans)
    cat("\nTukey (link):\n"); print(pt_link$tukey)
    cat("\nLSD (link):\n");   print(pt_link$lsd)
    cat("\nLSD letters (link):\n"); print(pt_link$cld_lsd)
    
    cat("\n--- emmeans on the RESPONSE (original) scale ---\n")
    pt_resp <- post_tests_tukey_lsd(best_model, type_resp = TRUE)
    print(pt_resp$emmeans)
    cat("\nTukey (response):\n"); print(pt_resp$tukey)
    cat("\nLSD (response):\n");   print(pt_resp$lsd)
    cat("\nLSD letters (response):\n"); print(pt_resp$cld_lsd)
    
    # Bar plot: means (original scale) + LSD letters
    print(barplot_means_with_lsd_letters(best_model, best_name, y))
    
    post_best[[y]] <- list(
      best_model = best_name,
      aic_table = sel$aic_table,
      post_link = pt_link,
      post_response = pt_resp
    )
    
  } else {
    
    pt <- post_tests_tukey_lsd(best_model, type_resp = FALSE)
    print(pt$emmeans)
    cat("\nTukey:\n"); print(pt$tukey)
    cat("\nLSD:\n");   print(pt$lsd)
    cat("\nLSD letters:\n"); print(pt$cld_lsd)
    
    # Bar plot: means (original scale) + LSD letters
    print(barplot_means_with_lsd_letters(best_model, best_name, y))
    
    post_best[[y]] <- list(
      best_model = best_name,
      aic_table = sel$aic_table,
      post = pt
    )
  }
}

# =========================
# 10) CLEAN DATA: Pearson residual boxplots for all variables and models
# =========================
for(y in vars){
  
  models <- fit_glm_set(dat_clean, y)
  
  diag_gauss <- diag_influence(models$gauss) %>% mutate(var = y, model = "gauss")
  diag_gamma <- diag_influence(models$gamma) %>% mutate(var = y, model = "gamma")
  diag_ig    <- diag_influence(models$invgauss) %>% mutate(var = y, model = "invgauss")
  
  diag_tbl <- bind_rows(diag_gauss, diag_gamma, diag_ig)
  
  base_join <- dat_clean %>%
    mutate(.id = row_number()) %>%
    select(.id, trat)
  
  diag_tbl2 <- diag_tbl %>%
    left_join(base_join, by = ".id")
  
  print(boxplot_resid(filter(diag_tbl2, model == "gauss"), y, "gauss (clean data)"))
  print(boxplot_resid(filter(diag_tbl2, model == "gamma"), y, "gamma (clean data)"))
  print(boxplot_resid(filter(diag_tbl2, model == "invgauss"), y, "inv-gauss (clean data)"))
}

# =========================
# 10.1) Bar plots (clean data): means and mean residual bars
# =========================
for(y in vars){
  sum_tbl_clean <- mean_residual_bar_data(dat_clean, y)
  print(barplot_means(sum_tbl_clean, y))
  print(barplot_mean_residuals(sum_tbl_clean, y, use = "se"))
}

# =========================
# IMPORTANT NOTES
# =========================
# 1) Gamma and Inverse-Gaussian require the response to be strictly > 0.
#    If you have zeros/negatives, those fits may fail.
# 2) ANOVA tests:
#    - Gaussian: F test (anova(..., test="F"))
#    - Gamma/Inverse-Gaussian: LR/Chi-square (anova(..., test="Chisq"))
# 3) LSD is implemented as adjust="none" (no multiple-comparison correction).
# 4) Post-tests are run ONLY for the AIC-selected model, using the parsimony rule:
#    if Gaussian is within ΔAIC < 2 of the best model, choose Gaussian.