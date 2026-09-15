########################################################################
# 1. Load required libraries
########################################################################
library(haven)
library(dplyr)
library(gtsummary)
library(flextable)
library(officer)
library(logistf)
library(broom)
library(purrr)

########################################################################
# 1b. Read the "variance_stats" SAS dataset (row-level data).
########################################################################
variance_file_path <- "C:\\Users\\hujoe\\Dropbox\\Carlisle\\R-diabetes\\SAS datafiles\\variance_stats.sas7bdat"
variance_data <- read_sas(variance_file_path)
names(variance_data)

########################################################################
# 2. Derive the zombieformanuscript equivalent from variance_data
########################################################################
sas_max <- function(x) if (all(is.na(x))) x[1] else max(x, na.rm = TRUE)
sas_min <- function(x) if (all(is.na(x))) x[1] else min(x, na.rm = TRUE)

sas_data <- variance_data %>%
  group_by(PMID) %>%
  summarise(
    Authors__judgement      = sas_max(Authors__judgement),
    Support_for_judgement   = sas_max(Support_for_judgement),
    funding1                = sas_max(funding1),
    First_Author            = sas_max(First_Author),
    max_month               = sas_max(month),
    max_Hba1c_effect1       = sas_min(effect1),
    max_Hba1c_effect        = sas_min(effect),
    lncvrsigma              = sas_max(lncvrsigma),
    effectsigma             = sas_max(effectsigma),
    sigmacorrelation        = sas_max(sigmacorrelation),
    tukeycorrelation        = sas_max(tukeycorrelation),
    sigmacorrelation_GLY    = sas_max(sigmacorrelation_GLY),
    tukeycorrelation_GLY    = sas_max(tukeycorrelation_GLY),
    pathological_table      = sas_max(pathological_table),
    non_genuine_data        = sas_max(non_genuine_data),
    Other_bias              = sas_max(Other_bias),
    random_quality          = sas_max(random_quality),
    blinding_quality        = sas_max(blinding_quality),
    ascertainment           = sas_max(ascertainment),
    nGroups                 = sas_max(nGroups),
    numb_authors            = sas_max(numb_authors),
    Publication_Year        = sas_max(Publication_Year),
    sample_size             = sas_max(sample_size),
    origin_data             = sas_max(origin_data),
    retraction_countries    = sas_max(retraction_countries),
    pubmed                  = sas_max(pubmed),
    crossover               = sas_max(crossover),
    zombie05                = sas_max(zombie05),
    zombie01                = sas_max(zombie01),
    zombie001               = sas_max(zombie001),
    flag1                   = sas_max(flag1),
    unknown_funding         = sas_max(unknown_funding),
    subjective              = sas_max(subjective),
    p_onesided              = sas_max(P_onesided),
    p_twosided              = sas_max(P_twosided),
    CochraneID              = sas_max(CochraneID),
    CochraneID1             = sas_max(CochraneID1),
    CochraneID2             = sas_max(CochraneID2),
    CochraneID3             = sas_max(CochraneID3),
    nMeasures               = sas_max(nMeasures),
    correlation_available   = sas_max(correlation_available),
    duration                = sas_max(duration),
    .groups = "drop"
  ) %>%
  arrange(PMID)

names(sas_data)

########################################################################
# 3. Variable labels and factor coding
########################################################################
continuous_labels <- list(
  random_quality             = "Randomization",
  blinding_quality           = "Blinding",
  ascertainment              = "Ascertainment",
  Publication_Year           = "Publication year",
  sample_size                = "Sample size",
  duration                   = "Duration of trial",
  pubmed                     = "Is trial Pubmed indexed?",
  retraction_countries                   = "Is trial from country with high retraction ranking?",
  unknown_funding            = "Does trial fail to report funding source?",
  zombie01                   = "Does trial have a Carlisle-Stouffer-Fisher p-value < 0.01?",
  tukeycorrelation_GLY       = "Does trial have a Tukey outlier in FPG correlation?",
  non_genuine_data           = "Does trial have a table with three or more \u03c1s > |0.99|?",
  pathological_table         = "Does trial have a table with three or more \u03c1s > |0.99| or a glycemic |\u03c1| >1?"
)

# Numeric variables with few unique values that gtsummary would otherwise
# auto-classify as categorical; forced to continuous in every tbl_summary().
force_continuous_vars <- c("random_quality", "blinding_quality", "ascertainment")

recode_factors <- function(df) {
  df %>%
    mutate(
      retraction_countries             = factor(retraction_countries,             levels = c(0, 1), labels = c("No", "Yes")),
      unknown_funding      = factor(unknown_funding,      levels = c(0, 1), labels = c("No", "Yes")),
      pubmed               = factor(pubmed,               levels = c(0, 1), labels = c("No", "Yes")),
      zombie01             = factor(zombie01,             levels = c(0, 1), labels = c("No", "Yes")),
      tukeycorrelation_GLY = factor(tukeycorrelation_GLY, levels = c(0, 1), labels = c("No", "Yes")),
      correlation_available = factor(
        correlation_available,
        levels = c(0, 1),
        labels = c("HbA1c corr not calculable", "HbA1c corr calculable")
      ),
      tukeycorrelation = factor(
        tukeycorrelation,
        levels = c(0, 1),
        labels = c("Non-anomalous \u03c1", "Anomalous \u03c1")
      ),
      # Parametric 3-sigma definition of an anomalous HbA1c rho (Tables S2-S3)
      sigmacorrelation = factor(
        sigmacorrelation,
        levels = c(0, 1),
        labels = c("Non-anomalous \u03c1", "Anomalous \u03c1")
      ),
      pathological_table   = factor(pathological_table,   levels = c(0, 1), labels = c("No", "Yes")),
      non_genuine_data     = factor(non_genuine_data,     levels = c(0, 1), labels = c("No", "Yes"))
    )
}

# Publication_Year_c: one unit = one decade, 0 = year 2000 (aids Firth
# convergence). sample_size_10: sample size per 10 subjects.
data_pmid <- sas_data %>%
  mutate(sample_size_10     = sample_size / 10,
         Publication_Year_c = (Publication_Year - 2000) / 10)
data_pmid <- recode_factors(data_pmid)

########################################################################
# 3b. Variable grouping shared by all three tables
#     Each table shows the same three bold group headings; a table only
#     displays those group members that are present in its variable set.
########################################################################
var_groups <- list(
  "Data-intrinsic integrity-concerns" = c(
    "zombie01", "tukeycorrelation_GLY", "non_genuine_data", "pathological_table"
  ),
  "Contextual and structural trial correlates " = c(
    "Publication_Year", "sample_size", "duration", "pubmed",
    "retraction_countries", "unknown_funding"
  ),
  "Methodological-quality scores" = c(
    "random_quality", "blinding_quality", "ascertainment"
  )
)

# Table 2 uses the centred/scaled versions of two variables
var_groups_tbl2 <- var_groups
var_groups_tbl2[[2]] <- c(
  "Publication_Year_c", "sample_size_10", "duration", "pubmed",
  "retraction_countries", "unknown_funding"
)

# Insert a bold group-heading row before each group in a gtsummary table
# (applied AFTER add_p()/add_n(), so no statistics are affected).
add_group_headers <- function(tbl, groups) {
  tbl %>%
    modify_table_body(function(body) {
      pieces <- list()
      for (i in seq_along(groups)) {
        vars <- intersect(groups[[i]], unique(body$variable))
        if (length(vars) == 0) next
        hdr <- body[NA_integer_, ]            # one all-NA row, same columns/types
        hdr$variable <- paste0("hdr_", i)
        hdr$row_type <- "label"
        hdr$label    <- names(groups)[i]
        rows <- body %>%
          filter(variable %in% vars) %>%
          arrange(match(variable, vars))
        pieces[[length(pieces) + 1]] <- bind_rows(hdr, rows)
      }
      bind_rows(pieces)
    }) %>%
    modify_table_styling(
      columns     = label,
      rows        = startsWith(variable, "hdr_"),
      text_format = "bold"
    )
}

# After as_flex_table(): indent variable labels under their heading and
# levels/missing rows one step further (header rows stay flush left).
indent_under_headers <- function(ft, tbl) {
  body   <- tbl$table_body
  is_hdr <- startsWith(body$variable, "hdr_")
  is_lab <- body$row_type == "label" & !is_hdr
  is_lvl <- !is_hdr & !is_lab
  if (any(is_lab)) ft <- ft %>% padding(i = which(is_lab), j = 1, padding.left = 12)
  if (any(is_lvl)) ft <- ft %>% padding(i = which(is_lvl), j = 1, padding.left = 24)
  ft
}

########################################################################
# 4. Table S1 (n = 305): by correlation_available
########################################################################
tableS1_vars <- c(
  "random_quality", "blinding_quality", "ascertainment",
  "Publication_Year", "sample_size", "duration", "pubmed",
  "zombie01", "retraction_countries", "unknown_funding"
)

binary_vars_S1 <- c("pubmed", "zombie01", "retraction_countries", "unknown_funding")

tblS1 <- data_pmid %>%
  select(all_of(tableS1_vars), correlation_available) %>%
  tbl_summary(
    by = correlation_available,
    label = continuous_labels,
    # Binary No/Yes factors forced to "categorical" so both levels are shown
    # rather than collapsing to a single "Yes" row.
    type = list(
      all_of(force_continuous_vars) ~ "continuous",
      all_of(binary_vars_S1)        ~ "categorical"
    ),
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous()  ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    missing = "ifany",
    missing_text = "Missing"
  ) %>%
  add_p(
    test = list(
      all_continuous()  ~ "kruskal.test",
      all_categorical() ~ "fisher.test"
    ),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  add_n() %>%
  add_group_headers(var_groups)

########################################################################
# 5. Table 1 (n = 90): correlation_available == 1, by tukeycorrelation
########################################################################
data_avail <- data_pmid %>%
  filter(correlation_available == "HbA1c corr calculable")

cat("Subset with correlation_available=1:", nrow(data_avail), "rows\n")

table1_vars <- c(
  "random_quality", "blinding_quality", "ascertainment",
  "Publication_Year", "sample_size", "duration", "pubmed",
  "retraction_countries", "unknown_funding", "zombie01",
  "tukeycorrelation_GLY", "non_genuine_data"
)

binary_vars_1 <- c(binary_vars_S1,
                   "tukeycorrelation_GLY", "non_genuine_data")

# Table 1 (by = tukeycorrelation) and Table S2 (by = sigmacorrelation) are
# identical apart from the stratifying variable, so one function builds both.
make_char_table <- function(by_var) {
  data_avail %>%
  select(all_of(table1_vars), all_of(by_var)) %>%
  tbl_summary(
    by = all_of(by_var),
    label = continuous_labels,
    type = list(
      all_of(force_continuous_vars) ~ "continuous",
      all_of(binary_vars_1)         ~ "categorical"
    ),
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous()  ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    missing = "ifany",
    missing_text = "Missing"
  ) %>%
  add_p(
    test = list(
      all_continuous()  ~ "kruskal.test",
      all_categorical() ~ "fisher.test"
    ),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  add_n() %>%
  add_group_headers(var_groups)
}

tbl1  <- make_char_table("tukeycorrelation")   # Table 1:  Tukey inner fences
tblS2 <- make_char_table("sigmacorrelation")   # Table S2: parametric 3-sigma

########################################################################
# 6. Table 2 / Table S3 (n = 90): Univariate Firth logistic regression
#    Table 2 outcome = tukeycorrelation; Table S3 outcome = sigmacorrelation
########################################################################
predictors_tbl2 <- c(
  "random_quality", "blinding_quality", "ascertainment",
  "Publication_Year_c", "sample_size_10", "duration", "pubmed",
  "retraction_countries", "unknown_funding", "zombie01",
  "tukeycorrelation_GLY", "non_genuine_data"
)

predictor_labels <- c(
  random_quality            = "Randomization",
  blinding_quality          = "Blinding",
  ascertainment             = "Ascertainment",
  Publication_Year_c        = "Publication year (per decade, centered at 2000)",
  sample_size_10            = "Sample size (per 10)",
  duration                  = "Duration of trial",
  pubmed                    = "Is trial Pubmed indexed?",
  retraction_countries                  = "Is trial from country with high retraction ranking?",
  unknown_funding           = "Does trial fail to report funding source?",
  zombie01                  = "Does trial have a Carlisle-Stouffer-Fisher p-value < 0.01?",
  tukeycorrelation_GLY      = "Does trial have a Tukey outlier in FPG correlation?",
  non_genuine_data          = "Does trial have a table with three or more \u03c1s > |0.99|?",
  pathological_table        = "Does trial have a table with three or more \u03c1s > |0.99| or a glycemic |\u03c1| >1?"
)

run_firth_single <- function(predictor, data, outcome_var = "outcome") {
  formula_str <- paste(outcome_var, "~", predictor)

  # Capture the fitted object AND any logistf warnings (e.g. iteration limit).
  warning_msg  <- NULL
  fit <- tryCatch(
    withCallingHandlers(
      logistf(as.formula(formula_str), data = data, pl = TRUE, firth = TRUE),
      warning = function(w) {
        warning_msg <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) NULL
  )

  if (!is.null(warning_msg)) {
    cat(sprintf("[logistf WARNING] predictor = %-35s | %s\n", predictor, warning_msg))
  }

  converged <- is.null(warning_msg)

  if (is.null(fit)) {
    return(tibble(Variable = predictor, Term = predictor, converged = FALSE,
                  OR = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_,
                  p_value = NA_real_, n_obs = NA_integer_))
  }
  coef_names <- names(coef(fit))[-1]
  map_dfr(coef_names, function(term) {
    idx <- which(names(coef(fit)) == term)
    tibble(
      Variable = predictor, Term = term, converged = converged,
      OR       = exp(coef(fit)[idx]),
      CI_lower = exp(fit$ci.lower[idx]),
      CI_upper = exp(fit$ci.upper[idx]),
      p_value  = fit$prob[idx],
      n_obs    = fit$n
    )
  })
}

format_or <- function(or, lo, hi) {
  ifelse(is.na(or), "\u2014",
         sprintf("%.2f (%.2f\u2013%.2f)", or, lo, hi))
}

firth_footnote <- "\u2020 Firth logistic regression did not converge (maximum iterations exceeded); estimate and CI may be unreliable."

# Fit the univariate Firth models for one outcome definition and return the
# formatted flextable. outcome_factor is the 0/1 factor to model
# ("tukeycorrelation" for Table 2, "sigmacorrelation" for Table S3).
make_firth_table <- function(outcome_factor, footnotes = firth_footnote) {

data_avail_lr <- data_avail %>%
  mutate(outcome = as.integer(.data[[outcome_factor]]) - 1L)

tbl2_results <- map_dfr(predictors_tbl2, ~ run_firth_single(.x, data_avail_lr))

cat(sprintf("\n=== Firth logistic regression convergence summary (outcome = %s) ===\n",
            outcome_factor))
conv_summary <- tbl2_results %>%
  distinct(Variable, converged) %>%
  arrange(converged, Variable)
for (i in seq_len(nrow(conv_summary))) {
  status <- if (conv_summary$converged[i]) "converged OK" else "*** DID NOT CONVERGE ***"
  cat(sprintf("  %-40s %s\n", conv_summary$Variable[i], status))
}
cat("=====================================================\n\n")

tbl2_rows <- tbl2_results %>%
  mutate(
    Label     = unname(predictor_labels[Variable]),
    OR_CI     = format_or(OR, CI_lower, CI_upper),
    OR_CI     = ifelse(!converged, paste0(OR_CI, " \u2020"), OR_CI),
    p_display = ifelse(is.na(p_value), "\u2014",
                       ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))),
    n_obs     = as.character(n_obs),
    is_hdr    = FALSE
  ) %>%
  select(Variable, Label, OR_CI, p_display, n_obs, is_hdr)

# Insert the same three bold group headings used in Tables S1 and 1
tbl2_display <- map_dfr(seq_along(var_groups_tbl2), function(i) {
  vars <- intersect(var_groups_tbl2[[i]], unique(tbl2_rows$Variable))
  if (length(vars) == 0) return(NULL)
  bind_rows(
    tibble(Variable = paste0("hdr_", i), Label = names(var_groups_tbl2)[i],
           OR_CI = "", p_display = "", n_obs = "", is_hdr = TRUE),
    tbl2_rows %>% filter(Variable %in% vars) %>% arrange(match(Variable, vars))
  )
})

hdr_rows_2 <- which(tbl2_display$is_hdr)
var_rows_2 <- which(!tbl2_display$is_hdr)

ft2 <- flextable(tbl2_display %>% select(Label, OR_CI, p_display, n_obs)) %>%
  set_header_labels(
    Label     = "Variable",
    OR_CI     = "Odds Ratio (95% CI)",
    p_display = "P-value",
    n_obs     = "N"
  ) %>%
  add_footer_lines(footnotes) %>%
  autofit() %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  bold(i = hdr_rows_2, j = "Label", part = "body") %>%
  padding(i = var_rows_2, j = "Label", padding.left = 12, part = "body") %>%
  fontsize(size = 9, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  align(j = c("OR_CI", "p_display", "n_obs"), align = "center", part = "all") %>%
  align(j = "Label", align = "left", part = "all") %>%
  align(part = "footer", align = "left") %>%
  width(j = "Label",     width = 3.5) %>%
  width(j = "OR_CI",     width = 2.2) %>%
  width(j = "p_display", width = 1.0) %>%
  width(j = "n_obs",     width = 0.6)

ft2
}

# Footnote shown under Tables 1 and 2 (Tukey definition)
tukey_footnote <- "Anomalous \u03c1 defined by Tukey inner fences. Parallel analyses using parametric 3-sigma limits are in Tables S2\u2013S3."

ft2  <- make_firth_table("tukeycorrelation",
                         footnotes = c(tukey_footnote, firth_footnote))   # Table 2
ftS3 <- make_firth_table("sigmacorrelation")                              # Table S3

########################################################################
# 7. Convert gtsummary Tables S1, 1 & S2 to flextable with margin-safe widths
#    (landscape page, usable width ~9.5 in)
########################################################################
apply_ft_format <- function(ft, label_width = 2.5, data_col_width = 1.1) {
  ft <- ft %>%
    autofit() %>%
    theme_booktabs() %>%
    bold(part = "header") %>%
    fontsize(size = 9, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    align(j = 1, align = "left", part = "all")

  ft <- ft %>% width(j = 1, width = label_width)

  ncols <- length(ft$col_keys)
  if (ncols > 1) {
    ft <- ft %>% width(j = 2:ncols, width = data_col_width)
  }
  ft
}

ftS1 <- tblS1 %>% as_flex_table() %>%
  apply_ft_format(label_width = 2.5, data_col_width = 1.3) %>%
  indent_under_headers(tblS1)
ft1  <- tbl1  %>% as_flex_table() %>%
  apply_ft_format(label_width = 2.5, data_col_width = 1.3) %>%
  indent_under_headers(tbl1) %>%
  add_footer_lines(tukey_footnote) %>%
  align(part = "footer", align = "left") %>%
  fontsize(size = 9, part = "footer") %>%
  font(fontname = "Times New Roman", part = "footer")
ftS2 <- tblS2 %>% as_flex_table() %>%
  apply_ft_format(label_width = 2.5, data_col_width = 1.3) %>%
  indent_under_headers(tblS2)

########################################################################
# 8. Export to Word (LANDSCAPE orientation)
########################################################################
output_dir <- "C:/Users/hujoe/Dropbox/Carlisle/R-diabetes/Submission_manuscripts/correlation/Clinical Trials"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_file <- file.path(
  output_dir,
  paste0("Systematic_Review_Tukey_NGD", format(Sys.Date(), "%Y-%m-%d"), ".docx")
)

landscape_props <- prop_section(
  page_size   = page_size(orient = "landscape", width = 11, height = 8.5),
  page_margins = page_mar(top = 0.75, bottom = 0.75, left = 0.75, right = 0.75)
)

doc <- read_docx() %>%

  # --- Table S1 ---
  body_add_par(
    "Table S1. Trial characteristics by whether an HbA1c ρ could be recovered",
    style = "heading 1"
  ) %>%
  body_add_flextable(ftS1) %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_landscape() %>%

  # --- Table 1 ---
  body_add_par(
    "Table 1. Trial characteristics by anomalous versus non-anomalous HbA1c ρ",
    style = "heading 1"
  ) %>%
  body_add_flextable(ft1) %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_landscape() %>%

  # --- Table 2 ---
  body_add_par(
    "Table 2. Univariate Firth-penalized logistic regression of anomalous HbA1c ρ",
    style = "heading 1"
  ) %>%
  body_add_flextable(ft2) %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_landscape() %>%

  # --- Table S2 ---
  body_add_par(
    "Table S2. Trial characteristics by anomalous versus non-anomalous HbA1c ρ, parametric 3-sigma definition",
    style = "heading 1"
  ) %>%
  body_add_flextable(ftS2) %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_landscape() %>%

  # --- Table S3 ---
  body_add_par(
    "Table S3. Univariate Firth-penalized logistic regression of anomalous HbA1c ρ, parametric 3-sigma definition",
    style = "heading 1"
  ) %>%
  body_add_flextable(ftS3) %>%
  body_add_par("", style = "Normal") %>%
  body_end_section_landscape()

print(doc, target = output_file)

cat("\n\u2713 Word document successfully written to:\n", output_file, "\n")
