#  MLS wages ~ performance & popularity
#  v2: plots pop up inline; four new additions


# Repro + options
set.seed(123)
options(contrasts = c("contr.treatment","contr.poly"))
options(modelsummary_factory_latex = "kableExtra")
options(modelsummary_format_numeric_latex = "plain")

## Packages
packs <- c(
  "tidyverse","fixest","modelsummary","readxl","clubSandwich",
  "broom","marginaleffects","sandwich","glue","knitr","ggrepel","patchwork"
)
to_install <- setdiff(packs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(packs, library, character.only = TRUE))

# data

df <- read_excel("/Users/umairs/Desktop/Work/Dissertation/DissertationFiles/2024mlsStatsFinal_clean.xlsx")

needed <- c(
  "player","club","rosterDesignation","position","nationality","age","mins",
  "salaryChange","salary2025",
  "attack_score","creation_score","progression_score","dribbling_score","defense_score","overall_score",
  "popularity","popularity_z"
)
stopifnot(all(needed %in% names(df)))

# clean

df <- df %>%
  mutate(
    rosterDesignation = if_else(is.na(rosterDesignation) | trimws(rosterDesignation)=="", "Regular", rosterDesignation),
    salary2024 = salary2025 - salaryChange,
    salary2024 = ifelse(salary2024 > 0, salary2024, NA_real_),
    log_w2025  = ifelse(salary2025 > 0, log(salary2025), NA_real_),
    log_w2024  = ifelse(salary2024 > 0, log(salary2024), NA_real_),
    dlog_wage  = log_w2025 - log_w2024,
    club = factor(club),
    position = factor(position),
    rosterDesignation = factor(rosterDesignation)
  ) %>%
  mutate(
    position = forcats::fct_relevel(position, "D"),
    rosterDesignation = forcats::fct_relevel(rosterDesignation, "Regular")
  )

# plot helpers

theme_set(theme_minimal(base_size = 12))
cap <- function(x, upper) pmin(x, upper)

YRANGE  <- range(df$log_w2025, na.rm = TRUE)
P99_POP <- quantile(df$popularity_z, 0.99, na.rm = TRUE)

# summary

datasummary_skim(
  df %>% select(salary2025, dlog_wage, overall_score, popularity_z, age, mins),
  output = "markdown"
)

# salary x designation
df <- df %>%
  mutate(rosterDesignation_plot =
           forcats::fct_reorder(rosterDesignation, log_w2025, .fun = median, na.rm = TRUE))

p_box <- ggplot(df, aes(y = rosterDesignation_plot, x = log_w2025)) +
  geom_boxplot(outlier.alpha = .5, width = .6) +
  labs(title = "Log Salary (2025) by Roster Designation",
       y = "Roster Designation", x = "Log Salary (2025)")
print(p_box)

# exploratory

perf_vec <- c("attack_score","creation_score","progression_score","dribbling_score","defense_score")

# Build N lookup from actual data counts, then map by label
rd_counts <- table(as.character(df$rosterDesignation))
add_n <- function(f){ paste0(f, " (N=", rd_counts[f], ")") }
df <- df %>% mutate(rosterDesignation_lab = forcats::fct_relabel(rosterDesignation, add_n))

# Components vs log salary
perf_long <- df %>% pivot_longer(cols = all_of(perf_vec), names_to = "component", values_to = "score")
p_comp <- ggplot(perf_long, aes(x = score, y = log_w2025)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ component, scales = "fixed") +
  coord_cartesian(ylim = YRANGE) +
  labs(title = "Performance Components vs Log Salary (2025)",
       x = "Component Score (z)", y = "Log Salary (2025)")
print(p_comp)

# Popularity vs salary by designation
p_pop <- ggplot(df, aes(x = cap(popularity_z, P99_POP), y = log_w2025)) +
  geom_point(alpha = 0.45, size = 1) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ rosterDesignation_lab, ncol = 3) +
  labs(title = "Popularity vs Log Salary (2025) by Roster Designation",
       x = "Popularity (z, clamped at 99th pct for plotting)", y = "Log Salary (2025)")
print(p_pop)

# Player-level salary vs performance scatter with labels
# Highlights the most interesting outliers — overperformers & underperformers relative to pay
simple_fit        <- lm(log_w2025 ~ overall_score, data = df)
df$pred_simple    <- fitted(simple_fit)
df$resid_simple   <- residuals(simple_fit)
df$outlier_label  <- with(df, case_when(
  resid_simple > quantile(resid_simple, 0.92, na.rm = TRUE) ~ player,
  resid_simple < quantile(resid_simple, 0.08, na.rm = TRUE) ~ player,
  TRUE ~ NA_character_
))

p_scatter <- ggplot(df, aes(x = overall_score, y = log_w2025)) +
  geom_point(aes(color = resid_simple), alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_text_repel(aes(label = outlier_label), size = 2.8, max.overlaps = 15, na.rm = TRUE) +
  scale_color_gradient2(low = "steelblue", mid = "grey80", high = "firebrick",
                        midpoint = 0, name = "Residual") +
  labs(title = "Overall Performance vs Log Salary (2025)",
       subtitle = "Red = overpaid relative to performance | Blue = underpaid",
       x = "Overall Performance Score (z)", y = "Log Salary (2025)") +
  theme(legend.position = "right")
print(p_scatter)

# corr matrix

vars <- df %>% select(
  overall_score, attack_score, creation_score, progression_score,
  dribbling_score, defense_score, popularity_z, age, mins
)
C <- cor(vars, use = "pairwise.complete.obs", method = "pearson")
knitr::kable(round(C, 2))

# specifications
controls <- "age + I(age^2) + mins + position + rosterDesignation"
f_levels_overall   <- as.formula(paste("log_w2025 ~ overall_score + popularity_z +", controls, "| club"))
f_levels_overall_x <- as.formula(paste("log_w2025 ~ overall_score * popularity_z +", controls, "| club"))
f_levels_vector    <- as.formula(paste("log_w2025 ~", paste(c(perf_vec,"popularity_z",controls), collapse=" + "), "| club"))
f_growth_overall   <- as.formula(paste("dlog_wage  ~ overall_score + popularity_z +", controls, "| club"))
f_growth_vector    <- as.formula(paste("dlog_wage  ~", paste(c(perf_vec,"popularity_z",controls), collapse=" + "), "| club"))

# estimation of core models

m1 <- feols(f_levels_overall,   data = df, cluster = ~ club)
m2 <- feols(f_levels_overall_x, data = df, cluster = ~ club)
m3 <- feols(f_levels_vector,    data = df, cluster = ~ club)
m4 <- feols(f_growth_overall,   data = df, cluster = ~ club)
m5 <- feols(f_growth_vector,    data = df, cluster = ~ club)

# baseliens regressions

gof_map2 <- tibble::tribble(
  ~raw,            ~clean,        ~fmt,
  "nobs",          "Observations", 0,
  "rmse",          "RMSE",         2,
  "r.squared",     "R²",           2,
  "adj.r.squared", "Adj. R²",      2
)

modelsummary(
  list(
    "Levels: overall"        = m1,
    "Levels: overall × pop"  = m2,
    "Levels: perf vector"    = m3,
    "Δlog wage: overall"     = m4,
    "Δlog wage: perf vector" = m5
  ),
  stars = TRUE,
  fmt = 2,
  coef_rename = c(
    "overall_score"              = "Performance (overall z)",
    "popularity_z"               = "Popularity (z)",
    "overall_score:popularity_z" = "Perf × Pop",
    "attack_score"      = "Attack (z)",
    "creation_score"    = "Creation (z)",
    "progression_score" = "Progression (z)",
    "dribbling_score"   = "Dribbling (z)",
    "defense_score"     = "Defense (z)"
  ),
  gof_map = gof_map2,
  gof_omit = "IC|Log|Within|Std.Errors",
  output = "markdown"
)
cat("\nNote: Omitted categories are Defender (position) and Regular (roster designation).\n\n")

# wald tests

wt_m3 <- wald(m3, keep = perf_vec)
wt_m5 <- wald(m5, keep = perf_vec)
wald_tbl <- tibble::tibble(
  Model  = c("M3: Levels, perf vector", "M5: Growth, perf vector"),
  F_stat = c(unname(wt_m3$stat), unname(wt_m5$stat)),
  df1    = c(unname(wt_m3$df1),  unname(wt_m5$df1)),
  df2    = c(unname(wt_m3$df2),  unname(wt_m5$df2)),
  p_value= c(unname(wt_m3$p),    unname(wt_m5$p))
)
print(wald_tbl)

# robustness

p99_val <- quantile(df$popularity_z, 0.99, na.rm = TRUE)
df_trim <- dplyr::filter(df, popularity_z <= p99_val)
m_trim  <- feols(f_levels_overall, data = df_trim, cluster = ~ club)

ranks <- rank(df$popularity, na.last = "keep", ties.method = "average")
df    <- df %>% mutate(popularity_rint = qnorm((ranks - 0.5) / sum(!is.na(ranks))))
f_levels_overall_rint <- as.formula(paste("log_w2025 ~ overall_score + popularity_rint +", controls, "| club"))
m_rint <- feols(f_levels_overall_rint, data = df, cluster = ~ club)

m_wls  <- feols(f_levels_overall, data = df, weights = ~ mins, cluster = ~ club)
m_2way <- feols(f_levels_overall, data = df, cluster = ~ club + position)
V_CR2_m1 <- clubSandwich::vcovCR(m1, cluster = df$club, type = "CR2")

robust_models <- list(
  "Baseline (club cluster)"  = m1,
  "CR2 (club)"               = m1,
  "WLS (weights = mins)"     = m_wls,
  "2-way cluster (club+pos)" = m_2way,
  "Trim top 1% pop_z"        = m_trim,
  "RINT popularity"          = m_rint
)
robust_vcov <- list(
  sandwich::vcovCL(m1, cluster = ~ club),
  V_CR2_m1,
  sandwich::vcovCL(m_wls, cluster = ~ club),
  sandwich::vcovCL(m_2way, cluster = ~ club + position),
  sandwich::vcovCL(m_trim, cluster = ~ club),
  sandwich::vcovCL(m_rint, cluster = ~ club)
)
modelsummary(
  robust_models,
  vcov = robust_vcov,
  fmt = 2,
  coef_map = c("overall_score" = "Performance (overall z)",
               "popularity_z"  = "Popularity (z)",
               "popularity_rint" = "Popularity (RINT)"),
  estimate  = "{estimate}",
  statistic = "({std.error})",
  gof_map   = tibble::tribble(
    ~raw,            ~clean,        ~fmt,
    "nobs",          "Observations", 0,
    "rmse",          "RMSE",         2,
    "r.squared",     "R²",           2,
    "adj.r.squared", "Adj. R²",      2
  ),
  gof_omit  = "IC|Log|Within|Std.Errors",
  output    = "markdown"
)

# interaction visual

V_club_m2 <- sandwich::vcovCL(m2, cluster = ~ club)
p_marg <- plot_slopes(
  m2,
  variables = "popularity_z",
  condition = list(overall_score = seq(-2, 2, by = 0.5)),
  vcov = V_club_m2
) + labs(title = "Marginal Effect of Popularity across Performance",
         x = "Overall performance (z)", y = "∂ log Salary / ∂ Popularity (z)")
print(p_marg)

# defense score deeper dive

p_def <- ggplot(df, aes(x = defense_score, y = log_w2025)) +
  geom_point(alpha = 0.4, size = 1.5, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
  facet_wrap(~ position) +
  labs(title = "Defense Score vs Log Salary — by Position",
       subtitle = "Is defensive contribution rewarded differently across positions?",
       x = "Defense Score (z)", y = "Log Salary (2025)")
print(p_def)

# heterogenity x designation
sub_controls <- "age + I(age^2) + mins + position"
f_sub_m1 <- as.formula(paste("log_w2025 ~ overall_score + popularity_z +", sub_controls, "| club"))
f_sub_m2 <- as.formula(paste("log_w2025 ~ overall_score * popularity_z +", sub_controls, "| club"))
f_sub_m3 <- as.formula(paste("log_w2025 ~", paste(c(perf_vec,"popularity_z",sub_controls), collapse=" + "), "| club"))
f_sub_m4 <- as.formula(paste("dlog_wage  ~ overall_score + popularity_z +", sub_controls, "| club"))
f_sub_m5 <- as.formula(paste("dlog_wage  ~", paste(c(perf_vec,"popularity_z",sub_controls), collapse=" + "), "| club"))

do_het_print <- function(d, rdes_label){
  models <- list(
    "M1: Levels overall"        = feols(f_sub_m1, data = d, cluster = ~ club),
    "M2: Levels overall×pop"    = feols(f_sub_m2, data = d, cluster = ~ club),
    "M3: Levels perf vector"    = feols(f_sub_m3, data = d, cluster = ~ club),
    "M4: Δlog wage overall"     = feols(f_sub_m4, data = d, cluster = ~ club),
    "M5: Δlog wage perf vector" = feols(f_sub_m5, data = d, cluster = ~ club)
  )
  cat("\n\n========== Heterogeneity:", rdes_label, "==========\n\n")
  modelsummary(
    models,
    stars = TRUE,
    fmt = 2,
    coef_rename = c(
      "overall_score"              = "Performance (overall z)",
      "popularity_z"               = "Popularity (z)",
      "overall_score:popularity_z" = "Perf × Pop",
      "attack_score"      = "Attack (z)",
      "creation_score"    = "Creation (z)",
      "progression_score" = "Progression (z)",
      "dribbling_score"   = "Dribbling (z)",
      "defense_score"     = "Defense (z)"
    ),
    gof_map = tibble::tribble(
      ~raw,            ~clean,        ~fmt,
      "nobs",          "Observations", 0,
      "rmse",          "RMSE",         2,
      "r.squared",     "R²",           2,
      "adj.r.squared", "Adj. R²",      2
    ),
    gof_omit = "IC|Log|Within|Std.Errors",
    output   = "markdown"
  )
  invisible(models)
}

for (rdes in levels(df$rosterDesignation)) {
  dsub <- dplyr::filter(df, rosterDesignation == rdes)
  if (nrow(dsub) >= 10) do_het_print(dsub, rdes)
}

# diagnostics
r2_m3   <- fitstat(m3, "r2"); rmse_m3 <- fitstat(m3, "rmse")
p_pred <- ggplot(df, aes(x = predict(m3), y = log_w2025)) +
  geom_point(alpha = 0.5, size = 1.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Predicted vs Actual Log Salary (M3)",
       subtitle = sprintf("R² = %.2f   RMSE = %.2f", r2_m3, rmse_m3),
       x = "Predicted", y = "Actual")
print(p_pred)

res_m3 <- resid(m3); fit_m3 <- fitted(m3)
df_res <- tibble(.fitted = as.numeric(fit_m3), .resid = as.numeric(res_m3)) %>%
  mutate(std_resid = scale(.resid)[,1])
p_resid <- ggplot(df_res, aes(x = .fitted, y = std_resid)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_smooth(method = "loess", se = FALSE) +
  geom_hline(yintercept = c(-3, 0, 3), linetype = c("dotted","solid","dotted")) +
  labs(title = "Standardized Residuals vs Fitted (M3)",
       x = "Fitted values", y = "Standardized residuals")
print(p_resid)

qq <- ggplot(df_res, aes(sample = std_resid)) +
  stat_qq(alpha = 0.6) + stat_qq_line() +
  labs(title = "Normal Q-Q Plot (M3 standardized residuals)", x = "Theoretical", y = "Sample")
print(qq)

# roster counts
df %>%
  count(rosterDesignation, name = "N") %>%
  mutate(Share = round(N / sum(N), 3)) %>%
  arrange(desc(N)) %>%
  print()

# economic magnitudes
mean_salary <- mean(df$salary2025, na.rm = TRUE)
coef_m3 <- broom::tidy(m3)

get_bite <- function(term_label, term_name) {
  est <- coef_m3 %>% dplyr::filter(term == term_name) %>% dplyr::pull(estimate)
  tibble::tibble(
    term = term_label,
    elasticity_pct = 100 * est,
    dollars = mean_salary * est
  )
}

bites <- dplyr::bind_rows(
  get_bite("Attack (z)",      "attack_score"),
  get_bite("Progression (z)", "progression_score"),
  get_bite("Popularity (z)",  "popularity_z")
) %>%
  dplyr::mutate(
    elasticity_pct = sprintf("%.1f%%", elasticity_pct),
    dollars = scales::dollar(dollars, accuracy = 1),
    mean_salary = scales::dollar(mean_salary, accuracy = 1)
  )

knitr::kable(bites,
             caption = "Economic magnitudes implied by M3: mean salary × coefficient (≈ $ change per 1 SD).")

# extended economic magnitudes — all five components
bites_full <- dplyr::bind_rows(
  get_bite("Attack (z)",      "attack_score"),
  get_bite("Creation (z)",    "creation_score"),
  get_bite("Progression (z)", "progression_score"),
  get_bite("Dribbling (z)",   "dribbling_score"),
  get_bite("Defense (z)",     "defense_score"),
  get_bite("Popularity (z)",  "popularity_z")
) %>%
  dplyr::mutate(
    elasticity_pct = sprintf("%.1f%%", elasticity_pct),
    dollars = scales::dollar(dollars, accuracy = 1),
    mean_salary = scales::dollar(mean_salary, accuracy = 1)
  )

knitr::kable(bites_full,
             caption = "Full economic magnitudes (M3): all performance components + popularity.")

# robustness specs
m3_base <- m3
top_id  <- which.max(df$popularity_z)
m3_excl <- feols(f_levels_vector, data = df[-top_id, ], cluster = ~ club)

df <- df %>% dplyr::mutate(mins_k = mins/1000)
f_levels_vector_nl <- as.formula(paste(
  "log_w2025 ~", paste(c(perf_vec,"popularity_z","age","I(age^2)","mins_k","I(mins_k^2)","position","rosterDesignation"), collapse=" + "), "| club"
))
m3_nl <- feols(f_levels_vector_nl, data = df, cluster = ~ club)

V_CR2_base <- clubSandwich::vcovCR(m3_base, cluster = df$club, type = "CR2")
V_CR2_excl <- clubSandwich::vcovCR(m3_excl,  cluster = df$club[-top_id], type = "CR2")
V_CR2_nl   <- clubSandwich::vcovCR(m3_nl,   cluster = df$club, type = "CR2")

if ("boottest" %in% ls("package:fixest")) {
  wb_attack <- boottest(m3, "attack_score", cluster = ~ club, B = 4999, type = "rademacher")
  wb_prog   <- boottest(m3, "progression_score", cluster = ~ club, B = 4999, type = "rademacher")
  wb_pop    <- boottest(m3, "popularity_z", cluster = ~ club, B = 4999, type = "rademacher")

  wild_tbl <- tibble::tibble(
    term     = c("Attack (z)", "Progression (z)", "Popularity (z)"),
    estimate = coef(m3)[c("attack_score","progression_score","popularity_z")],
    p_wild   = c(wb_attack$p.value, wb_prog$p.value, wb_pop$p.value)
  )
  knitr::kable(wild_tbl, digits = 3,
               caption = "Wild-cluster bootstrap p-values (fixest::boottest, Rademacher, B=4,999)")
}

# coeff plot robust
pack <- function(mod, lab) {
  broom::tidy(mod) %>% dplyr::filter(term %in% c("attack_score","progression_score","popularity_z")) %>%
    dplyr::mutate(model = lab)
}
df_coef <- dplyr::bind_rows(
  pack(m3_base, "Baseline"),
  pack(m3_excl, "Excl. top-pop"),
  pack(m3_nl,   "Minutes²")
) %>%
  dplyr::mutate(term = dplyr::recode(term,
    attack_score="Attack (z)", progression_score="Progression (z)", popularity_z="Popularity (z)")
  )

p_coef <- ggplot(df_coef, aes(x = term, y = estimate, color = model, group = model)) +
  geom_point(position = position_dodge(width = .5), size = 2.5) +
  geom_errorbar(aes(ymin = estimate - 1.96*std.error, ymax = estimate + 1.96*std.error),
                width = 0, position = position_dodge(width = .5)) +
  labs(title = "Key coefficients across robustness checks",
       x = NULL, y = "Estimate (log points)", color = "Spec") +
  theme_minimal(base_size = 11)
print(p_coef)

# Club-level fixed effects — who pays above/below market?
# Extracts the club FE from M3 and ranks clubs by how much they pay
# above or below what performance + popularity would predict.
club_fe <- fixef(m3)$club %>%
  as.data.frame() %>%
  rownames_to_column("club") %>%
  rename(fe = ".") %>%
  mutate(fe_centered = fe - mean(fe)) %>%   # recenter: 0 = league average
  arrange(desc(fe_centered)) %>%
  mutate(club = forcats::fct_reorder(club, fe_centered))

p_fe <- ggplot(club_fe, aes(x = fe_centered, y = club)) +
  geom_col(aes(fill = fe_centered > 0), show.legend = FALSE) +
  scale_fill_manual(values = c("steelblue","firebrick")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Club Fixed Effects (M3) — Relative to League Average",
       subtitle = "Red = pays above league avg | Blue = pays below league avg\n(conditional on performance, popularity, position, age, minutes)",
       x = "Fixed Effect relative to league mean (log points)", y = NULL)
print(p_fe)
