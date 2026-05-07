# Visible Value: Performance, Popularity & Pay in MLS
# Full analysis — charts, models, robustness, verification

set.seed(123)
options(contrasts = c("contr.treatment", "contr.poly"))


# Packages

packs <- c(
  "tidyverse", "fixest", "modelsummary", "readr", "clubSandwich",
  "broom", "marginaleffects", "sandwich", "glue", "knitr",
  "ggrepel", "scales", "nnet"
)
to_install <- setdiff(packs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(packs, library, character.only = TRUE))
  

# Data

# salary columns come in as comma-formatted strings — fix that on load
df <- read_csv(
  "/Users/umairs/Downloads/diss2/2024mlsStatsFinal.csv",
  show_col_types = FALSE
) %>%
  mutate(
    salary2025   = parse_number(as.character(salary2025)),
    salaryChange = parse_number(as.character(salaryChange))
  )


# Cleaning & feature construction

df <- df %>%
  mutate(
    # popularity_z recomputed every run so it never drifts from the raw column
    popularity_z = (popularity - mean(popularity, na.rm = TRUE)) /
      sd(popularity, na.rm = TRUE),
    
    # players with no designation filed are Regular
    rosterDesignation = if_else(
      is.na(rosterDesignation) | trimws(rosterDesignation) == "",
      "Regular", rosterDesignation
    ),
    
    salary2024 = salary2025 - salaryChange,
    salary2024 = ifelse(salary2024 > 0, salary2024, NA_real_),
    log_w2025  = ifelse(salary2025 > 0, log(salary2025), NA_real_),
    log_w2024  = ifelse(salary2024 > 0, log(salary2024), NA_real_),
    dlog_wage  = log_w2025 - log_w2024,
    
    club              = factor(club),
    position          = factor(position),
    rosterDesignation = factor(rosterDesignation)
  ) %>%
  mutate(
    position          = forcats::fct_relevel(position, "D"),
    rosterDesignation = forcats::fct_relevel(rosterDesignation, "Regular"),
    # Gabriel Pec is filed under his full legal name in the data
    player_label = if_else(player == "Gabriel Chaves", "Gabriel Pec", player)
  )


# Shared helpers

theme_set(theme_minimal(base_size = 12))

perf_vec <- c("attack_score", "creation_score", "progression_score",
              "dribbling_score", "defense_score")

controls      <- "age + I(age^2) + mins + position + rosterDesignation"
controls_nord <- "age + I(age^2) + mins + position"   # drops designation (for chart 2)

# model formulas
f_levels_overall   <- as.formula(paste("log_w2025 ~ overall_score + popularity_z +", controls, "| club"))
f_levels_overall_x <- as.formula(paste("log_w2025 ~ overall_score * popularity_z +", controls, "| club"))
f_levels_vector    <- as.formula(paste("log_w2025 ~", paste(c(perf_vec, "popularity_z", controls),    collapse = " + "), "| club"))
f_growth_overall   <- as.formula(paste("dlog_wage  ~ overall_score + popularity_z +", controls, "| club"))
f_growth_vector    <- as.formula(paste("dlog_wage  ~", paste(c(perf_vec, "popularity_z", controls),   collapse = " + "), "| club"))
f_before           <- as.formula(paste("log_w2025 ~", paste(c(perf_vec, "popularity_z", controls_nord), collapse = " + "), "| club"))

# shared rename map for modelsummary calls
coef_names <- c(
  "overall_score"              = "Performance (overall z)",
  "popularity_z"               = "Popularity (z)",
  "overall_score:popularity_z" = "Perf × Pop",
  "attack_score"               = "Attack (z)",
  "creation_score"             = "Creation (z)",
  "progression_score"          = "Progression (z)",
  "dribbling_score"            = "Dribbling (z)",
  "defense_score"              = "Defense (z)"
)

gof_map <- tibble::tribble(
  ~raw,            ~clean,         ~fmt,
  "nobs",          "Observations",  0,
  "rmse",          "RMSE",          2,
  "r.squared",     "R²",            2,
  "adj.r.squared", "Adj. R²",       2
)


# Core models

m1       <- feols(f_levels_overall,   data = df, cluster = ~ club)
m2       <- feols(f_levels_overall_x, data = df, cluster = ~ club)
m3       <- feols(f_levels_vector,    data = df, cluster = ~ club)   # main model
m4       <- feols(f_growth_overall,   data = df, cluster = ~ club)
m5       <- feols(f_growth_vector,    data = df, cluster = ~ club)
m_before <- feols(f_before,           data = df, cluster = ~ club)   # no designation control


# Chart 1 — What Explains MLS Wages?
# Shapley run twice: with and without roster as a group

compute_shapley <- function(data, outcome, groups) {
  n      <- length(groups)
  combos <- expand.grid(rep(list(c(FALSE, TRUE)), n))
  names(combos) <- names(groups)
  
  r2 <- numeric(nrow(combos))
  for (i in seq_len(nrow(combos))) {
    idx <- which(unlist(combos[i, ]))
    if (length(idx) == 0) { r2[i] <- 0; next }
    vars <- unlist(groups[idx])
    fit  <- lm(as.formula(paste(outcome, "~", paste(vars, collapse = " + "))), data = data)
    r2[i] <- summary(fit)$r.squared
  }
  
  shap <- numeric(n)
  for (j in seq_len(n)) {
    for (i in seq_len(nrow(combos))) {
      if (combos[i, j]) next
      s      <- sum(unlist(combos[i, ]))
      w      <- factorial(s) * factorial(n - s - 1) / factorial(n)
      with_j <- combos[i, ]; with_j[, j] <- TRUE
      k      <- which(apply(combos, 1, function(r) all(r == unlist(with_j))))
      shap[j] <- shap[j] + w * (r2[k] - r2[i])
    }
  }
  
  tibble(group = names(groups), shapley = shap, share = shap / sum(shap))
}

df_sha <- df %>% filter(!is.na(log_w2025))

sha_without <- compute_shapley(df_sha, "log_w2025", list(
  Performance  = perf_vec,
  Popularity   = "popularity_z",
  Demographics = c("age", "I(age^2)", "mins", "position"),
  Club         = "club"
)) %>% mutate(model = "Without roster controls")

sha_with <- compute_shapley(df_sha, "log_w2025", list(
  Roster       = "rosterDesignation",
  Performance  = perf_vec,
  Popularity   = "popularity_z",
  Demographics = c("age", "I(age^2)", "mins", "position"),
  Club         = "club"
)) %>% mutate(model = "With roster controls")

sha_combined <- bind_rows(sha_without, sha_with) %>%
  mutate(
    group = factor(group, levels = c("Roster", "Performance", "Popularity",
                                     "Demographics", "Club")),
    model = factor(model, levels = c("Without roster controls",
                                     "With roster controls"))
  )

p_shapley <- ggplot(sha_combined, aes(x = share, y = group, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    aes(label = scales::percent(share, accuracy = 0.1)),
    position = position_dodge(width = 0.7), hjust = -0.1, size = 3.5
  ) +
  scale_fill_manual(values = c(
    "Without roster controls" = "firebrick",
    "With roster controls"    = "steelblue"
  )) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "What Explains MLS Wages?",
    subtitle = "Roster designation accounts for more than performance, popularity, and club combined",
    x = "Share of salary variance explained", y = NULL, fill = NULL,
    caption  = "SHAPLEY DECOMPOSITION OF LOG SALARY VARIANCE — AVERAGED ACROSS ALL VARIABLE ORDERINGS"
  ) +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

ggsave("chart1_shapley.png", p_shapley, width = 7.5, height = 5, dpi = 300)
print(p_shapley)


# Chart 2 — What the Slot Does to Pay
# Same coefficients, before and after designation enters the model

terms_ba  <- c("attack_score", "creation_score", "progression_score",
               "dribbling_score", "defense_score", "popularity_z")
labels_ba <- c("Attack", "Creation", "Progression", "Dribbling", "Defense", "Popularity")
term_map  <- setNames(labels_ba, terms_ba)

ba_df <- bind_rows(
  broom::tidy(m_before) %>% filter(term %in% terms_ba) %>% mutate(model = "Before roster controls"),
  broom::tidy(m3)       %>% filter(term %in% terms_ba) %>% mutate(model = "After roster controls")
) %>%
  mutate(
    term_label = factor(recode(term, !!!term_map), levels = labels_ba),
    term_label = forcats::fct_reorder(term_label, estimate, .fun = mean),
    model      = factor(model, levels = c("Before roster controls", "After roster controls"))
  )

p_before_after <- ggplot(ba_df, aes(x = estimate, y = term_label, colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_point(size = 3, position = position_dodge(width = 0.55)) +
  geom_errorbarh(
    aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error),
    height = 0, position = position_dodge(width = 0.55)
  ) +
  scale_colour_manual(values = c(
    "Before roster controls" = "firebrick",
    "After roster controls"  = "steelblue"
  )) +
  labs(
    title    = "What the Slot Does to Pay",
    subtitle = "The same output earns less once the league decides the tier",
    x = "Wage return per 1 SD increase (log points)", y = NULL, colour = NULL,
    caption  = "WAGE RETURN PER STANDARD DEVIATION IMPROVEMENT — BEFORE (RED) AND AFTER (BLUE) CONTROLLING FOR ROSTER DESIGNATION"
  ) +
  theme(legend.position = "bottom")

ggsave("chart2_before_after.png", p_before_after, width = 7, height = 5, dpi = 300)
print(p_before_after)


# Chart 3 — What's Done on the Pitch
# Raw scatter per dimension — slope tells the story

perf_labels <- c(
  attack_score      = "Attack",
  creation_score    = "Creation",
  progression_score = "Progression",
  dribbling_score   = "Dribbling",
  defense_score     = "Defense"
)

perf_long <- df %>%
  pivot_longer(cols = all_of(perf_vec), names_to = "component", values_to = "score") %>%
  mutate(component = factor(
    recode(component, !!!perf_labels),
    levels = c("Attack", "Progression", "Creation", "Dribbling", "Defense")
  ))

p_pitch <- ggplot(perf_long, aes(x = score, y = log_w2025)) +
  geom_point(alpha = 0.2, size = 0.9, color = "grey40") +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue", linewidth = 0.9) +
  facet_wrap(~ component, nrow = 1) +
  labs(
    title    = "What's Done on the Pitch",
    subtitle = "Attack and progression earn wage premiums. Defense earns nothing.",
    x = "Contribution (z-score)", y = "Log salary (2025)",
    caption  = "WAGE RETURN PER STANDARD DEVIATION, BY PERFORMANCE DIMENSION — CONTROLLING FOR CLUB FIXED EFFECTS"
  )

ggsave("chart3_defense.png", p_pitch, width = 11, height = 4.5, dpi = 300)
print(p_pitch)


# Chart 4 — Performance vs Pay
# Residuals colour each dot; labels are hand-picked, not algorithmic

simple_fit      <- lm(log_w2025 ~ overall_score, data = df)
df$resid_simple <- residuals(simple_fit)

labeled_players <- c(
  "Sergio Busquets", "Lionel Messi", "Jordi Alba", "Lorenzo Insigne", "Emil Forsberg",
  "Christian Benteke", "Aleksei Miranchuk", "Hany Mukhtar", "Ryan Gauld", "Carles Gil",
  "Riqui Puig", "Walker Zimmerman", "Hugo Cuypers", "Léo Chú", "Evander", "Luciano Acosta",
  "Olivier Giroud", "Sam Surridge", "Dejan Joveljić", "Jonathan Rodríguez", "Joseph Paintsil",
  "Federico Bernardeschi", "Denis Bouanga", "Ezequiel Ponce", "Jordan Morris", "Petar Musa",
  "Albert Rusnák", "Diego Rossi", "Luis Muriel", "Gabriel Chaves",
  "Giacomo Vrioni", "Mikael Uhre", "Kevin Cabral", "Cristian Arango", "Dániel Gazdag",
  "Eduard Löwen", "Marco Reus",
  "Alhassan Yusuf", "David Martínez", "Julián Fernández", "Santiago Moreno", "Luca Orellano",
  "Joseph Rosales", "Maximilian Arfsten", "Edwin Mosquera", "Lassi Lappalainen",
  "Mauricio Cuevas", "Iuri Tavares", "Malachi Jones", "Stephen Afrifa"
)

df <- df %>%
  mutate(scatter_label = if_else(player %in% labeled_players, player_label, NA_character_))

resid_lim <- ceiling(max(abs(df$resid_simple), na.rm = TRUE) * 10) / 10

p_scatter <- ggplot(df, aes(x = overall_score, y = log_w2025)) +
  geom_point(aes(color = resid_simple), alpha = 0.6, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_text_repel(
    aes(label = scatter_label), size = 2.8, max.overlaps = 50,
    na.rm = TRUE, segment.color = "grey50", segment.size = 0.3
  ) +
  scale_color_gradient2(
    low = "steelblue", mid = "grey80", high = "firebrick",
    midpoint = 0, limits = c(-resid_lim, resid_lim), name = "Residual"
  ) +
  labs(
    title    = "Performance vs Pay",
    subtitle = "Red = paid above what their output predicts | Blue = paid below",
    x = "Overall performance (z-score)", y = "Log salary (2025)",
    caption  = "OVERALL PERFORMANCE VS SALARY. RED = PAID ABOVE WHAT THE FULL MODEL PREDICTS, BLUE = PAID BELOW"
  ) +
  theme(legend.position = "right")

ggsave("chart5_scatter.png", p_scatter, width = 9, height = 6.5, dpi = 300)
print(p_scatter)


# Chart 5 — What Gets the Slot?
# Multinomial logit: Regular / TAM / DP across performance × visibility

df_logit <- df %>%
  filter(!is.na(log_w2025)) %>%
  mutate(
    tier = case_when(
      rosterDesignation == "Designated Player" ~ "Designated Player",
      rosterDesignation == "TAM Player"        ~ "TAM Player",
      TRUE                                     ~ "Regular"
    ),
    tier = factor(tier, levels = c("Regular", "TAM Player", "Designated Player"))
  )

m_multi <- nnet::multinom(
  tier ~ overall_score * popularity_z + age + I(age^2) + position,
  data = df_logit, trace = FALSE
)

med_age <- median(df_logit$age, na.rm = TRUE)

grid_multi <- data.frame(
  overall_score = c(0, 2, 0, 2),
  popularity_z  = c(0, 0, 2, 2),
  age           = med_age,
  position      = factor("M", levels = levels(df$position)),
  row_label     = factor(
    c("Avg player", "Top performer, low visibility",
      "Avg performer, high visibility", "Top performer, high visibility"),
    levels = c("Top performer, high visibility", "Avg performer, high visibility",
               "Top performer, low visibility", "Avg player")
  )
)

probs_long <- predict(m_multi, newdata = grid_multi, type = "probs") %>%
  as.data.frame() %>%
  bind_cols(row_label = grid_multi$row_label) %>%
  pivot_longer(
    cols = c("Regular", "TAM Player", "Designated Player"),
    names_to = "tier", values_to = "prob"
  ) %>%
  mutate(tier = factor(tier, levels = c("Regular", "TAM Player", "Designated Player")))

p_logit <- ggplot(probs_long, aes(x = prob, y = row_label, fill = tier)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(
    aes(label = scales::percent(prob, accuracy = 1)),
    position = position_stack(vjust = 0.5), size = 3.2, color = "white"
  ) +
  scale_fill_manual(values = c(
    "Regular"           = "#C8B89A",
    "TAM Player"        = "#B05C38",
    "Designated Player" = "#7A1515"
  )) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "What Gets the Slot?",
    subtitle = "How performance and visibility combine to determine earning potential",
    x = "Predicted probability", y = NULL, fill = NULL,
    caption  = "PREDICTED PROBABILITY OF A DESIGNATED PLAYER CONTRACT — BY PERFORMANCE TIER AND PUBLIC VISIBILITY"
  ) +
  theme(legend.position = "top")

ggsave("chart4_archetypes.png", p_logit, width = 7, height = 5, dpi = 300)
print(p_logit)


# Chart 6 — Club Wage Premium
# Fixed effects from m3, recentred around the league mean

club_fe <- fixef(m3)$club %>%
  as.data.frame() %>%
  rownames_to_column("club") %>%
  rename(fe = ".") %>%
  mutate(
    fe_centered = fe - mean(fe),
    club        = forcats::fct_reorder(club, fe_centered)
  )

p_fe <- ggplot(club_fe, aes(x = fe_centered, y = club)) +
  geom_col(aes(fill = fe_centered > 0), show.legend = FALSE) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title    = "Club Wage Premium — Relative to League Average",
    subtitle = "Controlling for player profile: performance, popularity, position, age, minutes",
    x = "Wage premium relative to league average", y = NULL,
    caption  = "CLUB WAGE PREMIUM OR DISCOUNT, RELATIVE TO LEAGUE AVERAGE — CONTROLLING FOR SQUAD COMPOSITION, NOT TOTAL PAYROLL"
  )

ggsave("chart6_club_fe.png", p_fe, width = 7, height = 6, dpi = 300)
print(p_fe)


# Regression tables

modelsummary(
  list(
    "Levels: overall"       = m1,
    "Levels: overall × pop" = m2,
    "Levels: perf vector"   = m3,
    "Δlog wage: overall"    = m4,
    "Δlog wage: perf vector"= m5
  ),
  stars       = TRUE,
  fmt         = 2,
  coef_rename = coef_names,
  gof_map     = gof_map,
  gof_omit    = "IC|Log|Within|Std.Errors",
  output      = "markdown"
)
cat("\nOmitted: Defender (position), Regular (designation)\n")


# Wald tests — are the performance components jointly significant?
wt_m3 <- wald(m3, keep = perf_vec)
wt_m5 <- wald(m5, keep = perf_vec)
tibble(
  Model   = c("M3: Levels, perf vector", "M5: Growth, perf vector"),
  F_stat  = c(unname(wt_m3$stat), unname(wt_m5$stat)),
  df1     = c(unname(wt_m3$df1),  unname(wt_m5$df1)),
  p_value = c(unname(wt_m3$p),    unname(wt_m5$p))
) %>% print()


# Robustness

# trim Messi outlier in popularity
p99_val <- quantile(df$popularity_z, 0.99, na.rm = TRUE)
m_trim  <- feols(f_levels_overall, data = filter(df, popularity_z <= p99_val), cluster = ~ club)

# rank-inverse normal transform for popularity
ranks <- rank(df$popularity, na.last = "keep", ties.method = "average")
df    <- df %>% mutate(popularity_rint = qnorm((ranks - 0.5) / sum(!is.na(ranks))))
m_rint <- feols(
  as.formula(paste("log_w2025 ~ overall_score + popularity_rint +", controls, "| club")),
  data = df, cluster = ~ club
)

# weighted by minutes played
m_wls  <- feols(f_levels_overall, data = df, weights = ~ mins, cluster = ~ club)

# two-way cluster
m_2way <- feols(f_levels_overall, data = df, cluster = ~ club + position)

# CR2 sandwich
V_CR2  <- clubSandwich::vcovCR(m1, cluster = df$club, type = "CR2")

modelsummary(
  list(
    "Baseline"          = m1,
    "CR2"               = m1,
    "WLS (mins)"        = m_wls,
    "2-way cluster"     = m_2way,
    "Trim top 1% pop"   = m_trim,
    "RINT popularity"   = m_rint
  ),
  vcov = list(
    sandwich::vcovCL(m1,     cluster = ~ club),
    V_CR2,
    sandwich::vcovCL(m_wls,  cluster = ~ club),
    sandwich::vcovCL(m_2way, cluster = ~ club + position),
    sandwich::vcovCL(m_trim, cluster = ~ club),
    sandwich::vcovCL(m_rint, cluster = ~ club)
  ),
  fmt      = 2,
  coef_map = c(
    "overall_score"   = "Performance (overall z)",
    "popularity_z"    = "Popularity (z)",
    "popularity_rint" = "Popularity (RINT)"
  ),
  gof_map  = gof_map,
  gof_omit = "IC|Log|Within|Std.Errors",
  output   = "markdown"
)


# Marginal effect of popularity across performance levels (from interaction model)
V_m2   <- sandwich::vcovCL(m2, cluster = ~ club)
p_marg <- plot_slopes(
  m2, variables = "popularity_z",
  condition = list(overall_score = seq(-2, 2, by = 0.5)),
  vcov = V_m2
) + labs(
  title = "Marginal effect of popularity, by performance level",
  x = "Overall performance (z)", y = "∂ log salary / ∂ popularity"
)
print(p_marg)


# Heterogeneity by designation

sub_controls <- "age + I(age^2) + mins + position"

run_het <- function(d, label) {
  models <- list(
    "M1"          = feols(as.formula(paste("log_w2025 ~ overall_score + popularity_z +", sub_controls, "| club")), data = d, cluster = ~ club),
    "M2 (× pop)"  = feols(as.formula(paste("log_w2025 ~ overall_score * popularity_z +", sub_controls, "| club")), data = d, cluster = ~ club),
    "M3 (vector)" = feols(as.formula(paste("log_w2025 ~", paste(c(perf_vec, "popularity_z", sub_controls), collapse = " + "), "| club")), data = d, cluster = ~ club),
    "M4 (growth)" = feols(as.formula(paste("dlog_wage  ~ overall_score + popularity_z +", sub_controls, "| club")), data = d, cluster = ~ club),
    "M5 (growth)" = feols(as.formula(paste("dlog_wage  ~", paste(c(perf_vec, "popularity_z", sub_controls), collapse = " + "), "| club")), data = d, cluster = ~ club)
  )
  cat(sprintf("\n\n--- %s ---\n\n", label))
  modelsummary(models, stars = TRUE, fmt = 2, coef_rename = coef_names,
               gof_map = gof_map, gof_omit = "IC|Log|Within|Std.Errors", output = "markdown")
}

for (rd in levels(df$rosterDesignation)) {
  dsub <- filter(df, rosterDesignation == rd)
  if (nrow(dsub) >= 10) run_het(dsub, rd)
}


# Diagnostics (M3)

r2_m3   <- fitstat(m3, "r2")
rmse_m3 <- fitstat(m3, "rmse")

ggplot(df, aes(x = predict(m3), y = log_w2025)) +
  geom_point(alpha = 0.5, size = 1.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = sprintf("Predicted vs actual — M3 (R²=%.2f, RMSE=%.2f)", r2_m3, rmse_m3),
       x = "Predicted", y = "Actual")

df_res <- tibble(.fitted = as.numeric(fitted(m3)),
                 .resid  = as.numeric(resid(m3))) %>%
  mutate(std_resid = scale(.resid)[, 1])

ggplot(df_res, aes(x = .fitted, y = std_resid)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_smooth(method = "loess", se = FALSE) +
  geom_hline(yintercept = c(-3, 0, 3), linetype = c("dotted", "solid", "dotted")) +
  labs(title = "Residuals vs fitted — M3", x = "Fitted", y = "Std. residual")

ggplot(df_res, aes(sample = std_resid)) +
  stat_qq(alpha = 0.6) + stat_qq_line() +
  labs(title = "Q-Q plot — M3", x = "Theoretical", y = "Sample")


# Economic magnitudes

mean_sal <- mean(df$salary2025, na.rm = TRUE)

broom::tidy(m3) %>%
  filter(term %in% c(perf_vec, "popularity_z")) %>%
  mutate(
    label          = recode(term, !!!coef_names),
    wage_return_pct = sprintf("%.1f%%", estimate * 100),
    dollar_impact   = scales::dollar(mean_sal * estimate, accuracy = 1)
  ) %>%
  select(label, wage_return_pct, dollar_impact) %>%
  knitr::kable(caption = sprintf(
    "Wage return per 1 SD — evaluated at mean salary %s (M3)",
    scales::dollar(mean_sal, accuracy = 1)
  ))


# Backend verification
# Numbers behind every claim on the site

cat("\n\n--- Shapley (with roster controls) ---\n")
print(sha_with %>% mutate(pct = scales::percent(share, accuracy = 0.1)))

cat("\n\n--- Before/after coefficients (log points) ---\n")
ba_df %>%
  select(term, model, estimate, std.error) %>%
  mutate(estimate = round(estimate, 3), se = round(std.error, 3)) %>%
  pivot_wider(names_from = model, values_from = c(estimate, se)) %>%
  print()

cat("\n\n--- M3 coefficients ---\n")
broom::tidy(m3) %>%
  filter(term %in% c(perf_vec, "popularity_z")) %>%
  mutate(est = round(estimate, 3), se = round(std.error, 3), p = round(p.value, 3)) %>%
  select(term, est, se, p) %>% print()

cat("\n\n--- Logit predicted probabilities ---\n")
predict(m_multi, newdata = grid_multi, type = "probs") %>%
  as.data.frame() %>%
  mutate(scenario = as.character(grid_multi$row_label)) %>%
  mutate(across(where(is.numeric), ~ scales::percent(round(.x, 3), accuracy = 0.1))) %>%
  print()

cat("\n\n--- Logit: age peak and position effects ---\n")
coef_multi  <- coef(m_multi)
age_peak_dp <- -coef_multi["Designated Player", "age"] /
  (2 * coef_multi["Designated Player", "I(age^2)"])
cat(sprintf("DP probability peaks at age: %.1f\n\n", age_peak_dp))
print(round(coef_multi, 3))

cat("\n\n--- Key player lookup ---\n")
df %>%
  filter(player %in% c("Lionel Messi", "Sergio Busquets", "Gabriel Chaves",
                       "Stephen Afrifa", "Maximilian Arfsten", "Riqui Puig",
                       "Luciano Acosta", "Carles Gil")) %>%
  mutate(name = if_else(player == "Gabriel Chaves", "Gabriel Pec", player)) %>%
  select(name, salary2025, overall_score, popularity_z, rosterDesignation) %>%
  arrange(desc(salary2025)) %>%
  print()