

library(readr)
library(dplyr)
library(glmmTMB)
library(ggplot2)
library(tidyr)
library(sjPlot)

library(broom.mixed)
library(tidytext)
library(grid)
library(patchwork)

## ==============================================================
## pre-setups
## ==============================================================

to_prop_if_needed <- function(x) {
  x <- as.numeric(x)
  if (max(x, na.rm = TRUE) > 1.5) x / 100 else x
}

logit_safe <- function(p) {
  p <- to_prop_if_needed(p)
  p2 <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  log(p2 / (1 - p2))
}

keep_existing <- function(x, nms) x[x %in% nms]

pred_df <- function(model, data, xvar, group_var, x_seq, group_levels) {
  ref_rel <- levels(data$MainReligion)[1]
  
  newdat <- expand.grid(
    x_seq,
    group_levels,
    ref_rel,
    stringsAsFactors = FALSE
  )
  names(newdat) <- c(xvar, group_var, "MainReligion")
  
  newdat[[group_var]]    <- factor(newdat[[group_var]], levels = group_levels)
  newdat$MainReligion    <- factor(newdat$MainReligion, levels = levels(data$MainReligion))
  newdat$log_GDP_sc      <- 0
  newdat$Median_age_sc   <- 0
  newdat$popdens_log_sc  <- 0
  
  pr <- predict(model, newdata = newdat, type = "response", se.fit = TRUE, re.form = NA)
  
  newdat$fit   <- pr$fit
  newdat$lower <- pr$fit - 1.96 * pr$se.fit
  newdat$upper <- pr$fit + 1.96 * pr$se.fit
  newdat
}

find_int_coef <- function(fe_names, term, group_var, lev) {
  cand <- paste0(term, ":", group_var, lev)
  if (cand %in% fe_names) return(cand)
  
  hits <- fe_names[grepl(paste0("^", term, ":"), fe_names) &
                     grepl(group_var, fe_names, fixed = TRUE) &
                     grepl(lev, fe_names, fixed = TRUE)]
  if (length(hits) >= 1) return(hits[1])
  NA_character_
}

get_slopes_by_group <- function(model, term, group_var, label_var = "Group") {
  fe <- fixef(model)$cond
  fe_names <- names(fe)
  V <- vcov(model)$cond
  
  mf <- model.frame(model)
  levs <- levels(mf[[group_var]])
  ref  <- levs[1]
  
  if (!term %in% fe_names) stop("Term not in fixed effects: ", term)
  beta_ref <- unname(fe[[term]])
  
  # ---- random slopes (Country) with CI (from broom.mixed) ----
  re_tidy <- broom.mixed::tidy(model, effects = "ran_vals", component = "cond", conf.int = TRUE)
  re_slopes <- re_tidy %>%
    filter(group == "Country", term == !!term) %>%
    transmute(
      Country   = as.character(level),
      ranef_est = estimate,
      ranef_lo  = conf.low,
      ranef_hi  = conf.high
    )
  
  # ---- country-specific total slopes per group level ----
  out <- lapply(levs, function(lev) {
    if (lev == ref) {
      beta_g <- beta_ref
    } else {
      cn <- find_int_coef(fe_names, term, group_var, lev)
      if (is.na(cn)) stop("Cannot find interaction coef for level: ", lev)
      beta_g <- beta_ref + fe[[cn]]
    }
    
    re_slopes %>%
      mutate(
        !!label_var := lev,
        slope_total   = beta_g + ranef_est,
        slope_total_l = beta_g + ranef_lo,
        slope_total_u = beta_g + ranef_hi,
        is_overall    = FALSE
      )
  })
  
  slopes_countries <- bind_rows(out) %>%
    mutate(!!label_var := factor(.data[[label_var]], levels = levs))
  
  # ---- Overall (fixed) slope per group level: est ± 1.96*SE (fixed-effect only) ----
  overall_rows <- lapply(levs, function(lev) {
    L <- rep(0, length(fe)); names(L) <- fe_names
    L[term] <- 1
    
    if (lev != ref) {
      cn <- find_int_coef(fe_names, term, group_var, lev)
      if (is.na(cn)) stop("Cannot find interaction coef for level: ", lev)
      L[cn] <- 1
    }
    
    est <- sum(L * fe)
    se  <- sqrt(as.numeric(t(L) %*% V %*% L))
    
    data.frame(
      Country       = "Overall (fixed)",
      ranef_est     = NA_real_,
      ranef_lo      = NA_real_,
      ranef_hi      = NA_real_,
      slope_total   = est,
      slope_total_l = est - 1.96 * se,
      slope_total_u = est + 1.96 * se,
      is_overall    = TRUE,
      stringsAsFactors = FALSE
    ) %>%
      mutate(!!label_var := lev)
  })
  
  overall_rows <- bind_rows(overall_rows) %>%
    mutate(!!label_var := factor(.data[[label_var]], levels = levs))
  
  bind_rows(slopes_countries, overall_rows)
}

plot_forest_by_group <- function(slopes_df, group_col,
                                 xlab,
                                 show_y_title = TRUE,
                                 show_y_labels = FALSE) {
  
  dd <- slopes_df %>%
    group_by(.data[[group_col]]) %>%
    mutate(
      # 让 Overall 在每个 facet 顶部（只影响排序，不改变数值）
      slope_for_order = if_else(is_overall,
                                max(slope_total, na.rm = TRUE) + 1e6,
                                slope_total),
      Country_f = reorder_within(Country, slope_for_order, .data[[group_col]])
    ) %>%
    ungroup()
  
  ggplot(dd, aes(x = slope_total, y = Country_f, colour = .data[[group_col]])) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    
    # errorbar: countries (mapped colour)
    geom_errorbarh(
      data = dd %>% filter(!is_overall),
      aes(xmin = slope_total_l, xmax = slope_total_u),
      height = 0.2, alpha = 0.7
    ) +
    # errorbar: overall (black), same height/alpha/linewidth
    geom_errorbarh(
      data = dd %>% filter(is_overall),
      aes(xmin = slope_total_l, xmax = slope_total_u),
      height = 0.2, alpha = 0.7,
      colour = "black"
    ) +
    
    # points: countries
    geom_point(
      data = dd %>% filter(!is_overall),
      size = 1.5
    ) +
    # points: overall (black), same size
    geom_point(
      data = dd %>% filter(is_overall),
      size = 1.5,
      colour = "black"
    ) +
    
    facet_wrap(stats::as.formula(paste0("~", group_col)), ncol = 1, scales = "free_y") +
    tidytext::scale_y_reordered() +
    coord_cartesian(clip = "off") +
    labs(x = xlab, y = if (show_y_title) "Country" else NULL) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text       = element_blank(),
      strip.background = element_blank(),
      panel.spacing    = unit(0.8, "lines"),
      axis.text.y      = if (show_y_labels) element_text(size = 7) else element_blank(),
      axis.ticks.y     = if (show_y_labels) element_line() else element_blank(),
      plot.margin      = margin(5.5, 50, 5.5, 5.5),
      legend.position  = "none"
    )
}


ctrl <- glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))

## ==============================================================
## PART A. Age at first marriage
## ==============================================================

## ----------------------------
## A0) Read + filter
## ----------------------------
df_age_raw <- read_csv("marriage age sex ratio with covariates_2.csv", show_col_types = FALSE) %>%
  mutate(Year = as.integer(Year)) %>%
  filter(!is.na(Year), Year >= 1970)

## ----------------------------
## A1) Aggregate same-year SR into 15–29 / 30–49
## ----------------------------
sr_cols_age <- list(
  same_1529 = c("SR_same_[15-19]", "SR_same_[20-24]", "SR_same_[25-29]"),
  same_3049 = c("SR_same_[30-34]", "SR_same_[35-39]", "SR_same_[40-44]", "SR_same_[45-49]")
)
sr_cols_age <- lapply(sr_cols_age, keep_existing, nms = names(df_age_raw))
if (length(sr_cols_age$same_1529) == 0 && length(sr_cols_age$same_3049) == 0) {
  stop("No SR_same_* columns found in marriage dataset. Check column names.")
}

df_age <- df_age_raw %>%
  mutate(
    Country      = factor(Country),
    MainReligion = factor(MainReligion),
    
    # covariates (log + scale)
    log_GDP     = log(GDP_per_capita + 1),
    popdens_log = log(Population_density + 1),
    log_GDP_sc     = as.numeric(scale(log_GDP)),
    Median_age_sc  = as.numeric(scale(Median_age)),
    popdens_log_sc = as.numeric(scale(popdens_log)),
    
    # band SR (same-year)
    SR_same_1529 = if (length(sr_cols_age$same_1529) > 0) rowMeans(across(all_of(sr_cols_age$same_1529)), na.rm = TRUE) else NA_real_,
    SR_same_3049 = if (length(sr_cols_age$same_3049) > 0) rowMeans(across(all_of(sr_cols_age$same_3049)), na.rm = TRUE) else NA_real_,
    
    # logit then STANDARDISE
    sr_same_1529_log    = logit_safe(SR_same_1529),
    sr_same_3049_log    = logit_safe(SR_same_3049),
    sr_same_1529_log_sc = as.numeric(scale(sr_same_1529_log)),
    sr_same_3049_log_sc = as.numeric(scale(sr_same_3049_log))
  )

## ----------------------------
## A2) Build sex-specific datasets (choose SR band by marriage-age band)
## ----------------------------
prep_sex_same <- function(dat, age_col, sex_label) {
  dat %>%
    mutate(
      age_at_marriage = .data[[age_col]],
      Sex = sex_label,
      AgeBand2 = case_when(
        age_at_marriage >= 15 & age_at_marriage < 30 ~ "15–29",
        age_at_marriage >= 30 & age_at_marriage < 50 ~ "30–49",
        TRUE ~ NA_character_
      ),
      AgeBand2 = factor(AgeBand2, levels = c("15–29","30–49")),
      sr_same_log_sc = case_when(
        AgeBand2 == "15–29" ~ sr_same_1529_log_sc,
        AgeBand2 == "30–49" ~ sr_same_3049_log_sc,
        TRUE ~ NA_real_
      )
    ) %>%
    select(
      Country, Year, MainReligion, Sex, AgeBand2,
      age_at_marriage, sr_same_log_sc,
      log_GDP_sc, Median_age_sc, popdens_log_sc
    ) %>%
    drop_na(age_at_marriage, sr_same_log_sc, log_GDP_sc, Median_age_sc, popdens_log_sc, MainReligion) %>%
    mutate(Sex = factor(Sex, levels = c("Men","Women")))
}

df_men   <- prep_sex_same(df_age, "Marriage_Age_Men",   "Men")
df_women <- prep_sex_same(df_age, "Marriage_Age_Women", "Women")

## ----------------------------
## A3) Fit 2 models (Men & Women)
## ----------------------------
fit_age_same <- function(dat) {
  glmmTMB(
    age_at_marriage ~ sr_same_log_sc +
      log_GDP_sc + Median_age_sc + popdens_log_sc + MainReligion +
      (1 | Year) + (sr_same_log_sc || Country),
    family  = gaussian(),
    data    = dat,
    control = ctrl
  )
}

mod_men_age   <- fit_age_same(df_men)
mod_women_age <- fit_age_same(df_women)

summary(mod_men_age)
summary(mod_women_age)

## ----------------------------
## A4) Plot
## ----------------------------
# unify x-range for both sexes
x_all  <- c(df_men$sr_same_log_sc, df_women$sr_same_log_sc)
x_lims <- quantile(x_all, probs = c(0.025, 0.975), na.rm = TRUE)
x_seq  <- seq(x_lims[1], x_lims[2], length.out = 200)

pred_m <- pred_df(mod_men_age,   df_men,   "sr_same_log_sc", "Sex", x_seq, c("Men"))
pred_w <- pred_df(mod_women_age, df_women, "sr_same_log_sc", "Sex", x_seq, c("Women"))
pred_m$Sex <- factor("Men", levels = c("Men","Women"))
pred_w$Sex <- factor("Women", levels = c("Men","Women"))
pred_age <- bind_rows(pred_m, pred_w)

p_age_curve <- ggplot(pred_age, aes(x = sr_same_log_sc, y = fit, colour = Sex, fill = Sex)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1) +
  coord_cartesian(xlim = x_lims, ylim = c(20, 40)) +
  labs(
    x = "Sex ratio (log, standardised)",
    y = "Predicted age at first marriage",
    colour = NULL, fill = NULL
  ) +
  guides(fill = "none") +   # <- 关键：只保留一种图例（colour）
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.margin = margin(10,10,10,10))+ guides(fill = "none")

# slopes (by Sex, not AgeBand2)
sl_m <- broom.mixed::tidy(mod_men_age, effects = "ran_vals", component = "cond", conf.int = TRUE) %>%
  filter(group == "Country", term == "sr_same_log_sc") %>%
  transmute(
    Country = as.character(level),
    slope_total   = unname(fixef(mod_men_age)$cond["sr_same_log_sc"]) + estimate,
    slope_total_l = unname(fixef(mod_men_age)$cond["sr_same_log_sc"]) + conf.low,
    slope_total_u = unname(fixef(mod_men_age)$cond["sr_same_log_sc"]) + conf.high,
    Sex = factor("Men", levels = c("Men","Women")),
    is_overall = FALSE
  ) %>%
  bind_rows({
    b <- fixef(mod_men_age)$cond
    V <- vcov(mod_men_age)$cond
    est <- unname(b["sr_same_log_sc"])
    se  <- sqrt(unname(V["sr_same_log_sc","sr_same_log_sc"]))
    tibble(
      Country = "Overall (fixed)",
      slope_total   = est,
      slope_total_l = est - 1.96 * se,
      slope_total_u = est + 1.96 * se,
      Sex = factor("Men", levels = c("Men","Women")),
      is_overall = TRUE
    )
  })


sl_w <- broom.mixed::tidy(mod_women_age, effects = "ran_vals", component = "cond", conf.int = TRUE) %>%
  filter(group == "Country", term == "sr_same_log_sc") %>%
  transmute(
    Country = as.character(level),
    slope_total   = unname(fixef(mod_women_age)$cond["sr_same_log_sc"]) + estimate,
    slope_total_l = unname(fixef(mod_women_age)$cond["sr_same_log_sc"]) + conf.low,
    slope_total_u = unname(fixef(mod_women_age)$cond["sr_same_log_sc"]) + conf.high,
    Sex = factor("Women", levels = c("Men","Women")),
    is_overall = FALSE
  ) %>%
  bind_rows({
    b <- fixef(mod_women_age)$cond
    V <- vcov(mod_women_age)$cond
    est <- unname(b["sr_same_log_sc"])
    se  <- sqrt(unname(V["sr_same_log_sc","sr_same_log_sc"]))
    tibble(
      Country = "Overall (fixed)",
      slope_total   = est,
      slope_total_l = est - 1.96 * se,
      slope_total_u = est + 1.96 * se,
      Sex = factor("Women", levels = c("Men","Women")),
      is_overall = TRUE
    )
  })




summarise_country_slopes_age <- function(slopes_df, outcome_label, sex_label) {
  slopes_df %>%
    filter(!is_overall) %>%
    group_by(Sex) %>% 
    summarise(
      n_country  = n_distinct(Country),
      mean_slope = mean(slope_total, na.rm = TRUE),
      sd_slope   = sd(slope_total,   na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(
      Outcome = outcome_label,
      SexLab  = sex_label
    ) %>%
    select(Outcome, SexLab, Sex, n_country, mean_slope, sd_slope)
}

age_men_slopes   <- summarise_country_slopes_age(sl_m, "Age at marriage", "Men")
age_women_slopes <- summarise_country_slopes_age(sl_w, "Age at marriage", "Women")

bind_rows(age_men_slopes, age_women_slopes)


sl_age <- bind_rows(sl_m, sl_w)

p_age_forest <- plot_forest_by_group(
  sl_age, group_col = "Sex",
  xlab = "Country-specific slope",
  show_y_title = TRUE, show_y_labels = FALSE
) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )+ guides(colour = "none")


fig_age <- (p_age_curve / p_age_forest) +
  plot_layout(heights = c(1, 1.35), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "right", legend.title = element_blank())


fig_age

## ==============================================================
## PART B. Female fertility (ASFR) 
## ==============================================================

## ----------------------------
## B0) Read + filter
## ----------------------------
df_asfr_raw <- read_csv("ASFR_SR_with_covariates.csv", show_col_types = FALSE) %>%
  mutate(Year = as.integer(Year)) %>%
  filter(!is.na(Year), Year >= 1970)

req2 <- c("Country","Year","age_group","ASFR","SR",
          "GDP_per_capita","Population_density","Median_age","MainReligion")
miss2 <- setdiff(req2, names(df_asfr_raw))
if (length(miss2) > 0) stop("Missing required columns in ASFR dataset: ", paste(miss2, collapse = ", "))

age_levels_7 <- c("[15-19]","[20-24]","[25-29]","[30-34]","[35-39]","[40-44]","[45-49]")

## ----------------------------
## B1) Prep
## ----------------------------
df_asfr <- df_asfr_raw %>%
  mutate(
    Country      = factor(Country),
    Year         = as.integer(Year),
    age_group    = factor(age_group, levels = age_levels_7),
    MainReligion = factor(MainReligion),
    
    # outcome
    log_ASFR = log(as.numeric(ASFR) + 1),
    
    # SR logit then STANDARDISE
    sr_log    = logit_safe(SR),
    sr_log_sc = as.numeric(scale(sr_log)),
    
    # covariates (log + scale)
    log_GDP     = log(GDP_per_capita + 1),
    popdens_log = log(Population_density + 1),
    
    log_GDP_sc     = as.numeric(scale(log_GDP)),
    Median_age_sc  = as.numeric(scale(Median_age)),
    popdens_log_sc = as.numeric(scale(popdens_log)),
    
    # collapse to 2 bands
    AgeBand2 = case_when(
      age_group %in% c("[15-19]","[20-24]","[25-29]") ~ "15–29",
      age_group %in% c("[30-34]","[35-39]","[40-44]","[45-49]") ~ "30–49",
      TRUE ~ NA_character_
    ),
    AgeBand2 = factor(AgeBand2, levels = c("15–29","30–49"))
  ) %>%
  filter(!is.na(AgeBand2))

df_asfr_mod2 <- df_asfr %>%
  select(
    Country, Year, AgeBand2, MainReligion,
    log_ASFR, sr_log_sc,
    log_GDP_sc, Median_age_sc, popdens_log_sc
  ) %>%
  drop_na(
    log_ASFR, sr_log_sc, AgeBand2, MainReligion,
    log_GDP_sc, Median_age_sc, popdens_log_sc
  )

## ----------------------------
## B2) Fit model (2 bands only)
## ----------------------------
mod_asfr_2 <- glmmTMB(
  log_ASFR ~ sr_log_sc * AgeBand2 +
    log_GDP_sc + Median_age_sc + popdens_log_sc + MainReligion +
    (1 | Year) + (sr_log_sc || Country),
  family  = gaussian(),
  data    = df_asfr_mod2,
  control = ctrl
)

summary(mod_asfr_2)

## ----------------------------
## B3) Plot
## ----------------------------
# prediction grid on log_ASFR scale, then back-transform to ASFR
plot_asfr_curves_2band <- function(data, model) {
  x_lims <- quantile(data$sr_log_sc, probs = c(0.025, 0.975), na.rm = TRUE)
  x_seq  <- seq(x_lims[1], x_lims[2], length.out = 200)
  ref_rel <- levels(data$MainReligion)[1]
  
  newdat <- expand.grid(
    sr_log_sc     = x_seq,
    AgeBand2      = levels(data$AgeBand2),
    MainReligion  = ref_rel
  )
  newdat$AgeBand2     <- factor(newdat$AgeBand2, levels = levels(data$AgeBand2))
  newdat$MainReligion <- factor(newdat$MainReligion, levels = levels(data$MainReligion))
  
  newdat$log_GDP_sc          <- 0
  newdat$Median_age_sc       <- 0
  newdat$popdens_log_sc      <- 0
  
  pr <- predict(model, newdata = newdat, type = "response", se.fit = TRUE, re.form = NA)
  
  newdat$fit   <- exp(pr$fit) - 1
  newdat$lower <- exp(pr$fit - 1.96 * pr$se.fit) - 1
  newdat$upper <- exp(pr$fit + 1.96 * pr$se.fit) - 1
  
  ggplot(newdat, aes(x = sr_log_sc, y = fit, colour = AgeBand2, fill = AgeBand2)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.18, colour = NA, show.legend = FALSE) +
    geom_line(linewidth = 1) +
    coord_cartesian(xlim = x_lims) +
    labs(
      x = "Sex ratio (log, standardised)",
      y = "Predicted ASFR",
      colour = NULL,
      fill = NULL
    ) +
    scale_colour_discrete(
      labels = c(
        "15–29" = "Youth (15–29)",
        "30–49" = "Adult (30–49)"
      )
    ) +
    guides(fill = "none") +
    theme_bw(base_size = 11)
  
}

p_asfr_curve <- plot_asfr_curves_2band(df_asfr_mod2, mod_asfr_2)

sl_asfr2 <- get_slopes_by_group(mod_asfr_2, term = "sr_log_sc", group_var = "AgeBand2", label_var = "AgeBand2")
p_asfr_forest <- plot_forest_by_group(sl_asfr2, group_col = "AgeBand2", xlab = "Country-specific slope",
                                      show_y_title = TRUE, show_y_labels = FALSE)+ guides(colour = "none")

fig_asfr <- (p_asfr_curve / p_asfr_forest) +
  plot_layout(heights = c(1, 1.35), guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(legend.position = "right", legend.title = element_blank())


fig_asfr

summarise_country_slopes_asfr <- function(slopes_df, outcome_label) {
  slopes_df %>%
    filter(!is_overall) %>%
    group_by(AgeBand2) %>%   # sl_asfr2 里是 AgeBand2
    summarise(
      n_country  = n_distinct(Country),
      mean_slope = mean(slope_total, na.rm = TRUE),
      sd_slope   = sd(slope_total,   na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(
      Outcome = outcome_label,
      Sex     = "Women"
    ) %>%
    select(Outcome, Sex, AgeBand2, n_country, mean_slope, sd_slope)
}

asfr_slopes <- summarise_country_slopes_asfr(sl_asfr2, "Fertility (ASFR)")
asfr_slopes


final_2panel <- (fig_age | fig_asfr) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = "topleft"
  )


final_2panel



tab_model(
  mod_men_age, mod_women_age, mod_asfr_2,
  dv.labels  = c("Marriage age in men", "Marriage age in women",
                 "Women fertility"),
  show.re.var = FALSE,
  show.icc    = FALSE,
  show.aic    = TRUE,
  transform   = NULL   # no exp()
)
