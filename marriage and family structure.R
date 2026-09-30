
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("pak")
pak::pak("base64enc", dependencies = TRUE)
pak::pak("tidyverse")

pak::pak("tidyverse")


install.packages(c("tidyverse","ggridges","viridis","scales"))


library(tidyverse)
library(ggridges)
library(viridis)
library(scales)

#===============================================================================
# Descriptive information about sex ratio panel data
#===============================================================================

new_file <- "marriage status (full ASR) with covariates.csv"
df <- read.csv(new_file, stringsAsFactors = FALSE, check.names = FALSE) |> as_tibble()

required_cols <- c("Country","Year","BroadAgeGroup","SR")
if (!all(required_cols %in% names(df))) {
  stop(paste("Missing required columns:",
             paste(setdiff(required_cols, names(df)), collapse = ", ")))
}
age_levels <- c("Under 15","Youth (15-29)","Adult (30-49)","Older (50+)")

dat_plot <- df |>
  transmute(
    Country   = as.character(.data$Country),
    Year      = as.integer(.data$Year),
    age_group = factor(.data$BroadAgeGroup, levels = age_levels, ordered = TRUE),
    SR_raw    = as.numeric(.data$SR)
  ) |>
  filter(!is.na(Country), !is.na(Year), !is.na(age_group), !is.na(SR_raw))

if (max(dat_plot$SR_raw, na.rm = TRUE) > 1.5) {
  dat_plot <- dat_plot |> mutate(SR = SR_raw/100)
} else {
  dat_plot <- dat_plot |> mutate(SR = SR_raw)
}


eps <- 1e-6
dat_plot <- dat_plot |>
  mutate(
    SR_c  = pmin(pmax(SR, eps), 1 - eps),
    Year_f = factor(Year, levels = sort(unique(Year)))
  )

anno <- dat_plot |>
  group_by(age_group, Year, Year_f) |>
  summarise(
    n_countries = n_distinct(Country),
    n_obs       = n(),
    .groups = "drop"
  ) |>
  mutate(
    label = paste0("(", n_countries, ", ", n_obs, ")")  # e.g. (302, 302)
  )

out_dir <- "ridgeline_by_age_group"
if (!dir.exists(out_dir)) dir.create(out_dir)

fill_pal <- viridis::viridis(256, option = "C")

plot_one_age <- function(g) {
  d <- dat_plot |> filter(age_group == g)
  if (nrow(d) == 0) {
    message("Skip: no data for ", g)
    return(invisible(NULL))
  }
  a <- anno |> filter(age_group == g)
  
  p <- ggplot(d, aes(x = SR_c, y = Year_f, height = ..density.., group = Year_f, fill = ..x..)) +
    stat_density_ridges(
      geom = "density_ridges_gradient",
      calc_ecdf = FALSE,
      scale = 2.2,
      rel_min_height = 0.001,
      size = 0.2
    ) +
    # Right-side label: (countries, observations)
    geom_text(
      data = a,
      aes(x = 1.02, y = Year_f, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 0.5,
      size = 3
    ) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.05)) +
    scale_fill_gradientn(colors = fill_pal, name = "SR") +
    labs(
      title = paste0("SR Distribution by Year — Age Group: ", as.character(g)),
      subtitle = "Right-side labels show (number of countries, total data points) per year",
      x = "Sex ratio (SR, %)", y = "Year"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
  
  print(p)
  ggsave(file.path(out_dir, paste0("ridgeline_", g, ".png")),
         p, width = 9.5, height = 6.5, dpi = 300)
}

for (g in age_levels) {
  plot_one_age(g)
}

message("Done: ridgeline plots saved to folder: ", out_dir)

#===============================================================================
# 1.1. MARRIAGE AND DIVORCE
#===============================================================================

## --------------------------------------------------------------
## 1. Load required packages
## --------------------------------------------------------------
library(dplyr)
library(tibble)
library(glmmTMB)
library(ggplot2)
library(patchwork)
library(cowplot)
library(scales)
library(performance)
library(skimr)
library(sjPlot)
library(sjmisc)
library(broom.mixed)
library(tidytext)
library(grid)     

## ==============================================================
## 2. data loading
## ==============================================================
df <- read.csv("marriage status (full ASR) with covariates.csv",
               stringsAsFactors = FALSE, check.names = FALSE) |>
  as_tibble()

skim(df)

## ==============================================================
## 3. data cleaning
## ==============================================================
# Male Marriage
df1.1 <- df %>%
  filter(between(Year, 1970, 2018),
         BroadAgeGroup != "Under 15",
         Sex == "Male") %>%
  mutate(
    MarriagePrevalence  = pmin(pmax(marriage_prevalence / 100, 0.0001), 0.9999),
    BroadAgeGroup = factor(BroadAgeGroup,
                           levels = c("Youth (15-29)", "Adult (30-49)", "Older (50+)")),
    Country       = factor(Country),
    MainReligion  = factor(MainReligion),
    Year          = factor(Year),
    log_GDP       = log(GDP_per_capita + 1),
    popdens_log   = log(Population_density + 1),
    SR1           = pmin(pmax(SR, 1e-6), 1 - 1e-6),
    sr            = SR1 / (1 - SR1),
    sr_log        = log(sr),
    sr_log_sc     = as.numeric(scale(sr_log))
  ) %>%
  mutate(
    across(c(log_GDP, Median_age, popdens_log),
           ~ as.numeric(scale(.)), .names = "{.col}_sc")
  ) %>%
  select(
    MarriagePrevalence,
    sr, sr_log, sr_log_sc,
    BroadAgeGroup,
    log_GDP_sc, Median_age_sc, popdens_log_sc,
    MainReligion, Year, Country
  ) %>%
  drop_na()

cat("df1.1 (Male Marriage): ", nrow(df1.1), " rows\n")

# Female Marriage
df1.2 <- df %>%
  filter(between(Year, 1970, 2018),
         BroadAgeGroup != "Under 15",
         Sex == "Female") %>%
  mutate(
    MarriagePrevalence  = pmin(pmax(marriage_prevalence / 100, 0.0001), 0.9999),
    BroadAgeGroup = factor(BroadAgeGroup,
                           levels = c("Youth (15-29)", "Adult (30-49)", "Older (50+)")),
    Country       = factor(Country),
    MainReligion  = factor(MainReligion),
    Year          = factor(Year),
    log_GDP       = log(GDP_per_capita + 1),
    popdens_log   = log(Population_density + 1),
    SR1           = pmin(pmax(SR, 1e-6), 1 - 1e-6),
    sr            = SR1 / (1 - SR1),
    sr_log        = log(sr),
    sr_log_sc     = as.numeric(scale(sr_log))
  ) %>%
  mutate(
    across(c(log_GDP, Median_age, popdens_log),
           ~ as.numeric(scale(.)), .names = "{.col}_sc")
  ) %>%
  select(
    MarriagePrevalence,
    sr, sr_log, sr_log_sc,
    BroadAgeGroup,
    log_GDP_sc, Median_age_sc, popdens_log_sc,
    MainReligion, Year, Country
  ) %>%
  drop_na()

cat("df1.2 (Female Marriage): ", nrow(df1.2), " rows\n")


# Male Divorce
df2.1 <- df %>%
  filter(between(Year, 1970, 2018),
         BroadAgeGroup != "Under 15",
         Sex == "Male") %>%
  mutate(
    DivorcePrevalence   = pmin(pmax(divorce_prevalence / 100, 0.0001), 0.9999),
    BroadAgeGroup = factor(BroadAgeGroup,
                           levels = c("Youth (15-29)", "Adult (30-49)", "Older (50+)")),
    Country       = factor(Country),
    MainReligion  = factor(MainReligion),
    Year          = factor(Year),
    log_GDP       = log(GDP_per_capita + 1),
    popdens_log   = log(Population_density + 1),
    SR1           = pmin(pmax(SR, 1e-6), 1 - 1e-6),
    sr            = SR1 / (1 - SR1),
    sr_log        = log(sr),
    sr_log_sc     = as.numeric(scale(sr_log))
  ) %>%
  mutate(
    across(c(log_GDP, Median_age, popdens_log),
           ~ as.numeric(scale(.)), .names = "{.col}_sc")
  ) %>%
  select(
    DivorcePrevalence,
    sr, sr_log, sr_log_sc,
    BroadAgeGroup,
    log_GDP_sc, Median_age_sc, popdens_log_sc,
    MainReligion, Year, Country
  ) %>%
  drop_na()

cat("df2.1 (Male Divorce): ", nrow(df2.1), " rows\n")

# Female Divorce
df2.2 <- df %>%
  filter(between(Year, 1970, 2018),
         BroadAgeGroup != "Under 15",
         Sex == "Female") %>%
  mutate(
    DivorcePrevalence   = pmin(pmax(divorce_prevalence / 100, 0.0001), 0.9999),
    BroadAgeGroup = factor(BroadAgeGroup,
                           levels = c("Youth (15-29)", "Adult (30-49)", "Older (50+)")),
    Country       = factor(Country),
    MainReligion  = factor(MainReligion),
    Year          = factor(Year),
    log_GDP       = log(GDP_per_capita + 1),
    popdens_log   = log(Population_density + 1),
    SR1           = pmin(pmax(SR, 1e-6), 1 - 1e-6),
    sr            = SR1 / (1 - SR1),
    sr_log        = log(sr),
    sr_log_sc     = as.numeric(scale(sr_log))
  ) %>%
  mutate(
    across(c(log_GDP, Median_age, popdens_log),
           ~ as.numeric(scale(.)), .names = "{.col}_sc")
  ) %>%
  select(
    DivorcePrevalence,
    sr, sr_log, sr_log_sc,
    BroadAgeGroup,
    log_GDP_sc, Median_age_sc, popdens_log_sc,
    MainReligion, Year, Country
  ) %>%
  drop_na()

cat("df2.2 (Female Divorce): ", nrow(df2.2), " rows\n")

## ==============================================================
## 5. models
## ==============================================================
mod1.1 <- glmmTMB(MarriagePrevalence ~ sr_log_sc * BroadAgeGroup + log_GDP_sc + Median_age_sc + popdens_log_sc + MainReligion + (1 | Year) + (sr_log_sc || Country),
                  family = beta_family(link = "logit"),data = df1.1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(mod1.1)

mod1.2 <- glmmTMB( MarriagePrevalence ~ sr_log_sc * BroadAgeGroup + log_GDP_sc + Median_age_sc + popdens_log_sc + MainReligion + (1 | Year) + (sr_log_sc || Country),
                   family = beta_family(link = "logit"), data = df1.2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(mod1.2)

mod2.1 <- glmmTMB( DivorcePrevalence ~ sr_log_sc * BroadAgeGroup + log_GDP_sc + Median_age_sc + popdens_log_sc + MainReligion + (1 | Year) + (sr_log_sc || Country),
                   family = beta_family(link = "logit"), data = df2.1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(mod2.1)

mod2.2 <- glmmTMB( DivorcePrevalence ~ sr_log_sc * BroadAgeGroup + log_GDP_sc + Median_age_sc + popdens_log_sc + MainReligion + (1 | Year) + (sr_log_sc || Country),
                   family = beta_family(link = "logit"), data = df2.2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
summary(mod2.2)

## ==============================================================
## 6. result output
## ==============================================================
tab_model(
  mod1.1, mod1.2, mod2.1, mod2.2,
  dv.labels  = c("Marriage (df1.1)", "Marriage (df1.2)",
                 "Divorce (df2.1)",  "Divorce (df2.2)"),
  show.re.var = FALSE,
  show.icc    = FALSE,
  show.aic    = TRUE,
  transform   = NULL  
)

## ==============================================================
## 7. visualisation
## ==============================================================
plot_interaction <- function(data, model, response_var, y_max,
                             y_label = NULL, x_label = NULL) {
  x_lims <- quantile(data$sr_log_sc, probs = c(0.025, 0.975), na.rm = TRUE)
  
  ref_religion <- levels(data$MainReligion)[1]
  newdat <- expand.grid(
    sr_log_sc     = seq(x_lims[1], x_lims[2], length.out = 200),
    BroadAgeGroup = levels(data$BroadAgeGroup),
    MainReligion  = ref_religion
  )
  
  newdat$BroadAgeGroup <- factor(newdat$BroadAgeGroup,
                                 levels = levels(data$BroadAgeGroup))
  newdat$MainReligion  <- factor(newdat$MainReligion,
                                 levels = levels(data$MainReligion))
  
  newdat$log_GDP_sc     <- 0
  newdat$Median_age_sc  <- 0
  newdat$popdens_log_sc <- 0
  
  pred <- predict(
    model,
    newdata = newdat,
    type    = "link",
    se.fit  = TRUE,
    re.form = NA
  )
  
  newdat$fit   <- plogis(pred$fit)
  newdat$lower <- plogis(pred$fit - 1.96 * pred$se.fit)
  newdat$upper <- plogis(pred$fit + 1.96 * pred$se.fit)
  
  p <- ggplot() +
    geom_ribbon(
      data = newdat,
      aes(x = sr_log_sc,
          ymin = lower * 100,
          ymax = upper * 100,
          fill = BroadAgeGroup),
      alpha = 0.18, colour = NA
    ) +
    geom_line(
      data = newdat,
      aes(x = sr_log_sc,
          y = fit * 100,
          colour = BroadAgeGroup),
      linewidth = 1
    ) +
    coord_cartesian(xlim = x_lims, ylim = c(0, y_max)) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    scale_colour_discrete(name = "Age Group") +
    scale_fill_discrete(guide = "none") +
    labs(x = x_label, y = y_label) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(10, 10, 10, 10)
    )
  
  print(p)
  invisible(p)
}


p1 <- plot_interaction(
  df1.1, mod1.1, "MarriagePrevalence", 100,
  y_label = "Married/in-union prevalence  (%)",
  x_label = "Sex ratio (log, standardised)"
)

p2 <- plot_interaction(
  df1.2, mod1.2, "MarriagePrevalence", 100,
  y_label = NULL,
  x_label = "Sex ratio (log, standardised)"
) +
  theme(
    legend.position = "none",
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank()
  )

p3 <- plot_interaction(
  df2.1, mod2.1, "DivorcePrevalence", 20,
  y_label = "Divorced/separated prevalence  (%)",
  x_label = "Sex ratio (log, standardised)"
)

p4 <- plot_interaction(
  df2.2, mod2.2, "DivorcePrevalence", 20,
  y_label = NULL,
  x_label = "Sex ratio (log, standardised)"
) +
  theme(
    legend.position = "none",
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank()
  )


get_country_slopes_by_age <- function(model,
                                      term    = "sr_log_sc",
                                      n_label = 0) {
  

  fe <- fixef(model)$cond
  
  mf      <- model.frame(model)
  age_var <- "BroadAgeGroup"
  if (!age_var %in% names(mf)) stop("BroadAgeGroup not found in model frame.")
  age_levels <- levels(mf[[age_var]])
  ref_level  <- age_levels[1]
  
  if (!term %in% names(fe)) stop("Term ", term, " not found in fixed effects.")
  beta_sr <- unname(fe[term])
  

  re_tidy <- broom.mixed::tidy(
    model,
    effects   = "ran_vals",
    component = "cond",
    conf.int  = TRUE
  )
  
  re_slopes <- re_tidy %>%
    filter(group == "Country", term == !!term) %>%
    transmute(
      Country   = as.character(level),
      ranef_est = estimate,
      ranef_lo  = conf.low,
      ranef_hi  = conf.high
    )
  

  slopes_list <- lapply(age_levels, function(lev) {
    if (lev == ref_level) {
      beta_age <- beta_sr
    } else {
      coef_name <- paste0(term, ":", age_var, lev)
      if (!coef_name %in% names(fe)) stop("Interaction coefficient ", coef_name, " not found.")
      beta_age <- beta_sr + fe[[coef_name]]
    }
    
    re_slopes %>%
      mutate(
        AgeGroup      = lev,
        slope_total   = beta_age + ranef_est,
        slope_total_l = beta_age + ranef_lo,
        slope_total_u = beta_age + ranef_hi,
        is_overall    = FALSE
      )
  })
  
  slopes_all <- bind_rows(slopes_list) %>%
    mutate(
      AgeGroup      = factor(AgeGroup, levels = age_levels),
      BroadAgeGroup = AgeGroup
    )
  

  slopes_all <- slopes_all %>%
    group_by(AgeGroup) %>%
    mutate(
      Country_f   = reorder_within(Country, slope_total, AgeGroup),
      n_in_group  = n(),
      rank_slope  = rank(slope_total, ties.method = "first"),
      label_country = case_when(
        n_label > 0 & rank_slope <= n_label ~ Country,
        n_label > 0 & rank_slope >  n_in_group - n_label ~ Country,
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup()
  

  V <- vcov(model)$cond
  
  overall_rows <- lapply(age_levels, function(lev) {
    
    L <- rep(0, length(fe)); names(L) <- names(fe)
    L[term] <- 1
    
    if (lev != ref_level) {
      coef_name <- paste0(term, ":", age_var, lev)
      if (!coef_name %in% names(fe)) stop("Interaction coefficient ", coef_name, " not found.")
      L[coef_name] <- 1
    }
    
    est <- sum(L * fe)
    se  <- sqrt(as.numeric(t(L) %*% V %*% L))
    
    data.frame(
      Country       = "Overall (fixed)",
      AgeGroup      = lev,
      ranef_est     = NA_real_,
      ranef_lo      = NA_real_,
      ranef_hi      = NA_real_,
      slope_total   = est,
      slope_total_l = est - 1.96 * se,
      slope_total_u = est + 1.96 * se,
      Country_f     = NA_character_,
      n_in_group    = NA_integer_,
      rank_slope    = NA_real_,
      label_country = NA_character_,
      is_overall    = TRUE,
      stringsAsFactors = FALSE
    )
  })
  
  overall_rows <- bind_rows(overall_rows) %>%
    mutate(
      AgeGroup      = factor(AgeGroup, levels = age_levels),
      BroadAgeGroup = AgeGroup
    )
  
  slopes_all <- bind_rows(slopes_all, overall_rows)
  

  slopes_all <- slopes_all %>%
    group_by(AgeGroup) %>%
    mutate(
      slope_for_order = if_else(is_overall,
                                max(slope_total, na.rm = TRUE) + 1e6,
                                slope_total),
      Country_f = reorder_within(Country, slope_for_order, AgeGroup)
    ) %>%
    ungroup()
  
  return(slopes_all)
}


plot_country_forest_by_age <- function(slopes_by_age,
                                       xlab,
                                       show_y_title   = TRUE,
                                       show_y_labels  = FALSE) {
  
  if (!"label_country" %in% names(slopes_by_age)) {
    slopes_by_age$label_country <- NA_character_
  }
  
  label_data <- slopes_by_age %>%
    filter(!is.na(label_country))
  
  ggplot(
    slopes_by_age,
    aes(x = slope_total,
        y = Country_f,
        colour = BroadAgeGroup)
  ) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_errorbarh(
      data = slopes_by_age %>% dplyr::filter(!is_overall),
      aes(xmin = slope_total_l, xmax = slope_total_u),
      height = 0.2, alpha = 0.7
    ) +
    geom_errorbarh(
      data = slopes_by_age %>% dplyr::filter(is_overall),
      aes(xmin = slope_total_l, xmax = slope_total_u),
      height = 0.2, alpha = 0.7,
      colour = "black"
    ) +
    geom_point(
      data = slopes_by_age %>% dplyr::filter(!is_overall),
      size = 1.6
    ) +
    geom_point(
      data = slopes_by_age %>% dplyr::filter(is_overall),
      size = 1.6,
      colour = "black"
    ) +
    geom_text(
      data  = label_data,
      aes(label = label_country),
      hjust = -0.1,
      size  = 3,
      show.legend = FALSE
    ) +
    facet_wrap(~ AgeGroup, ncol = 1, scales = "free_y") +
    scale_y_reordered() +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 20)) +  # 23=实心菱形(可填充)
    scale_size_manual(values  = c(`FALSE` = 1.6, `TRUE` = 2.0)) +
    scale_linewidth_manual(values = c(`FALSE` = 0.6, `TRUE` = 1.2)) +
    coord_cartesian(clip = "off") +
    labs(
      x     = xlab,
      y     = if (show_y_title) "Country" else NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.text        = element_blank(),
      strip.background  = element_blank(),
      panel.spacing     = unit(0.8, "lines"),
      axis.text.y       = if (show_y_labels) element_text(size = 7) else element_blank(),
      axis.ticks.y      = if (show_y_labels) element_line() else element_blank(),
      plot.margin       = margin(5.5, 50, 5.5, 5.5),
      legend.position   = "none"
    )
}


slopes_m11_age <- get_country_slopes_by_age(mod1.1, term = "sr_log_sc", n_label = 0)
slopes_m12_age <- get_country_slopes_by_age(mod1.2, term = "sr_log_sc", n_label = 0)
slopes_d21_age <- get_country_slopes_by_age(mod2.1, term = "sr_log_sc", n_label = 0)
slopes_d22_age <- get_country_slopes_by_age(mod2.2, term = "sr_log_sc", n_label = 0)

summarise_country_slopes <- function(slopes_by_age,
                                     outcome_label,
                                     sex_label) {
  
  slopes_by_age %>%
    filter(!is_overall) %>%
    group_by(AgeGroup) %>%
    summarise(
      n_country  = n_distinct(Country),
      mean_logit = mean(slope_total, na.rm = TRUE),
      sd_logit   = sd(slope_total, na.rm = TRUE),
      mean_OR    = mean(exp(slope_total), na.rm = TRUE),
      sd_OR      = sd(exp(slope_total), na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(
      Outcome = outcome_label,
      Sex     = sex_label
    ) %>%
    select(
      Outcome, Sex, AgeGroup, n_country,
      mean_logit, sd_logit,
      mean_OR,   sd_OR
    )
}

range_m11 <- summarise_country_slopes(slopes_m11_age, "Marriage", "Male")
range_m12 <- summarise_country_slopes(slopes_m12_age, "Marriage", "Female")
range_d21 <- summarise_country_slopes(slopes_d21_age, "Divorce",  "Male")
range_d22 <- summarise_country_slopes(slopes_d22_age, "Divorce",  "Female")

slope_meansd_all <- bind_rows(range_m11, range_m12, range_d21, range_d22)
slope_meansd_all


xlab_marriage <- "Country-specific slope"
xlab_divorce  <- "Country-specific slope"

p_m11_age <- plot_country_forest_by_age(
  slopes_m11_age,
  xlab          = xlab_marriage,
  show_y_title  = TRUE,
  show_y_labels = FALSE
)

p_m12_age <- plot_country_forest_by_age(
  slopes_m12_age,
  xlab          = xlab_marriage,
  show_y_title  = FALSE,
  show_y_labels = FALSE
)

p_d21_age <- plot_country_forest_by_age(
  slopes_d21_age,
  xlab          = xlab_divorce,
  show_y_title  = TRUE,
  show_y_labels = FALSE
)

p_d22_age <- plot_country_forest_by_age(
  slopes_d22_age,
  xlab          = xlab_divorce,
  show_y_title  = FALSE,
  show_y_labels = FALSE
)


all_8panel <- (p1 | p2 | p3 | p4) /
  (p_m11_age | p_m12_age | p_d21_age | p_d22_age) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    text              = element_text(size = 11),
    axis.title        = element_text(size = 11),
    axis.text         = element_text(size = 9),
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 9),
    plot.tag          = element_text(size = 14, face = "bold"),
    plot.tag.position = "topleft"
  )

all_8panel
p_m11_age
p_m12_age
p_d21_age
p_d22_age



#===============================================================================
# 2. FAMILY STRUCTURE
#===============================================================================
## --------------------------------------------------------------
## 1. Read and filter data 
## --------------------------------------------------------------
df3 <- read.csv("family structure with covariates.csv",
               stringsAsFactors = FALSE, check.names = FALSE) |> 
  as_tibble()


target_ages <- c("18-34", "35-59", "60+")

df_hist <- df3 %>%
  filter(Age %in% target_ages) %>%
  mutate(
    Age = factor(Age, levels = target_ages),
    # Convert all to percentage
    Female_Single_Parent_pct = Single_Parent_female ,
    Male_Single_Parent_pct   = Single_Parent_male,
    Sex_Ratio_pct            = SR_2010_2018*100
  ) %>%
  select(Age, Female_Single_Parent_pct, Male_Single_Parent_pct, Sex_Ratio_pct) %>%
  pivot_longer(
    cols = ends_with("_pct"),
    names_to = "Variable",
    values_to = "Value_pct"
  ) %>%
  mutate(
    Variable = case_when(
      Variable == "Female_Single_Parent_pct" ~ "Female Single Parent Rate",
      Variable == "Male_Single_Parent_pct"   ~ "Male Single Parent Rate",
      Variable == "Sex_Ratio_pct"            ~ "Sex Ratio (2010–2018)"
    ),
    Variable = factor(Variable,
                      levels = c("Female Single Parent Rate",
                                 "Male Single Parent Rate",
                                 "Sex Ratio (2010–2018)"))
  )


ggplot(df_hist, aes(x = Value_pct, fill = Age)) +
  geom_histogram(
    alpha = 0.6,
    color = "white",
    position = "identity",   
    bins = 20                
  ) +
  facet_grid(Variable ~ Age, scales = "free_x") +
  
  # X-axis: all in %
  scale_x_continuous(
    labels = percent_format(accuracy = 1, scale = 1),
    name = "Percentage (%)"
  ) +
  
  scale_y_continuous(name = "Count") +
  
  # Colors for age groups
  scale_fill_manual(
    values = c("18-34" = "#66BB6A", "35-59" = "#FFA726", "60+" = "#EF5350"),
    name = "Age Group"
  ) +
  
  # Clean theme
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 10),
    panel.spacing = unit(1, "lines")
  ) +
  
  labs(
    title = "Histogram of Single Parenthood and Sex Ratio by Age Group",
    subtitle = "Raw data | All values in percentage (%)",
    x = "Percentage (%)",
    y = "Frequency (Count)",
    caption = "Data: df3 | Bins = 20"
  )


# single father
df3.1 <- df3 %>%
  filter(Age != "0-17") %>%
  mutate(
    SingleFather  = pmin(pmax(Single_Parent_male / 100, 0.0001), 0.9999),
    BroadAgeGroup = factor(Age, levels = c("18-34", "35-59", "60+")),
    Country       = factor(Country_Standardized),
    MainReligion  = factor(MainReligion_mode),
    log_GDP       = log(GDP_per_capita + 1),
    popdens_log   = log(Population_density + 1),
    SR1   = pmin(pmax(SR_2010_2018, 1e-6), 1 - 1e-6),
    sr  = SR1 / (1 - SR1),
    sr_log    = log(sr),
    sr_log_sc = as.numeric(scale(sr_log))
  ) %>%
  mutate(
    across(c(log_GDP, Median_age, popdens_log),
           ~ as.numeric(scale(.)), .names = "{.col}_sc")
  ) %>%
  select(SingleFather, sr, sr_log, sr_log_sc, BroadAgeGroup,
         log_GDP_sc, Median_age_sc, popdens_log_sc,
         MainReligion, Country) %>%
  tidyr::drop_na()

cat("Final dataset: ", nrow(df3.1), " rows\n")

# single mother
df3.2 <- df3 %>%
  filter(Age != "0-17") %>%
  mutate(
    SingleMother  = pmin(pmax(Single_Parent_female / 100, 0.0001), 0.9999),
    BroadAgeGroup = factor(Age, levels = c("18-34", "35-59", "60+")),
    Country       = factor(Country_Standardized),
    MainReligion  = factor(MainReligion_mode),
    log_GDP       = log(GDP_per_capita + 1),
    popdens_log   = log(Population_density + 1),
    SR1   = pmin(pmax(SR_2010_2018, 1e-6), 1 - 1e-6),
    sr  = SR1 / (1 - SR1),
    sr_log    = log(sr),
    sr_log_sc = as.numeric(scale(sr_log))
  ) %>%
  mutate(
    across(c(log_GDP, Median_age, popdens_log),
           ~ as.numeric(scale(.)), .names = "{.col}_sc")
  ) %>%
  select(SingleMother, sr, sr_log, sr_log_sc, BroadAgeGroup,
         log_GDP_sc, Median_age_sc, popdens_log_sc,
         MainReligion, Country) %>%
  tidyr::drop_na()

cat("Final dataset: ", nrow(df3.2), " rows\n")


## --------------------------------------------------------------
## 3. Models
## --------------------------------------------------------------

mod3.1 <- glmmTMB(
  SingleFather ~  sr_log_sc * BroadAgeGroup + log_GDP_sc + Median_age_sc + popdens_log_sc + MainReligion + (1| Country),
  family = beta_family(link = "logit"),
  data = df3.1, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))


mod3.2 <- glmmTMB(
  SingleMother ~ sr_log_sc  * BroadAgeGroup + log_GDP_sc + Median_age_sc + popdens_log_sc + MainReligion + (1| Country),
  family = beta_family(link = "logit"),
  data = df3.2, control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

summary(mod3.1)
summary(mod3.2)


VarCorr(mod3.1)
r2_nakagawa(mod3.1)
icc(mod3.1)


VarCorr(mod3.2)
r2_nakagawa(mod3.2)
icc(mod3.2)



tab_model(
  mod3.1, mod3.2,
  dv.labels = c("Single Father", "Single Mother"),
  show.re.var = FALSE,
  show.icc    = FALSE,
  show.aic    = TRUE,
  transform   = NULL  
)


## --------------------------------------------------------------
## 3. plotting
## --------------------------------------------------------------

plot_interaction_2 <- function(data, model, response_var, y_max,
                               y_label = NULL, x_label = NULL) {

  x_lims <- quantile(data$sr_log_sc, probs = c(0.025, 0.975), na.rm = TRUE)
  

  newdat <- expand.grid(
    sr_log_sc     = seq(x_lims[1], x_lims[2], length.out = 200),
    BroadAgeGroup = levels(data$BroadAgeGroup)
  )
  newdat$log_GDP_sc     <- 0
  newdat$Median_age_sc  <- 0
  newdat$popdens_log_sc <- 0
  

  if ("MainReligion" %in% names(model$frame)) {
    ref_religion <- levels(model$frame$MainReligion)[1] 
    newdat$MainReligion <- factor(ref_religion, levels = levels(model$frame$MainReligion))
  }
  

  pred <- predict(model, newdata = newdat, type = "link", se.fit = TRUE, re.form = NA)
  
  newdat <- newdat %>%
    mutate(
      fit   = plogis(pred$fit),
      lower = plogis(pred$fit - 1.96 * pred$se.fit),
      upper = plogis(pred$fit + 1.96 * pred$se.fit)
    )
  

  p <- ggplot() +
    geom_ribbon(
      data = newdat,
      aes(x = sr_log_sc, ymin = lower * 100, ymax = upper * 100, fill = BroadAgeGroup),
      alpha = 0.18, colour = NA
    ) +
    geom_line(
      data = newdat,
      aes(x = sr_log_sc, y = fit * 100, colour = BroadAgeGroup),
      linewidth = 1
    ) +
    coord_cartesian(xlim = x_lims, ylim = c(0, y_max)) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    scale_colour_discrete(name = "Age Group") +
    scale_fill_discrete(name = "Age Group") +
    labs(x = x_label, y = y_label) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  return(p)
}

p5 <- plot_interaction_2(df3.1, mod3.1, "SingleFather", 30,
                       y_label = "Single father (%)", x_label = "Sex Ratio (log, standardised)")
p6 <- plot_interaction_2(df3.2, mod3.2, "SingleMother", 30,
                       y_label = "Single mother (%)", x_label = "Sex Ratio (log, standardised)")


final_plot <- (p5 | p6)+
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "a",         
    tag_prefix = "(",
    tag_suffix = ")",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11),
    plot.tag        = element_text(size = 14, face = "bold"),
    plot.tag.position = "topleft"
  ) &
  guides(color = guide_legend(title = "Age Group", nrow = 1))

print(final_plot)

