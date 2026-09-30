

# ---- Packages ----
library(tidyverse)
library(ggridges)
library(viridis)
library(scales)
library(patchwork)

# ---- Read data ----
file_in <- "male_share_broad_cleaned_fixed.csv"
df <- read.csv(file_in, stringsAsFactors = FALSE, check.names = FALSE) |> as_tibble()

# ---- Required columns ----
req <- c("Country", "Year", "age_group", "SR")
if (!all(req %in% names(df))) {
  stop("Missing required columns: ",
       paste(setdiff(req, names(df)), collapse = ", "))
}

# ---- Constants ----
age_levels <- c("0-14","15-29","30-49","50+")
age_labels <- c("Under 15","Youth (15-29)","Adult (30-49)","Older (50+)")
fill_pal   <- viridis::viridis(256, option = "C")

# ---- Prepare data (stable /100 logic) ----
dat_plot <- df |>
  transmute(
    Country   = as.character(.data$Country),
    Year      = as.integer(.data$Year),
    age_group = as.character(.data$age_group),
    SR_raw    = as.numeric(.data$SR)
  ) |>
  filter(!is.na(Country), !is.na(Year), !is.na(age_group), !is.na(SR_raw)) |>
  filter(age_group %in% age_levels) |>
  mutate(
    age_group = factor(age_group, levels = age_levels, ordered = TRUE),
    age_lab   = factor(age_labels[match(age_group, age_levels)],
                       levels = age_labels, ordered = TRUE)
  )

if (max(dat_plot$SR_raw, na.rm = TRUE) > 1.5) {
  dat_plot <- dat_plot |> mutate(SR = SR_raw/100)
} else {
  dat_plot <- dat_plot |> mutate(SR = SR_raw)
}

dat_plot <- dat_plot |> mutate(SR_c = pmin(pmax(SR, 1e-6), 1 - 1e-6))


# ---- Shared axis range ----
sr_rng <- as.numeric(quantile(dat_plot$SR_c, probs = c(0.005, 0.995), na.rm = TRUE))
sr_rng <- c(max(0, sr_rng[1] - 0.01), min(1, sr_rng[2] + 0.01))
sr_rng[1] <- min(sr_rng[1], 0.5)
sr_rng[2] <- max(sr_rng[2], 0.5)
label_x <- sr_rng[2] + 0.02 * diff(sr_rng)

# ---- Common theme----
theme_common <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10, colour = "black"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    plot.margin = margin(6, 10, 6, 6)
  )

# ==========================================================
# (a) Global pooled distribution (ridgeline by age)
# ==========================================================
anno_a <- dat_plot |>
  group_by(age_lab) |>
  summarise(n_countries = n_distinct(Country), n_obs = n(), .groups = "drop") |>
  mutate(label = paste0("(", n_countries, ", ", n_obs, ")"))

p_a <- ggplot(dat_plot, aes(x = SR_c, y = age_lab, height = ..density..,
                            group = age_lab, fill = ..x..)) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = FALSE,
    scale = 1.8,
    rel_min_height = 0.001
  ) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_text(
    data = anno_a,
    aes(x = label_x, y = age_lab, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 0.5, size = 3.2
  ) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01), limits = sr_rng)+
  scale_fill_gradientn(colors = fill_pal, name = "SR", limits = sr_rng) +
  labs(title = "(a)", x = "Sex ratio", y = "Age group") +
  coord_cartesian(clip = "off") +
  theme_common +
  theme(plot.margin = margin(6, 55, 6, 6)) +
  guides(fill = "none")  # remove SR legend

# ==========================================================
# (b) Extreme 5% countries (within-country age ridgelines; all years pooled)
# ==========================================================
cty_mean <- dat_plot |>
  group_by(Country) |>
  summarise(mean_sr = mean(SR_c, na.rm = TRUE), .groups = "drop")

q_cut <- quantile(cty_mean$mean_sr, probs = c(0.025, 0.975), na.rm = TRUE)

extreme_cty <- cty_mean |>
  filter(mean_sr <= q_cut[1] | mean_sr >= q_cut[2]) |>
  arrange(mean_sr) |>
  pull(Country)

dat_ext <- dat_plot |>
  filter(Country %in% extreme_cty) |>
  mutate(Country_f = factor(Country, levels = extreme_cty))

p_b <- ggplot(dat_ext, aes(x = SR_c, y = age_lab, height = ..density..,
                           group = age_lab, fill = ..x..)) +
  stat_density_ridges(
    geom = "density_ridges_gradient",
    calc_ecdf = FALSE,
    scale = 1.35,
    rel_min_height = 0.01
  ) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  facet_wrap(~ Country_f, ncol = 4) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01),limits = sr_rng)+
  scale_fill_gradientn(colors = fill_pal, name = "SR", limits = sr_rng) +
  labs(title = "(c)", x = "Sex ratio", y = "Age group") +
  coord_cartesian(clip = "off") +
  theme_common +
  theme(strip.text = element_text(face = "bold", size = 9)) +
  guides(fill = "none")  # remove SR legend

# ==========================================================
# (c) 5-year boxplots + within-period age-gradient line (median)
# ==========================================================
bin5 <- function(y) floor(y/5) * 5
max_year <- max(dat_plot$Year, na.rm = TRUE)

dat5 <- dat_plot |>
  mutate(
    Year5_start = bin5(Year),
    Year5_lab = if_else(Year5_start + 4 >= max_year,
                        paste0(Year5_start, "-"),
                        paste0(Year5_start, "-", Year5_start + 4))
  )

x_breaks <- dat5 |> distinct(Year5_start) |> arrange(Year5_start) |> pull(Year5_start)
x_labels <- dat5 |> distinct(Year5_start, Year5_lab) |> arrange(Year5_start) |> pull(Year5_lab)

trend5 <- dat5 |>
  group_by(Year5_start, Year5_lab, age_lab) |>
  summarise(med = median(SR_c, na.rm = TRUE), .groups = "drop")

dodge_w <- 2.6
age_offset_tbl <- tibble(
  age_lab = factor(age_labels, levels = age_labels, ordered = TRUE),
  offset = c(-0.975, -0.325, 0.325, 0.975)
)

trend5 <- trend5 |>
  left_join(age_offset_tbl, by = "age_lab") |>
  mutate(x_dodged = Year5_start + offset)

p_c <- ggplot(dat5, aes(x = Year5_start, y = SR_c,
                        color = age_lab, fill = age_lab,
                        group = interaction(Year5_start, age_lab))) +
  geom_boxplot(
    position = position_dodge(width = dodge_w),
    outlier.size = 0.6,
    alpha = 0.15
  ) +
  geom_line(
    data = trend5,
    aes(x = x_dodged, y = med, group = Year5_start),
    inherit.aes = FALSE,
    color = "black", alpha = 0.45,
    linewidth = 0.6
  ) +
  geom_point(
    data = trend5,
    aes(x = x_dodged, y = med),
    inherit.aes = FALSE,
    color = "black", alpha = 0.6,
    size = 1.2
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01),limits = sr_rng)+
  labs(title = "(b)", x = "Years", y = "Sex ratio",
       color = "Age group", fill = "Age group") +
  theme_common +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

# ==========================================================
# Combine:
# ==========================================================
top_row <- (p_a | p_c) + plot_layout(widths = c(1, 4))

fig_all <- (top_row / p_b) +
  plot_layout(heights = c(1, 3)) &
  theme(legend.position = "right")

print(fig_all)


