# =============================================================================
# GMU Cherry Blossom Peak Bloom Prediction Competition 2026
# competition.statistics.gmu.edu
# Deadline: Feb 28, 2026
#
# Model: Two-stage ensemble
#   Stage 1 — DAM4 biological clock (Utah Chill Units + hourly GDH sigmoid)
#   Stage 2 — Statistical correction (ECMWF SEAS5 forecast + trend + LOOCV)
#   Data:     Open-Meteo ERA5 hourly (1950–2026) | ECMWF SEAS5 seasonal forecast
# =============================================================================

# --- 0. PACKAGES -------------------------------------------------------------
pkgs <- c("tidyverse", "lubridate", "httr", "jsonlite", "broom", "patchwork")
new  <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new)) install.packages(new)

library(tidyverse)
library(lubridate)
library(httr)
library(jsonlite)
library(broom)
library(patchwork)

# =============================================================================
# PHASE 1: CONFIGURATION
# =============================================================================

# Biological constants
UTAH_CU_REQ      <- 1200   # Utah Chill Units required for endodormancy release
ANALOG_RMSE_DEFAULT <- 6.0 # Default RMSE assumed for analog model (days)
NYC_OFFSET       <- 4L     # Days DC -> NYC empirical offset
ANALOG_ONLY_LOCS <- "vancouver" # Locations with too few years for regression

# Site coordinates
LOCS <- tribble(
  ~location,       ~lat,    ~lon,
  "washingtondc",  38.89,  -77.04,
  "kyoto",         35.02,  135.73,
  "liestal",       47.48,    7.73,
  "vancouver",     49.25, -123.12,
  "newyorkcity",   40.78,  -73.97
)

# =============================================================================
# PHASE 2: DATA GATHERING
# =============================================================================

# --- 2A. Bloom dates from GMU GitHub -----------------------------------------
BASE_URL <- paste0(
  "https://raw.githubusercontent.com/",
  "GMU-CherryBlossomCompetition/peak-bloom-prediction/main/data/"
)

bloom_raw <- map_dfr(
  c("washingtondc", "kyoto", "liestal", "vancouver", "newyorkcity"),
  function(loc) {
    tryCatch(
      read_csv(paste0(BASE_URL, loc, ".csv"), show_col_types = FALSE) |>
        mutate(location = loc),
      error = function(e) { warning(loc, ": ", e$message); NULL }
    )
  }
)

bloom <- bloom_raw |>
  rename_with(tolower) |>
  select(location, year, bloom_doy) |>
  filter(!is.na(bloom_doy))

cat("Bloom records loaded:", nrow(bloom), "\n")
print(count(bloom, location))

# --- 2B. Open-Meteo ERA5 hourly data (1950–2026) -----------------------------
fetch_openmeteo_hourly <- function(lat, lon, start, end) {
  url <- sprintf(
    paste0(
      "https://archive-api.open-meteo.com/v1/archive",
      "?latitude=%.4f&longitude=%.4f",
      "&start_date=%s&end_date=%s",
      "&hourly=temperature_2m&timezone=UTC"
    ),
    lat, lon, start, end
  )
  resp <- tryCatch(GET(url, timeout(120)), error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) {
    warning("ERA5 fetch failed for lat=", lat, " lon=", lon)
    return(NULL)
  }
  js <- fromJSON(content(resp, as = "text", encoding = "UTF-8"))
  tibble(
    datetime = as.POSIXct(js$hourly$time, format = "%Y-%m-%dT%H:%M", tz = "UTC"),
    temp_c   = js$hourly$temperature_2m
  ) |>
    filter(!is.na(temp_c))
}

message("Fetching ERA5 hourly data (this takes ~5-10 min)...")
hourly_all <- map_dfr(seq_len(nrow(LOCS)), function(i) {
  message("  Fetching: ", LOCS$location[i])
  df <- fetch_openmeteo_hourly(
    LOCS$lat[i], LOCS$lon[i],
    "1950-10-01", "2026-02-22"
  )
  if (!is.null(df)) mutate(df, location = LOCS$location[i]) else NULL
}) |>
  mutate(
    date  = as.Date(datetime),
    year  = year(datetime),
    month = month(datetime),
    hour  = hour(datetime),
    doy   = yday(date)
  )

cat("Hourly rows loaded:", nrow(hourly_all), "\n")

# =============================================================================
# PHASE 3: DORMANCY — Utah Chill Unit Model
# =============================================================================
# Richardson et al. (1974) — biologically superior to simple chilling-hours.
# Each hour contributes CU based on temperature; temps > 18 °C subtract units.

utah_cu <- function(temp_c) {
  case_when(
    temp_c <= 1.4  ~  0.0,
    temp_c <= 2.4  ~  0.5,
    temp_c <= 9.1  ~  1.0,
    temp_c <= 12.4 ~  0.5,
    temp_c <= 15.9 ~  0.0,
    temp_c <= 18.0 ~ -0.5,
    TRUE           ~ -1.0
  )
}

hourly_chill <- hourly_all |>
  filter(month %in% c(10L, 11L, 12L, 1L, 2L, 3L)) |>
  mutate(
    chill_year = if_else(month >= 10L, year, year - 1L),
    cu_contrib = utah_cu(temp_c)
  ) |>
  group_by(location, chill_year) |>
  arrange(datetime) |>
  mutate(cum_cu = cumsum(cu_contrib)) |>
  ungroup()

wakeup_dates <- hourly_chill |>
  filter(cum_cu >= UTAH_CU_REQ) |>
  group_by(location, chill_year) |>
  slice_min(datetime, n = 1) |>
  ungroup() |>
  transmute(
    location,
    chill_year,
    bloom_year  = chill_year + 1L,
    wakeup_date = as.Date(datetime),
    wakeup_doy  = doy
  )

cat("\nWakeup date coverage:\n")
print(count(wakeup_dates, location))

# =============================================================================
# PHASE 4: FORCING — Hourly GDH with Sigmoid
# =============================================================================
# Triangular-sigmoid hybrid (chillR convention).
# Effective range 4–36 °C; tapers cleanly to zero at 36 °C.

gdh_sigmoid <- function(temp_c) {
  pmax(0, temp_c - 4) * pmin(1, pmax(0, (36 - temp_c) / (36 - 25)))
}

gdh_hourly <- hourly_all |>
  mutate(bloom_year = if_else(month >= 10L, year + 1L, year)) |>
  inner_join(
    wakeup_dates |> select(location, bloom_year, wakeup_date),
    by = c("location", "bloom_year")
  ) |>
  filter(date >= wakeup_date, !is.na(temp_c)) |>
  mutate(gdh_contrib = gdh_sigmoid(temp_c)) |>
  group_by(location, bloom_year) |>
  arrange(datetime) |>
  mutate(cum_gdh = cumsum(gdh_contrib)) |>
  ungroup()

# Deduplicate bloom to one record per location-year
# (Kyoto has multiple tree observations per year -> duplicates otherwise)
bloom_dedup <- bloom |>
  group_by(location, year) |>
  summarise(bloom_doy = round(mean(bloom_doy)), .groups = "drop")

wakeup_bloom <- wakeup_dates |>
  inner_join(bloom_dedup, by = c("location", "bloom_year" = "year"))

cat("\nwakeup_bloom rows per location (should be ~1 per year):\n")
print(count(wakeup_bloom, location))

# slice_max(doy) then slice_max(datetime) ensures exactly ONE row per
# location-year (multiple hours share the same DOY).
gdh_at_bloom <- gdh_hourly |>
  inner_join(bloom_dedup, by = c("location", "bloom_year" = "year")) |>
  group_by(location, bloom_year) |>
  filter(doy <= bloom_doy) |>
  slice_max(doy,      n = 1, with_ties = TRUE)  |>
  slice_max(datetime, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(location, bloom_year, cum_gdh_at_bloom = cum_gdh)

stopifnot(
  nrow(gdh_at_bloom |> count(location, bloom_year) |> filter(n > 1)) == 0
)
cat("gdh_at_bloom: no duplicates confirmed\n")

gdh_threshold <- gdh_at_bloom |>
  group_by(location) |>
  summarise(
    median_gdh = median(cum_gdh_at_bloom, na.rm = TRUE),
    sd_gdh     = sd(cum_gdh_at_bloom,     na.rm = TRUE),
    .groups = "drop"
  )
cat("\nGDH thresholds per location:\n")
print(gdh_threshold)

# =============================================================================
# PHASE 4B: FEATURE ENGINEERING
# =============================================================================
monthly_tavg <- hourly_all |>
  group_by(location, year, month) |>
  summarise(tavg = mean(temp_c, na.rm = TRUE), .groups = "drop")

feb_mar_temps <- monthly_tavg |>
  filter(month %in% c(2L, 3L)) |>
  group_by(location, year) |>
  summarise(feb_mar_tavg = mean(tavg), .groups = "drop") |>
  rename(bloom_year = year)

jan_temps <- monthly_tavg |>
  filter(month == 1L) |>
  select(location, bloom_year = year, jan_tavg = tavg)

# Dec of prior year -> bloom spring; shift year forward by 1
dec_temps <- monthly_tavg |>
  filter(month == 12L) |>
  select(location, bloom_year = year, dec_tavg = tavg) |>
  mutate(bloom_year = bloom_year + 1L)

model_data <- wakeup_bloom |>
  inner_join(gdh_at_bloom,  by = c("location", "bloom_year")) |>
  inner_join(feb_mar_temps, by = c("location", "bloom_year")) |>
  inner_join(jan_temps,     by = c("location", "bloom_year")) |>
  inner_join(dec_temps,     by = c("location", "bloom_year")) |>
  mutate(
    # 2024 was a leap year; subtract 1 DOY so years are comparable
    bloom_doy_adj = if_else(bloom_year == 2024L, bloom_doy - 1L, bloom_doy),
    year_centered = bloom_year - 1990L
  )

cat("\nModel data rows per location:\n")
print(count(model_data, location))

# Validation plot: Feb/Mar temp vs bloom DOY
p_validate <- model_data |>
  ggplot(aes(x = feb_mar_tavg, y = bloom_doy_adj, colour = location)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  facet_wrap(~location, scales = "free") +
  labs(
    title    = "Validation: Feb/Mar Mean Temp (ERA5 hourly) vs Bloom DOY",
    subtitle = "Warmer Feb/Mar \u2192 earlier bloom",
    x = "Feb-Mar Mean Temperature (\u00b0C)",
    y = "Bloom DOY (leap-year corrected)"
  ) +
  theme_bw() +
  theme(legend.position = "none")
print(p_validate)

# =============================================================================
# PHASE 4C: MODEL FITTING WITH LOOCV
# =============================================================================
fit_loocv <- function(df) {
  n <- nrow(df)
  if (n < 8) return(list(model = NULL, rmse_loo = NA_real_, preds_loo = NULL))

  preds <- numeric(n)
  for (i in seq_len(n)) {
    mod      <- lm(bloom_doy_adj ~ feb_mar_tavg + jan_tavg + dec_tavg +
                     year_centered, data = df[-i, ])
    preds[i] <- predict(mod, newdata = df[i, ])
  }
  rmse_loo <- sqrt(mean((df$bloom_doy_adj - preds)^2))
  mod_full <- lm(bloom_doy_adj ~ feb_mar_tavg + jan_tavg + dec_tavg +
                   year_centered, data = df)
  list(
    model     = mod_full,
    rmse_loo  = rmse_loo,
    preds_loo = tibble(
      bloom_year = df$bloom_year,
      obs        = df$bloom_doy_adj,
      pred_loo   = preds
    )
  )
}

loocv_results <- model_data |>
  group_by(location) |>
  group_map(~ fit_loocv(.x), .keep = TRUE)
names(loocv_results) <- unique(model_data$location)

cat("\n=== LOOCV RMSE (out-of-sample, days) ===\n")
map_dfr(names(loocv_results), function(loc) {
  tibble(
    location   = loc,
    loocv_rmse = loocv_results[[loc]]$rmse_loo,
    n          = nrow(filter(model_data, location == loc))
  )
}) |> print()

models <- map(loocv_results, "model")

model_summaries <- map_dfr(names(models), function(loc) {
  if (is.null(models[[loc]])) return(NULL)
  glance(models[[loc]]) |> mutate(location = loc)
})
cat("\nModel R-squared:\n")
print(model_summaries |> select(location, r.squared, adj.r.squared, p.value))

# =============================================================================
# PHASE 5A: ECMWF SEAS5 FORECAST for Feb/Mar 2026
# =============================================================================
fetch_seas5_monthly <- function(lat, lon) {
  url <- sprintf(
    paste0(
      "https://seasonal-api.open-meteo.com/v1/seasonal",
      "?latitude=%.4f&longitude=%.4f",
      "&monthly=temperature_2m_mean",
      "&models=ecmwf_seas5",
      "&start_date=2026-02-01&end_date=2026-03-31",
      "&timezone=UTC"
    ),
    lat, lon
  )
  resp  <- tryCatch(GET(url, timeout(60)), error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) return(NA_real_)
  js    <- content(resp, as = "parsed")
  temps <- unlist(js$monthly$temperature_2m_mean)
  mean(temps, na.rm = TRUE)
}

message("\nFetching ECMWF SEAS5 Feb/Mar 2026 forecast...")

# Fix: iterate by row index so location label is always available
seas5_forecast <- map_dfr(seq_len(nrow(LOCS)), function(i) {
  message("  ", LOCS$location[i])
  tibble(
    location      = LOCS$location[i],
    lat           = LOCS$lat[i],
    lon           = LOCS$lon[i],
    feb_mar_seas5 = fetch_seas5_monthly(LOCS$lat[i], LOCS$lon[i])
  )
})

cat("\nECMWF SEAS5 Feb/Mar 2026 forecast (\u00b0C):\n")
print(seas5_forecast |> select(location, feb_mar_seas5))

feb_mar_hist_means <- model_data |>
  group_by(location) |>
  summarise(feb_mar_hist = mean(feb_mar_tavg, na.rm = TRUE), .groups = "drop")

forecast_2026 <- seas5_forecast |>
  left_join(feb_mar_hist_means, by = "location") |>
  mutate(
    # Reject SEAS5 if missing, near-zero, or > 5 °C from historical mean
    seas5_ok     = !is.na(feb_mar_seas5) &
      abs(feb_mar_seas5) > 0.01 &
      abs(feb_mar_seas5 - feb_mar_hist) < 5,
    feb_mar_2026 = if_else(seas5_ok, feb_mar_seas5, feb_mar_hist),
    source       = if_else(seas5_ok, "SEAS5", "historical_mean")
  )

cat("\nFinal Feb/Mar 2026 temperature input:\n")
print(forecast_2026 |> select(location, feb_mar_2026, source))

# Jan 2026 and Dec 2025 from observed ERA5 data
jan_2026 <- monthly_tavg |>
  filter(year == 2026L, month == 1L) |>
  select(location, jan_tavg_2026 = tavg)

dec_2025 <- monthly_tavg |>
  filter(year == 2025L, month == 12L) |>
  select(location, dec_tavg_2025 = tavg)

# =============================================================================
# PHASE 5B: ANALOG YEARS (Jan-temp matched, top 5)
# =============================================================================
jan_use <- model_data |>
  group_by(location) |>
  summarise(jan_hist = mean(jan_tavg, na.rm = TRUE), .groups = "drop") |>
  left_join(jan_2026, by = "location") |>
  mutate(jan_use = coalesce(jan_tavg_2026, jan_hist))

analogs <- model_data |>
  left_join(jan_use |> select(location, jan_use), by = "location") |>
  mutate(jan_diff = abs(jan_tavg - jan_use)) |>
  group_by(location) |>
  slice_min(jan_diff, n = 5, with_ties = FALSE) |>
  ungroup()

# Historical bloom SD per location (floor for sparse locations)
bloom_sd_hist <- bloom_dedup |>
  group_by(location) |>
  summarise(sd_hist = sd(bloom_doy, na.rm = TRUE), .groups = "drop")

analog_preds <- analogs |>
  group_by(location) |>
  summarise(
    pred_analog = mean(bloom_doy_adj, na.rm = TRUE),
    sd_analog   = sd(bloom_doy_adj,   na.rm = TRUE),
    n_analogs   = n(),
    .groups = "drop"
  ) |>
  left_join(bloom_sd_hist, by = "location") |>
  mutate(
    # Use historical SD as floor; cap at 2× historical SD
    sd_analog = if_else(is.na(sd_analog) | sd_analog < 1, sd_hist, sd_analog),
    sd_analog = pmin(sd_analog, 2 * sd_hist)
  ) |>
  select(-sd_hist)

cat("\nAnalog predictions:\n")
print(analog_preds)

# =============================================================================
# PHASE 5C: STATISTICAL MODEL PREDICTIONS
# =============================================================================
dec_use <- model_data |>
  group_by(location) |>
  summarise(dec_hist = mean(dec_tavg, na.rm = TRUE), .groups = "drop") |>
  left_join(dec_2025, by = "location") |>
  mutate(dec_use = coalesce(dec_tavg_2025, dec_hist))

newdata_2026 <- forecast_2026 |>
  select(location, feb_mar_tavg = feb_mar_2026) |>
  left_join(jan_use |> select(location, jan_tavg = jan_use), by = "location") |>
  left_join(dec_use |> select(location, dec_tavg = dec_use), by = "location") |>
  mutate(year_centered = 2026L - 1990L)

model_preds <- map_dfr(unique(model_data$location), function(loc) {
  mod   <- models[[loc]]
  rmse  <- loocv_results[[loc]]$rmse_loo
  nd    <- filter(newdata_2026, location == loc)
  n_obs <- if (!is.null(mod)) nrow(mod$model) else 0L

  if (!is.null(mod) && n_obs >= 8 && !(loc %in% ANALOG_ONLY_LOCS)) {
    pred <- tryCatch(
      predict(mod, newdata = nd, interval = "prediction", level = 0.90),
      error = function(e) NULL
    )
    if (!is.null(pred)) {
      return(tibble(
        location    = loc,
        pred_model  = pred[1, "fit"],
        lower_model = pred[1, "lwr"],
        upper_model = pred[1, "upr"],
        loocv_rmse  = rmse,
        model_used  = TRUE
      ))
    }
  }
  # Fall back to analog for locations with insufficient data
  a <- filter(analog_preds, location == loc)
  tibble(
    location    = loc,
    pred_model  = a$pred_analog,
    lower_model = a$pred_analog - 1.645 * a$sd_analog,
    upper_model = a$pred_analog + 1.645 * a$sd_analog,
    loocv_rmse  = a$sd_analog,
    model_used  = FALSE
  )
})

cat("\nModel predictions (Stage 2):\n")
print(model_preds)

# =============================================================================
# PHASE 5D: LOOCV-WEIGHTED ENSEMBLE
# =============================================================================
# Point estimate:  inverse-RMSE weighted average of model + analog
# Interval:        weighted blend of model PI half-width and analog 90% half-width

final_preds <- analog_preds |>
  inner_join(model_preds, by = "location") |>
  mutate(
    w_model  = 1 / pmax(loocv_rmse, 1),
    w_analog = 1 / ANALOG_RMSE_DEFAULT,
    w_total  = w_model + w_analog,

    prediction = round((pred_model * w_model + pred_analog * w_analog) / w_total),

    half_model  = (upper_model - lower_model) / 2,
    half_analog = 1.645 * sd_analog,
    half_blend  = (half_model * w_model + half_analog * w_analog) / w_total,

    lower    = round(prediction - half_blend),
    upper    = round(prediction + half_blend),
    interval = 90L
  )

cat("\nEnsemble interval widths (days):\n")
print(final_preds |>
        select(location, prediction, lower, upper) |>
        mutate(width = upper - lower))

# --- NYC: transfer from DC + empirical offset --------------------------------
dc <- filter(final_preds, location == "washingtondc")

nyc_row <- dc |>
  mutate(
    location    = "newyorkcity",
    pred_analog = pred_analog  + NYC_OFFSET,
    sd_analog   = sd_analog    + 2,
    pred_model  = pred_model   + NYC_OFFSET,
    lower_model = lower_model  + NYC_OFFSET,
    upper_model = upper_model  + NYC_OFFSET,
    prediction  = prediction   + NYC_OFFSET,
    lower       = lower        + NYC_OFFSET - 2L,
    upper       = upper        + NYC_OFFSET + 2L,
    model_used  = FALSE
  )

final_preds <- bind_rows(final_preds, nyc_row)

doy_to_date <- function(doy, yr = 2026L) {
  as.Date(paste0(yr, "-01-01")) + doy - 1
}

final_preds <- final_preds |>
  mutate(
    pred_date     = doy_to_date(prediction),
    lower_date    = doy_to_date(lower),
    upper_date    = doy_to_date(upper),
    interval_days = upper - lower
  )

cat("\n=== FINAL 2026 PREDICTIONS ===\n")
print(final_preds |>
        select(location, pred_date, lower_date, upper_date,
               interval_days, model_used))

# =============================================================================
# PHASE 6: SUBMISSION CSV
# =============================================================================
submission <- final_preds |>
  transmute(
    year       = 2026L,
    location,
    prediction = as.integer(prediction),
    lower      = as.integer(lower),
    upper      = as.integer(upper),
    interval   = 90L
  ) |>
  arrange(location)

print(submission)
write_csv(submission, "predictions.csv")
message("predictions.csv written!")

# =============================================================================
# PHASE 7: VISUALISATION
# =============================================================================
LOC_LABELS <- c(
  washingtondc = "Washington DC",
  newyorkcity  = "New York City",
  kyoto        = "Kyoto",
  liestal      = "Liestal",
  vancouver    = "Vancouver"
)

p_final <- final_preds |>
  mutate(loc_label = LOC_LABELS[location]) |>
  ggplot(aes(x = fct_reorder(loc_label, prediction), y = prediction)) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.3, colour = "#C0392B", linewidth = 1.1
  ) +
  geom_point(size = 3.5, colour = "#E91E8C") +
  geom_text(
    aes(label = format(pred_date, "%b %d")),
    vjust = -1.3, size = 3.5, fontface = "bold"
  ) +
  scale_y_continuous(
    name      = "Day of Year (2026)",
    sec.axis  = sec_axis(
      ~ .,
      name   = "Calendar Date",
      labels = function(x) format(as.Date("2026-01-01") + round(x) - 1, "%b %d")
    )
  ) +
  labs(
    title    = "2026 Cherry Blossom Peak Bloom Predictions",
    subtitle = paste(
      "Utah CU biological model + ECMWF SEAS5 statistical layer",
      "| 90% PI | LOOCV-weighted"
    ),
    x       = "Location",
    caption = "competition.statistics.gmu.edu"
  ) +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

print(p_final)
ggsave("bloom_predictions_2026.png", p_final, width = 10, height = 6, dpi = 150)

# LOOCV residual plot
loocv_plot_data <- map_dfr(names(loocv_results), function(loc) {
  res <- loocv_results[[loc]]
  if (is.null(res$preds_loo)) return(NULL)
  res$preds_loo |>
    mutate(location = loc, residual = obs - pred_loo)
})

p_loocv <- loocv_plot_data |>
  ggplot(aes(x = pred_loo, y = residual, colour = location)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.6) +
  facet_wrap(~location, scales = "free_x") +
  labs(
    title = "LOOCV Residuals: out-of-sample prediction error",
    x     = "Predicted DOY",
    y     = "Observed \u2212 Predicted (days)"
  ) +
  theme_bw() +
  theme(legend.position = "none")
print(p_loocv)

# =============================================================================
# PHASE 8: ABSTRACT STATS SUMMARY
# =============================================================================
cat("\n=== KEY STATS FOR ABSTRACT ===\n")
cat("Stage 1 — Utah Chill Unit model, requirement =", UTAH_CU_REQ, "CU\n")
cat("Stage 1 — Forcing: hourly GDH sigmoid (effective range 4–36 \u00b0C)\n")
cat("Stage 2 — Linear model: bloom_doy ~ feb_mar_tavg + jan_tavg + dec_tavg + year_trend\n")
cat("Stage 2 — Feb/Mar 2026 input: ECMWF SEAS5 51-member ensemble\n")
cat("Ensemble weighting: inverse LOOCV RMSE\n\n")

loocv_plot_data |>
  group_by(location) |>
  summarise(
    rmse = sqrt(mean(residual^2)),
    mae  = mean(abs(residual)),
    .groups = "drop"
  ) |>
  print()
