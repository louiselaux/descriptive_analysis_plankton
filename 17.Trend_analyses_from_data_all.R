##### Plankton analyses #####

# Load the libraries

library(readr)
library(purrr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(trend)      
library(strucchange) 
library(stringr)
library(stlplus)
library(nlme)
library(broom)
library(ggstats)
library(vegan)
library(worrms)
library(tidyverse)


# load the file
tablo <- read_tsv("data/tablo.tsv")


# Significant stars 
signif_stars <- function(x) case_when(
  x < 0.001 ~ "***",
  x < 0.01  ~ "**",
  x < 0.05  ~ "*",
  x < 0.1   ~ ".",
  TRUE      ~ ""
)

# GLS summary
glance.gls <- function(m){
  s <- summary(m)
  f <- predict(m)
  mss <- sum((f - mean(f))^2)
  rss <- sum(residuals(m)^2)
  rsq <- mss/(mss + rss)
  shap <- shapiro.test(m$residuals)
  a <- pacf(residuals(m, type="normalized"), plot = FALSE)
  tibble(
    gls_r.squared = rsq,
    gls_statistic = s$tTable[2, "t-value"],
    gls_p.value   = s$tTable[2, "p-value"],
    gls_intercept = m$coefficients[1],
    gls_slope     = m$coefficients[2],
    gls_cor.struct= class(m$modelStruct$corStruct)[1],
    gls_acf       = a$acf[1]
  )
}

# Main function 
run_trends_for <- function(df_src, group_col, top_k = NULL, log10p1 = TRUE,
                           stl_period = 26, stl_twindow = 10*26) {
  
  gsym <- rlang::ensym(group_col)
  
  # Aggregate per date
  agg <- df_src %>%
    transmute(object_date = as.Date(object_date),
              name = !!gsym,
              value = concentration) %>%
    mutate(name = if_else(is.na(name) | name=="", "Unknown", as.character(name))) %>%
    group_by(object_date, name) %>%
    summarise(value = sum(value, na.rm = TRUE), .groups = "drop")
  
  # Put 0 if NA for one date 
  all_dates <- sort(unique(agg$object_date))
  all_names <- sort(unique(agg$name))
  abundance_complete <- tidyr::crossing(object_date = all_dates, name = all_names) %>%
    left_join(agg, by = c("object_date","name")) %>%
    mutate(value = tidyr::replace_na(value, 0))
  
  # If we want to look at the most abundant groups only
  if (!is.null(top_k)) {
    keep <- abundance_complete %>%
      group_by(name) %>%
      summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total)) %>%
      slice_head(n = top_k) %>%
      pull(name)
    abundance_complete <- abundance_complete %>% filter(name %in% keep)
  }
  
  # Regularize on a grid of 14 days
  ref <- tibble(
    target_date = seq(as.Date("1967-01-05"), as.Date("2022-12-31"), by = 14),
    year = year(target_date)
  )
  # Remove 27 th
  pbs <- ref %>% count(year) %>% filter(n > 26)
  ref <- filter(ref, !(year %in% pbs$year & yday(target_date) > 360))
  
  avail_dates <- unique(abundance_complete$object_date)
  ref <- ref %>%
    mutate(
      closest_date = castr::closest(ref$target_date, avail_dates),
      date_diff    = abs(closest_date - target_date) %>% as.numeric()
    )
  
  # Merge
  table_reg <- left_join(ref, abundance_complete,
                         by = c("closest_date" = "object_date"),
                         relationship = "many-to-many")
  
  # Keep only real date if in tolerance interval
  tol_days_keep <- 3
  
  table_reg <- table_reg %>%
    mutate(
      value_obs_keep = if_else(date_diff <= tol_days_keep, value, NA_real_)
    )
  
  #Interpolation in using all available dates
  interp_grid <- table_reg %>%
    arrange(name, target_date, closest_date) %>%
    group_by(name) %>%
    group_modify(~{
      x <- .x
      
      # 
      obs_all_raw <- x %>%
        select(closest_date, value) %>%
        filter(!is.na(closest_date), !is.na(value)) %>%
        arrange(closest_date)
      
      if (nrow(obs_all_raw) < 2) {
        # not enough to interpolate
        x$value_final <- x$value_obs_keep
        return(x)
      }
      
      # linear interpolation 
      xx <- as.numeric(obs_all_raw$closest_date) 
      yy <- obs_all_raw$value
      xout <- as.numeric(x$target_date)
      
      y_lin <- approx(
        x = xx, y = yy, xout = xout,
        method = "linear", ties = "ordered"  
      )$y
      
      # do not extrapolate out of observed dates
      in_range <- (x$target_date >= min(obs_all_raw$closest_date) &
                     x$target_date <= max(obs_all_raw$closest_date))
      y_lin[!in_range] <- NA_real_
      
      # take intrpolated values only if real values not available
      x$value_final <- if_else(!is.na(x$value_obs_keep), x$value_obs_keep,
                               if_else(in_range, y_lin, NA_real_))
      x
    }) %>%
    ungroup()
  
  
  #
  interp_grid <- if (log10p1) mutate(interp_grid, value_final = log10(value_final + 1)) else interp_grid
  
  # STL on interpolated series
  dstl <- interp_grid %>%
    group_by(name) %>%
    group_modify(~{
      x <- .x
      if (all(is.na(x$value_final))) return(tibble())
      dec <- stlplus(x$value_final, x$target_date, n.p = stl_period,
                     s.window = "periodic", t.window = stl_twindow)
      out <- dec$data %>% select(raw, seasonal, trend, remainder)
      bind_cols(select(x, target_date:date_diff, value_final), out)
    }) %>%
    ungroup() %>%
    group_by(name) %>%
    group_modify(~{
      x <- arrange(.x, target_date)
      if (all(is.na(x$raw))) return(tibble())
      i0 <- min(which(!is.na(x$raw)))
      x[i0:nrow(x), ]
    }) %>%
    ungroup() %>%
    mutate(year = year(target_date),
           deseason = trend + remainder)
  
  
  # Trends
  stats <- dstl %>%
    group_by(name) %>%
    group_modify(~{
      x <- .x
      if (all(is.na(x$deseason))) return(tibble())
      # Fill in the deseason
      x$deseason_filled <- castr::interpolate(x = x$target_date,
                                              y = x$deseason,
                                              xout = x$target_date)
      # MK
      mkt <- trend::mk.test(x$deseason_filled)
      # GLS 
      m  <- gls(deseason_filled ~ target_date, data = x)
      a  <- pacf(residuals(m, type = "normalized"), plot = FALSE)
      if (abs(a$acf[1]) > 0.2) {
        m <- gls(deseason_filled ~ target_date, data = x,
                 correlation = corAR1(round(a$acf[1], 1)))
      }
      bind_cols(
        broom::glance(mkt) %>% select(mk_p.value = p.value),
        glance.gls(m)
      )
    }) %>%
    ungroup() %>%
    mutate(
      gls_acf    = abs(gls_acf),
      mk_signif  = signif_stars(mk_p.value),
      gls_signif = signif_stars(gls_p.value)
    )
  
  # Look at it
  p <- ggplot() +
    facet_wrap(name ~ ., scales = "free_y", ncol = 3) +
    geom_path(
      aes(target_date, deseason),
      data = dstl,
      colour = "grey20"
    ) +
    geom_abline(
      aes(slope = gls_slope, intercept = gls_intercept),
      data = stats %>% filter(gls_signif %in% c("*","**","***")),
      colour = "red", size = 1.2, alpha = 0.85
    ) +
    labs(x = "Date", y = "Deseasonalized",
         title = paste0("Gls trend ", rlang::as_name(gsym),
                        if (log10p1) " (log10+1)" else "")) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "#FFDAB9"),
      axis.text = element_text(size = 11),
      strip.text = element_text(size = 12)
    )
  
  #list(dstl = dstl, stats = stats, plot = p)
  list(dstl = dstl, stats = stats, plot = p, interp_grid = interp_grid)
  
}



##### Applicate it per group by order for example

res_order <- run_trends_for(
  df_src   = tablo %>% select(object_date, concentration, order_like),
  group_col= order_like,
  top_k    = 15
)
res_order$plot
res_order$stats %>% arrange(gls_p.value) %>% print(n=20)

##### Applicate if per functional group ######
res_fg <- run_trends_for(
  df_src   = tablo %>% select(object_date, concentration, functional_group),
  group_col= functional_group,
  top_k    = NULL         
)
res_fg$plot
res_fg$stats

# Label
res_fg <- run_trends_for(
  df_src   = tablo %>% select(object_date, concentration, label),
  group_col= label,
  top_k    = 20         
)
res_fg$plot
res_fg$stats


##### Plot to see if it seems consistent ######
plot_interp_with_raw <- function(res, names = NULL, start, end, log10p1 = TRUE) {
  stopifnot(!is.null(res$interp_grid))  # need to have the interpolation grid
  
  g <- res$interp_grid %>%
    mutate(kept = !is.na(value_obs_keep)) %>%
    filter(target_date >= as.Date(start), target_date <= as.Date(end))
  
  if (!is.null(names)) g <- g %>% filter(name %in% names)
  
  # Raw data at the real dates
  raw_pts <- g %>%
    select(name, closest_date, value) %>%
    filter(!is.na(closest_date), !is.na(value))
  
  # Same scaling as the line
  raw_pts <- raw_pts %>%
    mutate(value_plot = if (log10p1) log10(value + 1) else value)
  
  ggplot() +
    # 1) Interpolated data
    geom_line(data = g, aes(target_date, value_final)) +
    
    # 2) Target vs closest for kept points
    geom_segment(
      data = g %>% filter(kept),
      aes(x = pmin(target_date, closest_date),
          xend = pmax(target_date, closest_date),
          y = value_final, yend = value_final),
      inherit.aes = FALSE, alpha = 0.5
    ) +
    
    # 3) Red = kept points because in the tolerance interval
    geom_point(
      data = g %>% filter(kept),
      aes(target_date, value_final),
      shape = 16, size = 3, color = "red"
    ) +
    
    # 4) Points that were interpolated
    geom_point(
      data = g %>% filter(!kept & !is.na(value_final)),
      aes(target_date, value_final),
      shape = 4, size = 2, stroke = 1
    ) +
    
    # 5) Real points at the real available dates 
    geom_point(
      data = raw_pts,
      aes(closest_date, value_plot),
      shape = 1, size = 2, alpha = 0.6
    ) +
    
    facet_wrap(~ name, scales = "free_y") +
    labs(
      x = "Date",
      y = if (log10p1) "Valeur (log10(value+1))" else "Valeur",
      title = "Grille 14 j (—) + points bruts (○) + points gardés (●) + interpolés (×)",
      subtitle = paste(format(as.Date(start)), "→", format(as.Date(end)))
    ) +
    theme_bw()
}


# Applicate to something
plot_interp_with_raw(
  res_order,
  names = c("Calanoida","Cyclopoida"),
  start = "1995-01-01", end = "1996-12-31",
  log10p1 = TRUE
)



