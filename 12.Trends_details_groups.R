library(tidyverse)
library(castr)
library(stlplus)
library(nlme)
library(trend)
library(broom)
library(ggstats)

# Charger les données
zoo <- read_tsv("data/zoo_abundances.tsv.gz")

# Remove useless stuffs
zz <-zoo |> filter(!str_detect(label,"bubble"),
                   !str_detect(label,"artefact"),
                   !str_detect(label,"badfocus<artefact"),
                   !str_detect(label,"detritus"),
                   !str_detect(label,"^part"))


# Remove problematic date
zz <- zz %>% filter(!object_date=="1971-05-19")

# Lists of taxo groupings
copepoda_list <- c("Acartiidae", "Calanidae", "Calanoida", "Candaciidae","Centropagidae",
                   "Corycaeidae","Harpacticoida","Metridinidae","Oncaeidae","Oithonidae","Temoridae") 

gelatinouspred_list <- c("Aglaura", "Cnidaria", "Diphyidae", "Abylidae","Hydrozoa",
                         "Rhopalonema velatum") 
gelatinousfil_list <- c("Oikopleuridae","Fritillariidae","Doliolida","Salpida")
#rhizaria_list <- c("Phaeodaria","Foraminifera","Harosa")

# Add a taxo column
zz_group <- zz %>%
  mutate(taxoo = if_else(label %in% copepoda_list, "Copepoda", label),
         taxoo = if_else(taxoo %in% gelatinouspred_list, "Gelatpred", taxoo),
         taxoo = if_else(taxoo %in% gelatinousfil_list, "Gelatfilt", taxoo))
         #taxoo = if_else(taxoo %in% rhizaria_list, "Rhizaria", taxoo))
levels(as.factor(zz_group$taxoo))

# Calculation of total zoo conc per date
group_total <- zz_group %>% 
  group_by(object_date, taxoo) %>% 
  summarise(conc = sum(concentration, na.rm = TRUE)) %>%
  arrange(object_date)
# Change data table
abundance_ts <- group_total %>%
  rename(date = object_date, name = taxoo, value = conc)


##### Intermediate step : replace NA per 0 #####
#Replace plankton NA per 0
#crossing of all the possibilities

# Change data table
abundance_ts <- group_total %>%
  rename(date = object_date, name = taxoo, value = conc)

all_dates <- unique(abundance_ts$date)
all_taxa  <- unique(abundance_ts$name)

complete_grid <- crossing(date = all_dates, name = all_taxa)

abundance_complete <- complete_grid %>%
  left_join(abundance_ts, by = c("date", "name")) %>%
  mutate(value = replace_na(value, 0))  # replace NA by 0


abundance_ts<- abundance_complete %>% filter(name %in% c("Phaeodaria","Foraminifera", "Harosa"))
#First step: Regularization

# Define a sequence of dates 
ref <- tibble(
  target_date=seq(from=as.Date("1967-01-05"), to=as.Date("2022-12-31"), by=14),
  year=year(target_date)
)

# identify years in which the number of obs is larger than usual
pbs <- ref %>%
  count(year) %>%
  filter(n>26)

# ->this is often an extra date in very late decembre => just remove it
ref <- filter(ref, !(year %in% pbs$year & yday(target_date) > 360))
ref %>%
  count(year) %>%
  filter(n>26)
# -> all OK

# Match data based on these reference dates
avail_dates <- unique(abundance_ts$date)
ref <- ref %>%
  mutate(
    closest_date = castr::closest(ref$target_date, avail_dates),
    date_diff = abs(closest_date - target_date) %>% as.numeric()
  )

# Insert the data based on the matched dates
table_reg <- left_join(ref, abundance_ts, by=c("closest_date"="date"), relationship="many-to-many")

# erase data for matches that are too far from the target
table_reg_stl<- table_reg %>%
  mutate(value = if_else(date_diff > 7, NA, value))
#table_reg_stl<-table_reg_stl%>%filter(year>"1999-05-26")
table_reg_stl <- table_reg_stl%>%mutate(value=log10(value+1))
dstl <- table_reg_stl %>%
  group_by(name) %>%
  group_modify(.f=function(x,y) {
    # message(y)
    # if all is missing, do not do anything
    if ( all(is.na(x$value)) ) {
      out <- NULL
      # else perform stl
    } else {
      dec <- stlplus(x$value, x$target_date, n.p=26, s.window="periodic", t.window=10*26)
      out <- dec$data |> select(raw, seasonal, trend, remainder)
    }
    out <- bind_cols(select(x, target_date:date_diff), out)
  }) |>
  ungroup() |>
  
  # cut the part before the variable becomes available for the first time
  group_by(name) |>
  group_modify(.f=function(x, y) {
    if (all(is.na(x$raw))) {
      out <- NULL
    } else {
      x <- arrange(x, target_date)
      start_idx <- min(which(!is.na(x$raw)))
      out <- x[start_idx:nrow(x),]
    }
    return(out)
  }) |>
  ungroup()

# Plot the result
dstl %>%
  pivot_longer(raw:remainder, names_to="component") |>
  mutate(component=factor(component, levels=c("raw", "trend", "seasonal", "remainder"))) |>
  ggplot() + facet_wrap(~interaction(name, component), scale="free", nrow=4) +
  geom_path(aes(x=target_date, y=value))




## GLS regression ----

library("broom")
library("nlme")

glance.gls <- function(m) {
  s <- summary(m)
  
  # r.squared
  f <- predict(m)
  mss <- sum((f - mean(f))^2)
  rss <- sum(residuals(m)^2)
  rsq <- mss / (mss + rss)
  
  # residuals
  shap <- shapiro.test(m$residuals)
  a <- pacf(residuals(m, type="normalized"), plot=FALSE)
  
  tibble(
    r.squared = rsq,
    
    statistic = s$tTable[2, "t-value"],
    p.value = s$tTable[2, "p-value"],
    
    intercept = m$coefficients[1],
    slope = m$coefficients[2],
    
    shapiro.p.value = shap$p.value,
    cor.struct=class(m$modelStruct$corStruct)[1],
    acf1 = a$acf[1],
    acf2 = a$acf[2]
  )
}

# compute all trends
stats <- dstl %>%
  mutate(deseason = trend+remainder) |>
  group_by(name) |>
  group_modify(.f=function(x, y) {
    # message(y)
    if (all(is.na(x$raw))) {
      return(data.frame())
    }
    # 0. fill missing values through linear interpolation
    x$deseason_filled <- castr::interpolate(x=x$target_date, y=x$deseason, xout=x$target_date)
    
    # return(x)
    # 1. Mann-Kendall trend test
    mkt <- trend::mk.test(x$deseason_filled)
    
    # 2. GLS regression
    # simple model
    m <- gls(deseason_filled ~ target_date, data=x)
    a <- pacf(residuals(m, type="normalized"), plot=FALSE)
    # if autocorrelation is too strong
    if (abs(a$acf[1]) > 0.2) {
      # add AR1 model on residuals
      m <- gls(deseason_filled ~ target_date, data=x, cor=corAR1(round(a$acf[1], 1)))
      a <- pacf(residuals(m, type="normalized"), plot=FALSE)
      
      
    }
    
    # extract diagnostic information for both approaches
    bind_cols(
      glance(mkt) |> select(mk_p.value=p.value),
      glance.gls(m) |> select(r.squared, p.value=p.value, intercept, slope, cor.struct, acf=acf1) |> rename_with(function(n) {str_c("gls_", n)})
    )
  }) |>
  ungroup()



# display the result
dstl<- dstl%>%mutate(year=year(target_date))
# add date range
start_stop <- dstl %>%
  group_by(name) %>%
  summarise(start=min(year), end=max(year))

# add significance stars
# mk signif + gls non-signif may mean a non linear trend
signif_stars <- function(x) {
  case_when(
    x < 0.001 ~ "***",
    x < 0.01  ~ "**",
    x < 0.05  ~ "*",
    x < 0.1   ~ ".",
    TRUE      ~ ""
  )
}
stats <- stats %>%
  mutate(
    gls_acf = abs(gls_acf),
    mk_signif = signif_stars(mk_p.value),
    gls_signif = signif_stars(gls_p.value)
  ) %>%
  left_join(start_stop)

#save it
dstl<- dstl%>%mutate(deseason=trend+remainder)


##### Plot #####


ggplot() +
  facet_wrap(name ~ ., scales = "free_y", ncol = 3) +
  geom_path(
    aes(target_date, deseason),
    data = dstl,
    colour = "grey20"
  ) +
  geom_abline(
    aes(slope = gls_slope, intercept = gls_intercept),
    data = stats %>%
      filter(
               gls_signif %in% c("*", "**", "***")),
    colour = "red", size = 1.5, alpha = 0.8
  ) +
  xlab("Date") +
  ylab("Deseasonalized component") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.ticks.length = unit(0.3, "cm"),
    strip.background = element_rect(fill = "#FFDAB9")
  )
