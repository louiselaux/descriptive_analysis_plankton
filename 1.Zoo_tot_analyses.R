##### Plankton analyses #####

# Load the libraries

library(readr)
library(purrr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(trend)      
library(strucchange) 



# load the file
zoo <- read_tsv("data/zoo_abundances.tsv.gz")

# Calculation of total zoo conc per date
zoo_total <- zoo %>% 
  group_by(object_date) %>% 
  summarise(conc = sum(concentration, na.rm = TRUE)) %>%
  arrange(object_date)

# Look at it 
ggplot(zoo_total, aes(x = object_date, y = conc)) + 
  geom_line() + 
  scale_y_log10() + 
  labs(title = "Total concentration of zooplankton", y = "ind.m^-3") + 
  theme_minimal()

##### Step 1: Breaking Point analyses #####
pettitt_result <- pettitt.test(zoo_total$conc)
pettitt_result
change_point_date <- zoo_total$object_date[1172]

ggplot(zoo_total, aes(x = object_date, y = conc)) +
  geom_path(color = "black") +  # Série temporelle
  #geom_vline(xintercept = as.Date(change_point_date), linetype = "dashed", color = "red", size=1.5) +  # Point de rupture
  labs(
    title = "",
    x = "Date",
    y = "Concentrations totales ind.m^-3 "
  ) +
  theme_minimal()+
  theme(
    text = element_text(size = 16),  # Taille générale du texte
    axis.title = element_text(size = 16),  # Taille des titres des axes
    axis.text = element_text(size = 14),  # Taille des étiquettes des axes
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5), # Titre du graphique centré
    legend.text = element_text(size = 14),  # Taille du texte de la légende
    legend.title = element_text(size = 14)  # Taille du titre de la légende
  ) + scale_y_log10() + geom_vline(xintercept = as.Date(change_point_date), linetype = "dashed", color = "red", label = "Pettitt")

#Strucchange
library(strucchange)

breakpoints_result <- breakpoints(conc ~ 1,breaks=6, data = zoo_total)  # Constant model

summary(breakpoints_result)

# Dates
break_dates <- zoo_total$object_date[breakpoints_result$breakpoints]
print(break_dates)

ggplot(zoo_total, aes(x = object_date, y = conc)) +
  geom_path(color = "black") +  # Série temporelle
  geom_vline(xintercept = as.Date(change_point_date), linetype = "dashed", color = "red", label = "Pettitt") +  # Point de Pettitt
  geom_vline(xintercept = as.Date(break_dates), linetype = "dashed", color = "pink", label = "Strucchange") +  # Points de strucchange
  labs(
    title = "Comparison of breaking points (Pettitt vs Strucchange)",
    x = "Date",
    y = "Concentration"
  ) +
  theme_minimal()+scale_y_log10()


##### Step 2 : Trend analyses from this breaking point #####


#Biovolume calculation
#plankton_biovol<-plankton_tot%>%
  # convert the pixel variables in metric form
#  mutate(
#    area = object_area*(0.0106)^2,
#    area_excluded = object_area_exc*(0.0106)^2,
#    major = object_major*0.0106,
#    minor = object_minor*0.0106
#  ) |>
  # compute biovolume volume for each individual
#  mutate(
#    radius=sqrt(area/pi),
#    spherical_vol=(4/3)*pi*(radius)^3,
#    biovol=spherical_vol*concentration)

# Define a regular time sequence of 15 days
ref <- tibble(
  target_date=seq(from=as.Date("1967-01-05"), to=as.Date("2022-12-31"), by=14),
  year=year(target_date)
)
# identify years in which the number of obs is larger than usual
pbs <- ref |>
  count(year) |>
  filter(n>26)
# ->this is often an extra date in very late decembre => just remove it
ref <- filter(ref, !(year %in% pbs$year & yday(target_date) > 360))
ref %>%
  count(year) %>%
  filter(n>26)
# -> all OK

# change the name

plankton_group <- zoo_total
plankton_group$date <- zoo_total$object_date

# match data based on these reference dates
avail_dates_p <- unique(plankton_group$date)
avail_dates_p<-as.data.frame(avail_dates_p)
avail_dates_p<- avail_dates_p%>%mutate(year=year(avail_dates_p))
ref_p<-ref%>%
  mutate(
    closest_date = castr::closest(ref$target_date, avail_dates_p$avail_dates_p),
    date_diff = abs(closest_date - target_date) |> as.numeric()
  )

ref_p%>%group_by(year)%>%summarise(count=n())%>%filter(count!=26)

plankton_group<-zoo_total
dtest_p <- left_join(ref_p, plankton_group, by=c("closest_date"="date"), relationship="many-to-many")

# erase data for matches that are too far from the target
dtest_p <- dtest_p |>
  mutate(conc = if_else(date_diff > 14, NA, conc))

ggplot(dtest_p) +
  geom_point(aes(x=target_date, y=conc), size=0.2) + theme_bw()+labs(x="date", y="valeur")

dtest_p<-dtest_p%>%filter(target_date>"1980-04-09")
# -> OK
#Select point_B data whose dates are in closest date

#filter(target_date>"2000-02-02")


#Data visualization
y<- dtest_p
y%>%ggplot()+geom_line(aes(x=target_date,y=conc))+scale_y_log10()
#y<-y%>%filter(year>2000)
# -> OK

dec_p <- stlplus(y$conc,y$target_date, n.p=26, s.window="periodic", t.window=10)
plot(dec_p)
yy<- dec_p$data[,1:4]
yy$date<- y$target_date
yy<- yy%>%mutate(deseason=raw-seasonal)
yy <-yy %>%
  filter(!is.na(deseason))

#####Regression gls#####

statss_plankton <- yy %>%
  do({
    x <- .
    # 1. Mann-Kendall trend test
    mkt <- mk.test(yy$deseason)
    
    # 2. GLS regression
    # simple model
    m <- gls(deseason ~ date, data=yy)
    
    # Autocorrelation check
    a <- pacf(residuals(m, type="normalized"), plot=FALSE)
    if (abs(a$acf[1]) > 0.2) {
      # Add AR1 model on residuals
      m <- gls(deseason ~ date, data=yy, cor=corAR1(round(a$acf[1], 1)))
      a <- pacf(residuals(m, type="normalized"), plot=FALSE)
    }
    
    # Extract diagnostic information for both approaches
    mkt_results <- glance(mkt) |> select(p.value)
    
    # Extract summary information for GLS
    m_summary <- summary(m)
    
    # Create a tibble with the required statistics
    tibble(
      p.value.mk = mkt_results$p.value,
      r.squared.gls = 1 - (sum(m_summary$residuals^2) / sum((yy$deseason - mean(yy$deseason))^2)),
      p.value.gls = summary(m)$tTable["date", "p-value"],
      intercept.gls = m_summary$tTable[1, "Value"],
      slope.gls = m_summary$tTable[2, "Value"],
      cor.struct.gls = "AR1",  # or modify based on your cor structure
      acf.gls = a$acf[1]
    )
  }) |>
  ungroup() |>
  mutate(
    acf = abs(acf.gls),
    signif.mk = signif_stars(p.value.mk),
    signif.gls = signif_stars(p.value.gls)
  ) |>
  rename(
    p.value.mk = p.value.mk,
    r.squared.mk = r.squared.gls,
    p.value.gls = p.value.gls,
    intercept.gls = intercept.gls,
    slope.gls = slope.gls,
    cor.struct.gls = cor.struct.gls,
    acf.gls = acf.gls
  )


ggplot()  +
  geom_path(aes(date, deseason), data=yy)+
  geom_abline(aes(slope=slope.gls, intercept=intercept.gls), color="red",data=statss_plankton) + theme_minimal()+theme(
    text = element_text(size = 16),  # Taille générale du texte
    axis.title = element_text(size = 18),  # Taille des titres des axes
    axis.text = element_text(size = 14),  # Taille des étiquettes des axes
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Titre du graphique centré
    legend.text = element_text(size = 14),  # Taille du texte de la légende
    legend.title = element_text(size = 16)  # Taille du titre de la légende
  )+ scale_y_log10()
statss_plankton

