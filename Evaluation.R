##-------------------------- NECESSARY PACKAGES------------------------------------------------------
library(tidyverse)
library(tidymodels)
library(lubridate)
library(ranger)

##--------------------------BRING IN TARGETS AND REFORECAST DATA --------------------------------------
ttf <- seq(ymd("2024-01-01"), ymd("2024-12-31"), by = "day")
site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
  dplyr::filter(phenology == 1)

focal_sites <- site_data$field_site_id

url <- "https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D/phenology-targets.csv.gz"
targets <- readr::read_csv(url, show_col_types = FALSE)

targets <- targets |> 
  filter(site_id %in% focal_sites,
         variable == "gcc_90", 
         as.Date(datetime) %in% ttf) |> 
  group_by(site_id) |> 
  mutate(lag_gcc = lag(observation,1))

# Reget reforecast data
forecast_df_EFI <- read.csv("C:/Users/willi/Desktop/NEON-Forecast-Challenge-Phenology/reforecast_df_2024.csv") 
forecast_df_EFI |> 
  ggplot(aes(x=datetime, y=prediction, group = parameter)) +
  geom_line() +
  facet_wrap(~site_id) +
  labs(title = paste0('Forecast generated for ', forecast_df_EFI$variable[1], ' on ', forecast_df_EFI$reference_datetime[1]))+ 
  theme_bw(base_size = 10) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
##--------------------------------------------- EVALUATION -------------------------------------------

library(scoringRules)

forecasted_data <- forecast_df_EFI |>  
  mutate(horizon = as.numeric(mdy(datetime) - mdy(reference_datetime))) |> 
  group_by(datetime, horizon, site_id, reference_datetime) |> 
  summarize(mean = mean(prediction)) |> 
  mutate(datetime = mdy(datetime),
         reference_datetime = mdy(reference_datetime))

ref_dates = seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by="month")

model_crps <- data.frame(site_id = character(), horizon = integer(), crps = numeric(), datetime = as.Date(numeric()), reference_datetime = as.Date(numeric()))  
for(i in 1:length(focal_sites)){
  curr_site <- focal_sites[i]
  for (t in 1:length(ref_dates)){
    rel_datetimes <- seq(ref_dates[t], ref_dates[t] + months(1))
    for(d in 1:length(rel_datetimes)){
      targets_ref <- targets |> 
        filter(as.Date(datetime) == rel_datetimes[d],
               site_id == curr_site)
      model_ref <- forecasted_data |> 
        filter(datetime == rel_datetimes[d],
               site_id == curr_site)
      horizon <- model_ref$horizon
      if(is.na(targets_ref$observation[1]) || is.na(model_ref$mean[1])){
        model_crps <- model_crps |> add_row(horizon = horizon, crps = NA, datetime = rel_datetimes[d], site_id = curr_site, reference_datetime = ref_dates[t])
      }else{
        crps <- crps_sample(targets_ref$observation, model_ref$mean)
        model_crps<- model_crps |> add_row(horizon = horizon, crps = crps, datetime = rel_datetimes[d], site_id = curr_site, reference_datetime = ref_dates[t])
      }
      print(paste(curr_site, "-", ref_dates[t], "-", rel_datetimes[d]))
    }
  }
}

veg_t <- unique(site_data$phenocam_vegetation)
sites_by_vt <- cbind(site_id = site_data$field_site_id, veg_type = site_data$phenocam_vegetation)
vt_site_crps <- merge(model_crps, sites_by_vt, by = 'site_id') |> cbind(model_id = rep("(Random)_Forests_Turn_Green",nrow(model_crps)))


season_lineup <- data.frame(datetime = ymd(unique(vt_site_crps$datetime)), season = NA)
for(d in 1:length(season_lineup$datetime)){
  if(month(season_lineup$datetime[d]) %in% c(12, 1,2)){
    season_lineup$season[d] = "Winter"
  }else if(month(season_lineup$datetime[d]) %in% c(3,4,5)){
    season_lineup$season[d] = "Spring"
  }else if(month(season_lineup$datetime[d]) %in% c(6,7,8)){
    season_lineup$season[d] = "Summer"
  }else if(month(season_lineup$datetime[d]) %in% c(9, 10, 11)){
    season_lineup$season[d] = "Fall"
  }
}

#ALL SITES average

model_crps_avg <- model_crps |> 
  group_by(horizon) |> 
  summarize(crps = as.numeric(mean(crps, na.rm = TRUE))) |> 
  mutate(model = "(Random)_Forests_Turn_Green")

model_crps_szn <- model_crps |> 
  group_by(datetime) |> 
  summarize(crps = as.numeric(mean(crps, na.rm = TRUE))) |> 
  mutate(model_id = "(Random)_Forests_Turn_Green") |> 
  merge(season_lineup, by = 'datetime')

baseline_models <- arrow::open_dataset("s3://anonymous@bio230014-bucket01/challenges/scores/bundled-parquet/project_id=neon4cast/duration=P1D/variable=gcc_90?endpoint_override=sdsc.osn.xsede.org") |> 
  filter(site_id %in% focal_sites,
         reference_datetime %in% ref_dates,
         model_id %in% c("climatology", "persistenceRW")) |> 
  collect()

crps_bl <- baseline_models |> 
  mutate(horizon = as.numeric(datetime - reference_datetime)) |> 
  summarize(mean_crps = mean(crps, na.rm = TRUE), .by = c("model_id", "horizon")) |> 
  pivot_wider(names_from = model_id, values_from = mean_crps) |> 
  subset(!(horizon %in% c(32, 33, 34, 35))) |> 
  pivot_longer(cols = 2:3, names_to = "model", values_to = "crps") |> 
  rbind(model_crps_avg)

bl_szn_site_crps <- baseline_models |> 
  mutate(datetime = ymd(datetime)) |> 
  merge(season_lineup, by = 'datetime') |> 
  group_by(datetime, season, model_id) |> 
  summarise(crps = mean(crps, na.rm = TRUE)) |> 
  rbind(model_crps_szn)
  
bl_vt_site_crps <- merge(baseline_models, sites_by_vt, by = 'site_id')


crps_bl_vt <- bl_vt_site_crps |> 
  mutate(horizon = as.numeric(ymd(datetime) - ymd(reference_datetime))) |> 
  group_by(datetime, site_id, model_id, horizon, veg_type) |> 
  summarise(crps = mean(crps, na.rm = T)) |> 
  subset(!(horizon %in% c(32, 33, 34, 35))) |> 
  rbind(vt_site_crps) |> 
  group_by(model_id,veg_type) |> 
  summarise(crps = mean(crps, na.rm = T))

  

ggplot(crps_bl, aes(y=crps, x=horizon, color = model)) + geom_line()

ggplot(crps_bl_vt, aes(x=veg_type, y=crps, fill = model_id)) + geom_col(position = position_dodge())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(bl_szn_site_crps, aes(x=season, y=crps, fill = model_id)) + geom_col(position = position_dodge())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
