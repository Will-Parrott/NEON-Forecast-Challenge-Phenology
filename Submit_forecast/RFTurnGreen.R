#SCOPE
#This will be a 30-day forecast of phenology, specifically Daily Green Chromatic coordinate, in terrestrial NEON sites.
#Sites: 47 terrestrial sites: ABBY-YELL

#//////////////////////////////////////////////////////////////// DATA SETUP ////////////////////////////////////////////////////////////
#Model id

model_id <- "(Random)_Forests_Turn_Green"

#install docker non native packages
install.packages("tidymodels")
install.packages("ranger")

#Bring in necessary packages
library(tidyverse)
library(tidymodels)
library(lubridate)
library(ranger)


#Bring in targets

site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-targets/main/NEON_Field_Site_Metadata_20220412.csv") |> 
  dplyr::filter(phenology == 1)

focal_sites <- site_data$field_site_id

url <- "https://sdsc.osn.xsede.org/bio230014-bucket01/challenges/targets/project_id=neon4cast/duration=P1D/phenology-targets.csv.gz"
targets <- readr::read_csv(url, show_col_types = FALSE)

targets <- targets |> 
  filter(site_id %in% focal_sites,
         variable == "gcc_90") |> 
  group_by(site_id) |> 
  mutate(lag_gcc = lag(observation,1))


# Past weather
met_variables <- c("air_temperature", "surface_downwelling_shortwave_flux_in_air")

weather_past_s3 <- neon4cast::noaa_stage3()

weather_past <- weather_past_s3 |> 
  filter(site_id %in% focal_sites,
         datetime >= ymd('2017-01-01'),
         variable %in% met_variables) |> 
    dplyr::collect()

weather_past_daily <- weather_past |> 
  mutate(datetime = as_date(datetime)) |> 
  group_by(datetime, site_id, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  # convert air temperature to Celsius if it is included in the weather data
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction)

# Future weather
forecast_date = Sys.Date()
noaa_date = forecast_date - days(1)

#Forecast parameters
forecast_horizon <- 30
forecasted_dates <- seq(from = ymd(forecast_date), to = ymd(forecast_date) + forecast_horizon, by = "day")
summer <- seq(from = 6, to = 9)
winter <- c(12, 1, 2, 3)
n_members <- 200

weather_future_s3 <- neon4cast::noaa_stage2(start_date = as.character(noaa_date-1)) #-1 is a temp fix while stage 2 is being dumb

weather_future <- weather_future_s3 |> 
  dplyr::filter(datetime >= forecast_date,
                site_id %in% focal_sites,
                variable %in% met_variables) |> 
  collect()

weather_future_daily <- weather_future |> 
  mutate(datetime = as_date(datetime)) |> 
  # mean daily forecasts at each site per ensemble
  group_by(datetime, site_id, parameter, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  # convert air temperature to Celsius if it is included in the weather data
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction) |> 
  select(any_of(c('datetime', 'site_id', met_variables, 'parameter')))


#///////////////////////////////////////////////////////////UNCERTAINTY//////////////////////////////////////////////////////////////////
###NOTE TO FUTURE SELF: This fitting and evaluation should be moved to INSIDE a for loop - fit the model by site.
model_performance <- data.frame(
  site_id = character(),
  rsq = numeric()
)

forecast_df <- NULL
for(s in 1:length(focal_sites)){
  targets_m <- targets |> 
    pivot_wider(names_from = 'variable', values_from = 'observation') |> 
    left_join(weather_past_daily, 
              by = c("datetime","site_id")) |> 
    filter(datetime >= "2020-09-24",
           site_id == focal_sites[s]) |> 
    na.omit()
  
  weather_ensemble_names <- unique(weather_future_daily$parameter)
  #----------MODEL SETUP------------
  set.seed(228888) # 22 88 8 8, ring ring
  
  split <- initial_split(targets_m, prop = 0.80, strata = site_id) #this stratification doesn't work
  train_gcc<- training(split)
  test_gcc <- testing(split)
  
  #folds
  folds <- vfold_cv(train_gcc, v= 10)
  
  #Feature engineering w/recipe: cleaning up the data so that we can do regression with random forests
  cleanup <- train_gcc |> 
    recipe(gcc_90~.) |> 
    step_rm(datetime, site_id, project_id, duration) |> 
    step_naomit(air_temperature, surface_downwelling_shortwave_flux_in_air, lag_gcc)
  
  
  #Setting up object for hyperparameter tuning
  tune_spec <- 
    rand_forest(
      mtry=tune(),
      trees=tune(),
      min_n=tune()
    ) |> 
    set_engine("ranger", num.threads = parallel::detectCores()) |> 
    set_mode("regression")
  
  #Establishing workflow
  wF <- 
    workflow() |> 
    add_recipe(cleanup) |> 
    add_model(tune_spec)
  
  #parameter tuning
  tune_par <- 
    wF |> 
    tune_grid(
      resamples = folds, 
      grid = 27,
      control = control_grid(save_pred = T),
      metrics = metric_set(rmse)
    )
  tune_par |> collect_metrics() |> arrange(mean)
  best_param <- tune_par |> select_best(metric = "rmse") #save best parameters
  
  #Finalized workflow
  ud_wF <- 
    wF |> finalize_workflow(best_param)
  
  #Fit to training data
  gcc_fit <- ud_wF |> fit(data = train_gcc)
  
  gcc_predict <- predict(gcc_fit, new_data = test_gcc)
  gcc_pred_test <- bind_cols(test_gcc, gcc_predict)
  #---------MODEL EVALUATE-----------
  multi_metric <- metric_set(rmse, rsq)
  metric_table <- gcc_pred_test |> 
    multi_metric(truth = gcc_90, estimate = .pred) 
  print(paste(focal_sites[s], "-", s))
  print(metric_table)
  
  metric_table <- metric_table |> 
    pivot_wider(names_from = .metric, values_from = .estimate)
  model_performance <- model_performance |> add_row(site_id = focal_sites[s], rsq = metric_table$rsq)

  ggplot(gcc_pred_test, aes(x=gcc_90, y=.pred)) + geom_point()
  metric_table
  
  #plot data for expected trends
  dates_2025 <- seq(as.Date("2025-01-01"), as.Date("2025-12-31"), by="day")
  targets_2025 <- targets_m |> filter(as.Date(datetime) %in% dates_2025)
  ggplot(targets_2025, aes(x = datetime, y = gcc_90)) + geom_line()
  
  #set dataframe for stable period to determine observation uncertainty
  stable_per <- seq(as.Date("2025-05-15"), as.Date("2025-06-15"), by = "day")
  targets_unc_ref <- targets_2025 |> filter(as.Date(datetime) %in% stable_per)
  init_sd <- sd(targets_unc_ref$gcc_90)
  init_uc_df <- rnorm(n_members, mean = 0, sd = init_sd)
  
  for(d in 1:length(forecasted_dates)){
    
    met_ens_id <- 0
    for(ens in 1:n_members){
      print(paste0(forecasted_dates[d], "-", ens))
      
      if(met_ens_id <= 30){
        met_ens_id <- met_ens_id + 1
        ens_nm <- paste0(ens, "-", met_ens_id)
      }else{
        met_ens_id <- 1
      }
      
      met_ens <- weather_ensemble_names[met_ens_id]
      
      if(d == 1){ #set current lagged value - last available data for d=1, otherwise use last available prediction
        lag_curr <- targets_m$lag_gcc[nrow(targets_m)] + init_uc_df[met_ens+1]
      }else{
        if(met_ens <30){
          curr_ens_pred <- forecast_df |> filter(parameter == met_ens + 1)
        }else{
          curr_ens_pred <- forecast_df |> filter(parameter == met_ens)
        }
        lag_curr <- curr_ens_pred$prediction[nrow(curr_ens_pred)]
      }
      
      #Establish prediction dataframe with necessary components
      pred_df <- weather_future_daily |> 
        filter(datetime == as.Date(forecasted_dates[d]), site_id == focal_sites[s], parameter == met_ens) |>
        mutate(lag_gcc = lag_curr,
               project_id = unique(targets_m$project_id),
               duration = unique(targets_m$duration))
      
      forecasted_gcc <- predict(gcc_fit, new_data = pred_df)
      
      curr_site_df <- tibble(datetime = pred_df$datetime,
                             reference_datetime = forecast_date,
                             site_id = focal_sites[s],
                             parameter = ens, #make sure this uses ens, not met_ens. Plot breaks if done latter way.
                             prediction = forecasted_gcc[[1]],
                             variable = "gcc_90")
      forecast_df<- as.data.frame(rbind(forecast_df, curr_site_df)) 
      
      #the targets data only goes to 2024, remove the na.rm() thing from the targets.m call.
    }
  }
  print(paste(focal_sites[s], " forecast run"))
}
write.csv(forecast_df, file = "forecast_df.csv")

ggplot(forecast_df, aes(datetime, prediction, group=parameter, color=site_id)) + geom_line()
#---- Covert to EFI standard ----

# Make forecast fit the EFI standards
forecast_df_EFI <- forecast_df %>%
  filter(datetime > forecast_date) %>%
  mutate(model_id = model_id,
         family = 'ensemble',
         duration = 'P1D',
         parameter = as.character(parameter),
         project_id = 'neon4cast') %>%
  select(datetime, reference_datetime, duration, site_id, family, parameter, variable, prediction, model_id, project_id)
#---------------------------#



# ----- Submit forecast -----
# Write the forecast to file
theme <- 'terrestrial'
date <- forecast_df_EFI$reference_datetime[1]
forecast_name <- paste0(forecast_df_EFI$model_id[1], ".csv")
forecast_file <- paste(theme, date, forecast_name, sep = '-')

write_csv(forecast_df_EFI, forecast_file)

neon4cast::forecast_output_validator(forecast_file)


neon4cast::submit(forecast_file =  forecast_file, ask = FALSE) # if ask = T (default), it will produce a pop-up box asking if you want to submit

#--------------------------#

forecast_df_EFI |> 
  ggplot(aes(x=datetime, y=prediction, group = parameter)) +
  geom_line() +
  facet_wrap(~site_id) +
  labs(title = paste0('Forecast generated for ', forecast_df_EFI$variable[1], ' on ', forecast_df_EFI$reference_datetime[1]))+ 
  theme_bw(base_size = 10) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

plot_file_name <- paste0("Submit_forecast/", forecast_df_EFI$variable[1], '-', forecast_df_EFI$reference_datetime[1], ".png")
ggsave(plot_file_name)


#References:
#https://rpubs.com/GChirinos/Tuning_Random_Forest
#https://www.tmwr.org/grid-search.html
