source("functions.R")
source("../check_nowcast_submission/plot_functions.R")

library(zoo)
library(lubridate)
library(dplyr)
library(tidyr)

# read truth data:
observed0 <- read.csv("../../data-truth/COVID-19/COVID-19_hospitalizations_preprocessed.csv",
                      colClasses = c("date" = "Date"))


# dates for which to produce nowcasts:
forecast_dates <- seq(from = ymd("2021-11-22"),
                      to = ymd("2022-04-29"),
                      by = "day"
                      )


for(i in seq_along(forecast_dates)){
  all_nc <- NULL
  
  forecast_date <- forecast_dates[i]
  cat(as.character(forecast_dates[i]), "\n")
  
  # limited by number of observations (in the early part, not relevant anymore)
  n_history_dispersion <- min(60, as.numeric(forecast_date - as.Date("2021-04-06")) - 40)
  
  # generate nowcasts for age groups:
  for(ag in c("00+", "00-04", "05-14", "15-34", "35-59", "60-79", "80+")){
    
    observed_temp <- subset(observed0, location == "DE" & age_group == ag)
    # prepare for plotting:
    observed_for_plot <- truth_as_of(observed_temp, age_group = ag,
                                     location = "DE",
                                     date = Sys.Date())
    
    # truth data as of forecast_date for plot:
    observed_for_plot_old <- truth_as_of(observed_temp, age_group = ag,
                                         location = "DE",
                                         date = forecast_date)
    
    # generate truth data as of forecast_date
    observed_temp <- back_in_time_df(observed_temp, date = forecast_date)
    
    # compute nowcast:
    
    # undebug(compute_nowcast)
    nc <- compute_nowcast_revised(observed = observed_temp, 
                          location = "DE", 
                          age_group = ag, 
                          n_history_expectations = 60, 
                          n_history_dispersion = n_history_dispersion,
                          min_horizon = 0,
                          max_horizon = 28)
    
    # generate a plot:
    # undebug(plot_forecast)
    # plot_forecast(forecasts = nc,
    #               location = "DE", age_group = ag,
    #               truth = observed_for_plot, target_type = paste("inc hosp"),
    #               levels_coverage = c(0.5, 0.95),
    #               start = as.Date(forecast_date) - 35,
    #               end = as.Date(forecast_date) + 28,
    #               forecast_date = forecast_date)
    # axis(1)
    # title(paste(forecast_date, "-", ag))
    # lines(observed_for_plot_old$date, observed_for_plot_old$value, col = "darkgrey", lty  ="dashed")
    
    if(is.null(all_nc)){
      all_nc <- nc
    }else{
      all_nc <- rbind(all_nc, nc)
    }
  }
  
  # generate nowcasts for federal states:
  for(loc in c("DE-BB", "DE-BE", "DE-BW", "DE-BY",
               "DE-HB", "DE-HE", "DE-HH", "DE-MV", 
               "DE-NI", "DE-NW", "DE-RP", "DE-SH", 
               "DE-SL", "DE-SN", "DE-ST", "DE-TH")){
    
    observed_temp <- subset(observed0, location == loc & age_group == "00+")
    # prepare for plotting:
    observed_for_plot <- truth_as_of(observed_temp, age_group = "00+",
                                     location = loc,
                                     date = Sys.Date())
    
    # truth data as of forecast_date for plot:
    observed_for_plot_old <- truth_as_of(observed_temp, age_group = "00+",
                                         location = loc,
                                         date = forecast_date)
    
    # generate truth data as of forecast_date
    observed_temp <- back_in_time_df(observed_temp, date = forecast_date)
    
    # compute nowcast:
    # undebug(compute_nowcast)
    nc <- compute_nowcast_revised(observed = observed_temp, 
                          location = loc, 
                          age_group = "00+",
                          n_history_expectations = 60,
                          n_history_dispersion = n_history_dispersion,
                          min_horizon = 0,
                          max_horizon = 28)
    
    # generate a plot:
    # plot_forecast(forecasts = nc,
    #               location = loc, age_group = "00+",
    #               truth = observed_for_plot, target_type = paste("inc hosp"),
    #               levels_coverage = c(0.5, 0.95),
    #               start = as.Date(forecast_date) - 35,
    #               end = as.Date(forecast_date) + 28,
    #               forecast_date = forecast_date)
    # axis(1)
    # title(paste(forecast_date, "-", loc))
    # lines(observed_for_plot_old$date, observed_for_plot_old$value, col = "darkgrey", lty  ="dashed")
    
    if(is.null(all_nc)){
      all_nc <- nc
    }else{
      all_nc <- rbind(all_nc, nc)
    }
  }
  
  # write out:
  write.csv(all_nc, file = paste0("../../data-processed_retrospective/KIT-simple_nowcast_revised/", forecast_date, "-KIT-simple_nowcast.csv"), row.names = FALSE)
}

# Quick plot of a single nowcast 
df_wide <- nc |> 
  filter(type == "quantile") |>
  pivot_wider(id_cols  = c(forecast_date, target_end_date,
                           age_group, location),
              names_from = quantile,
              names_prefix = "q_",
              values_from = value) |>
  select(-forecast_date, -age_group, -location, 
         -q_0.1, - q_0.9) |>
  mutate(model = "KIT simple nowcast revised")

ggplot(df_wide) + 
  geom_line(aes( x = target_end_date, y = `q_0.5`)) +
  geom_ribbon(aes(x = target_end_date, ymin = `q_0.025`, ymax = `q_0.975`),
              alpha = 0.1)+
  geom_ribbon(aes(x = target_end_date, ymin = `q_0.25`, ymax = `q_0.75`),
              alpha = 0.1)


