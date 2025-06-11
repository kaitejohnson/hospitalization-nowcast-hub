source("functions.R")
source("../check_nowcast_submission/plot_functions.R")

library(zoo)

# read truth data:
observed0 <- read.csv("../../data-truth/COVID-19/COVID-19_hospitalizations_preprocessed.csv",
                      colClasses = c("date" = "Date"))




forecast_date <-"2022-04-01"
ag <- "00+"

observed_temp <- subset(observed0, location == "DE" & age_group == ag)
observed_temp <- back_in_time_df(observed_temp, date = forecast_date)


# Specify the function args
observed <- observed_temp
location  <-"DE" 
age_group <- ag 
n_history_expectations <- 60
n_history_dispersion <- 60
min_horizon <- 0
max_horizon <- 28
max_delay <- 40

# reporting triangle as matrix
matr_observed <- as.matrix(observed[, grepl("value", colnames(observed))])
#reduce to max delay:
matr_observed <- cbind(matr_observed[, 1:(max_delay + 1)],
                       matrix(rowSums(matr_observed[, -(1:(max_delay))],
                                      na.rm = TRUE), ncol = 1))

colnames(matr_observed)[max_delay + 2] <- paste0("value_", max_delay + 1, "d")
observed[nrow(matr_observed) - 0:max_delay, max_delay + 2] <- NA
rownames(matr_observed) <- observed$date

# # Make matr_observed really small so we are testing on a very restricted use 
# # case we can break down more easily 
# max_t <- nrow(matr_observed)
# matr_observed <- matr_observed[((max_t-10): max_t), 1:5]
# matr_observed[7:9, 3]<- 0
# max_delay <- 4

# max_horizon <-3
# n_history_dispersion <- 3
# n_history_expectations <- 5




nr <- nrow(matr_observed)
nc <- ncol(matr_observed)
# compute point forecasts
expectation_to_add <- # full expectations
  expectation_to_add_already_observed <- # expectations of the sum over already observable quantities
  to_add_already_observed <- # sums over the respective observed quantities
  matrix(NA, nrow = nr, ncol = max_horizon + 1)


# generate point forecasts for current date and n_history_dispersion preceding weeks
# these are necessary to estimate dispersion parameters
for(t in (nr - n_history_dispersion):nr){
  matr_observed_temp <- back_in_time(matr_observed, t)
  point_forecasts_temp <- compute_expectations(matr_observed_temp, n_history = n_history_expectations)
  
  for(d in min_horizon:max_horizon){
    inds_nowc <- indices_nowcast(matr_observed_temp, d = d, w=7, 
                                 n_history_expectations = n_history_expectations)
    inds_already_observed <- tail(!is.na(matr_observed[1:t, ]), n_history_expectations)
    
    expectation_to_add[t, d + 1] <- sum(point_forecasts_temp*inds_nowc, na.rm = TRUE)
    expectation_to_add_already_observed[t, d + 1] <- sum(point_forecasts_temp*inds_already_observed*inds_nowc, na.rm = TRUE)
    to_add_already_observed[t, d + 1] <- sum(tail(matr_observed[1:t, ], n_history_expectations)*
                                               inds_already_observed*inds_nowc, na.rm = TRUE)
  }
}

# remove last row to estimate dispersion
expectation_to_add_already_observed <- expectation_to_add_already_observed[-nrow(expectation_to_add_already_observed), ]
to_add_already_observed <- to_add_already_observed[-nrow(to_add_already_observed), ]

# estimate dispersion
size_params <- numeric(max_horizon +1)
for(i in min_horizon:max_horizon){
  size_params[i + 1] <- fit_nb(x = to_add_already_observed[, i + 1], 
                               mu = expectation_to_add_already_observed[, i + 1] + 0.1)
  # plus 0.1 to avoid ill-defined NB
}


# Estimate dispersion using baselinenowcast
triangle <- matr_observed[, 1:41]
colnames(triangle) <- NULL

# fix 0s 
triangle[321:323, 38]<- 0

pt_nowcast_mat <- generate_pt_nowcast_mat(triangle, n = n_history_expectations,
                                          max_delay = max_delay)
trunc_rep_tri_list <- truncate_triangles(triangle,
                                         n = n_history_dispersion)
reporting_triangle_list <- generate_triangles(trunc_rep_tri_list)
pt_nowcast_mat_list <- generate_pt_nowcast_mat_list(reporting_triangle_list,
                                                    max_delay = max_delay,
                                                    n = n_history_expectations)

disp_bnc <- estimate_dispersion(
  pt_nowcast_mat_list,
  trunc_rep_tri_list,
  reporting_triangle_list,
  fun_to_aggregate = sum,
  k = 7
)

df <- data.frame(
  horizon = 0:max_horizon,
  dispersion_kit = size_params,
  disp_bnc = disp_bnc[1:(max_horizon+1)]
)
p <- ggplot(df) +
  geom_line(aes(x = horizon, y = dispersion_kit), color = "orange4", size = 2) +
  geom_line(aes(x = horizon, y = disp_bnc), color = "darkgreen")
p

ggsave(p, filename = "compare_disp_params_full_data_7d_sum_w_last_col.png")
