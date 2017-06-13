library(reservoir)  # load reservoir package for bootcast function

load("example_inflow")  # load example inflow time series

# For this example we'll produce the forecasts for years 2005 through 2010...
startYr <- 2005
endYr <- 2010
Q_ <- window(Q, start = c(startYr, 1), end = c(endYr, 12))

# SET PARAMETERS
H = 12  # 12 month forecast horizon
error = 0.5  # "injected error"

# GENERATE FORECASTS FOR OBSERVED INFLOW USING KNN BOOTSTRAP METHOD
fcast <- bootcast(Q, H = H,
         start_yr = startYr,
         end_yr = endYr,
         k = 7, d = 2, sampling_mode = "all")

# MARTINGALE MODEL OF FORECAST EVOLUTION
source("MMFE.R")

## CALL THE MARTINGALE MODEL
fcast_ <- Martingale(Q_, H, error = error, fcast = fcast)


## PLOT RESULTS...

layout(1:2)

## PLOT THE INITIAL FORECAST
plot(Q_, col = "darkgrey", ylim = c(0, max(Q_)), ylab = "Q")
for (i in 1:length(Q_)) {
  lines(ts(fcast[i, ], start = time(Q_)[i], 
           frequency = frequency(Q_)), col = "red")
}

## PLOT THE MARTINGALE FORECAST
plot(Q_, col = "darkgrey", ylim = c(0, max(Q_)), ylab = "Q")
for (i in 1:length(Q_)) {
  lines(ts(fcast_[i, ], start = time(Q_)[i], 
           frequency = frequency(Q_)), col = "blue")
}
legend("top", paste0("Injected error = ", error))


