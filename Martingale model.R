library(reservoir)  # load reservoir package for bootcast function

load("example_inflow")  # load example inflow time series
plot(Q)  # view observed inflow time series

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
Martingale <- function(Q, H, error, fcast){
  
  cyc <- cycle(Q)
  len <- length(Q)
  frq <- frequency(Q)
  cycF <- matrix(rep(cyc, length.out = H * len + (H - 1))[-((1 : (H-1) * len + 1) + seq(0,1, length.out = H - 1))],ncol = H)
  fcast_error <- rep(error, H)
  fcast_matrix <- matrix(rep(Q, len),nrow = len, byrow = TRUE)
  fcast_update <- matrix(nrow = (len - H + 1), ncol = H)

  for (i in 1:(len)){
    if(len - i < H - 1){
      fcast_matrix[i, i:len] <- fcast[i,1 : (len - i + 1)]
    }else{
      fcast_matrix[i, i:(i + H - 1)] <- fcast[i,]
    }
  }

  for (i in 1 : (len - H + 1)){
    fcast_update[i,] <- fcast_matrix[i, i : (i + H - 1)] - fcast_matrix[i + 1, i : (i + H - 1)]
  }
  
  VCV_matrix <- cov(fcast_update)
  V <- chol(VCV_matrix)
  X <- matrix(rnorm((len + H - 1) * H), len + H - 1, H)
  
  for (i in 1:H){
    X[,i] <- rnorm(len + H - 1 , 0, fcast_error[i])
  }

  Y <- (t(apply(X, 1, "%*%",(V))))
  fcast_rep <- matrix(nrow = len, ncol = H)
  Yx <- rbind(fcast_update,fcast_update)

  for(i in 1 : len){
    fcast_rep[i,1] <- Q[i] + Y[i,1]
  }
  
  for(j in 2 : H){
    for(i in 1: (len - (j - 1))){
      fcast_rep[i, j] <- fcast_rep[i + 1,j - 1] + Y[i, j]
    }
  }
  
  fcast_rep[which(fcast_rep < 0)] <- min(Q)  # Prevents negative values in forecast
  fcast_rep[is.na(fcast_rep)] <- 0
  return(fcast_rep)
}

## CALL THE MARTINGALE MODEL
fcast_ <- Martingale(Q_, H, error = error, fcast = fcast)


## PLOT RESULTS...

layout(1:2)

## PLOT THE INITIAL FORECAST
plot(Q_, col = "darkgrey", ylim = c(0, max(Q_)), ylab = "Q")
for (i in 1:length(Q_)) {
  lines(ts(fcast[i, ], start = time(Q_)[i], 
           frequency = frq), col = "red")
}

## PLOT THE MARTINGALE FORECAST
plot(Q_, col = "darkgrey", ylim = c(0, max(Q_)), ylab = "Q")
for (i in 1:length(Q_)) {
  lines(ts(fcast_[i, ], start = time(Q_)[i], 
           frequency = frq), col = "blue")
}
legend("top", paste0("Injected error = ", error))


