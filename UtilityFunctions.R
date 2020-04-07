## Utility Functions
# Last edited: 022120

## (1) buchanan_log10N, (2) gompertz_log10N, (3) baranyi_log10N
# Source: All equations for the 3 models below are copied from the Nlms package in R (https://rdrr.io/cran/nlsMicrobio/src/R/growthmodels.R)

# Purpose: Calculate log10N using respective growth model (either buchanan, gompertz, or barayani)

# Parameters: Same for (1) buchanan, (2) gompertz, and (3) barayani
# (i) t: time in hours
# (ii) lag: length of lag phase
# (iii) mumax: growth rate (LOG10 CFU/mL per hour)
# (iv) LOG10N0: initial microbial concentration 
# (v) LOG10Nmax: carrying capacity

# Functions:
buchanan_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax){
  ans <- LOG10N0 + (t >= lag) * (t <= (lag + (LOG10Nmax - LOG10N0) * log(10)/mumax)) * mumax * (t - lag)/log(10) + (t >= lag) * (t > (lag + (LOG10Nmax - LOG10N0) * log(10)/mumax)) * (LOG10Nmax -     LOG10N0)
  return(ans)
}

gompertz_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10N0 + (LOG10Nmax - LOG10N0) * exp(-exp(mumax * exp(1) * (lag - t)/((LOG10Nmax - LOG10N0) * log(10)) + 1))
  return(ans)
}

baranyi_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10Nmax + log10((-1 + exp(mumax * lag) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax - LOG10N0)))
  return(ans)
}


## (4) log10N_func

#Purpose: Implement the appropriate growth model based on provided model_name, in order to calculate log10N 

#Parameters:
# (i) t, (ii) lag, (iii) mumax, (iv) LOG10N0, & (v) LOG10Nmax: Same as for above functions for the 3 growth models
# (vi) model_name: Model to use for calcuating log10N

# Function:
log10N_func <- function(t, lag, mumax, LOG10N0, LOG10Nmax, model_name="buchanan") {
  if (model_name == "buchanan") {
    return(buchanan_log10N(t, lag, mumax, LOG10N0, LOG10Nmax) )
  }
  else if(model_name == 'baranyi') {
    return(baranyi_log10N(t, lag, mumax, LOG10N0, LOG10Nmax) )
  }
  else if(model_name == 'gompertz') {
    return(gompertz_log10N(t, lag, mumax, LOG10N0, LOG10Nmax) )
  }
  else {
    stop(paste0(model_name, " is not a valid model name. Must be one of buchanan, baranyi, gompertz"))
  }
}


## (5) muAtNewTemp

# Source: this function uses Ratkowsky's square root model which describes the effect of temperature on the growth of microorganisms (https://www.ncbi.nlm.nih.gov/pubmed/22417595)

# Purpose: Calculate the new mu parameter at new temperature.

# Parameters: 
# (i) newTemp: New temperature for which we calculate mu
# (ii) oldMu: Previous mu value to adjust
# (iii) oldTemp: Temperature corresponding to previous mu; NOTE: If you don't specify "oldTemp", 
#then automatically oldTemp = 6C; this value (6C) was the temp at which growth curve experiments were originally performed at)
# (iv) T0: Parameter used to calculate new mu; NOTE: If you don't specify "T0", then automatically T0 = -3.62C; this value (-3.62C) 
#was determined using Ratkowsky's square root model and Paenibacillus ordorifer growth curves obtained at 4, 7, and 32C in BHI broth (N.H. Martin unpublished data)

# Function:
muAtNewTemp <- function(newTemp, oldMu, oldTemp = 6, T0 = -3.62) { 
  numerator <- newTemp - T0
  denom <- oldTemp - T0
  newMu <- ((numerator / denom)^2) * oldMu
  
  return(newMu)
}


## (6) adjustLag

# Purpose: Adjust the lag phase based on the Zwietering 1994 paper.

# Parameters:  
# (i) t: Current timestep (in days)
# (ii) oldLag: Lag time (in days) at the previous temperature 
# (iii) newLag: Lag time (in days) at the current temparature 
# (iv) restartExp: If true, then lag phase restarts even if already in exponential growth phase; NOTE: If you don't 
#specify "restartExp", then automatically restartExp = T;
# (v) adjustmentConstant: Amount to adjust lag; NOTE: If you don't provide "adjustmentConstant", 
#then automatically adjustmentConstant = 0.25; Zwietering 1994 paper suggests using adjustmentConstant = 0.25

# Function:
adjustLag <- function (t, oldLag, newLag, restartExp = T, adjustmentConstant = 0.25) {
  #determine the amount of lag phase completed
  remainingLag <- 1 - (t / oldLag)
  if(restartExp) {
    remainingLag <- ifelse(remainingLag < 0, 0, remainingLag)
  }
  else {
    adjustedLag <- ifelse(remainingLag <=0, oldLag,
                          t + remainingLag * newLag + adjustmentConstant*newLag)
  }
  
  adjustedLag <- t + remainingLag*newLag + adjustmentConstant*newLag
  return(adjustedLag)
}


## (7) lagAtNewTemp

# Purpose: Calculate the new lag parameter at new temperature.

# Parameters:
# (i) t: Current timestep (in days)
# (ii) newTemp: New temperature for which we calculate lag
# (iii) oldLag: Previous lag value to adjust
# (iv) oldTemp: Temperature corresponding to previous lag; NOTE: If you don't specify "oldTemp", 
#then automatically oldTemp = 6C; this value (6C) was the temp at which growth curve experiments were originally performed at)
# (v) T0: Parameter used to calculate new mu; NOTE: If you don't specify "T0", then automatically T0 = -3.62C; this value (-3.62C) 
#was determined using Ratkowsky's square root model and Paenibacillus ordorifer growth curves obtained at 4, 7, and 32C in BHI broth (N.H. Martin unpublished data)

# Function:
lagAtNewTemp <- function (t, newTemp, oldLag, oldTemp = 6, T0 = -3.62) {
  numerator <- oldTemp -T0
  denom <- newTemp - T0
  newLag <- ( (numerator / denom)^2) * oldLag
  return(newLag)
}

## (8) getPrevRow

getPrevRow <- function(df, sim_run, milk_unit, day) {
  old_temp <- df[df$lot_id == sim_run & df$milk_unit == milk_unit & df$day==day-1,] 
}

# Purpose: 

# Parameters:
# (i) df: dataframe
# (ii) sim_run:
# (iii) milk_unit:
# (iv) day:
