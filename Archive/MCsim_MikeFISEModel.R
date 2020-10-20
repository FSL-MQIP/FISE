##################################################
#  Monte Carlo Psychrotolerant Sporeformer Simulation v2.0
##################################################
## Project: Dairy Spoilage Model
## Script purpose: Simulate spore growth for half gallon samples of milk over a given
##                 number of days.
## Date:  June 01, 2019
## Original Author: Ariel Buehler, Cornell University, ajb466@cornell.edu
## V2.0 Modifications by: Mike Phillips, Cornell University, mdp38@cornell.edu
##################################################
## Notes: Version  2.0 adds a temperature distribution to the original model, 
## which assumed a constant 6.0 degree temperature for all samples.
##################################################


#use a seed value for reproducibility
seed_value = 42;
set.seed(seed_value)


# Utility Functions   -----
##muAtNewTemp
#Purpose: Calculate the new mu parameter at new temperature.
#Params:  newTemp: the new temperature for which we calculate mu
#         oldMu: the previous mu value to adjust
#         oldTemp: the temperature corresponding to previous mu
#         T0:    Parameter used to calculate new mu
muAtNewTemp <- function(newTemp, oldMu, oldTemp = 6, T0 = -3.62) {
  numerator <- newTemp - T0
  denom <- oldTemp - T0
  newMu <- ((numerator / denom)^2) * oldMu
  
  return(newMu)
}

##lagAtNewTemp
#Purpose: Calculate the new lag parameter at new temperature.
#Params:  newTemp: the new temperature for which we calculate lag
#         oldLag: the previous lag value to adjust
#         oldTemp: the temperature corresponding to previous lag
#         T0:    Parameter used to calculate new lag
lagAtNewTemp <- function (newTemp, oldLag, oldTemp = 6, T0 = -3.62) {
  numerator <- oldTemp -T0
  denom <- newTemp - T0
  newLag <- ( (numerator / denom)^2) * oldLag
  return(newLag)
}

#Function to calculate log10N
#Verified and tested correct, same as nlsMicrobio function
#Purpose: This implements the growth model
log10N_func = function(t,lag,mumax,LOG10N0,LOG10Nmax){
  ans <- LOG10N0 + (t >= lag) * (t <= (lag + (LOG10Nmax - LOG10N0) *     log(10)/mumax)) * mumax * (t - lag)/log(10) + (t >= lag) * (t > (lag + (LOG10Nmax - LOG10N0) * log(10)/mumax)) * (LOG10Nmax -     LOG10N0)
  return(ans)
}



# Data frame creation and setup   ----

#Set up data frame to store count at each day
#Size is for n_sim bulk tanks, n_half_gal half gallon lots, n_day days (14-24)
n_sim <-1000    #1000 is for testing and exploring, experiments require at least 10k
n_halfgal <-10
n_day <- 11

#Repeat each element of the sequence 1..100, 110 times 11,000
BT <- rep(seq(1, n_sim), each = n_halfgal * n_day)
#Repeat 1..10 1100 times  11,000
half_gal <- rep(seq(1, n_halfgal), times = n_day * n_sim)
#Vector of FALSE, 11,000
AT <- vector(mode="logical", n_sim * n_halfgal * n_day)
#Repeat   (Repeat 14..24 each number 10 times) 100 times
day <- rep(rep(seq(14, 24), each = n_halfgal), times = n_sim)
count <- vector(mode = "logical", n_sim * n_halfgal * n_day)

#matrix with columns:
#  BT   half_gal    AT    day   count
data <- data.frame(BT, half_gal, AT, day, count)

#Now import the data from our input files and begin filling in our data frames
#input files
frequency_file <- "Frequency.csv"
growth_file <- "GrowthParameters.csv"
init_file <- "InitialCountsMPN.csv"

#Import frequency data and get the rpoB allelic type
freq_import <- read.csv(frequency_file, stringsAsFactors = FALSE, header = TRUE)
freq_data = freq_import$rpoB.allelic.type

#Import growth parameter data
growth_import <-read.csv(growth_file, stringsAsFactors = FALSE)

#Import initial count logMPN data
initialcount_import <- read.csv(init_file, stringsAsFactors = FALSE)
#MPN Column
initialcount_data = initialcount_import[,3]
#LOG MPN Column
initialcountlog_data = initialcount_import[,4]

# Calculate samples used in the monte carlo   ----

#Now sample the MPN distributions and the temperature distribution
#Sample logMPN from normal distribution
logMPN_mean <- c(-0.7226627)
logMPN_sd <- c(.9901429)
logMPN_samp = rnorm(n_sim, logMPN_mean, logMPN_sd)
MPN_samp = 10^logMPN_samp
MPN_samp_halfgal = MPN_samp * 1900 #MPN per half gallon (~1900 ml in half gallon)


#Sample temperature from normal distribution
#data from Consumer Phase Risk Assessment for Listeria monocytogenes in Deli Meats (2006)
temp_mean <- 4.096
temp_sd <-   2.381
temp_sample <- rnorm((n_sim), temp_mean, temp_sd)


#Generate initial MPN for each half gallon from Poisson distribution
#Also sample AT for each half gallon
MPN_init<-vector()
allele <- vector()
for (i in 1:n_sim){
  MPN_init_samp <-rep(rpois(n_halfgal, MPN_samp_halfgal[i]), times = n_day)
  MPN_init<-c(MPN_init, MPN_init_samp)
  allele_samp <- rep(sample(freq_data, n_halfgal, replace = T), times = n_day)
  allele <- c(allele, allele_samp)
}

#Convert MPN_init from half-gallon to mLs
MPN_init_mL <- MPN_init / 1900
#remove 0's from the data and replace with detection limit
MPN_init_mL[MPN_init_mL == 0] <- 0.01;

#Now we add in those calculations to our original dataframe
data$logMPN_init <- log10(MPN_init_mL) #Add initial logMPN to data frame
data$AT<-allele #Add in AT data

##Now we will calculate the log10N for each row in the data frame
##Get the AT and day from the data frame, get growth parameters depending on the AT
# Simulation   ----
for (i in 1:(n_sim *n_halfgal * n_day)){
  #Find row in growth parameter data that corresponds to allele sample
  allele_index <- which(growth_import$rpoBAT == data$AT[i]) 
  
  #calculate the new growth parameters using the square root model and our
  #sampled temperature
  newT <- temp_sample[data$BT[i]]
  newLag <- lagAtNewTemp(newT, growth_import$lag[allele_index])
  newMu <-  muAtNewTemp(newT, growth_import$mumax[allele_index])
  
  #Calculate the log10N count using our new growth parameters
  data$count[i] <- log10N_func(data$day[i], newLag, newMu,data$logMPN_init[i],growth_import$LOG10Nmax[allele_index])

}






