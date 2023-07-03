# parameter estimation for the basic SIR model and smooth data

# use libraries
library(readxl)
library(zoo)
library(data.table)
library(deSolve)
library(ggplot2)

# get the file with the data
file <- "Fallzahlen_Gesamtuebersicht.xlsx"  # has to be locally available

# read the excel file with the data in it - the table starts in line 3, hence we 
# skip the first 2 rows. "datf" stands for data frame. 
datf <- read_excel(file, skip = 2)

setDT(datf)

# rename the columns in datf to get rid of " "
colnames(datf)[2] <- "AnzahlCOVID19Faelle"
colnames(datf)[3] <- "DifferenzVortagFaelle"

# data from the 4th wave + the 2 weeks prior (important for the calculation of 
# the currently infected individuals)
datf_wave4plus <- subset(datf, Berichtsdatum >= as.POSIXct('2021-10-01 00:00') &
                           Berichtsdatum <= as.POSIXct('2022-01-04 00:00'))

# get only the data from the 4th wave
datf_wave4 <- subset(datf, Berichtsdatum >= as.POSIXct('2021-10-13 00:00') & 
                       Berichtsdatum <= as.POSIXct('2022-01-04 00:00'))

# smooth the data for the larger collection
DifferenzVortagFaelleSmooth7 <- zoo::rollmean(datf_wave4plus$DifferenzVortagFaelle, k=7)
time <- 1:length(DifferenzVortagFaelleSmooth7)
smoothdataPLUS <- data.table(time, DifferenzVortagFaelleSmooth7)

# smooth the data for the daily new detected cases - careful, the newly created 
# vector does not have the same length as the old one, hence we create a new 
# data table for this one
DifferenzVortagFaelleSmooth7 <- zoo::rollmean(datf_wave4$DifferenzVortagFaelle, k=7)
time <- 1:length(DifferenzVortagFaelleSmooth7)
smoothdata <- data.table(time, DifferenzVortagFaelleSmooth7)

# this function will help us calculate the number of currently infected people
currently_infected <- function(recovery_time, start_value, new_cases){
  # pre-allocate vector (less time consuming that appending values to empty vector)
  realI <- numeric(length(new_cases))  # realI is a vector filled with 0s
  realI[1] <- start_value + new_cases[1]
  for (i in 2:length(new_cases)){
    # if no one has recovered yet, we can just add the daily new cases
    if (i <= recovery_time){
      realI[i] = realI[i-1] + new_cases[i]
      # once people start to recover, we need to subtract them from the total I
    } else {
      realI[i] = realI[i-1] + new_cases[i] - new_cases[i-recovery_time]
    }
  }
  return(realI)
}

# we choose a recovery time of 10 days as Stefanie did
# and assume that on 30.09.2021, there were a total of 0 infected individuals (I = 0).
curr_inf <- currently_infected(10, 0, smoothdataPLUS$DifferenzVortagFaelleSmooth7)

# we add this new vector to our data table
smoothdata$AktuellInfizierte <- curr_inf[-(1:12)]

I_start <- smoothdata$DifferenzVortagFaelleSmooth7[1]
N <- 8e7

intervention_period <- c(43)  # stricter restrictions starting from Nov 24
contact_reduction <- c(0.5)  # after intervention: beta = contact_reduction*beta

# D e fi n e the SIR model
model_SIR <- function(t, x, params){
  S <- x[1]
  I <- x[2]
  R <- x[3]
  # parameters to be fitted
  beta <- params["beta0"]
  # vectorized if-else: condition, true, false
  beta_t <- if_else(t <= intervention_period[1],
                    beta0,
                    beta0 * contact_reduction[1])
  # code the model equations; we already estimated the recovery rate alpha to be
  # equal to 1/10 (10 days needed for recovery)
  dSdt <- -beta * x[1] * x[2] / N
  dIdt <- beta * x[1] * x[2] / N - (1/10) * x[2]
  dRdt <- (1/10) * x[2]
  # return result as a list
  list(c(dSdt, dIdt, dRdt))
}

# calculate the sum of squares error
SSE <- function(parameters){
  names(parameters) <- c("beta0") 
  times <- seq(from=1, to=83, by=1)
  xstart <- c(S=N-I_start, I=I_start, R=0)
  out <- as.data.table(ode(func = model_SIR,
                           y = xstart,
                           times = times,
                           parms = parameters))
  # calculate the number of daily detected cases from the model
  # out$casesperday <- out[, abs(S - c(N, S[1:length(S) - 1]))]
  data_model <- out[, c("time", "I")]
  data_real <- smoothdata[, c("time", "AktuellInfizierte")]
  data_validation <- merge(data_model, data_real, by="time", all=FALSE)
  SSE <- sum((data_validation[, I] - data_validation[, AktuellInfizierte])^2)
  return(SSE)
}

# optimize the parameter beta with the L-BFGS-B method
opt <- optim(c(beta=0.1), SSE, method="L-BFGS-B", lower=c(0), upper=c(3))

opt_par <- setNames(opt$par, c("beta"))

# compute the AIC
SSE_val <- opt$value
n <- length(datf_wave4$Berichtsdatum)
k <- 1  # one parameter to be fitted, the other one is estimated
AIC <- n * log(SSE_val / n) + 2 * (k + 1)

# output of the SIR model with fitted parameters
times <- seq(from=1, to=83, by=1)
x_start <- c(S=N-I_start, I=I_start, R=0)

out <- as.data.table(ode(func=model_SIR, 
                         y=x_start, 
                         times=times, 
                         parms=opt_par))

# plot the fitted model together with the real data 
# (PLOT FOR THE 7-DAY MOVING AVERAGE ON THE ACTIVE CASES)
ggplot() +
  geom_point(data=smoothdata, aes(x=time, y=AktuellInfizierte, colour="data")) +
  geom_line(data=out, aes(x=time, y=I, colour="model")) +
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="currently infected people in Germany", 
       title="Active COVID-19 cases in Germany (using a 7-day-mean value)") + 
  theme(legend.title=element_blank())

# plot of daily new infections
day1_new <- smoothdata$DifferenzVortagFaelleSmooth7[1]
out$casesperday <- out[, abs(S - c((S[1] - day1_new), S[1:length(S) - 1]))]

ggplot() +
  geom_point(data=smoothdata, aes(x=time, y=DifferenzVortagFaelleSmooth7, colour="data")) +
  geom_line(data=out, aes(x=time, y=casesperday, colour="model")) +
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="newly infected people in Germany", 
       title="Daily detected COVID-19 cases in Germany (using a 7-day-mean value)") + 
  theme(legend.title=element_blank())