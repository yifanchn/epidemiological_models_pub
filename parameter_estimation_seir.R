# parameter estimation for the SEIR model 
# (as seen in Stefanie Fuderer's Master's Thesis)

# use libraries
library(readxl)
library(data.table)
library(dplyr)
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

# generate new column with the row number (important for the parameter estimation)
datf_wave4$time <- seq.int(nrow(datf_wave4))

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
curr_inf <- currently_infected(10, 0, datf_wave4plus$DifferenzVortagFaelle)

# we add this new vector to our data table
datf_wave4$AktuellInfizierte <- curr_inf[-(1:12)]

# the number of exposed individuals at 13.10. assuming gamma = 1/5.5
E_start <- sum(datf_wave4plus$DifferenzVortagFaelle[8:12]) + 0.5*datf_wave4plus$DifferenzVortagFaelle[7]
I_start <- datf_wave4$AktuellInfizierte[1]
N <- 8e7

intervention_period <- c(43)  # stricter restrictions starting from Nov 24
# contact_reduction <- c(0.5)  # after intervention: beta = contact_reduction*beta

# Define the SIR model
model_SEIR <- function(t, x, params){
  S <- x[1]
  E <- x[2]
  I <- x[3]
  R <- x[4]
  # parameters to be fitted
  lambda0 <- params["lambda0"]
  contact_red <- params["contact_red"]
  # vectorized if-else: condition, true, false
  lambda_t <- if_else(t <= intervention_period[1],
                    lambda0,
                    lambda0 * contact_red)
                    #lambda0 * contact_reduction[1])
  # code the model equations; we already estimated the recovery rate alpha to be
  # equal to 1/10 (10 days needed for recovery), and the duration of the latent
  # phase is estimated to be about 5.5 days, hence gamma = 1/5.5
  dSdt <- -lambda_t * S * I / N
  dEdt <- lambda_t * S * I / N - (1/5.5) * E
  dIdt <- (1/5.5) * E - (1/10) * I
  dRdt <- (1/10)*I
  # return result as a list
  list(c(dSdt, dEdt, dIdt, dRdt))
}

# calculate the sum of squares error
SSE <- function(parameters){
  names(parameters) <- c("lambda0", "contact_red") 
  times <- seq(from=1, to=83, by=1)
  xstart <- c(S=N-I_start-E_start, E=E_start, I=I_start, R=0)
  out <- as.data.table(ode(func = model_SEIR,
                           y = xstart,
                           times = times,
                           parms = parameters))
  data_model <- out[, c("time", "I")]
  data_real <- datf_wave4[, c("time", "AktuellInfizierte")]
  data_validation <- merge(data_model, data_real, by="time", all=FALSE)
  SSE <- sum((data_validation[, I] - data_validation[, AktuellInfizierte])^2)
  return(SSE)
}

# optimize the parameter beta with the L-BFGS-B method
opt <- optim(c(lambda0=0.1, contact_red=0.1), SSE, method="L-BFGS-B", lower=c(0, 0), upper=c(3, 1))

opt_par <- setNames(opt$par, c("lambda0", "contact_red"))

# compute the AIC
SSE_val <- opt$value
n <- length(datf_wave4$time)
# k <- 1  # one parameter to be fitted, the other ones are estimated
k <- 2
AIC <- n * log(SSE_val / n) + 2 * (k + 1)

# output of the SIR model with fitted parameters
times <- seq(from=1, to=83, by=1)
x_start <- c(S=N-I_start-E_start, E=E_start, I=I_start, R=0)

out <- as.data.table(ode(func=model_SEIR, 
                         y=x_start, 
                         times=times, 
                         #parms=setNames(c(0.0005), c("lambda0"))))
                         parms=opt_par))

# plot of currently active infections (PLOT OF I)
ggplot() +
  geom_point(data=datf_wave4, aes(x=time, y=AktuellInfizierte, colour="data")) +
  geom_line(data=out, aes(x=time, y=I, colour="model (SEIR)")) +
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="currently infected people in Germany", 
       title="Active COVID-19 cases in Germany") + 
  theme(legend.title=element_blank())

# plot of daily new infections
day1_new <- datf_wave4$DifferenzVortagFaelle[1]
out$casesperday <- out[, abs(S - c((S[1] -day1_new), S[1:length(S) - 1]))]

ggplot() +
  geom_point(data=datf_wave4, aes(x=time, y=DifferenzVortagFaelle, colour="data")) +
  geom_line(data=out, aes(x=time, y=casesperday, colour="model (SEIR)")) +
  labs(x="days between 2021/10/13 and 2022/01/04",
       y="newly infected people in Germany",
       title="Daily detected COVID-19 cases in Germany") +
  theme(legend.title=element_blank())

# plot of cumulated cases
out$cumulatedcases <- cumsum(out$casesperday)
datf_wave4$KumuliertNurVierteWelle <- cumsum(datf_wave4$DifferenzVortagFaelle)

ggplot() +
  geom_point(data=datf_wave4, aes(x=time, y=KumuliertNurVierteWelle, colour="data")) +
  geom_line(data=out, aes(x=time, y=cumulatedcases, colour="model (SEIR)")) +
  labs(x="days between 2021/10/13 and 2022/01/04",
       y="total detected cases in Germany",
       title="Total detected COVID-19 cases in Germany") +
  theme(legend.title=element_blank())