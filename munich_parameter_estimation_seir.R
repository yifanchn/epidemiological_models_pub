# parameter estimation for the SIRS model using data from Munich
# (as seen in Stefanie Fuderer's Master's Thesis)

# use libraries
library(readxl)
library(data.table)
library(dplyr)
library(deSolve)
library(ggplot2)

# define the path of the csv file
file <- "cases_and_recovery.csv"  # has to be locally available

# load data from the csv file
data_from_csv <- read.csv(file, sep=",", quote="")

# find out how R saved the column names
# colnames(data_from_csv)

# converting the dates to POSIXct format (= "real" R date format)
data_from_csv$X..Meldedatum.. <- as.POSIXct(data_from_csv$X..Meldedatum..)
data_from_csv$X..Refdatum.. <- as.POSIXct(data_from_csv$X..Refdatum..)


# extract only the rows of the data frame, where people recovered and the date
# of recovery is known (date format for POSIXct-queries: YYYY-mm-dd hh:mm)
data_rec_and_inf_known <- subset(data_from_csv, X..AnzahlGenesen... > 0 & X..IstErkrankungsbeginn.. == 1)

# data from the 4th wave + the 2 weeks prior (important for the calculation of 
# the currently infected individuals)
data_wave4plus <- subset(data_rec_and_inf_known, 
                         X..Refdatum.. >= as.POSIXct('2021-10-01 00:00') & 
                           X..Refdatum.. <= as.POSIXct('2022-01-04 00:00'))

# filter for the data from the 4th COVID wave
data_wave4 <- subset(data_rec_and_inf_known,
                     X..Refdatum.. >= as.POSIXct('2021-10-13 00:00') &
                       X..Refdatum.. <= as.POSIXct('2022-01-04 00:00'))

# sort the entries of the 4th wave data by start of symptoms and add the numbers
# of cases where the symptoms started on the same day.
new_cases_munich_plus <- data_wave4plus %>%
  group_by(X..Refdatum..) %>%
  summarise(across(X..AnzahlFall.., sum))

new_cases_munich <- data_wave4 %>%
  group_by(X..Refdatum..) %>%
  summarise(across(X..AnzahlFall.., sum))

new_cases_munich$time <- seq.int(nrow(new_cases_munich))

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
curr_inf <- currently_infected(10, 0, new_cases_munich_plus$X..AnzahlFall..)

# we add this new vector to our data table
new_cases_munich$AktuellInfizierte <- curr_inf[-(1:12)]

I_start <- new_cases_munich$AktuellInfizierte[1]
N <- 1500000

intervention_period <- c(43)  # stricter restrictions starting from Nov 24; note 
# that there was another tightening on Dec 04, which would be day 53

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
  data_real <- new_cases_munich[, c("time", "AktuellInfizierte")]
  data_validation <- merge(data_model, data_real, by="time", all=FALSE)
  SSE <- sum((data_validation[, I] - data_validation[, AktuellInfizierte])^2)
  return(SSE)
}

# optimize the parameter beta with the L-BFGS-B method
opt <- optim(c(lambda0=0.1, contact_red=0.1), SSE, method="L-BFGS-B", lower=c(0, 0), upper=c(3, 1))

opt_par <- setNames(opt$par, c("lambda0", "contact_red"))

# compute the AIC
SSE_val <- opt$value
n <- length(new_cases_munich$time)
k <- 2  # two parameters to be fitted
AIC <- n * log(SSE_val / n) + 2 * (k + 1)

# output of the SIR model with fitted parameters
times <- seq(from=1, to=83, by=1)
x_start <- c(S=N-I_start-E_start, E=E_start, I=I_start, R=0)

out <- as.data.table(ode(func=model_SEIR, 
                         y=x_start, 
                         times=times, 
                         parms=opt_par))

# plot of currently active infections (PLOT OF I)
ggplot() +
  geom_point(data=new_cases_munich, aes(x=time, y=AktuellInfizierte, colour="data")) +
  geom_line(data=out, aes(x=time, y=I, colour="model (SEIR)")) +
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="currently infected people in Munich", 
       title="Active COVID-19 cases in Munich") + 
  theme(legend.title=element_blank())

# plot of daily new infections
day1_new <- new_cases_munich$X..AnzahlFall..[1]
out$casesperday <- out[, abs(S - c((S[1] -day1_new), S[1:length(S) - 1]))]

ggplot() +
  geom_point(data=new_cases_munich, aes(x=time, y=X..AnzahlFall.., colour="data")) +
  geom_line(data=out, aes(x=time, y=casesperday, colour="model (SEIR)")) +
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="newly infected people in Munich", 
       title="Daily detected COVID-19 cases in Munich") + 
  theme(legend.title=element_blank())

# plot of cumulated cases
out$cumulatedcases <- cumsum(out$casesperday)
new_cases_munich$KumuliertNurVierteWelle <- cumsum(new_cases_munich$X..AnzahlFall..)

ggplot() +
  geom_point(data=new_cases_munich, aes(x=time, y=KumuliertNurVierteWelle, colour="data")) +
  geom_line(data=out, aes(x=time, y=cumulatedcases, colour="model (SEIR)")) +
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="total detected cases in Germany", 
       title="Total detected COVID-19 cases in Germany") + 
  theme(legend.title=element_blank())
