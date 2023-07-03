library(dplyr)
library(ggplot2)
library(zoo)

# define the path of the csv file
file <- "cases_and_recovery.csv"  # has to be saved in the same folder

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

# filter for the data from the 4th COVID wave
data_wave4 <- subset(data_rec_and_inf_known, 
                     X..Refdatum.. >= as.POSIXct('2021-10-13 00:00') & 
                       X..Refdatum.. <= as.POSIXct('2022-01-04 00:00'))

# compute the time difference between the begin of the symptoms and the report 
# to the RKI. For this, we need our own function as we need to use difftime() 
# and as.numeric() because we want to obtain the number of days instead of a
# difftime-object.
MyDifftime <- function(var1, var2){
  x <- as.numeric(difftime(var1, var2))
  return(x)
}
time_difference <- mapply(MyDifftime, data_wave4$X..Meldedatum.., data_wave4$X..Refdatum..)

# add the time differences as a column to our data frame for wave 4
data_wave4_diff <- data_wave4
data_wave4_diff$TimeDelta <- time_difference

# ISSUE: the column TimeDelta does not show how long certain individuals have 
# been sick. This information is not available from the csv file.

# anyway, we can already do a lot by just assuming that the time until recovery
# is equal to 14 days.
# sort the entries of the 4th wave data by start of symptoms and add the numbers 
# of cases where the symptoms started on the same day.
new_cases_munich <- data_wave4_diff %>%
  group_by(X..Refdatum..) %>%
  summarise(across(X..AnzahlFall.., sum))

# smooth the data for the daily new detected cases - careful, the newly created 
# vector does not have the same length as the old one, hence we create a new 
# data table for this one
DifferenzVortagFaelleSmooth7 <- zoo::rollmean(new_cases_munich$X..AnzahlFall.., k=7)
timeline <- 1:length(DifferenzVortagFaelleSmooth7)
smoothdata <- data.table(timeline, DifferenzVortagFaelleSmooth7)

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

# we choose a recovery time of 14 days as proposed by the RKI 
# (see https://github.com/robert-koch-institut/SARS-CoV-2-Infektionen_in_Deutschland#hinweis-zu-genesenen)
# and assume that on 12.10.2021, there were a total of 0 infected individuals (I = 0).
# Of course this is a simplifying assumption!
curr_inf_munich <- currently_infected(14, 0, new_cases_munich$X..AnzahlFall..)
new_cases_munich$AktuellInfizierte <- curr_inf_munich

# compute the cumulative cases
cumulative_infections <- cumsum(new_cases_munich$X..AnzahlFall..)
new_cases_munich$KumulierteFaelle <- cumulative_infections

# PLOT FOR CURRENTLY ACTIVE INFECTIONS (I IN THE SIR MODEL)
ggplot(new_cases_munich, aes(x=X..Refdatum.., y=AktuellInfizierte, colour="data")) + 
  geom_point() +
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="currently infected people in Munich", 
       title="Active COVID-19 cases in Munich")

# PLOT FOR DAILY NEW CASES
ggplot(new_cases_munich, aes(x=X..Refdatum.., y=X..AnzahlFall.., colour="data")) + 
  geom_point() + 
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="newly infected people in Munich", 
       title="Daily detected COVID-19 cases in Munich")

# PLOT FOR CUMULATIVE CASES
ggplot(new_cases_munich, aes(x=X..Refdatum.., y=KumulierteFaelle, colour="data")) + 
  geom_point() + 
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="total detected cases in Munich", 
       title="Cumulative detected COVID-19 cases in Munich")

# PLOT FOR THE 7-DAY MOVING AVERAGE ON THE NEW CASES
ggplot(smoothdata, aes(x=timeline, y=DifferenzVortagFaelleSmooth7, colour="data")) + 
  geom_point() + 
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="newly infected people in Munich", 
       title="Daily detected COVID-19 cases in Munich (using a 7-day-mean value)")