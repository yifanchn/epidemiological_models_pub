# install library if not already done
# install.packages("readxl")

# use libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(data.table)
library(zoo)

# define the file we want to read (as.character has to be used because otherwise 
# R is confused)
file <- "Fallzahlen_Gesamtuebersicht.xlsx"  # has to be saved in the same directory

# read the excel file with the data in it - the table starts in line 3, hence we 
# skip the first 2 rows. "datf" stands for data frame. 
datf <- read_excel(file, skip = 2)

setDT(datf)

# rename the columns in datf to get rid of " "
colnames(datf)[2] <- "AnzahlCOVID19Faelle"
colnames(datf)[3] <- "DifferenzVortagFaelle"

# get only the data from the 4th wave
datf_wave4 <- subset(datf, Berichtsdatum >= as.POSIXct('2021-10-13 00:00') & 
                       Berichtsdatum <= as.POSIXct('2022-01-04 00:00'))

# smooth the data for the daily new detected cases - careful, the newly created 
# vector does not have the same length as the old one, hence we create a new 
# data table for this one
DifferenzVortagFaelleSmooth7 <- zoo::rollmean(datf_wave4$DifferenzVortagFaelle, k=7)
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
curr_inf <- currently_infected(14, 0, datf_wave4$DifferenzVortagFaelle)

# we add this new vector to our data table
datf_wave4$AktuellInfizierte <- curr_inf

# PLOT FOR CURRENTLY ACTIVE INFECTIONS (I IN THE SIR MODEL)
# instead of Berichtsdatum, one can use 1:83
ggplot(datf_wave4, aes(x=Berichtsdatum, y=AktuellInfizierte, colour="data")) + 
  geom_point() + 
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="currently infected people in Germany", 
       title="Active COVID-19 cases in Germany")

# PLOT FOR DAILY NEW CASES
ggplot(datf_wave4, aes(x=Berichtsdatum, y=DifferenzVortagFaelle, colour="data")) + 
  geom_point() + 
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="newly infected people in Germany", 
       title="Daily detected COVID-19 cases in Germany")

# PLOT FOR CUMULATIVE CASES
ggplot(datf_wave4, aes(x=Berichtsdatum, y=AnzahlCOVID19Faelle, colour="data")) + 
  geom_point() + 
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="total detected cases in Germany", 
       title="Cumulative detected COVID-19 cases in Germany")

# PLOT FOR THE 7-DAY MOVING AVERAGE ON THE NEW CASES
ggplot(smoothdata, aes(x=timeline, y=DifferenzVortagFaelleSmooth7, colour="data")) + 
  geom_point() + 
  labs(x="days between 2021/10/13 and 2022/01/04", 
       y="newly infected people in Germany", 
       title="Daily detected COVID-19 cases in Germany (using a 7-day-mean value)")
