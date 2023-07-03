# epidemiological_models_pub

What to discover in this repository:

## Folder sir
This folder contains an RShiny app which allows you to play around with the parameters of the SIR, SEIR and SIRS model. You can also compare the number of actively infected individuals I(t) from the models with the actual number from the 4th COVID wave in Munich.
- data: folder in which a csv-file with the daily new infections in Munich can be found. This folder needs to be saved in the same directory as app.R and helpers.R for the app to run without problems.
- app.R: contains the UI related features of the app. If you want to run the app, you need to open app.R in RStudio and click "Run App" in the upper right corner of the editor window.
- helpers.R: the functions that are necessary for solving and plotting the ODEs.

## Folder data_plots
This folder contains two data files and two scripts that give you some insights into the data that is used in our data fitting examples.
- Fallzahlen_Gesamtuebersicht.xlsx and cases_and_recovery.csv: data of COVID cases starting in 2020, for Germany and Munich, respectively. cases_and_recovery.csv is much more detailed and differs between the date of report of an infection and its actual outbreak. If you would like to learn more about this data, please refer to https://github.com/robert-koch-institut/SARS-CoV-2-Infektionen_in_Deutschland.
- plot_incidence.R: creates plots with the data from Germany, 4th wave. You obtain a plot of the currently active infections (= I(t) in the infection models), the daily new cases, the cumulated cases and a plot of a 7-day moving average on the daily new cases.
- work_with_munich_data.R: does the same as plot_incidence.R, but with the data from Munich.

## The ODE models SIR, SEIR and SIRS
If you want to get a basic understanding of how the ODE models are implemented and solved using R, you should have a look at the following files:
- sir_closed.R: contains a function that solves the ODE of the basic SIR model and plots it. Specify your parameter values, the timeline and the initial values using the variables `parms`, `times` and `xstart`, run the script, and enter `plot_sir(parms, times, xstart)` in the console to obtain the plot.
- seir.R: same as sir_closed.R, but for the SEIR model.
- sirs.R: analogous.

## The data fitting files
For all of the above models, data fitting is conducted using the SSE method for both the data from Germany and Munich (4th COVID wave).
- parameter_estimation_sir.R, ...seir.R, ...sirs.R: uses the data from Germany, fits the parameters of the respective model and generates some plots. Assumes that on November 24th, new contact restrictions were implemented by the gouvernment. Also computes the AIC for each model.
- parameter_estimation_smooth_data_sir.R: uses smoothed data from Germany and fits `beta` from the basic SIR model.
- munich_parameter_estimation_sir.R, ...seir.R, ...sirs.R: same as the parameter estimators for Germany, but with data from Munich. munich_parameter_estimation_sir.R also contains a plot with a log-scaled y-axis for illustration.
