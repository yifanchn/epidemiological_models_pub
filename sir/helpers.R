## SIR Model without population dynamics
## https://kinglab.eeb.lsa.umich.edu/480/nls/de.html

library(deSolve)
library(data.table)
library(tidyr)
library(ggplot2)

# parameter values
#N <- 1484226
N <- 1500000
dt <- fread("data/new_cases_munich.csv")
dt[, V1 := NULL]

plot_model <- function(parms, model, real){
  if (real == FALSE){
    #times <- seq(from=0,to=300,by=0.1)
    if (model == 1){
      parms <- parms[1:2]
      plot_sir(parms)
    } else if (model == 2){
      plot_sirs(parms)
    } else {
      plot_seir(parms)
    }
  } else {
    #times <- seq(from=1,to=84,by=0.1)
    if (model == 1){
      parms <- parms[1:2]
      plot_sir_real(parms)
    } else if (model == 2){
      plot_sirs_real(parms)
    } else {
      plot_seir_real(parms)
    }
  }
}

## Helper functions
## With real data for I
plot_sir_real <- function(parms){
  times <- seq(from=1,to=84,by=0.1)
  sir_model <- function(t, x, params){
    ## define state variables
    S <- x[1]
    I <- x[2]
    R <- x[3]
    ## extract parameter values
    beta <- params["beta"]
    alpha <- params["alpha"]
    N <- S+I+R
    ## model equations
    dSdt <- -beta*S*I/N
    dIdt <- beta*S*I/N-alpha*I
    dRdt <- alpha*I
    dxdt <- c(dSdt,dIdt,dRdt)
    list(dxdt)
  }
  
  I0 <- dt[time == 1, AktuellInfizierte]
  xstart = c(S=N-I0,I=I0,R=0)
  
  result <- as.data.table(ode(
    func=sir_model,
    y=xstart,
    times=times,
    parms=parms
  ))
  
  ggplot()+
    geom_line(data = result, aes(x=times, y=I, color="model I(t)"), size=1.5)+
    geom_point(data = dt, aes(x=time, y=AktuellInfizierte, color="real I(t)"))+
    theme_classic()+
    labs(x='time (in days)',y='number of individuals', colour="")+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=15), #change legend title font size
          legend.text = element_text(size=15)) #change legend text font size
}

plot_sirs_real <- function(parms){
  times <- seq(from=1,to=84,by=0.1)
  sirs_model <- function(t, x, params){
    ## define state variables
    S <- x[1]
    I <- x[2]
    R <- x[3]
    ## extract the parameters
    beta <- params["beta"]
    alpha <- params["alpha"]
    gamma <- params["gamma"]
    #N <- S+I+R
    ## model equations
    dSdt <- -beta*S*I/N +gamma*R
    dIdt <- beta*S*I/N-alpha*I
    dRdt <- alpha*I-gamma*R
    dxdt <- c(dSdt,dIdt,dRdt)
    
    list(dxdt)
  }
  
  # parameter values
  #parms <- c(beta=2e-3, alpha=1/10, gamma=1/30)
  #times <- seq(from=0,to=30,by=0.1)
  I0 <- dt[time == 1, AktuellInfizierte]
  xstart = c(S=N-I0,I=I0,R=0)
  
  result <- as.data.table(ode(
    func=sirs_model,
    y=xstart,
    times=times,
    parms=parms
  ))
  
  ggplot()+
    geom_line(data = result, aes(x=times, y=I, color="model I(t)"), size=1.5)+
    geom_point(data = dt, aes(x=time, y=AktuellInfizierte, color="real I(t)"))+
    theme_classic()+
    labs(x='time (in days)',y='number of individuals', colour="")+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=15), #change legend title font size
          legend.text = element_text(size=15)) #change legend text font size
}

plot_seir_real <- function(parms){
  times <- seq(from=1,to=84,by=0.1)
  seir_model <- function (t, x, params) {
    ## define state variables
    S <- x[1]
    E <- x[2]
    I <- x[3]
    R <- x[4]
    
    ## extract parameter values
    beta <- params["beta"]
    gamma <- params["gamma"]
    alpha <- params["alpha"]
    #N <- S+E+I+R
    
    ## model equations
    dSdt <- -beta*I*S/N
    dEdt <- beta*I*S/N-gamma*E
    dIdt <- gamma*E-alpha*I
    dRdt <- alpha*I
    dxdt <- c(dSdt, dEdt, dIdt, dRdt)
    
    list(dxdt)
  }
  
  #parms <- c(lambda=4, alpha=1/13, gamma=1/5)
  #times <- seq(from=0,to=30,by=0.1)
  I0 <- dt[time == 1, AktuellInfizierte]
  E0 <- dt[time == 1, AnzahlFall]
  xstart = c(S=N-I0-E0,E=E0,I=I0, R=0)
  
  result <- as.data.table(ode(
    func=seir_model,
    y=xstart,
    times=times,
    parms=parms
  ))
  
  ggplot()+
    geom_line(data = result, aes(x=times, y=I, color="model I(t)"), size=1.5)+
    geom_point(data = dt, aes(x=time, y=AktuellInfizierte, color="real I(t)"))+
    theme_classic()+
    labs(x='time (in days)',y='number of individuals')+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=15), #change legend title font size
          legend.text = element_text(size=15)) #change legend text font size
}




## Theoretical models only
plot_sir <- function(parms){
  times <- seq(from=0,to=300,by=0.1)
  sir_model <- function(t, x, params){
    ## define state variables
    S <- x[1]
    I <- x[2]
    R <- x[3]
    ## extract parameter values
    beta <- params["beta"]
    alpha <- params["alpha"]
    N <- S+I+R
    ## model equations
    dSdt <- -beta*S*I/N
    dIdt <- beta*S*I/N-alpha*I
    dRdt <- alpha*I
    dxdt <- c(dSdt,dIdt,dRdt)
    list(dxdt)
  }
  
  xstart = c(S=N-1,I=1,R=0)
  
  result <- as.data.table(ode(
    func=sir_model,
    y=xstart,
    times=times,
    parms=parms
  ))
  
  result %>%
    gather(variable,value,-time) %>%
    ggplot(aes(x=time,y=value,color=variable))+
    geom_line(size=1.5)+
    theme_classic()+
    labs(x='time (in days)',y='number of individuals', color="state")+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=15), #change legend title font size
          legend.text = element_text(size=15)) #change legend text font size
  
}
  
plot_sirs <- function(parms){
  times <- seq(from=0,to=3500,by=0.1)
  sirs_model <- function(t, x, params){
    ## define state variables
    S <- x[1]
    I <- x[2]
    R <- x[3]
    ## extract the parameters
    beta <- params["beta"]
    alpha <- params["alpha"]
    gamma <- params["gamma"]
    #N <- S+I+R
    ## model equations
    dSdt <- -beta*S*I/N +gamma*R
    dIdt <- beta*S*I/N-alpha*I
    dRdt <- alpha*I-gamma*R
    dxdt <- c(dSdt,dIdt,dRdt)
    
    list(dxdt)
  }
  
  # parameter values
  #parms <- c(beta=2e-3, alpha=1/10, gamma=1/30)
  #times <- seq(from=0,to=30,by=0.1)
  xstart = c(S=N-1,I=1,R=0)
  
  result <- as.data.table(ode(
    func=sirs_model,
    y=xstart,
    times=times,
    parms=parms
  ))
  
  result %>%
    gather(variable,value,-time) %>%
    ggplot(aes(x=time,y=value,color=variable))+
    geom_line(size=1.5)+
    theme_classic()+
    labs(x='time (in days)',y='number of individuals', color="state")+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=15), #change legend title font size
          legend.text = element_text(size=15)) #change legend text font size
}

plot_seir <- function(parms){
  times <- seq(from=0,to=3500,by=0.1)
  seir_model <- function (t, x, params) {
    ## define state variables
    S <- x[1]
    E <- x[2]
    I <- x[3]
    R <- x[4]
    
    ## extract parameter values
    beta <- params["beta"]
    gamma <- params["gamma"]
    alpha <- params["alpha"]
    #N <- S+E+I+R
    
    ## model equations
    dSdt <- -beta*I*S/N
    dEdt <- beta*I*S/N-gamma*E
    dIdt <- gamma*E-alpha*I
    dRdt <- alpha*I
    dxdt <- c(dSdt, dEdt, dIdt, dRdt)
    
    list(dxdt)
  }
  
  #parms <- c(lambda=4, alpha=1/13, gamma=1/5)
  #times <- seq(from=0,to=30,by=0.1)
  xstart = c(S=N-1,E=0,I=1, R=0)
  
  result <- as.data.table(ode(
    func=seir_model,
    y=xstart,
    times=times,
    parms=parms
  ))
  
  result %>%
    gather(variable,value,-time) %>%
    ggplot(aes(x=time,y=value,color=variable))+
    geom_line(size=1.5)+
    theme_classic()+
    labs(x='time (in days)',y='number of individuals', color="state")+
    theme(legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=15), #change legend title font size
          legend.text = element_text(size=15)) #change legend text font size
}




