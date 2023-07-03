## SEIR model without population dynamics
## Further reading can be found here: https://kinglab.eeb.lsa.umich.edu/480/nls/de.html

library(deSolve)
library(data.table)
library(tidyr)
library(ggplot2)

# parameter values
parms <- c(beta=0.5, alpha=1/13, gamma=1/5)
times <- seq(from=0,to=30,by=0.1)
xstart <- c(S=999,E=0,I=1,R=0)

plot_seir <- function(parms, times, xstart){
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
    N <- S+E+I+R
    
    ## model equations
    dSdt <- -beta*S*I/N
    dEdt <- beta*S*I/N-gamma*E
    dIdt <- gamma*E-alpha*I
    dRdt <- alpha*I
    dxdt <- c(dSdt, dEdt, dIdt, dRdt)
    
    list(dxdt)
  }
  
  result <- as.data.table(ode(
    func=seir_model,
    y=xstart,
    times=times,
    parms=parms
  ))
  
  result %>%
    gather(variable,value,-time) %>%
    ggplot(aes(x=time,y=value,color=variable))+
    geom_line(size=1)+
    theme_classic()+
    labs(x='time (in days)',y='number of individuals', color="state")
}

