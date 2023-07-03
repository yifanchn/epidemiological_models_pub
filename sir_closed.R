## SIR Model without population dynamics
## Further reading can be found here: https://kinglab.eeb.lsa.umich.edu/480/nls/de.html

library(deSolve)
library(data.table)
library(tidyr)
library(ggplot2)

## parameter values
parms <- c(beta=0.5, alpha=1/13)
times <- seq(from=0,to=30,by=0.1)
xstart <- c(S=999,I=1,R=0)

plot_sir <- function(parms, times, xstart){
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
  
  result <- as.data.table(ode(
    func=sir_model,
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
  


