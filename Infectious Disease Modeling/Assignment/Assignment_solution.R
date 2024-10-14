###############################################################################################

#Exercise solution

###############################################################################################
library(deSolve)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)

## Answers to the question 1 and 2 ##
## SEIRS model (a model with latent state and waning immunity) ##
#1. Define model function
# SEIRS model with births and deaths
OpenSEIRS<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S + E + I + R
    sigma = 1/t_lat # 1/latent period
    
    #SIR w/ demography equations from lecture
    dS <- -beta*S*I/N + birth*N - death*S + omega*R
    dE <- beta*S*I/N - sigma*E - death*E
    dI <- sigma*E - death*I - gamma*I
    dR <- gamma*I - death*R - omega*R
    
    dC <- beta*S*I/N
    
    # return the rates of change as a list
    list(c(dS, dE, dI, dR, dC))
  })
}

#2. Define parameters and starting compartment sizes
parameters <- c(beta = 0.2, # place holder for effective contact rate (aka transmission rate)
                gamma = 1/14, #recovery rate (1/duration infection)
                birth = 120/10000/365,#0.00003, #birth rate (per capita)
                death = 120/10000/365, #all-cause mortality rate
                omega = 1/(30.5*3),#0.01, # waning immunity
                t_lat = 5 # latent period from E
                
)

beta_low <- 3*(parameters[["death"]] + parameters[["gamma"]])*(parameters[["death"]] + 1/parameters[["t_lat"]])*parameters[["t_lat"]]
beta_high <- 10*(parameters[["death"]] + parameters[["gamma"]])*(parameters[["death"]] + 1/parameters[["t_lat"]])*parameters[["t_lat"]]

state <- c(S = 10000-1, #population of 10,000, 1 person starts of infected
           E = 0, 
           I = 1, 
           R = 0,
           C = 0 #track cumulative number of infections
)


T_end <- 300 #run model for 300 time steps (e.g. days)
times <- seq(0, T_end, by = 1) #runs the model for 300 time steps (e.g. days), and computes output at each time step 

#Run the base-case with lower beta (R0=3)
parameters['beta'] <- beta_low
out_lowbeta <- data.frame(ode(y = state, times = times, func = OpenSEIRS, parms = parameters))
out_lowbeta$totI <- out_lowbeta$E + out_lowbeta$I
#Run the base-case with higher beta (R0=10)
parameters['beta'] <- beta_high
out_highbeta <- data.frame(ode(y = state, times = times, func = OpenSEIRS, parms = parameters))
out_highbeta$totI <- out_highbeta$E + out_highbeta$I

# Plot S,E,I,R over time
out_lowbeta_t <- melt(out_lowbeta, id.vars='time')
ggplot(out_lowbeta_t)+
  geom_line(aes(time,value,color=variable))+
  ggtitle("R0 = 3")+
  theme_bw()
max(out_lowbeta$C) #cumulative infections

out_highbeta_t <- melt(out_highbeta, id.vars='time')
ggplot(out_highbeta_t)+
  geom_line(aes(time,value,color=variable))+
  ggtitle("R0 = 10")+
  theme_bw()
max(out_highbeta$C) #cumulative infections

# Plot E + I over time & peak of epidemics
totI <- data.frame(time = out_lowbeta$time, lowbeta = out_lowbeta$totI, highbeta = out_highbeta$totI)
ggplot(totI)+
  geom_line(aes(time,lowbeta, color='low beta'))+
  geom_line(aes(time,highbeta, color='high beta'))+  
  ylab("E + I")+
  geom_vline(xintercept = c(which.max(totI$lowbeta),which.max(totI$highbeta)),
             linetype = c("dashed","dashed"),
             color = c("#0072B2","#D55E00"))+
  annotate(geom="text",
           label = c(as.character(which.max(totI$lowbeta)),as.character(which.max(totI$highbeta))),
           x = c(which.max(totI$lowbeta),which.max(totI$highbeta)),
           y = c(max(totI$lowbeta), max(totI$highbeta)),
           color = c("#0072B2","#D55E00"),
           hjust = -0.5)+
  ggtitle("Total number of infections (E+I)")+
  theme_bw()+
  scale_color_manual(values = c("#D55E00","#0072B2"))
# Plot only I
onlyI <- data.frame(time = out_lowbeta$time, lowbeta = out_lowbeta$I, highbeta = out_highbeta$I)
ggplot(onlyI)+
  geom_line(aes(time,lowbeta, color='low beta'))+
  geom_line(aes(time,highbeta, color='high beta'))+  
  ylab("Infectious population")+
  geom_vline(xintercept = c(which.max(onlyI$lowbeta),which.max(onlyI$highbeta)),
             linetype = c("dashed","dashed"),
             color = c("#0072B2","#D55E00"))+
  annotate(geom="text",
           label = c(as.character(which.max(onlyI$lowbeta)),as.character(which.max(onlyI$highbeta))),
           x = c(which.max(onlyI$lowbeta),which.max(onlyI$highbeta)),
           y = c(max(onlyI$lowbeta), max(onlyI$highbeta)),
           color = c("#0072B2","#D55E00"),
           hjust = -0.5)+
  ggtitle("Total number of infectious individuals (I)")+
  theme_bw()+
  scale_color_manual(values = c("#D55E00","#0072B2"))