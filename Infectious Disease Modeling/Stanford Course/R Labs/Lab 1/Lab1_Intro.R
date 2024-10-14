###############################################################################################

#This file contains code corresponding to parts 1-5 of Lab #1/Homework #1
#Stanford University HUMBIO 154D/HRP 204
#Spring Quarter 2020

###############################################################################################

#1. Load packages
#if you are installing packages for the first time (only run once):
#install.packages("deSolve")
#install.packages("ggplot2")

#if you've already installed the packages:
library(deSolve) #differential equation solver
library(ggplot2) #plotting package

#2. Start by defining:
  # A vector of model parameters (e.g. effective contact rate, recovery rate)
  # A vector of initial compartment sizes (often 1 infected person and the rest susceptible)
  # A vector of time steps corresponding to how long we want to run the model.

parameters <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                gamma = 0.3 #recovery rate (1/duration infection)
)

state <- c(S = 99999, #population of 100,000, 1 person starts of infected
           I = 1, 
           R = 0
)

T_end <- 500 #run model for 500 time steps
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps, and computes output at each time step 

#3. Define function for a basic SIR model without demography
# This will be used with the deSolve package to simulate how your population moves between compartments over time

BasicSIR<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{ #this tells R that "S" refers to the "S" in the "state" vector, "beta" refers to the "beta" in "parameters", etc.
    
    N = S + I + R #define N (total population size)
    
    #SIR model equations from lecture - rates of change in and out of each compartment 
    dS <- -beta*S*I/N
    dI <- beta*S*I/N - gamma*I
    dR <- gamma*I
    
    #return the rates of change as a list
    list(c(dS, dI, dR)) 
  })
}

#4. Run the SIR model using the ode() function that is part of the deSolve package.
output <- ode(y = state, times = times, func = BasicSIR, parms = parameters)

#Note: ode() is a function in the deSolve package. 
# You must fill it in using the syntax: 
#    "ode(y= [vector of initial compartment sizes], 
#         times = [vector of time steps], 
#         func = [name of SIR model function], 
#         parms = [vector or list of parameter values])"  

#5. Examine the model output by:
#   typing "View(output)" in the console, 
#   clicking on "output" in the Environment (top right of RStudio), OR
#   printing model output to the console.
# It is often easier to examine model output by plotting. 
# We define a function that will plot the epidemic curve of our model output from step #4. 

View(output)

print(output[200:204,])
print(head(output))
print(tail(output))

show_SIR_model_results<-function(output) {
  df1 <- data.frame(output)
  ggplot() + 
    geom_line(data = df1, aes(x = time, y = S), color = "blue") +
    geom_line(data = df1, aes(x = time, y = I), color = "red") +
    geom_line(data = df1, aes(x = time, y = R), color = "green") +
    xlab('Time') +
    ylab('Count') 
}

show_SIR_model_results(output)
