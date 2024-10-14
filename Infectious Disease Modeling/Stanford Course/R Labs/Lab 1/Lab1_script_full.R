###############################################################################################

#This file contains code corresponding to all of Lab #1/Homework #1
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

#6. 
#calculate R0
R0 <- parameters["beta"]/parameters["gamma"]
print(R0)

#model-predicted prevalence at t=50
output_df <- data.frame(output)
output_df$N <- output_df$S + output_df$I + output_df$R
print(output_df$I[50]/output_df$N[50])

S_star <- output_df$S[T_end+1]/output_df$N[T_end+1] #S prevalence at equilibrium
Rt_star <- R0*S_star #effective reproductive number
I_star <- output_df$I[T_end+1]/output_df$N[T_end+1] #I prevalence at equilibrium

print(Rt_star)
print(I_star)

#7.
#p=20%
parameters <- c(beta = 2, #2 is 10*20%
                gamma = 0.3
)
R0 <- parameters["beta"]/parameters["gamma"]
print(R0)
output <- ode(y = state, times = times, func = BasicSIR, parms = parameters)
show_SIR_model_results(output)


#what if p was 2% instead?
parameters <- c(beta = 0.2, #0.2 is 10*2%
                gamma = 0.3
)
R0 <- parameters["beta"]/parameters["gamma"]
print(R0)
output <- ode(y = state, times = times, func = BasicSIR, parms = parameters)
show_SIR_model_results(output)

#8.
parameters <- c(beta = 0.5,
                gamma = 0.0
)
output <- ode(y = state, times = times, func = BasicSIR, parms = parameters)
show_SIR_model_results(output)

#9. 
OpenSIR<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    N = S + I + R
    
    #SIR w/ demography equations from lecture
    dS <- -beta*S*I/N + birth*N - death*S + omega*R
    dI <- beta*S*I/N - death*I - gamma*I
    dR <- gamma*I - death*R - omega*R
    
    # return the rates of change as a list
    list(c(dS, dI, dR))
  })
}

parameters <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                gamma = 0.3, #recovery rate (1/duration infection)
                birth = 0.03, #birth rate (per capita)
                death = 0.03, #all-cause mortality rate
                omega = 0.0
)

#10.
R0 <- parameters["beta"]/(parameters["gamma"] + parameters["death"])
print(R0)

S_star <- 1/R0
I_star <- (parameters["death"]/parameters["beta"])*(R0-1)
R_star <- 1-(S_star + I_star)
print(c(S_star, I_star, R_star)) #prevalences at equilibrium

#11. 
output <- ode(y = state, times = times, func = OpenSIR, parms = parameters)
show_SIR_model_results(output)

output_df <- data.frame(output)
output_df$N <- output_df$S + output_df$I + output_df$R
print(output_df$I[T_end+1]/output_df$N[T_end+1])

#to "zoom in" on the I compartment (to observe oscillations)
show_SIR_I_output <- function(output) {
  df1 <- data.frame(output)
  ggplot() + 
    geom_line(data = df1, aes(x = time, y = I), color="red") +
    xlab('Time') +
    ylab('Count') 
}

show_SIR_I_output (output)

#12.
phase_diagram <-function(output) {
  #convert model output to a data frame (type of object in R) and reshape it so that each row represents a time-compartment combination (e.g. how many ppl in S compartment at time 10) - the latter just makes it easier to plot your results
  df1 <- data.frame(output)
  df1$N <- df1$S + df1$I + df1$R
  df1[,2:ncol(df1)] <- df1[,2:ncol(df1)]/df1$N #divide by pop size to calculate pop prevalence
  ggplot() + 
    geom_point(data = df1, aes(x = S*100, y = I*100, color=time)) + #plot S_prev against I_prev
    xlab('Susceptible Prevalence (%)') + ylab('Infected Prevalence (%)') + #label x and y-axes
    theme_bw() #this isn't necessary but makes the graph look a bit nicer
} 
phase_diagram(output)

#13. 
A <- 1/(parameters["death"]*(R0-1))
print(A)

G <- 1/(parameters["death"]+parameters["gamma"]) #duration of infectiousness
T <- 2*pi*sqrt(A*G)
print(T)

#14. 
#plotting function that is technically more complex but can be applied to all model output (not just SIR) - you could also just use show_SIR_model_results() for this as before.
show_model_results<-function(output) {
  #convert model output to a data frame (type of object in R) and reshape it so that each row represents a time-compartment combination (e.g. how many ppl in S compartment at time 10) - the latter just makes it easier to plot your results
  df1 <- data.frame(output)
  df2 <- reshape(df1, varying=colnames(df1)[2:ncol(df1)], 
                 v.name="sizes",  timevar="Compartment",
                 times=colnames(df1)[2:ncol(df1)], idvar="time", direction="long")
  ggplot() + 
    geom_line(data = df2, aes(x = time, y = sizes, color=Compartment)) + #plot the number of ppl in each compartment over time
    xlab('Time') + ylab('Count') + #label x and y-axes
    theme_bw() #this isn't necessary but makes the graph look a bit nicer
} 

#original (for comparison)
parameters <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                gamma = 0.3, #recovery rate (1/duration infection)
                birth = 0.03, #birth rate (per capita)
                death = 0.03, #all-cause mortality rate
                omega = 0.0
)
output <- ode(y = state, times = times, func = OpenSIR, parms = parameters)
show_model_results(output)
phase_diagram(output)

parameters["beta"] <- 2.0
output <- ode(y = state, times = times, func = OpenSIR, parms = parameters)
show_model_results(output)
phase_diagram(output)

parameters["beta"] <- 0.5
parameters["gamma"] <- 0.2
output <- ode(y = state, times = times, func = OpenSIR, parms = parameters)
show_model_results(output)
phase_diagram(output)

parameters["gamma"] <- 0.3
parameters["birth"] <- 0.1
parameters["death"] <- 0.1
output <- ode(y = state, times = times, func = OpenSIR, parms = parameters)
show_model_results(output)
phase_diagram(output)

#16. 
parameters <- c(beta = 0.5, #effective contact rate (aka transmission rate)
                gamma = 0.3, #recovery rate (1/duration infection)
                birth = 0.03, #birth rate (per capita)
                death = 0.03, #all-cause mortality rate
                omega = 0.1
)
output <- ode(y = state, times = times, func = OpenSIR, parms = parameters)
show_model_results(output)

output_df <- data.frame(output)
output_df$N <- output_df$S + output_df$I + output_df$R
print(output_df$I[T_end+1]/output_df$N[T_end+1])
