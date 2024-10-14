###############################################################################################

#Dynamic Compartmental Models & Related Code - you will find these helpful for Lab 1/Homework 1
#Stanford University HUMBIO 154D/HRP 204
#Spring Quarter 2020

###############################################################################################

#ORGANIZATION OF THIS FILE
#Part 1: install/load packages (starts line 17)
#Part 2: define SIR and similar model functions (starts line 26)
#Part 3: define functions to plot model output (starts line 106)
#Part 4: code demonstration: defines model parameters, runs model, plots output (starts line 181)

###############################################################################################

#1. PACKAGES - you will need 2 packages: deSolve and ggplot2
#to install for the first time (you should only need to run this once):
install.packages("deSolve")
install.packages("ggplot2")

#to load after installing:
library(deSolve) #differential equation solver
library(ggplot2) #plotting package

#2. DEFINE DYNAMIC COMPARTMENTAL MODEL FUNCTIONS

#SIR model without demography
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

#SIR model with demography and optional waning immunity
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

#SIR model with demography and infection-induced mortality
OpenSIR_mort <-function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
        N = S + I + R
        
        #SIR w/ demography equations and infection-induced mortality from lecture
        dS <- -beta*S*I/N + birth - death*S + omega*R #note: births are no longer dependent on pop. size N 
        dI <- beta*S*I/N - death*I - gamma*I - (rho/(1-rho))*(gamma+death)*I
        dR <- gamma*I - death*R - omega*R
        
        # return the rates of change as a list
        list(c(dS, dI, dR))
    })
}

#SEIR model with demography
OpenSEIR<-function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
        N = S + E + I + R
        
        #SIR w/ demography equations from lecture
        dS <- -beta*S*I/N + birth*N - death*S + omega*R
        dE <- beta*S*I/N - sigma*E - death*E
        dI <- sigma*E - death*I - gamma*I
        dR <- gamma*I - death*R - omega*R
        
        # return the rates of change as a list
        list(c(dS, dE, dI, dR))
    })
}

#SICR model (C=carriers) with demography
OpenSICR<-function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
        N = S + I + C + R
        
        #SIR w/ demography equations from lecture
        dS <- -beta*S*I/N - epsilon*beta*S*C/N + birth*N - death*S + omega*R
        dI <- beta*S*I/N + epsilon*beta*S*C/N - death*I - gamma*I 
        dC <- gamma*q*I - gamma_c*C - death*C
        dR <- gamma*(1-q)*I +gamma_c*C - death*R - omega*R
        
        # return the rates of change as a list
        list(c(dS, dI, dC, dR))
    })
}

#3. DEFINE PLOTTING FUNCTIONS

#plot all compartments over time

#SIR model only
show_SIR_model_results<-function(output) {
    df1 <- data.frame(output)
    ggplot() + 
        geom_line(data = df1, aes(x = time, y = S), color = "blue") +
        geom_line(data = df1, aes(x = time, y = I), color = "red") +
        geom_line(data = df1, aes(x = time, y = R), color = "green") +
        xlab('Time') +
        ylab('Count') 
}

#SEIR model only
show_SEIR_model_results<-function(output) {
    df1 <- data.frame(output)
    ggplot() + 
        geom_line(data = df1, aes(x = time, y = S), color = "blue") +
        geom_line(data = df1, aes(x = time, y = I), color = "red") +
        geom_line(data = df1, aes(x = time, y = R), color = "green") +
        geom_line(data = df1, aes(x = time, y = E), color = "purple")
        xlab('Time') +
        ylab('Count') 
}

#All types of models with nice legend (but a bit more technical)
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

#plot the I compartment only from the SIR model (to zoom in on oscillatory dynamics, for example)
show_SIR_I_output <- function(output) {
    df1 <- data.frame(output)
    ggplot() + 
        geom_line(data = df1, aes(x = time, y = I), color="red") +
        xlab('Time') +
        ylab('Count') 
}

#plot a single compartment (of your choice) over time from any type of model
plot_1_compartment<-function(output, compartment) {
    #convert model output to a data frame (type of object in R) and reshape it so that each row represents a time-compartment combination (e.g. how many ppl in S compartment at time 10) - the latter just makes it easier to plot your results
    df1 <- data.frame(output)
    df2 <- reshape(df1, varying=colnames(df1)[2:ncol(df1)], 
                   v.name="sizes",  timevar="Compartment",
                   times=colnames(df1)[2:ncol(df1)], idvar="time", direction="long")
    ggplot() + 
        geom_line(data = subset(df2, Compartment==compartment), 
                  aes(x = time, y = sizes, color=Compartment)) + #plot the number of ppl in each compartment over time
        xlab('Time') + ylab('Count') + #label x and y-axes
        theme_bw() #this isn't necessary but makes the graph look a bit nicer
} 

#phase diagram (S prevalence vs. I prevalence over time)
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

#4. BASIC DEMONSTRATION (run the above code first) using the basic SIR model without demography

#define inputs to the ODE solver (parameters, initial compartment sizes, time steps)
parameters <- c(beta = 0.5, #effective contact rate
                gamma = 0.3 #recovery rate (1/duration infection)
)

state <- c(S = 99999, #population of 100,000, 1 person starts of infected
           I = 1, 
           R = 0
)

T_end <- 500 #run model for 500 time steps
times <- seq(0, T_end, by = 1) #runs the model for 500 time steps, and computes output at each time step 

#run ODE solver - you must use this syntax within the ode() function:
output <- ode(y = state, times = times, func = BasicSIR, parms = parameters)

#plot output
show_SIR_model_results(output)
show_model_results(output)
show_SIR_I_output(output)
plot_1_compartment(output, "I")
phase_diagram(output)
