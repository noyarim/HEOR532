###### Section 0: Libraries to load; If you do not have them install from Tools menu
library(deSolve)
library(ggplot2)
library(plsgenomics)

#### Helper Functions
show_SIR_model_results<-function(output) {
## blue (not shown) is susceptibles
## red is infectious
## green is recovered
## purple is vaccinated
  df1 <- data.frame(output)
  ggplot() + 
    #    geom_line(data = df1, aes(x = time, y = X), color = "blue") +
    geom_line(data = df1, aes(x = time, y = Y), color = "red") +
    geom_line(data = df1, aes(x = time, y = Z), color = "green") +
    geom_line(data = df1, aes(x = time, y = V), color = "purple") +
    xlab('Time') +
    ylab('Count') 
}

###### Section 1: Examine the function for the Open SIR Model With Imperfect/Waning Vaccination
OpenSIRVax<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    N = X + Y + Z + V
    dX <- (1-(p_vax*p_imm_vax))*birth*N + waning*V - beta*X*Y/N - death*X
    dY <- beta*X*Y/N - death*Y - recover*Y
    dZ <- recover*Y - death*Z
    dV <- ((p_vax*p_imm_vax))*birth*N - waning*V - death*V
    
    # return the rate of change
    list(c(dX, dY, dZ, dV))
  }) # end with(as.list ...
}

###### Section 2: Initializing parameters, initial state, and First model run w/o vaccination

birth = 0.05
death = birth
parameters_novax <- c(birth     = birth, # birth rate
                      death     = death, # death rate = birth rate so stable population
                      beta      = 0.4,   # trasmission rate
                      recover   = 0.2,   # recovery rate
                      p_vax     = 0.0,   # proportion of infants vaccinated
                      p_imm_vax = 0.0,   # probability of immunity conditional on vaccination
                      waning    = 0.0    # waning rate for vaccine induced immunity
)

state <- c(X = 99999,
           Y = 1,
           Z = 0,
           V = 0
)

times <- seq(0, 500, by = 1)
output_novax <- ode(y = state, times = times, func = OpenSIRVax, parms = parameters_novax)
show_SIR_model_results(output_novax)

###### Section 3: Two more helper functions

compute_cumulative_infection_time<-function(output) {
  df1 <- data.frame(output)
  return(sum(df1$Y))
}

compute_cumulative_infection_time(output_novax)

###### Section 4: Running the model with perfect vaccination but incomplete coverage

parameters_vaxperf <- c(birth     = birth, # birth rate
                        death     = death, # death rate = birth rate so stable population
                        beta      = 0.4,   # trasmission rate
                        recover   = 0.2,   # recovery rate
                        p_vax     = 0.1,   # proportion of infants vaccinated
                        p_imm_vax = 1.0,   # probability of immunity conditional on vaccination
                        waning    = 0.0    # waning rate for vaccine induced immunity
)

output_vaxperf <- ode(y = state, times = times, func = OpenSIRVax, parms = parameters_vaxperf)
show_SIR_model_results(output_vaxperf)

compute_cumulative_infection_time(output_vaxperf)

###### Section 5: Simple threshold analysis

call_ode<-function(params) {
  output <- ode(y = state, time = times, func = OpenSIRVax, parms = params)
  return(compute_cumulative_infection_time(output))
}

params_to_try <- parameters_novax
for(i in 1:90) {
  
  new_params    <- parameters_vaxperf
  new_params[5] <- i/100
    
  params_to_try <- rbind(params_to_try, new_params)
  rownames(params_to_try)[i+1] <- paste0("vax",i)
}

head(params_to_try)
tail(params_to_try)

cum_infections <- t(apply(params_to_try, 1, call_ode))
plot(x = seq(0:90)/100, y=cum_infections)

###### Section 6: Imperfect vaccination
  
parameters_vaximperf <- c(birth     = birth, # birth rate
                        death     = death, # death rate = birth rate so stable population
                        beta      = 0.4,   # trasmission rate
                        recover   = 0.2,   # recovery rate
                        p_vax     = 0.1,   # proportion of infants vaccinated
                        p_imm_vax = 0.4,   # probability of immunity conditional on vaccination
                        waning    = 0.0    # waning rate for vaccine induced immunity
)

output_vaximperf <- ode(y = state, times = times, func = OpenSIRVax, parms = parameters_vaximperf)
show_SIR_model_results(output_vaximperf)

###### Section 7: Multi-way sensitivity analysis

call_ode(params = parameters_vaxperf)
call_ode(params = parameters_vaximperf)
call_ode(params = parameters_novax)

results_matrix = matrix(0, nrow=20, ncol=20)
for(i in 0:19) {
  for(j in 0:19) {
    curr_params <- parameters_novax
    curr_params[5] = i/40
    curr_params[6] = 1 - j/40
    results_matrix[i+1,j+1] = call_ode(curr_params)
  }
}

matrix.heatmap(results_matrix)

###### Section 8: Waning Immunity and Challenge opportunity

parameters_vaximperfwane <- c(birth     = birth, # birth rate
                          death     = death, # death rate = birth rate so stable population
                          beta      = 0.4,   # trasmission rate
                          recover   = 0.2,   # recovery rate
                          p_vax     = 0.1,   # proportion of infants vaccinated
                          p_imm_vax = 0.4,   # probability of immunity conditional on vaccination
                          waning    = 0.5    # waning rate for vaccine induced immunity
)

output_vaximperfwane <- ode(y = state, times = times, func = OpenSIRVax, parms = parameters_vaximperfwane)
show_SIR_model_results(output_vaximperf)

call_ode(params = parameters_vaxperf)
call_ode(params = parameters_vaximperf)
call_ode(params = parameters_vaximperfwane)
call_ode(params = parameters_novax)

#### If you want a CHALLEGE, try coding and running a 2-way sensitivity
#### analysis of waning vs. p_imm_vax. What do you notice

