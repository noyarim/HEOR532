#######################################################################
## This code explains how to run an individual-level simulation model##
## This code is a  modified version of DARTH's microsimulation code 
## You can find the original code here:https://github.com/DARTH-git/Microsimulation-tutorial
## Reference: Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
##            Microsimulation modeling for health decision sciences using R: A tutorial. Med Decis Making. 2018;38(3):400–22.

## Writer: Kyueun Lee 
## Date created: Mar 28, 2023
## Date modified: Apr 11, 2023

## In the base model, following is assumed:
## (1) Individual is characterized by age, sex, and risk level
## (2) Individuals in 'high risk' are at a higher risk of developing illness
## (3) Everyone has the same mortality (regardless of sex and age) unless they develop illness
## (4) Sick individuals have higher mortality than the others
## (5) In 'Treatment' scenario, everyone with disease are on the treatment
## (6) Treatment improves utility during sickness but it does NOT recover illness

## ASSIGNMNET ##
## For assignment #1, generate cost and QALY outcomes by risk groups and discuss who would benefit from treatment more or less.

## For assignment #2, revise the assumption of equal mortality and model age- and sex-specific mortality using the data ".xlsx"
## Discuss how the results (Cost, QALY of two scenarios) changes and speculate why


#######################################################################





## Install required packages (Need to run only once)
packages <- c("reshape2","dplyr","ggplot2")
for (package in packages){
  if(! package %in% installed.packages()){
    install.packages(package)
  }
}

## Load required packages
library(reshape2)
library(dplyr)
library(ggplot2)

## 0. Model parameters

seed <- 1234 # seed for generating population
n_pop <- 100 # pop size
t_end <- 10 # total simulation time
cSex <- c("female","male") 
dSex <- c(female=0.5, male=0.5) # distribution of sex
cRisk <- c("low","high")
dRisk <- c(low=0.8, high=0.2) # initial distribution of low/high risk group
v_state <- c("H","S","D")
dHealth <- c(healthy=1, sick=0, dead=0) # initial distribution of health states
n_state <- length(v_state)
mu_age <- 50 # mean age
sd_age <- 4 # sd of age

# Transition probabilities
pHS_low <- 0.1 # Probability of getting sick if low risk 
pHS_high <- 0.6 # Probability of getting sick if high risk
pHS <- list(low=pHS_low, high=pHS_high)
pDie_H <- 0.01 # Mortality if healthy (Baseline mortality)
pDie_S <- 0.05 # Mortality if sick (Baseline + disease specific mortality)
pDie <- list(H=pDie_H, S = pDie_S)

# Use this dataset for assignment (Add age-specific mortality)
#dt_mort <- read.csv("data/us_mortality_2019.csv")

# Utilities
u.H <- 1 # Utility if healthy
u.S <- 0.4 # Utility if Sick
u.Trt <- 0.7 # Utility if the sick individuals are on treatment 

# Costs
c.H <- 100
c.S <- 5000
c.Trt <- 500
c.Scrn <- 300

# Discount factors
d.c <- 0.03 # Discount rate for cost
d.e <- 0.03 # Discount rate for utility


## 1. Generate population (i.e. generate initial distribution of attributes to follow up)
# Age
set.seed(seed)
v_age <- rnorm(n_pop, mean=mu_age, sd=sd_age) # age vector
m_age <- matrix(rep(v_age, t_end+1),nrow = n_pop)

# Sex
rn_sex <- runif(n_pop)
v_sex <- as.numeric(rn_sex < dSex['female']) # 1=female, 0=male 
v_sex <- sample(cSex, size = n_pop, prob = dSex, replace = TRUE)
table(v_sex)
m_sex <- matrix(rep(v_sex, t_end+1),nrow = n_pop)

# Risk of disease
v_risk <- sample(cRisk, size=n_pop, prob=dRisk, replace=TRUE)
table(v_risk)
m_risk <- matrix(rep(v_risk, t_end+1),nrow = n_pop)

# Health states
v_health <- rep("H",n_pop)
m_health <- matrix(rep(v_health, t_end+1),nrow = n_pop)

# Utility
m_E <- matrix(0, ncol=t_end+1, nrow=n_pop)

# Cost
m_C <- matrix(0, ncol=t_end+1, nrow=n_pop)

# A vector of discount factors
v.dwc <- 1 / (1 + d.c) ^ (0:t_end)   # calculate the cost discount weight based on the discount rate d.c 
v.dwe <- 1 / (1 + d.e) ^ (0:t_end)   # calculate the QALY discount weight based on the discount rate d.e


## 3. Functions 
# A function to run a microsimulation model
MicroSim <- function(m_health, n_pop, t_end, v_state, d.c, d.e, TR.out = TRUE, TS.out = TRUE, Trt = FALSE, seed = 1) {
  
  for(i in (1:n_pop)){
    set.seed(seed + i)                  # set the seed for every individual for the random number generator
    m_C[i, 1] <- Costs(m_health[i, 1], Trt)  # estimate costs per individual for the initial health state conditional on treatment
    m_E[i, 1] <- Effs (m_health[i, 1], Trt)  # estimate QALYs per individual for the initial health state conditional on treatment
  
    # Transition probability for individual i
    for (t in (1:t_end)){
      # Calculate transition probability for individual i given time t
      pTrans <- cal_p(parameters,health=m_health[i,t],sex=m_sex[i,t],risk=m_risk[i,t],age=m_age[i,t])
      
      # update health state
      m_health[i,t + 1] <- sample(v_state,size=1,prob=pTrans,replace=TRUE)
      
      # update age
      m_age[i,t + 1] <- m_age[i,t] + 1
      
      # Calculate cost and QALYs for individuals i given time t
      m_C[i, t + 1] <- Costs(m_health[i, t + 1], Trt)   # estimate costs per individual during cycle t + 1 conditional on treatment
      m_E[i, t + 1] <- Effs( m_health[i, t + 1], Trt)   # estimate QALYs per individual during cycle t + 1 conditional on treatment
    }
  }
  tc <- m_C %*% v.dwc       # total (discounted) cost per individual
  te <- m_E %*% v.dwe       # total (discounted) QALYs per individual 
  
  tc_hat <- mean(tc)        # average (discounted) cost 
  te_hat <- mean(te)        # average (discounted) QALYs
  
  if (TS.out == TRUE) {  # create a  matrix of transitions across states
    TS <- paste(m_health, cbind(m_health[, -1], NA), sep = "->") # transitions from one state to the other
    TS <- matrix(TS, nrow = n_pop)
    rownames(TS) <- paste("Ind",   1:n_pop, sep = " ")   # name the rows 
    colnames(TS) <- paste("Cycle", 0:t_end, sep = " ")   # name the columns 
  } else {
    TS <- NULL
  }
  
  if (TR.out == TRUE) { # create a trace from the individual trajectories
    TR <- t(apply(m_health, 2, function(x) table(factor(x, levels = v_state, ordered = TRUE))))
    TR <- TR / n_pop                                       # create a distribution trace
    rownames(TR) <- paste("Cycle", 0:t_end, sep = " ")     # name the rows 
    colnames(TR) <- v_state                                  # name the columns 
  } else {
    TR <- NULL
  }
  
  results <- list(m.M = m_health, m.C = m_C, m.E = m_E, tc = tc, te = te, tc_hat = tc_hat, te_hat = te_hat, TS = TS, TR = TR) # store the results from the simulation in a list  
  return(results)  # return the results

}

# A function to calculate transition probabilities for individual i given time t
cal_p <- function(parameters, health, sex, risk, age){
  # Transition probability from Healthy to Sick given risk group
  this_pHS <- pHS[[risk]]
  # Transition probability of dying given health state (Healthy or Sick)
  this_pDie <- pDie[[health]]
  
  pTrans <- rep(NA,n_state)
  pTrans[health == 'H'] <- c(1-this_pHS-this_pDie, this_pHS, this_pDie)
  pTrans[health == 'S'] <- c(0, 1-this_pDie, this_pDie)
  pTrans[health == 'D'] <- c(0,0,1)
  
  # Check if transition probabilities sum up to 1 
  ifelse(sum(pTrans)!=1,print("Transition probabilities do not sum up to 1"),return(pTrans))
}

# A function to calculate cost for individual i given time t
Costs <- function (health, Trt = FALSE) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  is the individual being treated? (default is FALSE) 
  
  c.it <- 0                                  # by default the cost for everyone is zero 
  c.it[health == "H"]  <- c.H                  # update the cost if healthy
  c.it[health == "S"] <- c.S + c.Trt * Trt   # update the cost if sick conditional on treatment
  return(c.it)        		                   # return the costs
}

# A function to calculate utility for individual i given time t
Effs <- function (health, Trt = FALSE) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Trt:  is the individual treated? (default is FALSE) 
  # cl:   cycle length (default is 1)
  
  u.it <- 0                      # by default the utility for everyone is zero
  u.it[health == "H"]  <- u.H      # update the utility if healthy
  u.it[health == "S"] <- Trt * u.Trt + (1 - Trt) * u.S  # update the utility if sick conditional on treatment
  QALYs <-  u.it             # calculate the QALYs during cycle t
  return(QALYs)              # return the QALYs
}


## 4. Run microsimulation with and without treatment 
# No intervention
sim_no_trt  <- MicroSim(m_health, n_pop, t_end, v_state, d.c, d.e, Trt = FALSE) # run for no treatment
# Treatment
sim_trt     <- MicroSim(m_health, n_pop, t_end, v_state, d.c, d.e, Trt = TRUE)  # run for treatment


##5. Generate graphs of intermediate outcomes (e.g. without treatment)
# A. Distribution of the time spent in 'Healthy' before decoming 'Sick' (i.e. when H->S happened)
TS_notrt <- sim_no_trt$TS
TS_HS <- apply(TS_notrt, 1, function(x) which(x=='H->S')[1])
TS_HD <- apply(TS_notrt, 1, function(x) which(x=='H->D')[1])

dt_TS <- data.frame(t_HS = TS_HS, t_HD = TS_HD)
dt_TS <- dt_TS %>%
  mutate(
    t_HS = ifelse(is.na(TS_HS) & is.na(TS_HD), t_end, t_HS),
    risk = v_risk,
    sex = v_sex
  )

ggplot(dt_TS, aes(x=t_HS,fill=risk)) +
  geom_histogram(binwidth = 0.6)+
  facet_wrap(~risk)+
  theme_bw()


# B. Prevalence of illness
# no treatment
TR_notrt <- sim_no_trt$TR # state distribution by cycle
n_alive_notrt <- TR_notrt[,"H"] + TR_notrt[,"S"] # number of alive individuals by cycle
prev_notrt <- TR_notrt[,"S"]/n_alive_notrt # Prevalence of illness by cycle
# treatment
TR_trt <- sim_trt$TR # state distribution by cycle
n_alive_trt <- TR_trt[,"H"] + TR_trt[,"S"] # number of alive individuals by cycle
prev_trt <- TR_trt[,"S"]/n_alive_trt # Prevalence of illness by cycle

dt_prev <- data.frame(t = seq(0,t_end), No_trt = prev_notrt, Trt= prev_trt)
# overall prevalence
ggplot(data=dt_prev)+
  geom_line(aes(x=t,y=No_trt), color='black')+
  geom_line(aes(x=t,y=Trt), color='red')+
  theme_bw()


## 5. Calculate QALYs and Costs
# Screening (complete this in assignment)
# store the mean costs (and the MCSE) of each strategy in a new variable v.C (vector costs)
v.C  <- c(sim_no_trt$tc_hat, sim_trt$tc_hat) 
se.C <- c(sd(sim_no_trt$tc), sd(sim_trt$tc)) / sqrt(n_pop)
# store the mean QALYs (and the MCSE) of each strategy in a new variable v.E (vector health outcomes)
v.E  <- c(sim_no_trt$te_hat, sim_trt$te_hat)
se.E <- c(sd(sim_no_trt$te), sd(sim_trt$te)) / sqrt(n_pop)

delta.C <- v.C[2] - v.C[1]                   # calculate incremental costs
delta.E <- v.E[2] - v.E[1]                   # calculate incremental QALYs
se.delta.E <- sd(sim_trt$te - sim_no_trt$te) / sqrt(n_pop) # Monte Carlo squared error (MCSE) of incremental costs
se.delta.C <- sd(sim_trt$tc - sim_no_trt$tc) / sqrt(n_pop) # Monte Carlo squared error (MCSE) of incremental QALYs
ICER    <- delta.C / delta.E                 # calculate the ICER
results <- c(delta.C, delta.E, ICER)         # store the values in a new variable

# Create full incremental cost-effectiveness analysis table
table_micro <- data.frame(
  c(round(v.C, 0),  ""),           # costs per arm
  c(round(se.C, 0), ""),           # MCSE for costs
  c(round(v.E, 3),  ""),           # health outcomes per arm
  c(round(se.E, 3), ""),           # MCSE for health outcomes
  c("", round(delta.C, 0),   ""),  # incremental costs
  c("", round(se.delta.C, 0),""),  # MCSE for incremental costs
  c("", round(delta.E, 3),   ""),  # incremental QALYs 
  c("", round(se.delta.E, 3),""),  # MCSE for health outcomes (QALYs) gained
  c("", round(ICER, 0),      "")   # ICER
)
rownames(table_micro) <- c(c("No Treatment","Treatment"), "* are MCSE values")  # name the rows
colnames(table_micro) <- c("Costs", "*",  "QALYs", "*", "Incremental Costs", "*", "QALYs Gained", "*", "ICER") # name the columns
table_micro  # print the table 

# Create the distribution of outcomes
dt_tc <- data.frame(No_treatment = sim_no_trt$tc, Treatment = sim_trt$tc)
ggplot(data=dt_tc)+
  geom_histogram(aes(No_treatment), fill = 'dark grey', alpha = 0.4, bins=20)+
  geom_histogram(aes(Treatment), fill = 'tomato', alpha = 0.4, bins=20)+
  geom_vline(aes(xintercept=sim_no_trt$tc_hat),color='dark grey')+
  geom_vline(aes(xintercept=sim_trt$tc_hat), color='tomato')+
  xlab("Total cost")+
  theme_bw()

# If we want to see the results by risk group..










