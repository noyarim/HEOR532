#######################################################################
## This code explains how to run an individual-level simulation model##
## Writer: Kyueun Lee 
## Date created: Mar 28, 2023
## Date modified: Mar 28, 2023
#######################################################################
## Install required packages (Need to run only once)


## Load required packages


## 0. setup
seed <- 1234 # seed for random number generation
n_pop <- 100 # pop size
t_end <- 10 # total simulation time
dSex <- c(female=0.5, male=0.5) # distribution of sex
dRisk <- c(low=0.8, high=0.2) # initial distribution of low/high risk group
dHealth <- c(healthy=1, sick=0, dead=0) # initial distribution of health states
mu_age <- 50 # mean age
sd_age <- 4 # sd of age
H <- 1 
S <- 2
D <- 3
## 1. Generate population (i.e. generate initial distribution of attributes to follow up)
# Age
v_age <- rnorm(n_pop, mean=mu_age, sd=sd_age) # age vector
m_age <- matrix(rep(v_age, t_end),nrow = t_end, byrow = TRUE)
# Sex
set.seed(seed)
rn_sex <- runif(n_pop)
v_sex <- as.numeric(rn_sex < dSex['female']) # 0=male, 1=female 
mean(v_sex)
m_sex <- matrix(rep(v_sex, t_end),nrow = t_end, byrow = TRUE)

# Risk of disease
rn_risk <- runif(n_pop)
v_risk <-ifelse(rn_risk<dRisk["low"],0,1) # 0=low, 1=high
mean(v_risk)
m_risk <- matrix(rep(v_risk, t_end),nrow = t_end, byrow = TRUE)

# Health state
v_health <- rep(1,n_pop)
m_health <- matrix(rep(v_health, t_end),nrow = t_end, byrow = TRUE)

# Compile initial attributes in a list 
l_pop <- list( age = m_age,
               sex = m_sex,
               risk = m_risk,
               health = m_health)

## 2. Transition probabilities
pHS_low <- 0.1
pHS_high <- 0.6
l_pHS <- list(L=pHS_low, H=pHS_high)

dt_mort <- read.csv("data/us_mortality_2019.csv")
l_mort <- list(F=dt_mort$f_mort, M=dt_mort$m_mort)

## 3. Simulate health transitions over time
## Functions
# A function to update health state
update_health <- function(this_v_health){
  # Individual index
  i_H <- which(this_v_health == H & this_v_risk == 0)
  i_S <- which(this_v_health == S & this_v_risk == 0)
  i_D <- which(this_v_health == D & this_v_risk == 0)

  # Healthy to Sick
  v_pHS[i_H] <- l_pHS[this_v_risk[i_H]]
  l_pHD <- l_mort[v_sex[i_H]]
  for (i in 1:length(l_pHD)){
    temp[i] <-l_pHD[i][age]
  }
  v_pHD[i_H] <- v_
  v_pHD <- c()
  v_pSD <- c()
  v_pSH <- 0
  rn_pHS <- runif(n_pop)
  
}



# Simulate over time
for (t in range(1,t_end)){
  
  # update health state
  m_health[t+1,] <- update_health(m_health[t,])
  # update risk group
  
  # update age
  

}



##4. Generate intermediate outcomes




## 5. Calculate QALYs and Costs



## 6. Model screening and treatment



## 7. 
