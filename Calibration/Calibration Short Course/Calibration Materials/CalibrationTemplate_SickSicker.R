###########################################################################
# This code was created by the DARTH workgroup (www.darthworkgroup.com). 
# When using or modifying this code, please do so with attribution and 
# cite our publications:

# - Alarid-Escudero F, Maclehose RF, Peralta Y, Kuntz KM, Enns EA. 
#   Non-identifiability in model calibration and implications for 
#   medical decision making. Med Decis Making. 2018; 38(7):810-821.

# - Jalal H, Pechlivanoglou P, Krijkamp E, Alarid-Escudero F, Enns E, 
#   Hunink MG. An Overview of R in Health Decision Sciences. 
#   Med Decis Making. 2017; 37(3): 735-746. 
###########################################################################


###################  Calibration Specifications  ###################

# Model: Sick-Sicker 4-state Markov Model
# Inputs to be calibrated: p.S1S2, hr.S1, hr.S2
# Targets: Surv, Prev, PropSick

# Search method: Random search using Latin-Hypercube Sampling
# Goodness-of-fit measure: Sum of Log-Likelihood

####################################################################


####################################################################
######  Load packages and function files  ####
####################################################################
# calibration functionality
library(lhs)

# visualization
library(plotrix)
library(psych)
library(scatterplot3d) # now that we have three inputs to estimate, we'll need higher dimension visualization


####################################################################
######  Load target data  ######
####################################################################
load("SickSickerTargets.RData")

#Plot the targets#
# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = SickSicker.targets$Surv$Time, y = SickSicker.targets$Surv$value, 
                ui = SickSicker.targets$Surv$ub,
                li = SickSicker.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")

# TARGET 2: (Prev")
plotrix::plotCI(x = SickSicker.targets$Prev$Time, y = SickSicker.targets$Prev$value, 
                ui = SickSicker.targets$Prev$ub,
                li = SickSicker.targets$Prev$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Prevalence")

# TARGET 3: ("PropSick")
plotrix::plotCI(x = SickSicker.targets$PropSick$Time, y = SickSicker.targets$PropSick$value, 
                ui = SickSicker.targets$PropSick$ub,
                li = SickSicker.targets$PropSick$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Proportion Sick (S1)")


####################################################################
######  Load model as a function  ######
####################################################################
# - inputs are parameters to be estimated through calibration
# - outputs correspond to the target data

source("Markov_Sick-Sicker - Function.R") # creates the function markov_sick_sicker()
#this function outputs survival, prev, and prop_sick over time at given probabilities#

# Check that it works
v.params.test = c(p.S1S2=.5, hr.S1=.1, hr.S2=.05)
markov_sick_sicker(v.params.test)


####################################################################
######  Specify calibration parameters  ######
####################################################################

# Specify seed (for reproducible sequence of random numbers)
set.seed(072218)

# number of random samples
n.samp = 1000

# names and number of input parameters to be calibrated
param.names = c("p.S1S2","hr.S1", "hr.S2")
n.param = length(param.names)

# range on input search space (plausible range)
lb <- c(p.S1S2 = 0.01, hr.S1 = 1.0, hr.S2 = 5) # lower bound
ub <- c(p.S1S2 = 0.50, hr.S1 = 4.5, hr.S2 = 15) # upper bound

# number of calibration targets
target.names = c("Surv", "Prev", "Prop")
n.target = length(target.names)


####################################################################
######  Calibrate!  ######
####################################################################

m.lhs.unit <- randomLHS(n.samp, n.param)

# Rescale to min/max of each parameter - transforming 0-1 unit LHS to any multivariate distribution you want
# this uses the inverse CDF
# in this case we've chosen to use a uniform distribution (on 0.04 to 0.16) - could also use qbeta, etc.
m.param.samp = matrix(nrow=n.samp,ncol=n.param)
for (i in 1:n.param){
  m.param.samp[,i] = qunif(m.lhs.unit[,i],
                           min = lb[i],
                           max = ub[i])
}


colnames(m.param.samp) = param.names

# view resulting parameter set samples
pairs.panels(m.param.samp)

###  Run the model for each set of input values ###

# initialize goodness-of-fit vector
m.GOF = matrix(nrow = n.samp, ncol = n.target) #each row corresponds to 1 parameter sample, each col is fit to a target (3 targets)
colnames(m.GOF) = paste0(target.names, "_fit")

# loop through sampled sets of input values
for (j in 1:n.samp){
  
  ###  Run model for a given parameter set  ###
  model.res = markov_sick_sicker(v.params = m.param.samp[j, ])
  
  
  ###  Calculate goodness-of-fit of model outputs to targets  ###
  
  # TARGET 1: Survival ("Surv")
  # log likelihood  
  m.GOF[j,1] = sum(dnorm(x = SickSicker.targets$Surv$value,
                         mean = model.res$Surv,
                         sd = SickSicker.targets$Surv$se,
                         log = T)) #dnorm is our likelihood function - normal pdf in this case. sum over dnorm for each time point
  #maxing normal likelihood is same as minimizing unweighted sum of squared errors
  #but using MLE gives us some nice stuff (like Hessian for Cov matrix) if we want it
  
  # weighted sum of squared errors (alternative to log likelihood)
  # w = 1/(CRS.targets$Surv$se^2)
  # m.GOF[j,1] = -sum(w*(CRS.targets$Surv$value - v.res)^2)
  
  
  # TARGET 2: (Prev)
  # log likelihood
  m.GOF[j,2] = sum(dnorm(x = SickSicker.targets$Prev$value,
                          mean = model.res$Prev,
                          sd = SickSicker.targets$Prev$se,
                          log = T))
  
  # TARGET 3: (PropSick)
  # log likelihood
  m.GOF[j,3] = sum(dnorm(x = SickSicker.targets$PropSick$value,
                         mean = model.res$PropSick,
                         sd = SickSicker.targets$PropSick$se,
                         log = T))
  
  
} # End loop over sampled parameter sets

###  Combine fits to the different targets into single GOF (bc we have 3 targets)  ###
# can give different targets different weights - if you just add, you assume same weights
v.weights = matrix(1, nrow = n.target, ncol = 1)
# matrix multiplication to calculate weight sum of each GOF matrix row
v.GOF.overall = c(m.GOF%*%v.weights)
# Store in GOF matrix with column name "Overall"
m.GOF = cbind(m.GOF,Overall_fit=v.GOF.overall)

####################################################################
######  Exploring best-fitting input sets  ######
####################################################################

# Arrange parameter sets in order of fit
m.calib.res = cbind(m.param.samp,m.GOF)
m.calib.res = m.calib.res[order(-m.calib.res[,"Overall_fit"]),]

# Examine the top 10 best-fitting sets
m.calib.res[1:10,]

#Plot the top 100 (top 10%)
scatterplot3d::scatterplot3d(x=m.calib.res[1:100,1], y=m.calib.res[1:100,2], z=m.calib.res[1:100,3],
                             xlab=colnames(m.calib.res)[1], ylab=colnames(m.calib.res)[2], zlab=colnames(m.calib.res)[3])

### Plot model-predicted output at mean vs targets ###
v.out.best <- markov_sick_sicker(m.calib.res[1,])

# TARGET 1: Survival ("Surv")
plotrix::plotCI(x = SickSicker.targets$Surv$Time, y = SickSicker.targets$Surv$value, 
                ui = SickSicker.targets$Surv$ub,
                li = SickSicker.targets$Surv$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Pr Survive")
points(x = SickSicker.targets$Surv$Time, 
       y = v.out.best$Surv, 
       pch = 8, col = "red")
legend("topright", 
       legend = c("Target", "Model-predicted output at mean"),
       col = c("black", "red"), pch = c(1, 8))

# TARGET 2: Prevalence ("Prev")
plotrix::plotCI(x = SickSicker.targets$Prev$Time, y = SickSicker.targets$Prev$value, 
                ui = SickSicker.targets$Prev$ub,
                li = SickSicker.targets$Prev$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Prevalence")
points(x = SickSicker.targets$Prev$Time, 
       y = v.out.best$Prev, 
       pch = 8, col = "red")
legend("topright", 
       legend = c("Target", "Model-predicted output at mean"),
       col = c("black", "red"), pch = c(1, 8))

# TARGET 3: Proportion Sick ("PropSick")
plotrix::plotCI(x = SickSicker.targets$PropSick$Time, y = SickSicker.targets$PropSick$value, 
                ui = SickSicker.targets$PropSick$ub,
                li = SickSicker.targets$PropSick$lb,
                ylim = c(0, 1), 
                xlab = "Time", ylab = "Proportion Sick")
points(x = SickSicker.targets$PropSick$Time, 
       y = v.out.best$PropSick, 
       pch = 8, col = "red")
legend("topright", 
       legend = c("Target", "Model-predicted output at mean"),
       col = c("black", "red"), pch = c(1, 8))
