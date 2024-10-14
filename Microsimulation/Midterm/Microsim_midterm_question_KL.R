library(tidyverse)

setwd("/Users/kyueunlee/Library/CloudStorage/GoogleDrive-leex5499@umn.edu/My Drive/Z Drive/Teaching/HEOR532/Microsimulation")
set.seed(1234)
inc_c1 <- rnorm(20,2500,500)
inc_c2 <- rnorm(20,2500,50)

set.seed(1234)
inc_e1 <- rnorm(20,0.7,0.1)
inc_e2 <- rnorm(20,0.7,0.01)

set.seed(1234)
no_c <- rnorm(20,5000,10)
no_e <- rnorm(20,3.2,0.1)

mirosim_result <- data.frame(no_c = no_c,
           no_e = no_e,
           ta_c = no_c + inc_c1,
           ta_e = no_e + inc_e1,
           tb_c = no_c + inc_c2,
           tb_e = no_e + inc_e2,
           inc_c1 = inc_c1,
           inc_e1 = inc_e1,
           inc_c2 = inc_c2,
           inc_e2 = inc_e2)
write.csv(mirosim_result, "Midterm_microsim_result.csv")
microsim_summary <- mirosim_result %>%
  summarize(mean_inc_c1 = mean(inc_c1),
            sd_inc_c1 = sd(inc_c1),
            mean_inc_e1 = mean(inc_e1),
            sd_inc_e1 = sd(inc_e1),
            mean_inc_c2 = mean(inc_c2),
            sd_inc_c2 = sd(inc_c2),
            mean_inc_e2 = mean(inc_e2),
            sd_inc_e2 = sd(inc_e2),
            icer_ta = mean_inc_c1/mean_inc_e1,
            icer_tb = mean_inc_c2/mean_inc_e2)
write.csv(microsim_summary, "Midterm_microsim_summary.csv")
