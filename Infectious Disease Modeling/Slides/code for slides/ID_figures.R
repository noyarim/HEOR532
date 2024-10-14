library(tidyverse)
library(scales)
library(deSolve)

gammas <- c(0.5, 1, 2)
times <- (0:100)/10

exp <- 1- pexp(times, 0.5) #1 minus CDF
exp2 <- exp(-times*0.5) #1 minus CDF

exp_all <- data.frame()
for(g in gammas) {
  exp <- 1-pexp(times, g)
  exp <- data.frame("gamma"=g, "times"=times, "inf"=exp)
  exp_all <- rbind(exp_all, exp)
}

exp_all$gamma <- as.factor(exp_all$gamma)
colors <- hue_pal()(3)
names(colors) <- levels(exp_all$gamma)

ggplot(exp_all, aes(x=times, y=inf, color=gamma)) +
  geom_line() + 
  scale_x_continuous(expand=expansion(mult=c(0,0)),
                     breaks=c(2,4,6,8)) + 
  scale_y_continuous(expand=expansion(mult=c(0, NA))) + 
  labs(x="Time", y="% infected", color="gamma") +
  theme_bw() + theme(panel.grid=element_blank())
ggsave("exp_dist.jpg", dpi=500, height=4, width=7)

ggplot(exp_all[exp_all$gamma==0.5,], 
       aes(x=times, y=inf, color=gamma, fill=gamma)) +
  geom_line() + 
  geom_area(alpha=0.5) +
  scale_x_continuous(expand=expansion(mult=c(0,0)),
                     breaks=c(2,4,6,8)) + 
  scale_y_continuous(expand=expansion(mult=c(0, NA))) + 
  scale_fill_manual(values=colors) +
  scale_color_manual(values=colors) + 
  labs(x="Time", y="% infected", color="gamma", fill="gamma") +
  theme_bw() + theme(panel.grid=element_blank())
ggsave("exp_dist2.jpg", dpi=500, height=4, width=7)

#births/deaths
life_exp <- c(60, 70, 80)
gammas <- 1/life_exp
times <- 0:110
life_exp_all <- data.frame()
for(i in 1:3) {
  exp <- 1-pexp(times, gammas[i])
  exp <- data.frame("gamma"=gammas[i], "life_exp"=life_exp[i], 
                    "times"=times, "alive"=exp)
  life_exp_all <- rbind(life_exp_all, exp)
}

ggplot(life_exp_all, aes(x=times, y=alive, color=as.factor(life_exp))) +
  geom_line() + 
  scale_x_continuous(expand=expansion(mult=c(0,0))) + 
  labs(x="Time", y="% alive", color="Life expectancy") +
  theme_bw() + theme(panel.grid=element_blank())
ggsave("life_exp_dist.jpg", dpi=500, height=4, width=7)

#model of births and deaths
births_deaths <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{ 
    
    N = state 
    
    #births and deaths
    dN <- -death*N + birth*N
    
    return(list(dN))
  })
}

state <- 1
T_end <- 1000 
times <- seq(0, T_end, by = 1)

output_grow <- ode(y = state, times = times, 
                   func = births_deaths, parms = c("death"=1/70, "birth"=1/65))
output_steady <- ode(y = state, times = times, 
                     func = births_deaths, parms = c("death"=1/70, "birth"=1/70))
output_shrink <- ode(y = state, times = times, 
                    func = births_deaths, parms = c("death"=1/70, "birth"=1/80))

out_all <- rbind(as.data.frame(output_grow) %>% mutate("scenario"="b > mu"),
                 as.data.frame(output_steady) %>% mutate("scenario"="b = mu"),
                 as.data.frame(output_shrink) %>% mutate("scenario"="b < mu"))
names(out_all) <- c("time", "pop", "scenario")

ggplot(out_all, aes(x=time, y=pop*1000, color=scenario)) + 
  geom_line() + 
  scale_x_continuous(expand=expansion(mult=c(0,0))) + 
  labs(x="Years", y="Population Size", color="") +
  theme_bw() + theme(panel.grid=element_blank())
ggsave("births_deaths.jpg", dpi=500, height=4, width=7)
