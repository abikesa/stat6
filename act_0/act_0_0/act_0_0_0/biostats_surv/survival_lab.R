rm(list = ls()) # good practice to start clean
# require(survival)
# install.packages("survival")
library(survival)

#### load data ####
dat <- read.csv("recid.csv")
colnames(dat)
head(dat)
View(dat)
summary(dat)
attach(dat)
mysurv <- Surv(week, arrest)
mysurv

##### K-M by group ####
# overall K-M
myKM1 <- survfit(mysurv ~ 1)
myKM1 # no median survival?
myKM1sum <- summary(myKM1)
myKM1sum
plot(myKM1, ylim = c(0.6, 1)) # basic plot
plot(myKM1, 
     main=("Kaplan-Meier estimate with 95% confidence bounds"), 
     xlab = "week", 
     ylab = "Survival function",
     ylim = c(0.6, 1)) # plot with axes and title
plot(myKM1sum$time, myKM1sum$surv, type = "s", ylim = c(0.6, 1)) 
# all information is in the "summary" object

# K-M by group
myKM2 <- survfit(mysurv ~ race)
myKM2sum <- summary(myKM2) 
myKM2 # not printed, but information is there
plot(myKM2, ylim = c(0.6, 1)) # basic plot
plot(myKM2, 
     conf.int = T, 
     col = c("blue", "red"),
     main=("Kaplan-Meier estimate with 95% confidence bounds by race"), 
     xlab = "week", 
     ylab = "Survival function", 
     ylim = c(0.6, 1))
# legend("bottomleft", c("Race = 1: black","Race = 0: others"), 
#        col = c("red", "blue"), lty = 1)

# log-rank test
survdiff(mysurv ~ race, rho = 0)

## advanced materials on log-rank test, check "cross-over"!
levels(myKM2sum$strata)
myKM2race0 <- data.frame(time = myKM2sum$time[myKM2sum$strata == "race=0"],
                         surv = myKM2sum$surv[myKM2sum$strata == "race=0"])
myKM2race1 <- data.frame(time = myKM2sum$time[myKM2sum$strata == "race=1"],
                         surv = myKM2sum$surv[myKM2sum$strata == "race=1"])
hz0 <- diff(-log(myKM2race0$surv)) / diff(myKM2race0$time)
hz1 <- diff(-log(myKM2race1$surv)) / diff(myKM2race1$time)
t0 <- (myKM2race0$time[-1] + myKM2race0$time[-nrow(myKM2race0)])/2
t1 <- (myKM2race1$time[-1] + myKM2race1$time[-nrow(myKM2race1)])/2
myhz <- data.frame(time = c(t0, t1),
                   hazard = c(hz0, hz1),
                   group = c(rep(0, length(t0)),
                             rep(1, length(t1))))
library(ggplot2)
ggplot(data = myhz, aes(x = time, y = hazard, color = factor(group))) + 
      geom_step() +
      geom_smooth(formula = "y~x", method = "loess")
# smooth0 <- loess.smooth(t0, hz0)
# smooth1 <- loess.smooth(t1, hz1)
# plot(smooth0, type = "s")
# lines(smooth1, type = "s")

#### Cox PHM, time-independent covariates ####
# beta coefficient
dat$age <- (dat$age - mean(dat$age))/sd(dat$age)
dat$educ <- dat$educ - 2
myfit1 <- coxph(mysurv ~ fin + age + race + wexp + mar + paro + prio,
         ties = "breslow")
summary(myfit1)
confint(myfit1)
# baseline (cumulative) hazard
bch <- basehaz(myfit1, centered = F)
plot(bch$time, bch$hazard, type = "l",
     xlab = "week", ylab = "baseline cumulative hazard")
plot(bch$time, exp(-bch$hazard), type = "l",
     xlab = "week", ylab = "baseline survival function")



#### Cox PHM, time-dependent covariates ####
library(reshape)
dat$id <- 1:nrow(dat)
newdat <- melt(dat, id.vars = colnames(dat)[c(1:10, ncol(dat))])
newdat$end <- sapply(newdat$variable, function(x) {
      return(as.numeric(strsplit(as.character(x), "emp")[[1]][2]))
})
newdat$start <- newdat$end - 1
newdat$emp <- newdat$value
head(newdat)
newdat <- newdat[newdat$end <= newdat$week, ]
newdat$arrest[newdat$end < newdat$week] <- 0
mytvcphm <- coxph(Surv(start, end, arrest) ~ cluster(id) + fin + age + race +
            wexp + mar + paro + prio + emp, data = newdat, 
      ties = "breslow")
### using the "employment status for previous week" as the covariate
newdat$empprev = c(NA, newdat$emp[-nrow(newdat)])
for (i in unique(newdat$id)) {
      newdat$empprev[newdat$id == i][1] <- NA
}
coxph(Surv(start, end, arrest) ~ cluster(id) + fin + age + race +
            wexp + mar + paro + prio + empprev, data = newdat, ties = "breslow")
