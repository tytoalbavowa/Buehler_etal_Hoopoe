#####
# R code for Bühler et al. 
#"Sexual segregation when foraging is relaxed under good environmental conditions in the European hoopoe Upupa epops"# 
# Urs G. Kormann 15.12.2020 
# urs.kormann@vogelwarte.ch

## questions to check: do we need quality == 2 only? 

#####  Testing whether sexes differed in foraging distance = loop max

setwd("~/Documents/Manuskripte/2018/Hoophoe/Urs_Analysen")
tmp3= read.csv("dfb3.csv")
library(ggplot2)
library(ggpubr)
library(nlme)

# How many animals
length(table(tmp3$id))

# How many Territorien
length(table(tmp3$territory))

# how many loops were done by 
#FEMALES
sum(table(tmp3$sex, tmp3$id)[1,])

#MALES
sum(table(tmp3$sex, tmp3$id)[2,])

##########
# 1A) Test whether sexes differ in their foraging distance 
# fit lme with the individual bird nested within territory as random intercept
sex_mean <-lme(log10(meandist.loop) ~ sex , random = ~1|territory/id, na.action = "na.omit",data = tmp3[tmp3$loop_nr != "0",])

# try out different variance structures
sex_mean_powr =update(sex_mean, weights=varPower())
sex_mean_ide =update(sex_mean, weights=varIdent(form=~1|sex))
sex_mean_exp =update(sex_mean, weights=varExp())

# compare using AICc
library(MuMIn)
AICc(sex_mean,sex_mean_ide,  sex_mean_powr, sex_mean_exp)

# varIde works best, deltaSICc app. 50

# check diagnostic plots
plot(sex_mean)
qqnorm(resid(sex_mean), main="qq-plot residuals")
qqline(resid(sex_mean))

plot(sex_mean_exp)
qqnorm(resid(sex_mean_exp), main="qq-plot residuals")
qqline(resid(sex_mean_exp))

plot(sex_mean_powr)
qqnorm(resid(sex_mean_powr), main="qq-plot residuals")
qqline(resid(sex_mean_powr))

plot(sex_mean_ide)
qqnorm(resid(sex_mean_ide), main="qq-plot residuals")
qqline(resid(sex_mean_ide))

# diagnostic plots not pretty, but good enough

# add temporal correlation structure, such that residuals of subsequent loops within an individual 
# show more similar residuals
sex_mean_ar =update(sex_mean_ide, correlation=corAR1(form=~loop_nr|territory/id))

# update using restricted likelyhood
sex_mean_reml =update(sex_mean_ar, method="REML")

# anova 
anova(sex_mean_reml)
# parameter estimates
summary(sex_mean_reml)

#### predict the differences in foraging distance between sexes
newdat <- data.frame(sex = levels(tmp3$sex)) 

pred1 <- predict(sex_mean_reml , newdat, level = 0, se = T)   
newdat$fit <- pred1$fit
newdat$upperCI<-  pred1$fit+1.96*pred1$se.fit
newdat$lowerCI<-  pred1$fit-1.96*pred1$se.fit

# backtransform to m scale
newdat[5:7] <- 10^newdat[2:4]

## how much further do 
newdat$fit.1[2]/newdat$fit.1[1]
## males forage 51% further away than females


## plot it
plot_mean_dist_all <-ggplot(newdat, aes(sex, fit.1)) +
  geom_jitter(data = tmp3,height = 0, width = 0.1, color = "grey", alpha = .5,aes (y = tmp3$meandist.loop, x = tmp3$sex ))+
  geom_point(size =3)+
  geom_linerange(data = newdat,aes(ymin=lowerCI.1,ymax=upperCI.1), size=1.5)+
  ylab("mean distance to nestbox [m]")+
  xlab("")+
  theme_classic()+
  ylim(0,3000)
plot_mean_dist


## alternatively, show only up to 600m for visualization
plot_mean_dist <-ggplot(newdat, aes(sex, fit.1)) +
  geom_jitter(data = tmp3,height = 0, width = 0.1, color = "grey", alpha = .5,aes (y = tmp3$meandist.loop, x = tmp3$sex ))+
  geom_point(size =3)+
  geom_linerange(data = newdat,aes(ymin=lowerCI.1,ymax=upperCI.1), size=1.5)+
  ylab("mean distance to nestbox [m]")+
  xlab("")+
  theme_classic()+
  ylim(0,600)
plot_mean_dist


#################################################
# 1B. Test whether sexes differ in their home range size
##################################################

setwd("~/Documents/Manuskripte/2018/Hoophoe/Urs_Analysen")
mydat  <- read.csv("Territory_occupancy.csv", sep=",")
str(mydat)
require(REEMtree)
require(vegan)  
require(ggpubr)
require(nlme)
require(MuMIn)

# how many dataset do we have
table(mydat$sex)
table( mydat$territory,mydat$sex, mydat$isopleth)

## 22 homeranges from 12 nests

######## first, analyze the 50% isopleth homerange
hist(subset(mydat, mydat$isopleth == 50)$area)
# looks ugly, log10 looks much better
hist(log10(subset(mydat, mydat$isopleth == 50)$area))

homerange<-lme(log10(area) ~ sex, random = ~1|territory,data = subset(mydat, mydat$isopleth == 50))
plot(homerange)
summary(homerange)

#try out variance functions
homerange_p =update(homerange, weights=varPower())
# no convergence
homerange_id =update(homerange, weights=varIdent(form=~1|sex))
homerange_exp =update(homerange, weights=varExp())
# no convergence
AICc(homerange, homerange_id)
# simply log10 transformed looks best

plot(homerange)
qqnorm(resid(homerange), main="qq-plot residuals")
qqline(resid(homerange))


# update with restricted likelyhood
homerange =update(homerange, method = "REML")
summary(homerange)
anova(homerange)

# no significant sex difference. Thus, show a boxplot for the illustration
# get predictions
newdat <- data.frame(sex = levels(mydat$sex)) 
predict(homerange , newdat, level = 0, se = T)   

#pred_dist <- predictmeans(homerange, "sex", trans=function(x) (10^x))
#pred_dist$mean_table

female
10^(0.5387963)
10^(0.5387963-2*0.1784837)
10^(0.5387963+2*0.1784837)

male
10^(0.6794979)
10^(0.6794979-2*0.1784837)
10^(0.6794979+2*0.1784837)


#######
# Look whether sexes differ in their  95% KDE 
hist(subset(mydat, mydat$isopleth == 95)$area)

##-> no transformation is best for 95
dat95 <-subset(mydat, mydat$isopleth == 95)

homerange<-lme(area ~ sex, random = ~1|territory,data = dat95)

homerange_p =update(homerange, weights=varPower())
homerange_id =update(homerange, weights=varIdent(form=~1|sex))
homerange_exp =update(homerange, weights=varExp())

AICc(homerange, homerange_id)
# again, untranformed is best
# Update model with restricted likelyhood
homerange =update(homerange, method = "REML")
summary(homerange)
anova(homerange)
### --> Again, no evidence that sexes differ in their homerange size, neither at 50 not 95%

##-> no transformation is best for 95
plot(homerange)
qqnorm(resid(homerange), main="qq-plot residuals")
qqline(resid(homerange))

newdat <- data.frame(sex = levels(mydat$sex)) 
predict(homerange , newdat, level = 0, se = T)   

female
51.83153
51.83153-2*23.92386
51.83153+2*23.92386

male
46.06088
#46.06088-2*23.92386
0
46.06088+2*23.92386


p <- ggplot(dat95, aes(sex, area))
a95 <-p + geom_boxplot()+
  ylab("homerange size [ha]")+ggtitle("95% HR kernel")+
  xlab("Sex")+
  theme_classic()

dat50 <-subset(mydat, mydat$isopleth == 50)
p <- ggplot(dat50, aes(sex, area))
a50 <-p + geom_boxplot()+
  ylab("homerange size [ha]")+ggtitle("50% HR kernel")+
  xlab("Sex")+
  theme_classic()

ggarrange(a50,a95,
          nrow = 1, ncol = 2,
          labels = c("A", "B"))

############################
#
# 2. Do sexes differ in their probability to return with a mole cricket as a function of foraging distance? 
#
##########################
tmp4$mole <- 1
tmp4$mole[tmp4$prey == "larvae"] <- 0

####### 
### plot it with Fränzis method
dat1 <- subset(tmp4, tmp5$meandist.loop <500)
dat1$log10di <- log10(dat1$meandist.loop)
tmp4$log10di <- log10(tmp4$meandist.loop)
mod <- glmer(mole ~  log10di +(1|territory/id) , na.action = "na.fail",tmp4, family=binomial)

summary(mod)

par(mfrow=c(2,2)) # divide the graphic window in 4 subregions
qqnorm(resid(mod)) # qq-plot of residuals
qqline(resid(propmod_mole2))
qqnorm(ranef(mod)$id[,1]) # qq-plot of the random effects
qqline(ranef(mod)$id[,1])
qqnorm(ranef(mod)$territory[,1]) # qq-plot of the random effects
qqline(ranef(mod)$territory[,1])

plot(fitted(mod), resid(mod)) # residuals vs fitted values
abline(h=0)

dat1$fitted <- fitted(mod) # fitted vs observed values
library(blmeco)
library(arm)
dispersion_glmer(mod)
[1]

?sim()
nsim <- 2000
bsim <- sim(mod, n.sim=nsim)
colnames(bsim@fixef) <- names(fixef(mod)) # for arm 1.6-10
fixef(mod)
apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))
dat1
summary(dat1)
newdat <- data.frame(log10di=log10(seq(1,2500,1)))
Xmat <- model.matrix(~log10di, newdat)

Xmat <- model.matrix(~log10di, data=newdat)

fitmat <- matrix(nrow=nrow(newdat), ncol=nsim)
newdat$pred <- plogis(Xmat %*% fixef(mod))

# get the credible intervals for these predictions:
predmat <- matrix(nrow=nrow(newdat),ncol=nsim)
for(i in 1:nsim) predmat[,i] <- plogis(Xmat %*% bsim@fixef[i,])
newdat$lwr <- apply(predmat,1,quantile,probs=0.025)
newdat$upr <- apply(predmat,1,quantile,probs=0.975)
newdat$dist <- seq(1,2500,1)

plot_pmole_dist <-
  ggplot(newdat, aes(dist, pred)) +
  geom_line(size=1.5)+
  geom_line(aes(dist, lwr), linetype = "dashed")+
  geom_line(aes(dist, upr), linetype = "dashed")+
  labs(x = "foraging distance [m]", y= "p (prey = mole cricket)")+
  theme_classic()+
  geom_point(data = dat1, alpha = 1/5, 
             mapping = aes( x = meandist.loop, y = mole))

############################
#
# 3. Do sexes differ in the foraging distance associated with different prey types?
#
############################

table(tmp3$prey, tmp3$loop_quality)

tmp4 <-tmp3
tmp4 <- subset(tmp3, tmp3$prey == "larvae" | tmp3$prey == "mole")
dim(tmp4)
dim(tmp3)
sum(table(tmp4$prey))
table(tmp4$prey)
levels(tmp4$prey)
# it sees like the model is only converging using the "good" loops, i.e. loops with 
# loop-quality = 2


### evt diesen hier einfügen 
##tmp5 <- tmp4
###

table(tmp4$prey)
tmp5 <- subset(tmp4, tmp4$loop_quality == "2")
tmp5 <- droplevels(tmp5)
tmp4 <- droplevels(tmp4)
m1mean_preyr <-lme(log10(meandist.loop) ~ sex*prey , random = ~1|territory/id, na.action = "na.omit",data = tmp4)
m1mean_prey_powr =update(m1mean_preyr, weights=varPower())
m1mean_preyr_ide =update(m1mean_preyr, weights=varIdent(form=~1|sex))
m1mean_preyr_exp =update(m1mean_preyr, weights=varExp())

AICc(m1mean_preyr,m1mean_prey_powr,  m1mean_preyr_ide, m1mean_preyr_exp)

## Varpower seems to be the best
plot(m1mean_prey_powr)
qqnorm(resid(m1mean_prey_powr), main="qq-plot residuals")
qqline(resid(m1mean_prey_powr))

#  add temporal correlation structure, such that residuals of subsequent loops within an individual 
m1mean_prey_powr_ar =update(m1mean_prey_powr, correlation=corAR1(form=~loop_nr|territory/id))

mean_prey_best =update(m1mean_prey_powr_ar, method="REML")

anova(mean_prey_best)
summary(mean_prey_best)


#pred_dist$mean_table

# significant sex*mole cricket interaction
require(effects)


# predict values
newdat <- expand.grid(sex = levels(tmp5$sex), prey = levels(tmp5$prey)) 
predict(mean_prey_best, newdat, level = 0,  se = T)   

newdat$estimate <- 10^(predict(mean_prey_best, newdat, level = 0,  se = T)$fit)  
newdat$lower <- 10^(log10(newdat$estimate) -2*predict(mean_prey_best, newdat, level = 0,  se = T)$se.fit )
newdat$upper <- 10^(log10(newdat$estimate) +2*predict(mean_prey_best, newdat, level = 0,  se = T)$se.fit)


## plot the results
newdat$prey_plot <- c("Larvae", "Mole cricket","Larvae", "Mole cricket" )
plot_prey_mean_dist <-ggplot(newdat, aes(prey, estimate)) +
  geom_jitter(data = tmp5,height = 0, width = 0.1, color = "grey", alpha = .5,aes (y = tmp5$meandist.loop, x = tmp5$prey ))+
  geom_point(size =3)+
  geom_linerange(data = newdat,aes(ymin=lower,ymax=upper), size=1.5)+
  ylab("mean distance to nestbox [m]")+
  xlab("")+
  facet_wrap(~sex)+
  ylim(0,600)+
  theme_classic()

ggarrange(plot_mean_dist,plot_prey_mean_dist,
          nrow = 1, ncol = 2,
          labels = c("A", "B"))

################
#
# Overvier of prey categories 
#
#################
dim(tmp3)
#711 feeding bouts 

sum(table(tmp3$prey))
# 571 bouts with prey

table(tmp3$prey)
# 571 prey items
571 - 2 
# 569 prey item ided
table(tmp3$prey, tmp3$sex)

###########
#
# 4. Does foraging distance change with habitat quality? (i.e. occupancy duration), 
#
##########

## do larger homeranges correlate with higher territory quality
## 
dat95
occ_1 <-lme(log10(area)~sex+ Occupancy, random = ~1|territory,data = dat95)
occ_2 <-lme(log10(area)~ Occupancy, random = ~1|territory,data = dat95)

AICc(occ_1, occ_2)
occ_1_p =update(occ_1, weights=varPower())
occ_1_id =update(occ_1, weights=varIdent(form=~1|sex))
occ_1_exp =update(occ_1, weights=varExp())

AICc(occ_1, occ_1_p, occ_1_id)
## -<id is best, but the untransformed has the best residuals
##-> no transformation is best for 95
plot(occ_1)
qqnorm(resid(occ_1), main="qq-plot residuals")
qqline(resid(occ_1))

occ_1_best =update(occ_1, method = "REML")
occ_2_best =update(occ_2, method = "REML")

summary(occ_1_best)
anova(occ_1_best)

summary(occ_2_best)
anova(occ_2_best)

anova(occ_2, occ_1)

AICc(occ_1, occ_2)
require(effects)
plot(allEffects(occ_1_best), type = "response", ci.style="lines", colors="black")

newdata <-data.frame(Occupancy = c(1:14),
                     territory = " B9")


newdata$lower<- (predict(occ_2, newdata =newdata,  level =0, se =T)$fit - 1.96*predict(newdata = newdata,occ_2, level =0, se =T)$se.fit)
newdata$fit<- (predict(occ_2, newdata =newdata,  level =0, se =T)$fit)  
newdata$upper<- (predict(occ_2, newdata =newdata,  level =0, se =T)$fit + 1.96*predict(newdata = newdata,occ_2, level =0, se =T)$se.fit) 

ggplot(newdata, aes(Occupancy, fit)) +
  geom_line(size=1.5)+
  geom_line(aes(Occupancy, lower), linetype = "dashed")+
  geom_line(aes(Occupancy, upper), linetype = "dashed")+
  geom_point(data = dat95, alpha = 1/2, 
             mapping = aes( x = Occupancy, y = log10(area)))+
  labs(x = "Territory occupancy [years]", y= expression(paste("log10 (home range size [ha])")))+
  theme_classic()

###########
#
# look at the picture from a differenc angel. Namely does the distance of the foiraging trip prediczt whether they bring in a caterpillar or a mole cricket? 
#
####
### add sex_prey_occupacy.R
##########
tmp5$mole <- 1
tmp5$mole[tmp5$prey == "larvae"] <- 0
names(tmp5)
dim(tmp5)
sum(tmp5$mole[tmp5$sex == "f"])

170 mole for males out of 281
110 mole for males out of 244
tmpx <-subset(tmp3, tmp3$loop_quality == "2")
with(tmpx, table(sex, prey))


require(lme4)

prop_mole0 <- glmer(mole ~  sex*log10(meandist.loop) + (1|territory/id) ,
                    na.action = "na.fail", subset(tmp5, tmp5$meandist.loop <500), family=binomial)

prop_mole1 <- glmer(mole ~  sex+log10(meandist.loop) + (1|territory/id) ,
                    na.action = "na.fail", subset(tmp5, tmp5$meandist.loop <500), family=binomial)
prop_mole2 <- glmer(mole ~  log10(meandist.loop) + (1|territory/id) ,
                    na.action = "na.fail", subset(tmp5, tmp5$meandist.loop <500), family=binomial)
prop_mole3 <- glmer(mole ~ 1 + (1|territory/id) ,
                    na.action = "na.fail", subset(tmp5, tmp5$meandist.loop <500), family=binomial)

prop_mole4 <- glmer(mole ~ sex+ (1|territory/id) ,
                    na.action = "na.fail", tmp5, family=binomial)
prop_mole5 <- glmer(mole ~ 1+ (1|territory/id) ,
                    na.action = "na.fail",tmp5, family=binomial)

summary(prop_mole4)

require("MuMIn")
require("nlme")
summary(prop_mole2)

library(car)
Anova(prop_mole2)


anova(prop_mole2)

anova(prop_mole4,prop_mole5)


anova(prop_mole4)

require(MuMIn)
AICc(prop_mole0,prop_mole1,prop_mole2,prop_mole3)
drop1(prop_mole0,prop_mole1)
drop1(prop_mole1)
anova(prop_mole0,prop_mole1)

anova(prop_mole1,prop_mole2)
anova(prop_mole3,prop_mole2)

str(allEffects(prop_mole1, partial.residuals=T), type = "response")
plot(allEffects(prop_mole0, partial.residuals=T), type = "response")

require(arm)

summary(allEffects(prop_mole1, partial.residuals=T), type = "response")

summary(allEffects(prop_mole, partial.residuals=T), type = "response")

allEffects(prop_mole0)


summary(prop_mole0)

allEffects(prop_mole2)[[1]]$lower

str(allEffects(prop_mole2))
invlogit(-1.4933)
invlogit(-3.74)


plot(allEffects(prop_mole2, partial.residuals=T), type = "response", ylab = "probability to return with mole",
     xlab = "mean loop distance [m]", 
     ci.style="lines", colors="black",rug=FALSE, main = "")

plot(allEffects(prop_mole0, partial.residuals=T), type = "response", ylab = "probability to return with mole",
     xlab = "mean loop distance [m]", 
     ci.style="lines", colors="black",rug=FALSE, main = "")



str(allEffects(prop_mole2), type = "response")


allEffects(prop_mole2)["log10(meandist.loop)"]
allEffects(prop_mole2)

[[“scores”]]
dist <- c(6, 100,200, 400, 500)
lower <- exp(c(-3.818, -0.525, 0.169, 0.815, 1.014)
             upper <- exp(c(0.702, 3.191, 3.92, 4.698, 4.959))
             pred <- exp(c(0.1739328, 0.7912602, 0.8854099, 0.9402985, 0.9519423  ))
             
             preddat<- data.frame(dist = dist, lower = lower, upper = upper, pred = pred)
             
             ggplot(preddat, aes(x =dist, y =pred))+
               geom_smooth()+
               geom_smooth(aes(x =dist, y =upper), color = "red")+
               geom_smooth(aes(x =dist, y =lower), color = "red")
             
             -3.818 -0.525 0.169 0.815 1.014
             
             
             coef(prop_mole2)
             require(lattice)
             plot(subset(tmp5, tmp5$meandist.loop <500)$meandist.loop, jitter(subset(tmp5, tmp5$meandist.loop <500)$mole))
             ??effects
             plot(allEffects(prop_mole2, partial.residuals=T), type = "response")
             
             
             str(allEffects(prop_mole2, partial.residuals=T))
             
             
             ###
             
             tmp6 <- subset(tmp5, tmp5$sex == "f")
             
             m1mean_preyr <-lme(log10(meandist.loop) ~ prey , random = ~1|id, na.action = "na.omit",data = tmp6)
             
             m1mean_prey_powr =update(m1mean_preyr, weights=varPower())
             m1mean_preyr_ide =update(m1mean_preyr, weights=varIdent(form=~1|sex))
             m1mean_preyr_exp =update(m1mean_preyr, weights=varExp())
             
             AICc(m1mean_preyr,m1mean_prey_powr,  m1mean_preyr_ide, m1mean_preyr_exp)
             
             ## powr seems to be the best
             
             plot(m1mean_preyr)
             plot(m1mean_preyr_exp)
             plot(m1mean_preyr_ide)
             plot(m1mean_prey_powr)
             
             
             qqnorm(resid(m1mean_preyr_exp), main="qq-plot residuals")
             qqline(resid(m1mean_preyr_exp))
             
             
             median_best_ar1 =update(m1mean_prey_powr, correlation=corAR1(form=~loop_nr|id))
             
             anova(median_best_ar1,m1mean_preyr_exp)
             
             mean_prey_best =update(m1mean_prey_powr, method="REML")
             
             anova(mean_prey_best)
             
             plot(allEffects(mean_prey_best))
             
             
             
             #########
             
             dat1 <- subset(tmp5, tmp5$meandist.loop <500)
             dat1$log10di <- log10(dat1$meandist.loop)
             mod <- glmer(mole ~  log10di + (1|territory/id) ,
                          na.action = "na.fail",dat1, family=binomial)
             
             summary(mod)
             prop_mole2
             par(mfrow=c(2,2)) # divide the graphic window in 4 subregions
             qqnorm(resid(mod)) # qq-plot of residuals
             qqline(resid(propmod_mole2))
             qqnorm(ranef(mod)$id[,1]) # qq-plot of the random effects
             qqline(ranef(mod)$id[,1])
             qqnorm(ranef(mod)$territory[,1]) # qq-plot of the random effects
             qqline(ranef(mod)$territory[,1])
             
             plot(fitted(mod), resid(mod)) # residuals vs fitted values
             abline(h=0)
             
             dat1$fitted <- fitted(mod) # fitted vs observed values
             library(blmeco)
             dispersion_glmer(prop_mole2)
             [1]
             mod <- prop_mole2
             nsim <- 2000
             bsim <- sim(mod, n.sim=nsim)
             colnames(bsim@fixef) <- names(fixef(mod)) # for arm 1.6-10
             fixef(mod)
             apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))
             dat1
             summary(dat1)
             newdat <- data.frame(log10di=log10(seq(1,500,1)))
             Xmat <- model.matrix(~log10di, newdat)
             
             Xmat <- model.matrix(~log10di, data=newdat)
             
             fitmat <- matrix(nrow=nrow(newdat), ncol=nsim)
             newdat$pred <- plogis(Xmat %*% fixef(mod))
             
             # get the credible intervals for these predictions:
             predmat <- matrix(nrow=nrow(newdat),ncol=nsim)
             for(i in 1:nsim) predmat[,i] <- plogis(Xmat %*% bsim@fixef[i,])
             newdat$lwr <- apply(predmat,1,quantile,probs=0.025)
             newdat$upr <- apply(predmat,1,quantile,probs=0.975)
             newdat$dist <- seq(1,500,1)
             
             pdf("prob_mole_dist.pdf",width=5,height=3,paper='special')
             
             plot_pmole_dist <-
               ggplot(newdat, aes(dist, pred)) +
               geom_line(size=1.5)+
               geom_line(aes(dist, lwr), linetype = "dashed")+
               geom_line(aes(dist, upr), linetype = "dashed")+
               labs(x = "foraging bout distance [m]", y= "p (prey = mole cricket)")+
               theme_classic()+
               geom_point(data = dat1, alpha = 1/5, 
                          mapping = aes( x = meandist.loop, y = mole))
             
             plot_pmole_dist
             dev.off()
             