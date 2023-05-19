##reported code shows all code that was used in manuscript

.libPaths("C:/Users/tp350/Google Drive/PhD/R/Rpackages")

install.packages("lme4")
install.packages("car")
install.packages("tidyverse")
install.packages("ggpubr")
install.packages("lmerTest")
install.packages("installr")
install.packages("aods3")
install.packages("DHARMa")
install.packages("rcompanion")
install.packages("performance")
install.packages("see")
install.packages("jtools")
install.packages("scales")
install.packages("performance")
install.packages("see")
install.packages("ggeffects")
install.packages("interactions")
install.packages("blme")
install.packages("plotrix")
install.packages("insight")
install.packages("lsmeans")
install.packages("rstatix")
library(car)
library(lme4)
library(tidyverse)
library(ggpubr)
library(lmerTest)
library(installr)
library(aods3)
library(DHARMa)
library(rcompanion)
library(performance)
library(see)
library(jtools)
library(scales)
library(performance)
library(see)
library(ggeffects)
library(interactions)
library(blme)
library(data.table)
library(plotrix)
library(insight)
library(lsmeans)
library(multcompView)
library(plyr)
library(rstatix)

###Total Offspring production

Prod2 = read.csv("Foraging_trips_Prod.csv", header = T)
###data on day largest group size seen 

tapply(Prod2$Total.Brood,Prod$Nest.type, range)
tapply(Prod2$Active,Prod$Nest.type, mean)
length(which(Prod2$Active == 0))

Prod2$Nest<- as.factor(Prod2$Nest)
Prod2$Date<- as.factor(Prod2$Date)
Prod2$Nest.type<- as.factor(Prod2$Nest.type)
Prod2$Active <- as.factor(Prod2$Active)
Prod2$Active <- as.numeric(as.character(Prod2$Active))
Prod2$TotalDays <- Prod2$Total.Brood/Prod2$Brood_per_active_day


###tapply and wilcox, are things the same within groups
tapply(Prod2$Group.size,Prod2$Nest.type,mean)
wilcox.test(Group.size~Nest.type, data=Prod2, exact = FALSE)

tapply(Prod2$Av_w_size,Prod2$Nest.type,mean)
wilcox.test(Av_w_size~Nest.type, data=Prod2, exact = FALSE)

tapply(Prod2$TotalDays,Prod2$Nest.type,mean)
wilcox.test(TotalDays~Nest.type, data=Prod2, exact = FALSE)



###labels for facet wrap graphs
nest.names <- list('Ql'= "Worker-queen", 'QR' = "Foundress-queen")
nest_labeller <- function(variable,value){
  return(nest.names[value])
}

##exploratory graph
summary(Prod2)
Prod2plot <- ggplot(Prod2, aes(x=Group.size, y=Total.Brood)) +
  geom_jitter(shape=18, width = 0.35) +
  facet_wrap(~Nest.type, labeller = nest_labeller) +
  theme_bw()
Prod2plot

summary(Prod2)
###reported models are in their post simplification format used to produce graphs
BroodCountTan3 <- glm(Total.Brood ~ Nest.type+Group.size +TotalDays, family=quasipoisson, data=Prod2)
BroodCountTan5 <- glm(Total.Brood ~ Nest.type, family=quasipoisson, data=Prod2)

summary(BroodCountTan3)
plot(BroodCountTan3)
res <- simulateResiduals(BroodCountTan3)
plot(res)
testDispersion(res)
plotNormalHistogram(Prod$Total.Brood)
effect_plot(BroodCountTan3, pred = Nest.type, interval = TRUE, plot.points = TRUE)
effect_plot(BroodCountTan3, pred = Group.size, interval = TRUE, plot.points = TRUE)
effect_plot(BroodCountTan3, pred = TotalDays, interval = TRUE, plot.points = TRUE)

summary(BroodCountTan5)
plot(BroodCountTan5)
res <- simulateResiduals(BroodCountTan5)
plot(res)
testDispersion(res)
plotNormalHistogram(Prod$Total.Brood)
effect_plot(BroodCountTan5, pred = Nest.type, interval = TRUE, plot.points = TRUE)

drop1(BroodCountTan3, test = "Chisq")

Brood.rg1 <- ref.grid(BroodCountTan3)
summary(Brood.rg1)
Brood.rg1$prediction <- exp(Brood.rg1$prediction)
Brood.rg1 <- as.data.frame(Brood.rg1)
Broodmeans <- lsmeans(BroodCountTan3, "Group.size", at = list(Group.size = c(1,2,3,4,5,6,7,8)))
Broodmeans2 <- lsmeans(BroodCountTan3, "TotalDays", at = list(TotalDays = c(6,7,8,9,10,11,12,13,14,15,16,17,18)))
Broodmeans <- as.data.frame(Broodmeans)
Broodmeans$lsmean <- exp(Broodmeans$lsmean)
Broodmeans$asymp.LCL <- exp(Broodmeans$asymp.LCL)
Broodmeans$asymp.UCL <- exp(Broodmeans$asymp.UCL)
Broodmeans2 <- as.data.frame(Broodmeans2)
Broodmeans2$lsmean <- exp(Broodmeans2$lsmean)
Broodmeans2$asymp.LCL <- exp(Broodmeans2$asymp.LCL)
Broodmeans2$asymp.UCL <- exp(Broodmeans2$asymp.UCL)

###Graph

newdata <- expand.grid(Nest.type=unique(Prod2$Nest.type), 
                       Group.size=seq(1,8,by=0.1),
                       TotalDays=unique(Prod2$TotalDays),
                       KEEP.OUT.ATTRS = TRUE, stringsAsFactors = TRUE)

newdata2 <- newdata %>% 
            filter(Nest.type == "QL" & Group.size <=6)
newdata3 <- newdata %>% 
  filter(Nest.type == "QR")

newdata <- rbind(newdata3,newdata2)


## prediction function to use in bootstrap routine
mypredictdf <- function (model, newdata, level=0.95){
  pred <- stats::predict(model, newdata = newdata, se =TRUE, type = "link")
  std <- qnorm(level/2 + 0.5)
  data.frame(newdata,
             y = model$family$linkinv(as.vector(pred$fit)),
             ymin = model$family$linkinv(as.vector(pred$fit - std * pred$se)),
             ymax = model$family$linkinv(as.vector(pred$fit + std * pred$se)), 
             se = as.vector(pred$se))
}
####making graph####
pdf <- mypredictdf(BroodCountTan3, newdata=newdata)

names(pdf)[4] <- "Total.Brood"
head(pdf)
gplot_pred1 <- ggplot(Prod2, aes(x=Group.size, y=Total.Brood)) +
  geom_jitter(shape=18, width = 0.1) +
  facet_wrap(~Nest.type, labeller = nest_labeller, scales = "free_x") +
  geom_smooth(data=pdf, se=FALSE, linetype = 1) +
  geom_smooth(data=pdf, aes(y = ymax), color = 'red', linetype = 2, se=FALSE) +
  geom_smooth(data=pdf, aes(y = ymin), color = 'red', linetype = 2, se=FALSE) +
  labs(title="Figure 1", x="Number of Foragers", y="Total Offspring") +
  theme_bw()
gplot_pred1

ggsave("Figure 1 Free.png", width = 5, height = 5)


###Foraging effort

Prod3 = read.csv("Foraging_trips_perday.csv", header = T)
###data broken down into each day of each nest
Prod3$Nest<- as.factor(Prod3$Nest)
Prod3$Date<- as.factor(Prod3$Date)
Prod3$Nest.type<- as.factor(Prod3$Nest.type)

Prod3$TripBrood <- Prod3$No..arrivals/Prod3$Brood_per_active_day
Prod3$TotalDays <- Prod3$Total.Brood/Prod3$Brood_per_active_day
Prod3$TotalTrips <- Prod3$TotalDays*Prod3$No..arrivals
Prod3$TripPerWorker <- Prod3$TotalTrips/Prod3$Group.size
Prod3$TripperWorkerDay <- Prod3$No..arrivals/Prod3$Group.size

summary(Prod3)

###Trip per brood
##exploratory graph
visual3 <- ggplot(Prod3, aes(x=Group.size, y=Brood_per_active_day, color=Nest.type)) +
  geom_point(shape=18) +
  labs(title="Cost of brood production", x="Group Size", y="Brood per day") +
  theme_classic() +
  scale_color_discrete(name="Queen Status",
                       breaks=c("QL", "QR"),
                       labels=c("Queenless", "Queenright"))
visual3

###reported models are in their post simplification format used to produce graphs
Cost <- glmer(TripBrood ~ Nest.type+Group.size + (1|Nest), family=gaussian(link = "log"), data=Prod3)
summary(Cost)
resCost <- simulateResiduals(Cost)
plot(resCost)
testDispersion(resCost)
effect_plot(Cost, pred = Group.size, interval = TRUE, plot.points = TRUE)
effect_plot(Cost, pred = Nest.type, interval = TRUE, plot.points = TRUE)
check_model(Cost, check = "reqq")
ydrop1(Cost, test="Chisq")
check_model
Cost.rg1 <- ref.grid(Cost)
summary(Cost.rg1)
Costmeans <- lsmeans(Cost, "Group.size", at = list(Group.size = c(1,2,3,4,5,6,7,8)))
Costmeans <- as.data.frame(Costmeans)
Costmeans$lsmean <- exp(Costmeans$lsmean)
Costmeans$asymp.LCL <- exp(Costmeans$asymp.LCL)
Costmeans$asymp.UCL <- exp(Costmeans$asymp.UCL)

###graph### 

####confidence intervals for graphs####

easyPredCI <- function(model,newdata=NULL,alpha=0.05) {
  ## baseline prediction, on the linear predictor (logit) scale:
  pred0 <- predict(model,re.form=NA,newdata=newdata)
  ## fixed-effects model matrix for new data
  X <- model.matrix(formula(model,fixed.only=TRUE)[-2],newdata)
  beta <- fixef(model) ## fixed-effects coefficients
  V <- vcov(model)     ## variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  ## inverse-link function
  linkinv <- family(model)$linkinv
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(conf.low=pred0-crit*pred.se,
                conf.high=pred0+crit*pred.se))
}

costframe <- expand.grid(Nest.type=unique(Prod3$Nest.type), 
                         Group.size=seq(1,8,by=0.1),
                         KEEP.OUT.ATTRS = TRUE, stringsAsFactors = TRUE)

costframe2 <- costframe %>% 
  filter(Nest.type == "QL" & Group.size <=5)
costframe3 <- costframe %>% 
  filter(Nest.type == "QR")

costframe <- rbind(costframe2,costframe3)

costframe$output<-predict(Cost, newdata=costframe, type="response", se.fit=TRUE, re.form=NA)

Cost_predCI <- easyPredCI(Cost,newdata=costframe)
Cost_predCI2<- as.data.frame(Cost_predCI)
costframe$conf.low <- Cost_predCI2$conf.low
costframe$conf.high <- Cost_predCI2$conf.high

summary(costframe)
costframe$TripBrood <- costframe$output
summary(costframe)
summary(Prod)


gplot_cost <- ggplot(Prod3, aes(x=Group.size, y=TripBrood)) +
  geom_jitter(shape=18, width = 0.35) +
  facet_wrap(~Nest.type, labeller = nest_labeller, scales = "free_x") +
  geom_smooth(data=costframe, se=FALSE, linetype = 1) +
  geom_smooth(data=costframe, aes(y = conf.high), color = 'red', linetype = 2, se=FALSE) +
  geom_smooth(data=costframe, aes(y = conf.low), color = 'red', linetype = 2, se=FALSE) +
  labs(title="Figure 2", x="Number of Foragers", y="Trips per offspring") +
  theme_bw()
gplot_cost

ggsave("Figure 2 Free.png", width = 5, height = 5)

####Trips per worker####

##exploratory graph
visual4 <- ggplot(Prod, aes(x=Group.size, y=(No..arrivals/Group.size), color=Nest.type)) +
  geom_point(shape=18) +
  labs(title="Are workers less efficent at larger group sizes?", x="No. of Workers", y="Trips per worker") +
  theme_classic() +
  scale_color_discrete(name="Queen Status",
                       breaks=c("QL", "QR"),
                       labels=c("Queenless", "Queenright"))
visual4


###reported models are in their post simplification format used to produce graphs
TripWorker<- glmer((No..arrivals/Group.size) ~ Nest.type+Group.size + (1|Nest), family=gaussian(link = "log"), data=Prod3)
summary(TripWorker)
effect_plot(TripWorker, pred = Group.size, interval = TRUE, plot.points = TRUE)
effect_plot(TripWorker, pred = Nest.type, interval = TRUE, plot.points = TRUE)

drop1(TripWorker, test = "Chisq")

resWorker <- simulateResiduals(TripWorker)
plot(resWorker)
testDispersion(resWorker)
plotNormalHistogram(Prod$TripPerWorker)
check_model(TripWorker, check = "reqq")

TripWorker.rg1 <- ref.grid(TripWorker)
summary(TripWorker.rg1)
TripWorkerMeans <- lsmeans(TripWorker, "Group.size", at = list(Group.size = c(1,2,3,4,5,6,7,8)))
TripWorkerMeans <- as.data.frame(TripWorkerMeans)
TripWorkerMeans$lsmean <- exp(TripWorkerMeans$lsmean)
TripWorkerMeans$asymp.LCL <- exp(TripWorkerMeans$asymp.LCL)
TripWorkerMeans$asymp.UCL <- exp(TripWorkerMeans$asymp.UCL)
TripWorkerMeans2 <- lsmeans(TripWorker, "Nest.type")
TripWorkerMeans2 <- as.data.frame(TripWorkerMeans2)
TripWorkerMeans2$lsmean <- exp(TripWorkerMeans2$lsmean)
TripWorkerMeans2$asymp.LCL <- exp(TripWorkerMeans2$asymp.LCL)
TripWorkerMeans2$asymp.UCL <- exp(TripWorkerMeans2$asymp.UCL)

####graph###

workframe <- expand.grid(Nest.type=unique(Prod3$Nest.type), 
                         Group.size=seq(1,8,by=0.1),
                         KEEP.OUT.ATTRS = TRUE, stringsAsFactors = TRUE)

workframe2 <- workframe %>% 
  filter(Nest.type == "QL" & Group.size <=5)
workframe3 <- workframe %>% 
  filter(Nest.type == "QR")

workframe <- rbind(workframe2,workframe3)

workframe$output<-predict(TripWorker, newdata=workframe, type="response", se.fit=TRUE, re.form=NA)

Work_predCI <- easyPredCI(TripWorker,newdata=workframe)
Work_predCI2<- as.data.frame(Work_predCI)
workframe$conf.low <- Work_predCI2$conf.low
workframe$conf.high <- Work_predCI2$conf.high

summary(workframe)
workframe$TripperWorkerDay <- workframe$output
summary(workframe)

gplot_work <- ggplot(Prod3, aes(x=Group.size, y=TripperWorkerDay)) +
  geom_jitter(shape=18, width = 0.35) +
  facet_wrap(~Nest.type, labeller = nest_labeller, scales = "free_x") +
  geom_smooth(data=workframe, se=FALSE, linetype = 1) +
  geom_smooth(data=workframe, aes(y = conf.high), color = 'red', linetype = 2, se=FALSE) +
  geom_smooth(data=workframe, aes(y = conf.low), color = 'red', linetype = 2, se=FALSE) +
  labs(title="Figure 4", x="Number of Foragers", y="Trips per Worker") +
  theme_bw()
gplot_work

ggsave("Figure 4 Free.png", width = 5, height = 5)



###Total foraging time 
##exploratory graph
visual2 <- ggplot(Prod3, aes(x=Group.size, y=Foraging.Time.Total, color=Nest.type)) +
  geom_point(shape=18) +
  labs(title="Time off nest", x="No. of Workers", y="Time Foraging") +
  theme_classic() +
  scale_color_discrete(name="Queen Status",
                       breaks=c("QL", "QR"),
                       labels=c("Queenless", "Queenright"))
visual2

###reported models are in their post simplification format used to produce graphs
GroupFor <- glmer(Foraging.Time.Total ~ Nest.type + Group.size + (1|Nest), family=gaussian(link ="log"), data=Prod3)
summary(GroupFor)
effect_plot(GroupFor, pred = Nest.type, interval = TRUE, plot.points = TRUE)
effect_plot(GroupFor, pred = Group.size, interval = TRUE, plot.points = TRUE)

check_model(GroupFor, check = "reqq")
resGroup <- simulateResiduals(GroupFor)
plot(resGroup)
plotNormalHistogram(Prod$Foraging.Time.Total)
plotNormalHistogram(Prod$Group.size)
testDispersion(resProp)
drop1(GroupFor, test = "Chisq")

GroupFor.rg1 <- ref.grid(GroupFor)
summary(GroupFor.rg1)
GroupFormeans <- lsmeans(GroupFor, "Group.size", at = list(Group.size = c(1,2,3,4,5,6,7,8)))
GroupFormeans <- as.data.frame(GroupFormeans)
GroupFormeans$lsmean <- exp(GroupFormeans$lsmean)
GroupFormeans$asymp.LCL <- exp(GroupFormeans$asymp.LCL)
GroupFormeans$asymp.UCL <- exp(GroupFormeans$asymp.UCL)

###graph
Forframe <- expand.grid(Nest.type=unique(Prod3$Nest.type), 
                        Group.size=seq(1,8,by=0.1),
                        KEEP.OUT.ATTRS = TRUE, stringsAsFactors = TRUE)

Forframe$output<-predict(GroupFor, newdata=Forframe, type="response", se.fit=TRUE, re.form=NA)

For_predCI <- easyPredCI(GroupFor,newdata=Forframe)
For_predCI2<- as.data.frame(For_predCI)
Forframe$conf.low <- For_predCI2$conf.low
Forframe$conf.high <- For_predCI2$conf.high

Forframe2 <- Forframe %>% 
  filter(Nest.type == "QL" & Group.size <=5)
Forframe3 <- Forframe %>% 
  filter(Nest.type == "QR")

Forframe <- rbind(Forframe2,Forframe3)


summary(Forframe)
Forframe$Foraging.Time.Total <- Forframe$output
summary(Forframe)

gplot_pred2 <- ggplot(Prod3, aes(x=Group.size, y=Foraging.Time.Total)) +
  geom_jitter(shape=18, width = 0.35) +
  facet_wrap(~Nest.type, labeller = nest_labeller, scales = "free_x") +
  geom_smooth(data=Forframe, se=FALSE, linetype = 1) +
  geom_smooth(data=Forframe, aes(y = conf.high), color = 'red', linetype = 2, se=FALSE) +
  geom_smooth(data=Forframe, aes(y = conf.low), color = 'red', linetype = 2, se=FALSE) +
  labs(title="Figure 3", x="Number of Foragers", y="Total time spent foraging (minutes)") +
  theme_bw()
gplot_pred2

ggsave("Figure 3 Free.png", width = 5, height = 5)
