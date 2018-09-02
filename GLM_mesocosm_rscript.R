#####Blackman et al 2018 GLM R script##########

#set working directory
#load data
#rename file

drb <- drb_mesocosm

#Check for confounding explanatory variables###
###############################################

plot(drb$Density, drb$TotalBio)
abline(lm(drb$TotalBio ~ drb$Density), lty=1, col="red")

#perform a Pearson's product moment correlation
cor.test(drb$Density,drb$TotalBio, method="pearson") 

#Pearson's product-moment correlation

#data:  drb$Density and drb$TotalBio
#t = 79.071, df = 52, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  0.9928556 0.9976109
#sample estimates:
#  cor 
#0.9958673
##massive correlation only include one of these in the model

#are density and average biomass confounding variables? 
plot(drb$Density, drb$AvBio)
abline(lm(drb$AvBio ~ drb$Density), lty=1, col="red")

#perform a Pearson's product moment correlation
cor.test(drb$Density,drb$AvBio, method="pearson")  

#Pearson's product-moment correlation
#data:  drb$Density and drb$AvBio
#t = -7.8899, df = 52, p-value = 1.901e-10
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#  -0.8399018 -0.5862547
#sample estimates:
#  cor 
#-0.7381444  
##negative correlation only include one of these in the model

#Test individual models 
#we have binomial data so doing a binomial glm with default logit link

##############################################################################################################
#GLM1: Hours + TotalBio
#start with full models (hours + TotalBio) then do model simplification

glmdrb1 <-glm(drb$stPCR ~ drb$Hours + drb$TotalBio, family=binomial())
summary(glmdrb1) 

#model checking
#so we can see all 4 model checking plots at once
par(mfrow=c(2,2)) 
plot(glmdrb1)
qqnorm(resid(glmdrb1,type="deviance"))
qqline(resid(glmdrb1,type="deviance"),col="red",lwd=2) 
hist(resid(glmdrb1,type="deviance"),freq=F)  

#check the if the residual deviance follows a chisq distribution. If it does the model is OK!
1-pchisq(49.368, 51) ##0.5384768 = not significant so the model is OK!!

#check the analysis of deviance table
anova(glmdrb1, test = "Chisq")

#check dispersal statistic of model, must be <2
overdispersion1 <- glmdh1$deviance/ glmdh1$df.residual
overdispersion1

##############################################################################################################
#GLM2: Hours + AvBio

glmdrb2 <-glm(drb$stPCR ~ drb$Hours + drb$AvBio, family=binomial())
summary(glmdrb2) 

#model checking
par(mfrow=c(2,2)) #so we can see all 4 model checking plots at once
plot(glmdrb2)
qqnorm(resid(glmdrb2,type="deviance"))
qqline(resid(glmdrb2,type="deviance"),col="red",lwd=2) 
hist(resid(glmdrb2,type="deviance"),freq=F)  

#check the if the residual deviance is significant
1-pchisq(52.305,51) ##0.4230644 = not significant so the model is OK

#check the analysis of deviance table
anova(glmdrb2, test = "Chisq")

#check dispersal statistic of model, must be <2
overdispersion2 <- glmdh2$deviance/ glmdh2$df.residual
overdispersion2

##############################################################################################################
#glm 3 : Hours + density

glmdrb3 <-glm(drb$stPCR ~ drb$Hours + drb$Density, family=binomial())
summary(glmdrb3) 

#####model checking#####
par(mfrow=c(2,2)) #so we can see all 4 model checking plots at once
plot(glmdrb3)
qqnorm(resid(glmdrb3,type="deviance"))
qqline(resid(glmdrb3,type="deviance"),col="red",lwd=2) 
hist(resid(glmdrb3,type="deviance"),freq=F)  

#check the if the residual deviance is significant
1-pchisq(49656,51) ##0 = significant so the model is not OK

####check the analysis of deviance table####
anova(glmdrb3, test = "Chisq")

#check dispersal statistic of model, must be <2
overdispersion3 <- glmdh3$deviance/ glmdh3$df.residual
overdispersion3

##############################################################################################################
#glmdh4: TotalBio

glmdrb4 <-glm(drb$stPCR ~ drb$TotalBio, family=binomial())
summary(glmdrb4) 

#model checking
par(mfrow=c(2,2)) #so we can see all 4 model checking plots at once
plot(glmdhrb4)
qqnorm(resid(glmdrb4,type="deviance"))
qqline(resid(glmdrb4,type="deviance"),col="red",lwd=2) 
hist(resid(glmdrb4,type="deviance"),freq=F)  

#check the if the residual deviance is significant
1-pchisq(51.718,52) ##0.4849402 = not significant so the model is OK!!

####check the analysis of deviance table####
anova(glmdrb4, test = "Chisq")

#check dispersal statistic of model, must be <2
overdispersion4 <- glmdh4$deviance/ glmdh4$df.residual
overdispersion4

##############################################################################################################
#glm 5:Hours

glmdrb5 <-glm(drb$stPCR ~ drb$Hours, family=binomial())
summary(glmdrb5)

#model checking
par(mfrow=c(2,2)) #so we can see all 4 model checking plots at once
plot(glmdrb5)
qqnorm(resid(glmdrb5,type="deviance"))
qqline(resid(glmdrb5,type="deviance"),col="red",lwd=2) 
hist(resid(glmdrb5,type="deviance"),freq=F)  

#check the if the residual deviance is significant
1-pchisq(57.584,52) # 0.2762671 = not significant model ok!

####check the analysis of deviance table####
anova(glmdrb5, test = "Chisq")

#check dispersal statistic of model, must be <2
overdispersion5 <- glmdh5$deviance/ glmdh5$df.residual
overdispersion5

##############################################################################################################
#glm 6:Density

glmdrb6 <-glm(drb$stPCR ~ drb$Density, family=binomial())
summary(glmdrb6)

#check the if the residual deviance is significant
1-pchisq(51.988,52) ##0.4743832 - not significant so model is ok!

#model checking glm6
par(mfrow=c(2,2)) #so we can see all 4 model checking plots at once

plot(glmdrb6)
qqnorm(resid(glmdrb6,type="deviance"))
qqline(resid(glmdrb6,type="deviance"),col="red",lwd=2) 
hist(resid(glmdrb6,type="deviance"),freq=F)  

####check the analysis of deviance table####
anova(glmdrb6, test = "Chisq")

#check dispersal statistic of model, must be <2
overdispersion6 <- glmdh6$deviance/ glmdh6$df.residual
overdispersion6

##############################################################################################################
#glm 7: AvBio

glmdrb7 <-glm(drb$stPCR ~ drb$AvBio, family=binomial())
summary(glmdrb7) 

#model checking glmdrb7
par(mfrow=c(2,2)) #so we can see all 4 model checking plots at once
plot(glmdrb7)
qqnorm(resid(glmdrb7,type="deviance"))
qqline(resid(glmdrb7,type="deviance"),col="red",lwd=2) 
hist(resid(glmdrb7,type="deviance"),freq=F)  

#check the if the residual deviance is significant
1-pchisq(54.573, 52) ## 0.3769788 = not significant so the model is OK

####check the analysis of deviance table####
anova(glmdrb7, test = "Chisq")

#check dispersal statistic of model, must be <2
overdispersion7 <- glmdh7$deviance/ glmdh7$df.residual
overdispersion7

####model checking and output table####
install.packages("AICcmodavg")
install.packages("bbmle")
install.packages("broom")

Modnames <- c("glmdrb1","glmdrb2","glmdrb3","glmdrb4", "glmdrb5","glmdrb6","glmdrb7")
summary.table <- do.call(rbind, lapply(list(glmdrb1,glmdrb2,glmdrb3,glmdrb4,glmdrb5,glmdrb6,glmdrb7), broom::glance))

table.cols <- c("df.residual", "deviance", "AIC")
reported.table <- summary.table[table.cols]
names(reported.table) <- c("Resid. Df", "Resid. Dev", "AIC")

reported.table[['dAIC']] <-  with(reported.table, AIC - min(AIC))
reported.table[['weight']] <- with(reported.table, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table$AIC <- NULL
reported.table$weight <- round(reported.table$weight, 2)
reported.table$dAIC <- round(reported.table$dAIC, 1)
row.names(reported.table) <- Modnames

######anova##
anova(glmdrb1, glmdrb2, glmdrb3, glmdrb4, glmdrb5, glmdrb6, glmdrb7, test = "Chisq")

##repeat for next IAS
###############################################################################################################


#Reference
#Burnham and Anderson 2002, delta AIC value <2 from the best model - important
#Dougherty et al 2016 AIC, uses delta AIC to compare models
#http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/

