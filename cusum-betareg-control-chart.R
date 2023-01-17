#############################################################################
#                                                                           #
#       PROGRAM: cusum-betareg-control-chart.R								              #
#                                                                           #
#       USAGE: Computation of the CUSUM beta regression control chart       #
#              using the quantile residual                                  #                                             #                                                                           #
#       AUTHOR: Cristine Rauber Oliveira	                                  #
#							                                                              #
#############################################################################



# getting the dataset from the rattle package
data(weatherAUS)
head(weatherAUS)
tail(weatherAUS)

# checking the dimension of the dataset
dim(weatherAUS)

#checking the locations available in this dataset
table(weatherAUS$Location)

# obtaining the data for Sydney, which is the city we are interested in 
# analysing the relative humidity
df_rh <- na.omit(weatherAUS[weatherAUS$Location == "Sydney",])
attach(df_rh)

# phase I data: we use these observations to fit the model and get 
# the parameter estimates
df_rh_p1 <- df_rh[1:845, c(3,4,5,6,7,15,17,19)]

# phase II data: we monitor these observations after fitting the model to 
# the phase I data and getting the parameter estimates
df_rh_p2 <- df_rh[846:1690, c(3,4,5,6,7,15,17,19)]

# feature names
names(df_rh_p1)

# summary of each feature
summary(df_rh_p1)

# transforming the target to a uniform scale (the beta regression model 
# only applies to the target in the unit interval (0,1))
df_rh_p1$Humidity3pm <- df_rh_p1$Humidity3pm/100
df_rh_p2$Humidity3pm <- df_rh_p2$Humidity3pm/100

# fitting the beta regression model to the data in phase I
fit <- betareg(Humidity3pm ~ MinTemp + MaxTemp + Rainfall + Evaporation + 
                 Pressure3pm + Cloud3pm | MinTemp + Sunshine + 
                 Pressure3pm, data = df_rh_p1)
summary(fit)

# obtaining the design matrices for both submodels
X_mu <- cbind(1, df_rh_p1[,c(1,2,3,4,7,8)])
X_phi <- cbind(1, df_rh_p1[,c(1,5,7)])

# linear predictor for each submodel
eta_mu <- as.matrix(X_mu)%*%fit$coefficients$mean
eta_phi <- as.matrix(X_phi)%*%fit$coefficients$precision

shape1 <- exp(eta_mu)/(1+exp(eta_mu))
shape2 <- exp(eta_phi)

# parameters of the beta regression distribution
p <- shape1*shape2
q <- shape2-(shape1*shape2)

rh <- df_rh_p1$Humidity3pm

# quantile residual
residual <- qnorm(pbeta(rh, p, q))

# plot of the residuals
plot(residual, caption = NULL, sub.caption = NULL, ylim = c(-4,4), 
     ylab = "Quantile residual", xlab = "Index")
abline(h = c(-3,3), lty = 2)
abline(h = 0, col = "red")

# qqplot of the residuals
qqnorm(residual, caption = NULL, main = "", ylab = "Empirical quantiles",  
       xlab = "Theoretical quantiles")
qqline(residual)

# design matrices for the data in phase II
X_mu2 <- cbind(1, df_rh_p2[,c(1,2,3,4,7,8)])
X_phi2 <- cbind(1, df_rh_p2[,c(1,5,7)])

# linear predictor for each submodel
eta_mu2 <- as.matrix(X_mu2)%*%fit$coefficients$mean
eta_phi2 <- as.matrix(X_phi2)%*%fit$coefficients$precision

shape1_2 <- exp(eta_mu2)/(1+exp(eta_mu2))
shape2_2 <- exp(eta_phi2)

p2 <- shape1_2*shape2_2
q2 <- shape2_2-(shape1_2*shape2_2)

# rh of phase II data 
rh2 <- df_rh_p2$Humidity3pm

set.seed(1)

# residuals for the phase II data
res2 <- qnorm(pbeta(rh2, p2, q2))

# control limits
control_limits <- cusum(res2, center = mean(residual), 
                        std.dev = sd(residual),
                        decision.interval = 3*1.842887, 
                        se.shift = 1, plot = F)

DI <- control_limits$decision.interval
cusum_pos <- as.vector(control_limits$pos)
cusum_neg <- as.vector(control_limits$neg)*(-1)

out_model1 <- c(which(cusum_neg > DI))
out_model2 <- c(which(cusum_pos > DI))

obs1 <- cusum_neg[out_model1]
obs2 <- cusum_pos[out_model2]

# red dots are observations out of control
plot(cusum_pos, xaxs = "r", type = "o", lwd = 1, pch = 1, 
     ylab = expression("Cumulative Sum"),
     xlab = expression("Observations"), 
     ylim = c(-0.4, 15))
points(cusum_neg, type = "o", xaxs = "r", lwd = 1, pch = 1) 
abline(h = DI, lwd = 1, lty = 1)
points(out_model1, obs1, pch = 16, col = "red")
points(out_model2, obs2, pch = 16, col = "red")