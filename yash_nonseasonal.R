# Following BOX-JENKINS MODELS
# data - https://www.federalreserve.gov/releases/h8/current/default.htm

# Importing the data
library(forecast)
library(xts)
library(ggplot2)
library(rugarch)
library(timeSeries)
library(tseries)
library(timetk)
library(TSA)
library(tsibble)
library(dplyr)
data <- read.csv("/Users/yashadmuthe/Desktop/FRB_H8.1.csv")
data

# Removing NA from the data

data2<- na.omit(data)

# Converting data into TS data and checking dimentions 
datats <- ts(data2$Bank.credit, frequency = 4, start = c(1963, 1))
class(datats)

autoplot(datats)
head(datats)
tail(datats)
dim(data)
str(datats)
cycle(datats) ##This will print the cycle across years.
frequency(datats)
mean(datats)
var(datats)

# ploting Time Series data 
plot.ts(datats)
abline(reg=lm(datats~time(datats)))  # This will fit in a line
plot(aggregate(datats,FUN=mean)) #This will aggregate the cycles and display a year on year trend
boxplot(datats~cycle(datats))

## checking for non-stationarity in data
# 1. Augmented dickey-fuller test

adf.test(datats)
adf.test(datats, k=2)
adf.test(datats, k=1)
# As the p value here is 0.01 which is less than 0.05, we reject null hypothesis 
# and the alternative hypothesis suggests the data is stationary.
# As the data is stationary we do not need to do differencing.


#2.KPSS test

kpss.test(datats)
#In kpss if the p-value is < significance level (say 0.05), then the series is non-stationary
# Here as the p-value is greater than 0.05, the data is stationary.

# checking ACF AND PACF
acf <- acf(as.vector(datats),plot = FALSE, lag.max = 100)
plot(acf)

# In ACF plot we see 11 spikes are significant, out of which 5 are more significant
# Also note that we ignore the spike at lag 0.
# We also check how well the present value of the series is related with its past values. 
# A time series can have components like trend, seasonality, cyclic and residual. 
# ACF considers all these components while finding correlations hence it's a 'complete auto-correlation plot'.
# From ACF plot we can say MA can be 1,2 or 5 

pacf <- pacf(as.vector(datats),lag.max = 100)
plot(pacf)
# In PACF plot we see 2 spikes are significant, out of which 1 are more significant
# From PACF, we can say AR can be 1

eacf(as.vector(datats))
# Using timetk
datats %>%
  tk_tbl(rename_index = "year") %>%
  plot_acf_diagnostics(
    .date_var = year,
    .value    = value,
    .lags     = 20,
    .show_white_noise_bars = TRUE, 
    .interactive = FALSE
  )


#Checking best model fit by comparing their AIC,BIC and log likelihood values.
mod<- arma(datats, order = c(1, 1),lag = NULL, coef = NULL,
           include.intercept = TRUE)
summary(mod)
#Fit:sigma^2 estimated as 9.978,  Conditional Sum-of-Squares = 1456.74,  AIC = 766.46

mod2<- arma(datats, order = c(1, 2))
summary(mod2)
#Fit:sigma^2 estimated as 9.989,  Conditional Sum-of-Squares = 1448.45,  AIC = 768.63

mod3<- arma(datats, order = c(1, 4))
summary(mod3)
#Fit:sigma^2 estimated as 9.9,  Conditional Sum-of-Squares = 1415.66,  AIC = 771.3

mod4<- arma(datats, order = c(1, 5))
summary(mod4)
#Fit:sigma^2 estimated as 9.896,  Conditional Sum-of-Squares = 1405.21,  AIC = 773.24

tsmod <-Arima(datats, order=c(1,0,1))
summary(tsmod)
#Fit : log likelihood = -379.42 AIC=766.84   AICc=767.12   BIC=778.83

tsmod2 <-Arima(datats, order=c(1,0,2))
show(tsmod2)
#Fit : log likelihood = -379.32 AIC=768.63   AICc=769.05   BIC=783.62

tsmod3 <-Arima(datats, order=c(1,0,5))
summary(tsmod3)
#Fit : log likelihood = -377.62 AIC=771.25   AICc=772.28   BIC=795.23

tsmod4 <-Arima(datats, order=c(1,0,20))
show(tsmod4)
#Fit : log likelihood = -367.26 AIC=780.53   AICc=789.43   BIC=849.46

#Comparing all the models we can confidently say that ARIMA (1,0,1) fits very well

#Also checking with Garch model
# install.packages("rugarch", repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com"))


s <- ugarchspec(mean.model=list(armaOrder = c(1,1)),
                variance.model = list(model= 'sGARCH'),
                distribution.model ="norm")


m <- ugarchfit(data = datats, spec = s, solver.control = list(trace=0))

#stargazer(list(m), title="Regression Results", type="text", keep.stat=c("n","ll","aic","bic"), out="kindaresults.doc")
show(m)
#By comparing loglikelihood of ARIMA and Garch model, we conclude ARIMA(1,0,1) is still the best fit.

# Plot of non-seasonal 
o1 <- datats %>%
  as_tsibble() %>%
  autoplot(value)
# Lag plot of non-seasonal 
o2 <- datats %>%
  as_tsibble() %>%
  gg_lag(y=value, geom = "point") 
# Plot it
plot_grid(o1, o2, ncol=1, rel_heights = c(1,2))


#checking residuals

r<- residuals(tsmod)

plot(r,ylab ='Standardized Residuals',type='o'); 
abline(h=0)

# Normality Test on Residuals
qqnorm(r,start=c(1963,1))
qqline(r,start=c(1963,1))

r
r2<- na.omit(r)
acf(as.vector(r2), lag.max = 100)

# Shapiro-wilk test
# The purpose of this test is to see if the data is normally distributed or not
# Null Hypothesis : The data is normally distributed
# Alternative Hypothesis : The data is not normally distributed
shapiro.test(datats)

# From the output, the p-value(0.3504) > 0.05 implying that the distribution of the data are not 
# significantly different from normal distribution. In other words, we can assume the normality.

# Ljung-Box test
# It is a statistical test of whether any group of autocorrelations of a time series 
# are different from 0. Instead of testing randomness of each distinct lag,
# it tests overall randomness based on number of lags.

# H0 : The series is i.i.d
# H1 : The series exhibits serial correlation

Box.test(datats, lag = 10, fitdf = 0, type = 'Lj')

#If p-value < 0.051: You can reject the null hypothesis assuming a 5% chance of making a mistake. So you can assume that your values are showing dependence on each other.
# As here we get a very small p-value, reject H0. The series is not white noise.

ggAcf(datats)
gghistogram(r) + ggtitle("Histogram of residuals")

# Forecasting

ts.plot(datats, xlim=c(1963,2003),main = "Prediction")
fit = datats - r
points(fit, type = 'l', col='red', lty =2)
prediction = predict(tsmod)
prediction$pred[1]

predict(tsmod, n.ahead =12)
forecast <- predict(tsmod, n.ahead = 12)$pred
forecast_se <- predict(tsmod, n.ahead = 12)$se
points(forecast, type = "l", col = 2)
points(forecast - 2*forecast_se, type = "l", col = 2, lty = 2)
points(forecast + 2*forecast_se, type = "l", col = 2, lty = 2)

fcst <- forecast(tsmod, h=6)
autoplot(fcst, include= 60)
print(summary(fcst))

# Forecast for next six Quarters are as follows
# Forecasts:
#   Point Forecast    Lo 80    Hi 80     Lo 95    Hi 95
# 2000 Q1       7.838334 3.775077 11.90159 1.6241192 14.05255
# 2000 Q2       7.875995 3.499985 12.25201 1.1834657 14.56852
# 2000 Q3       7.902917 3.375425 12.43041 0.9787158 14.82712
# 2000 Q4       7.922161 3.319188 12.52514 0.8825206 14.96180
# 2001 Q1       7.935918 3.294848 12.57699 0.8380132 15.03382
# 2001 Q2       7.945752 3.285334 12.60617 0.8182584 15.07325


