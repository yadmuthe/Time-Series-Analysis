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
data <- read.csv("/Users/yashadmuthe/Desktop/Delhi_temp.csv")
data


# Converting data into TS data and checking dimentions 
datats <- ts(data$Delhi_Temperature, frequency = 12, start = c(2001, 1))
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
# Here as the p-value(i.e 0.1) is greater than 0.05, the data is stationary.

# checking ACF AND PACF
acf <- acf(as.vector(datats),plot = FALSE, lag.max = 200)
plot(acf)

pacf <- pacf(as.vector(datats),lag.max = 100)
plot(pacf)

eacf(as.vector(datats))

# Using timetk to plot ACF and PACF

datats %>%
  tk_tbl(rename_index = "year") %>%
  plot_acf_diagnostics(
    .date_var = year,
    .value    = value,
    .lags     = 20,
    .show_white_noise_bars = TRUE, 
    .interactive = FALSE
  )

#Checking for models from ACF and PACF plots.
#arima-(2/3/8,0,9/10) (0,0,9/10)s=12


#Checking best model fit
mod1 <- arima(datats,order=c(2,0,9), seasonal = list(order = c(0,0,9), period = 12),method='CSS')
summary(mod1)
# sigma^2 estimated as 2.221:  part log likelihood = -261.79

mod2<- arima(datats,order=c(8,0,9), seasonal = list(order = c(0,0,9), period = 12),method='CSS')
summary(mod2)
# sigma^2 estimated as 0.8826:  part log likelihood = -195.34

mod3<- arima(datats,order=c(3,0,9), seasonal = list(order = c(0,0,9), period = 12),method='CSS')
summary(mod3)
# sigma^2 estimated as 1.843:  part log likelihood = -248.35

mod4<- arima(datats,order=c(2,0,9), seasonal = list(order = c(0,0,10), period = 12),method='CSS')
summary(mod4)
#sigma^2 estimated as 2.153:  part log likelihood = -259.54

tsmod <-arima(datats,order=c(8,0,9), seasonal = list(order = c(0,0,10), period = 12),method='CSS')
show(tsmod)
# sigma^2 estimated as 0.8902:  part log likelihood = -195.95

tsmod2 <-arima(datats,order=c(3,0,9), seasonal = list(order = c(0,0,10), period = 12),method='CSS')
show(tsmod2)
# sigma^2 estimated as 1.519:  part log likelihood = -234.41

tsmod3 <-arima(datats,order=c(2,0,10), seasonal = list(order = c(0,0,9), period = 12),method='CSS')
summary(tsmod3)
# sigma^2 estimated as 2.285:  part log likelihood = -263.83

tsmod4 <-arima(datats,order=c(8,0,10), seasonal = list(order = c(0,0,9), period = 12),method='CSS')
show(tsmod4)
# sigma^2 estimated as 1.509:  part log likelihood = -233.97

tsmod5 <-arima(datats,order=c(3,0,10), seasonal = list(order = c(0,0,9), period = 12),method='CSS')
show(tsmod5)
# sigma^2 estimated as 2.132:  part log likelihood = -258.85
# Comparing all the models we can confidently say that ARIMA(8,0,9)(0,0,9)s=12 fits very well


# install.packages("rugarch", repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com"))


s <- ugarchspec(mean.model=list(armaOrder = c(8,9)),
                variance.model = list(model= 'sGARCH'),
                distribution.model ="norm")

m <- ugarchfit(data = datats, spec = s, solver.control = list(trace=0))
show(m)
# Checked fitting different Garch models(8,1)(8,9)(2,2)(1,9) and after comparing with ARIMA
# found still the best model to be ARIMA(8,0,9)(0,0,9)s=12.

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
qqnorm(r,start=c(2001,1))
qqline(r,start=c(2001,1))

r
r2<- na.omit(r)
acf(as.vector(r2), lag.max = 100)

# Shapiro-wilk test
# The purpose of this test is to see if the data is normally distributed or not
# Null Hypothesis : The data is normally distributed
# Alternative Hyposthesis : The data is not normally distributed

shapiro.test(datats)
# From the output, the p-value(2.816e-07) > 0.05 implying that the distribution of the data are  
# significantly different from normal distribution. In other words, we can assume the non-normality.

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


ts.plot(datats, xlim=c(2000,2015),main = "Prediction")
fit = datats - r
points(fit, type = 'l', col='red', lty =2)
prediction = predict(mod2)
prediction$pred[1]

predict(mod2, n.ahead =12)
forecast <- predict(mod2, n.ahead = 12)$pred
forecast_se <- predict(mod2, n.ahead = 12)$se
points(forecast, type = "l", col = 2)
points(forecast - 2*forecast_se, type = "l", col = 2, lty = 2)
points(forecast + 2*forecast_se, type = "l", col = 2, lty = 2)





