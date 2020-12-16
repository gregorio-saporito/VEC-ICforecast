# load the libraries
library(readr)
library(devtools)
library(tidyverse)
library(zoo)
library(forecast)

# check the countries available in the dataset
countries <- unique(read_csv("data.csv")$country)

# choose Denmark and clean dataset
data <- read_csv("data.csv") %>%
  mutate(date = as.Date(date,format="%d/%m/%Y")) %>%
  dplyr::select(-year_week,-source,-url) %>%
  filter(indicator %in% c("Daily hospital occupancy","Daily ICU occupancy")) %>%
  filter(country == "Denmark") %>%
  spread(indicator,value)

# check min and max dates and missing dates
min(data$date)
max(data$date)

# generate the full daily dates which we would have without missing values
full_dates <- seq(from=min(data$date), to=max(data$date),
    by="day")

# there are 4 missing dates
setdiff(full_dates,data$date)

## Make a data frame with a full series of dates from the min date to the max date
## in the incomplete data frame
full_dates <- data.frame(date = full_dates)

## Merge the complete data frame with the incomplete to fill in the dates and add 
## NAs for missing values
my_complete_data <- merge(full_dates, data, by = "date", 
                          all.x = TRUE)

# we need to interpolate the missing  values for hospitalised
s_hosp <- ts(my_complete_data$`Daily hospital occupancy`)

for(i in c(1:length(s_hosp))){
  if(is.na(s_hosp[i])){
    a <- ts(s_hosp[1:(i-1)])
    b <- forecast(auto.arima(a))
    s_hosp[i] <- round(b$mean[1])
  }else{
    print("skip")
  }
}
  
# we need to interpolate missing values for IC
s_IC <- ts(my_complete_data$`Daily ICU occupancy`)

for(i in c(1:length(s_IC))){
  if(is.na(s_IC[i])){
    a <- ts(s_IC[1:(i-1)])
    b <- forecast(auto.arima(a))
    s_IC[i] <- round(b$mean[1])
  }else{
    print("skip")
  }
}

# integrate the interpolated missing values
data_clean <- tibble(my_complete_data) %>%
  mutate(fit_hosp=s_hosp,
         fit_IC=s_IC) %>%
  mutate(
    `Daily hospital occupancy` = ifelse(is.na(`Daily hospital occupancy`),fit_hosp,`Daily hospital occupancy`),
    `Daily ICU occupancy` = ifelse(is.na(`Daily ICU occupancy`),fit_IC,`Daily ICU occupancy`)
  ) %>%
  dplyr::select(-fit_hosp,-fit_IC) %>%
  gather(key="indicator",value="value",-date,-country) %>%
  # putting an arbitrarily small value to be able to take logs later
  mutate(value=ifelse(value<=0,0.001,value))

dates <- filter(data_clean,indicator=="Daily hospital occupancy")$date

HwS <- filter(data_clean,indicator=="Daily hospital occupancy")$value
IC <- filter(data_clean,indicator=="Daily ICU occupancy")$value

HwS <- zoo(log(HwS), order.by = dates)
IC <- zoo(log(IC), order.by = dates)

# being too close to zero the logarithm still returns an extreme negative spike for a value
# it is better to interpolate it
IC_fix <- ts(IC)

for(i in c(1:length(IC_fix))){
  if(IC_fix[i]<0){
    a <- ts(IC_fix[1:(i-1)])
    b <- forecast(auto.arima(a))
    IC_fix[i] <- round(b$mean[1])
  }else{
    print("skip")
  }
}

IC <- zoo(IC_fix, order.by = dates)

# (a) Create a single time series plot with two log prices
par(mfrow=c(1,1))
plot(HwS, type='l', main="WTI vs. Gas", col="blue", 
     ylim = c(0, 7)
     )
lines(IC, col="red")


#Q.b : ERS unit root test
library(urca)

# The results below show that we can reject the null of no unit root.
# test statistic has to be larger than the critical value
HwS.urers1 <- ur.ers(HwS, type="P-test")
summary(HwS.urers1)
IC.urers2 <- ur.ers(IC, type="P-test")
summary(IC.urers2)

# reject Null of no unit root -> exist unit root -> difference
dHwS <- diff(HwS)
dIC <- diff(IC)

dHwS.urers1 <- ur.ers(dHwS, type="P-test")
summary(dHwS.urers1)
dIC.urers2 <- ur.ers(dIC, type="P-test")
summary(dIC.urers2)

# Based on this result, we can say that both time series are I(1).
par(mfrow=c(1,1))
plot(dHwS, type='l', main="(log difference)WTI vs. (log difference)Gas", col="blue")
lines(dIC, col="red")

# ols test unit root of residuals
ols <- lm(IC ~ HwS)
plot(ols$residuals)

# no unit root
# there is a linear combination that is I(0) + previous evidence the series are cointegrated
ur_test <- ur.ers(ols$residuals, type="P-test")
s_urtest <- summary(ur_test)
s_urtest@teststat
s_urtest@cval[2]


#Q.c
# (c) Determine the number of lags to include in cointegration analysis.
# Run the Johansen’s trace and maximum eigenvalue cointegration tests
# The results of Johansen’s Trace test show that we can reject the null hypothesis. 
# Therefore we can conclude that there is cointegration.
library(vars)
y <- cbind(HwS, IC)
colnames(y) <- c("HwS","IC")
y <- na.trim(y)

# determine number of lags to be included in cointegration test and in VEC model
y.VAR.IC <- VARselect(y, type="const")
nlags <- y.VAR.IC$selection["AIC(n)"]
nlags

y <- window(y, start=min(dates), end=(max(dates)-7))

# perform cointegration test
library(urca)
# r=0 tests for the presence of cointegration, if the test statistic exceeds
# the %significance levle reject the null of no cointegration
y.CA <- ca.jo(y, ecdet="const", type="eigen", K=nlags, spec="transitory")
summary(y.CA)

# Conducts a likelihood ratio test for no inclusion of a linear trend in a VAR.
lttest(y.CA, r=1)

# Q.d
# estimate unrestriced VEC model
y.VEC <- cajorls(y.CA, r=1)
y.VEC

#Q.E
# to see t-statistics and p-values
# (e) Are the adjustment parameters α1 and α2 in the estimated VEC model
# statistically significant?
library(tsDyn)
coefA(y.VEC)
# For long run relationship to be stable, we need α1 < 0 or =0, α2 > 0 
# or =0 and at least one of them can not be equal 0.
summary(y.VEC$rlm)

# (f) Reestimate the VEC model with a restriction α1=0.
# test for restricted adjustment parameters alpha
rest.alpha <- matrix(c(0,1), c(2,1))
# estimate the restricted VAR
y.CA.ralpha <- alrtest(y.CA, A=rest.alpha, r=1)
summary(y.CA.ralpha)

# one month ahead forecast
# Create and plot sequence of one month ahead forecasts for the period

yall <- cbind(HwS, IC)
colnames(yall) <- c("HwS","IC")
yall <- na.trim(yall)

yall <- window(yall, start=min(dates), end=max(dates))
first.m <- min(dates)
last.m <- (max(dates)-7)

y.p1 <- window(yall, end=last.m)
y.p2 <- window(yall, start=last.m+1)

new <- as.character(seq(from=last.m, to=max(dates)-1, by="day"))

y.VAR.f1 <-data.frame()
y.VAR.f2 <-data.frame()

for(i in new){
  
  y <- window(yall, end = i)
  y.CA <- ca.jo(y, ecdet="const", type="eigen", K=nlags, spec="transitory")
  y.VAR <- vec2var(y.CA, r=1)
  y.VAR.updt <- predict(y.VAR, n.ahead=1)
  y.VAR.f1 <-rbind(y.VAR.f1, as.ts(y.VAR.updt$fcst$HwS))
  y.VAR.f2 <-rbind(y.VAR.f2, as.ts(y.VAR.updt$fcst$IC))
  
}


# forecast
y.VAR.f2

# actual
y.p2[,2]

# old
y.p1[,2]

old = tibble(
  date = dates[1:(length(dates)-7)],
  value = as.vector(y.p1[,2])
)

forecast = tibble(
  date = dates[(length(dates)-6):length(dates)],
  value = y.VAR.f2$fcst
)

actual = tibble(
  date = dates[(length(dates)-6):length(dates)],
  value = as.vector(y.p2[,2])
)

ggplot() +
  geom_line(data = old,aes(x=date,y=value),col='black') +
  geom_point(data = old,aes(x=date,y=value),col='blue') +
  geom_line(data = actual,aes(x=date,y=value),col='black') +
  geom_point(data=forecast,aes(x=date,y=value),col='red') +
  theme_classic()

# try to compare it with the performance of a random walk

#Q.I
dIC.p1 <- window(dIC, end=last.m)
dIC.p2 <- window(dIC, start=last.m+1)

par(mfrow=c(1,2))
Acf(dIC.p1, type="correlation", lag=48, main="ACF for IC")
Acf(dIC.p1, type="partial", lag=48, main="PACF for IC")

# apply random walk
arma01 <- Arima(dIC.p1, order=c(0, 1, 0))
tsdiag(arma01, gof.lag=36)

#Q.J
rol.f <- c()
for(i in new){
  y <- window(dIC, end = i)
  rol.updt <- arima(y, order=c(0,1,0))
  rol.f <- c(rol.f, forecast(rol.updt, 1)$mean)
}

dates_rw <- index(dIC)

old_rw = tibble(
  date = dates_rw[1:(length(dates_rw)-7)],
  value = as.vector(dIC[1:(length(dates_rw)-7)])
)

forecast_rw = tibble(
  date = dates_rw[(length(dates_rw)-6):length(dates_rw)],
  value = rol.f[1:(length(rol.f))]
)

actual_rw = tibble(
  date = dates_rw[(length(dates_rw)-6):length(dates_rw)],
  value = as.vector(dIC[(length(dates_rw)-6):length(dates_rw)])
)

ggplot() +
  geom_line(data = old_rw,aes(x=date,y=value),col='black') +
  geom_point(data = old_rw,aes(x=date,y=value),col='blue') +
  geom_line(data = actual_rw,aes(x=date,y=value),col='black') +
  geom_point(data=forecast_rw,aes(x=date,y=value),col='red') +
  theme_classic()


#Q.k
# (k) Compare the RMSE of the forecast based on the VEC model
# and the forecast based on the random walk

# look at MAE and MAPE
# VEC model
err_VEC <- accuracy(as.vector(y.p2[,2]), y.VAR.f2$fcst)
err_VEC

# from the random walk
err_rw <- accuracy(as.vector(dIC[(length(dates_rw)-6):length(dates_rw)]), rol.f[1:(length(rol.f))])
err_rw

