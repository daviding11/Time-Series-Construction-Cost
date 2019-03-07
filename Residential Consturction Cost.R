library (astsa )
library (tseries)
library (MASS)
library (forecast)
library (TSA)
library (GeneCycle)
library (ggplot2)
library (qpcR)


#Data exploration 
const <- read.csv("C:/Users/David Ing/Desktop/UCSB Classes/Pstat Class/PSTAT 174/total-construction-spending-residential.csv") 
const_res <- ts(const[,2],frequency = 1) 
plot.ts(const_res, main = "Residential Construction Spending from Jan 2002 - Sep 2018", xlab = "Time in Months", ylab = "Construction Cost") 
seasonplot(const_res,12, col = rainbow(16), year.labels = T, main = "Seasonal Plot of Residential Construction Costs")

#Box Cox 
require(MASS) 
bc <- boxcox(const_res ~ as.numeric(1:length(const_res))) 
lambda <- bc$x[which(bc$y == max(bc$y))] #lambda = 0.4646 
trans_const <- const_res ^ lambda 
trans_const1 <- ts(trans_const [1:(length(const_res) - 10)]) # removing 10 observations 
plot(trans_const1, main = "Transformed Residential Construction Spending       from Jan 2002 - Nov 2017", xlab = "Time in Months", ylab = "Transformed Construction Cost") 
op <- par(mfrow = c(1,2)) 
acf(trans_const1, lag.max = 50, main = "ACF of Transformed Data") 
pacf(trans_const1, lag.max = 50, main = "PACF of Transformed Data") 
par(op)
const.diff <- diff(trans_const1, lag = 12) #differencing for deseasonalization 
plot(const.diff, main=expression(nabla[12] ~ "Transformed Data") , xlab = "Time in Months", ylab = "")
const.diff.nt <- diff(const.diff, differences = 1, lag = 1)# differencing for detrending 
var(const.diff.nt) 
op <- par(mfrow = c(1,2)) 

#Differencing Second Time because lower variance
plot(const.diff.nt, main=expression(nabla ~ nabla[12] ~ "Transformed Data") , xlab = "Time in Months",                 ylab = "") 
plot(diff(const.diff,differences = 2, lag = 1), main=expression(nabla ~ nabla ~ nabla[12] ~ "Transformed Data") , xlab = "Time in Months", ylab = "") 
par(op) 
op <- par(mfrow = c(1,2)) 
acf(diff(const.diff,differences = 2, lag = 1), lag.max = 50, main = expression("ACF of " ~ nabla ~ nabla ~ nabla[12] ~ "Transformed Data")) 
pacf(diff(const.diff,differences = 2, lag = 1), lag.max = 50, main = expression("PACF of " ~ nabla ~ nabla ~ nabla[12] ~ "Transformed Data")) 
par(op) 

#Finding AIC values to select model 
AICc <- numeric() 
for (p in 0:4) {   
  for (q in 0:4) {     
    AICc <- c(AICc , sarima(trans_const2, p, 2, q, 1, 1, 0, 12, details=FALSE)$AICc)   } }
AICc <- matrix(AICc ,nrow=5,byrow=TRUE) 
rownames(AICc) <- c("p=0" , "p=1" , "p=2", "p=3", "p=4") 
colnames (AICc) <- c("q=0" , "q=1" , "q=2", "q=3", "q=4" ) 
AICc <- data.frame(AICc) 
AICc #gives dataframe of various p and q values to test  

#models with good AIC values and fit law of parsimony 
fit2 <- arima(trans_const1, order = c(2,1,0), seasonal = list(order = c(1,1,0), period = 12), method = "ML") 
fit3 <- arima(trans_const1, order = c(4,1,3), seasonal = list(order = c(1,1,0), period = 12), method = "ML") 
fit4 <- arima(trans_const1, order = c(2,1,1), seasonal = list(order = c(1,1,0), period = 12), method = "ML") fit2 
fit3 
fit4 

#Fit2 stationary/invertibility plots 
op <- par(mfrow = c(2,3)) 
plot.roots(NULL, polyroot(c(1,0.4476,0.1723)), main = "roots for ar part") 
plot.roots(NULL, polyroot(c(1,-0.1600)), main = "roots for sar part") 
par(op) 

#Fit3 stationarity/invertibility plots 
op <- par(mfrow = c(2,3)) 
plot.roots(NULL, polyroot(c(1,0.3296,-0.0689,-0.6734,0.5177)), main = "roots for ar part") 
plot.roots(NULL, polyroot(c(1,0.1278,0.3141,0.8728)), main ="roots for ma part") 
plot.roots(NULL, polyroot(c(1,-0.1722)), main = "roots for sar part") 
par(op) 

#Fit4 stationary/invertibility plots 
op <- par(mfrow = c(2,3)) 
plot.roots(NULL, polyroot(c(1,-0.3553,0.5599)), main = "roots for ar part") 
plot.roots(NULL, polyroot(c(1,0.8748)), main ="roots for ma part") 
plot.roots(NULL, polyroot(c(1,-0.1921)), main = "roots for sar part") 
par(op) 

#diagnostic checking using residuals
res2 <- residuals(fit2) 
res4 <- residuals(fit4)

#Normality 
op <- par(mfrow = c(1,2)) 
hist(res2, main = "Histogram of Residuals of Model 2", breaks = 15) 
qqnorm(res2, main = "QQ plot of Model 2") 
qqline(res2) 
par(op) 
op <- par(mfrow = c(1,2)) 
hist(res4, main = "Histogram of Residuals of Model 4", breaks = 15) 
qqnorm(res4, main = "QQ plot of Model 4") 
qqline(res4) 
par(op)
#Shapiro test for normality s
hapiro.test(res2) 
shapiro.test(res4) 

#Box pierce test for serial correlation 
box2 <- Box.test(res2, lag = 12, type = "Box-Pierce", fitdf = 2) 
box4 <- Box.test(res2, lag = 12, type = "Box-Pierce", fitdf = 3) 
box2 
box4 

#ACF PACF plots of residuals for  model 2 and 4 
op <- par(mfrow = c(1,2)) 
acf(res2, main = "ACF of Model 2") 
pacf(res2, main = "PACF of Model 2") 
par(op) 
op <- par(mfrow = c(1,2)) 
acf(res4, main ="ACF of Model 4") 
pacf(res4, main = "PACF of Model 4") par(op) 

#forecasting- transformed data 
pred.tr <- predict(fit2, n.ahead = 10) 
upper <- pred.tr$pred + 2*pred.tr$se 
lower <- pred.tr$pred - 2*pred.tr$se 
ts.plot(trans_const1, xlim = c(1, length(trans_const1)+10), main ="Forecasting Transformed Data", ylab = "Transformed Construction Cost") 
lines(upper, col = "BLUE", lty = "dashed") 
lines(lower, col = "BLUE", lty = "dashed") 
points((length(trans_const1) +1):(length(trans_const1)+10), pred.tr$pred, col = "red", pch = 16, cex = 0.5) 

#forecasting- original data 
pred.origin <- pred.tr$pred ^ (1/lambda) 
uorig <- upper ^ (1/lambda) 
lorig <- lower ^ (1/lambda) 
const_res2 <- const_res 
ts.plot(const_res2, xlim = c(1, length(const_res2)), main ="Forecasting Original Data", ylab = "Construction Spending") 
lines(uorig, col = "BLUE", lty = "dashed") 
lines(lorig, col = "BLUE", lty = "dashed") 
points((length(trans_const1) +1):(length(trans_const1)+10), pred.origin, col = "red", pch = 16, cex = 0.5) 

#comparison of original and forecasted points 
ts.plot((const_res2), xlim = c(length(const_res2) - 20, length(const_res2)), main = "Comparison of Forecasted Values and Original Values",         
        ylab = "Construction Costs") 
points(192:201, const_res2[192:201], col = "GREEN") 
points(192:201, pred.origin, col = "RED") 
points(192:201, uorig, lty = 1, col = "BLUE") 
points(192:201, lorig, lty = 1, col = "BLUE")


