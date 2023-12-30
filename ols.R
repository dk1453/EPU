library(vars)
library(readxl)
library(ggplot2)
epu <- read_excel("~/Desktop/typ/data/data.xlsx")
epu <- epu[, 2:25]
epu <- as.matrix(epu)


### generate the ols coefficients ###
model_ols <- VAR(epu, ic = 'AIC')
coef_ols <- Bcoef(model_ols)[,-25] # the last column is the const

### connectedness ###
res_ols <- resid(model_ols) #each country a column
mean(res_ols ^ 2) #mse of the static analysis
Q_ols <- (t(res_ols)%*%res_ols)/nrow(res_ols)
dim(Q_ols) # dim of Q has to be consistant with the dim of coef(withou const)

cnct_ols <- GFEVD(coef_ols, Q_ols, n.ahead=6, standardize = T)$GFEVD * 100 # std or not, no difference, * 100 to let row sum to be 100
View(cnct_ols)
rowSums(cnct_ols)
colSums(cnct_ols)
## static analysis ##
cnct_ols_st <- cnct_analysis(cnct_ols)
cnct_ols <- cbind(cnct_ols, cnct_ols_st$From)
cnct_ols <- rbind(cnct_ols, append(cnct_ols_st$To, 0), append(cnct_ols_st$Net, 0))
# write.csv(cnct_ols_st$To,"desktop/cnct_ols_TO.csv")
# write.csv(cnct_ols_st$From,"desktop/cnct_ols_FROM.csv")
# write.csv(cnct_ols,"desktop/typ/cnct_ols.csv")
cnct_ols_st$Sum

# quadrants
cnct_ols_quadrants <- read_excel("~/Desktop/cnct_ols_quadrants.xlsx")
png("desktop/var(2)/cnct_ols_quadrants_var(1).png", width = 3570, height = 1865, res = 350)
ggplot(cnct_ols_quadrants, aes(y=From, x=To, label=country)) +
  geom_text_repel(aes(y=From, x=To, label=country)) +
  geom_point() +
  scale_y_continuous(breaks=seq(0,100,10),limits=c(0,100)) + 
  scale_x_continuous(breaks=seq(0,200,10), limits=c(0,200)) +
  # #lims(y=c(70,110),x=c(-50,250)) +
  #coord_fixed() +  # 
  geom_vline(xintercept = mean(cnct_ols_quadrants$To)) + geom_hline(yintercept = mean(cnct_ols_quadrants$From)) + 
  geom_abline(intercept = 0)
dev.off()

## dynamic analysis ##


### generate the rolling window aic
VARselect(epu, lag.max = 5)
test_criteria <- VARselect(epu, lag.max = 1)
aic1 <- test_criteria$criteria
aic_rolling <- matrix(NA, 173, 1)
for (t in 1:173){
  criteria_ <- VARselect(epu[(0:49)+t,], lag.max = 1)
  aic_tem <- criteria_$criteria[1]
  aic_rolling[t,] <-aic_tem
}
View(aic_rolling)


### generate the rolling window cnectedness tables
cnct_ols_rw <- array(NA, dim = c(24, 24, 173))
res_ols_dy_collection <- array(NA, dim = c(49, 24, 173))
cnct_ols_fromrw <- array(NA, dim =c(24,1,173))
cnct_ols_torw <- array(NA, dim =c(24,1,173))
for (t in 1:173)
{
  model_ols_dy <- VAR(epu[(0:49)+t,], ic = 'AIC')
  coef_ols_dy <- Bcoef(model_ols_dy)[,-25]
  res_ols_dy <- resid(model_ols_dy) 
  Q_ols_dy <- (t(res_ols_dy)%*%res_ols_dy)/nrow(res_ols_dy)
  cnct_ols_mid <- GFEVD(coef_ols_dy, Q_ols_dy, n.ahead=6, standardize = T)$GFEVD * 100
  res_ols_dy_collection[,,t] <- as.matrix(res_ols_dy)
  cnct_ols_rw[,,t] <- as.matrix(cnct_ols_mid)
  
  cnct_ols_mid_dy <- cnct_analysis(cnct_ols_mid)
  cnct_ols_fromrw[,,t] <- as.matrix(cnct_ols_mid_dy$From)
  cnct_ols_torw[,,t] <- as.matrix(cnct_ols_mid_dy$To)
}
# var(cnct_ols_fromrw)
# var(cnct_ols_torw)
cnct_ols_rw[,,173]
mean(apply(res_ols_dy_collection, 3 , function(x) mean(x ^ 2))) #mse of the rolling window
### total, from, to, net spillover
cnct_ols_dy <- apply(cnct_ols_rw, 3, cnct_analysis)
##total
sum_spillover_ols <- matrix(NA, 173, 1)
for (j in 1:173){
  sum_spillover1 <- cnct_ols_dy[[j]]$Sum
  sum_spillover_ols[j,1] = as.numeric(sum_spillover1)
} 
sum_spillover_ols <- ts(sum_spillover_ols, frequency = 12, start = c(2007, 2))

##TO
##USA to
to_spillover_ols_usa <- matrix(NA, 173, 1)
for (j in 1:173){
  to_spillover1 <- cnct_ols_dy[[j]]$To[1]
  to_spillover_ols_usa[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_usa <- ts(to_spillover_ols_usa, frequency = 12, start = c(2007, 2))
##China to
to_spillover_ols_China <- matrix(NA, 173, 1)
for (j in 1:173){
  to_spillover1 <- cnct_ols_dy[[j]]$To[7]
  to_spillover_ols_China[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_China <- ts(to_spillover_ols_China, frequency = 12, start = c(2007, 2))
##Germany to
to_spillover_ols_Germany <- matrix(NA, 173, 1)
for (j in 1:173){
  to_spillover1 <- cnct_ols_dy[[j]]$To[12]
  to_spillover_ols_Germany[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_Germany <- ts(to_spillover_ols_Germany, frequency = 12, start = c(2007, 2))
##Greece to 
to_spillover_ols_Greece <- matrix(NA, 173, 1)
for (j in 1:173){
  to_spillover1 <- cnct_ols_dy[[j]]$To[13]
  to_spillover_ols_Greece[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_Greece <- ts(to_spillover_ols_Greece, frequency = 12, start = c(2007, 2))
##SK to 
to_spillover_ols_SK <- matrix(NA, 173, 1)
for (j in 1:173){
  to_spillover1 <- cnct_ols_dy[[j]]$To[19]
  to_spillover_ols_SK[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_SK <- ts(to_spillover_ols_SK, frequency = 12, start = c(2007, 2))
##Mexico to
to_spillover_ols_mx <- matrix(NA, 173, 1)
for (j in 1:173){
  to_spillover1 <- cnct_ols_dy[[j]]$To[20]
  to_spillover_ols_mx[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_mx <- ts(to_spillover_ols_mx, frequency = 12, start = c(2007, 2))
##Russia to
to_spillover_ols_Russia <- matrix(NA, 173, 1)
for (j in 1:173){
  to_spillover1 <- cnct_ols_dy[[j]]$To[21]
  to_spillover_ols_Russia[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_Russia <- ts(to_spillover_ols_Russia, frequency = 12, start = c(2007, 2))
##UK to
to_spillover_ols_UK <- matrix(NA, 173, 1)
for (j in 1:173){
  to_spillover1 <- cnct_ols_dy[[j]]$To[24]
  to_spillover_ols_UK[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_UK <- ts(to_spillover_ols_UK, frequency = 12, start = c(2007, 2))


##From
##USA from
from_spillover_ols_usa <- matrix(NA, 173, 1)
for (j in 1:173){
  from_spillover1 <- cnct_ols_dy[[j]]$From[1]
  from_spillover_ols_usa[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_usa <- ts(from_spillover_ols_usa, frequency = 12, start = c(2007, 2))
##China from
from_spillover_ols_China <- matrix(NA, 173, 1)
for (j in 1:173){
  from_spillover1 <- cnct_ols_dy[[j]]$From[7]
  from_spillover_ols_China[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_China <- ts(from_spillover_ols_China, frequency = 12, start = c(2007, 2))
##Germany from
from_spillover_ols_Germany <- matrix(NA, 173, 1)
for (j in 1:173){
  from_spillover1 <- cnct_ols_dy[[j]]$From[12]
  from_spillover_ols_Germany[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_Germany <- ts(from_spillover_ols_Germany, frequency = 12, start = c(2007, 2))
##Greece from
from_spillover_ols_Greece <- matrix(NA, 173, 1)
for (j in 1:173){
  from_spillover1 <- cnct_ols_dy[[j]]$From[13]
  from_spillover_ols_Greece[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_Greece <- ts(from_spillover_ols_Greece, frequency = 12, start = c(2007, 2))
##SK from
from_spillover_ols_SK <- matrix(NA, 173, 1)
for (j in 1:173){
  from_spillover1 <- cnct_ols_dy[[j]]$From[19]
  from_spillover_ols_SK[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_SK <- ts(from_spillover_ols_SK, frequency = 12, start = c(2007, 2))
##Mexico from
from_spillover_ols_mx <- matrix(NA, 173, 1)
for (j in 1:173){
  from_spillover1 <- cnct_ols_dy[[j]]$From[20]
  from_spillover_ols_mx[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_mx <- ts(from_spillover_ols_mx, frequency = 12, start = c(2007, 2))
##Russia from
from_spillover_ols_Russia <- matrix(NA, 173, 1)
for (j in 1:173){
  from_spillover1 <- cnct_ols_dy[[j]]$From[21]
  from_spillover_ols_Russia[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_Russia <- ts(from_spillover_ols_Russia, frequency = 12, start = c(2007, 2))
##UK from
from_spillover_ols_UK <- matrix(NA, 173, 1)
for (j in 1:173){
  from_spillover1 <- cnct_ols_dy[[j]]$From[24]
  from_spillover_ols_UK[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_UK <- ts(from_spillover_ols_UK, frequency = 12, start = c(2007, 2))


##Net
zero <- ts(matrix(1,200) - 1, frequency = 12, start = c(2006, 1))
##USA net
net_spillover_ols_usa <- matrix(NA, 173, 1)
for (j in 1:173){
  net_spillover1 <- cnct_ols_dy[[j]]$Net[1]
  net_spillover_ols_usa[j,1] = as.numeric(net_spillover1)
} 
net_spillover_ols_usa <- ts(net_spillover_ols_usa, frequency = 12, start = c(2007, 2))
##China net
net_spillover_ols_China <- matrix(NA, 173, 1)
for (j in 1:173){
  net_spillover1 <- cnct_ols_dy[[j]]$Net[7]
  net_spillover_ols_China[j,1] = as.numeric(net_spillover1)
} 
net_spillover_ols_China <- ts(net_spillover_ols_China, frequency = 12, start = c(2007, 2))
##Germany net
net_spillover_ols_Germany <- matrix(NA, 173, 1)
for (j in 1:173){
  net_spillover1 <- cnct_ols_dy[[j]]$Net[12]
  net_spillover_ols_Germany[j,1] = as.numeric(net_spillover1)
} 
net_spillover_ols_Germany <- ts(net_spillover_ols_Germany, frequency = 12, start = c(2007, 2))
##Greece net
net_spillover_ols_Greece <- matrix(NA, 173, 1)
for (j in 1:173){
  net_spillover1 <- cnct_ols_dy[[j]]$Net[13]
  net_spillover_ols_Greece[j,1] = as.numeric(net_spillover1)
} 
net_spillover_ols_Greece <- ts(net_spillover_ols_Greece, frequency = 12, start = c(2007, 2))
##SK net
net_spillover_ols_SK <- matrix(NA, 173, 1)
for (j in 1:173){
  net_spillover1 <- cnct_ols_dy[[j]]$Net[19]
  net_spillover_ols_SK[j,1] = as.numeric(net_spillover1)
} 
net_spillover_ols_SK <- ts(net_spillover_ols_SK, frequency = 12, start = c(2007, 2))
##Mexico net
net_spillover_ols_mx <- matrix(NA, 173, 1)
for (j in 1:173){
  net_spillover1 <- cnct_ols_dy[[j]]$Net[20]
  net_spillover_ols_mx[j,1] = as.numeric(net_spillover1)
} 
net_spillover_ols_mx <- ts(net_spillover_ols_mx, frequency = 12, start = c(2007, 2))
##Russia net
net_spillover_ols_Russia <- matrix(NA, 173, 1)
for (j in 1:173){
  net_spillover1 <- cnct_ols_dy[[j]]$Net[21]
  net_spillover_ols_Russia[j,1] = as.numeric(net_spillover1)
} 
net_spillover_ols_Russia <- ts(net_spillover_ols_Russia, frequency = 12, start = c(2007, 2))
##UK net
net_spillover_ols_UK <- matrix(NA, 173, 1)
for (j in 1:173){
  net_spillover1 <- cnct_ols_dy[[j]]$Net[24]
  net_spillover_ols_UK[j,1] = as.numeric(net_spillover1)
} 
net_spillover_ols_UK <- ts(net_spillover_ols_UK, frequency = 12, start = c(2007, 2))


### pairwise spillover
##us-cn, us-mexico
cnct_ols_US_CN <- apply(cnct_ols_rw, 3, cnct_pairwise_analysis, i = 1, j = 7)
cnct_ols_US_CN <- ts(cnct_ols_US_CN, frequency = 12, start = c(2007, 2))
cnct_ols_US_MX <- apply(cnct_ols_rw, 3, cnct_pairwise_analysis, i = 1, j = 20)
cnct_ols_US_MX <- ts(cnct_ols_US_MX, frequency = 12, start = c(2007, 2))

##us-russia, cn-hk
cnct_ols_US_RU <- apply(cnct_ols_rw, 3, cnct_pairwise_analysis, i = 1, j = 21)
cnct_ols_US_RU <- ts(cnct_ols_US_RU, frequency = 12, start = c(2007, 2))
cnct_ols_CN_HK <- apply(cnct_ols_rw, 3, cnct_pairwise_analysis, i = 7, j = 14)
cnct_ols_CN_HK <- ts(cnct_ols_CN_HK, frequency = 12, start = c(2007, 2))

##Germany-Greece, cn-sk
cnct_ols_GR_GE <- apply(cnct_ols_rw, 3, cnct_pairwise_analysis, i = 12, j = 13)
cnct_ols_GR_GE <- ts(cnct_ols_GR_GE, frequency = 12, start = c(2007, 2))
cnct_ols_CN_SK <- apply(cnct_ols_rw, 3, cnct_pairwise_analysis, i = 7, j = 19)
cnct_ols_CN_SK <- ts(cnct_ols_CN_SK, frequency = 12, start = c(2007, 2))








###plot####################################
##total##
png("desktop/typ/overall_spillovers_ols.png", width = 1200, height = 800, res = 200)
plot(sum_spillover_ols, type = "l", las = 2, ylab = 'Total Spillover', main = 'OLS')
axis(1, at = c(2007:2022), las = 2)
dev.off()
##representative countries##
## TO
png("desktop/typ/to_ols.png", width = 2824, height = 3010, res = 300)
par(mfrow=c(4,2))
plot(to_spillover_ols_usa, type = "l", las = 2,  ylab = 'To Spillover', main = 'USA')
axis(1, at = c(2007:2022), las = 2)
plot(to_spillover_ols_China, type = "l", las = 2,  ylab = 'To Spillover', main = 'China')
axis(1, at = c(2007:2022), las = 2)
plot(to_spillover_ols_Germany, type = "l", las = 2,  ylab = 'To Spillover', main = 'Germany')
axis(1, at = c(2007:2022), las = 2)
plot(to_spillover_ols_Greece, type = "l", las = 2,  ylab = 'To Spillover', main = 'Greece')
axis(1, at = c(2007:2022), las = 2)
plot(to_spillover_ols_SK, type = "l", las = 2,  ylab = 'To Spillover', main = 'South Korea')
axis(1, at = c(2007:2022), las = 2)
plot(to_spillover_ols_mx, type = "l", las = 2, ylab = 'To Spillover', main = 'Mexico')
axis(1, at = c(2007:2022), las = 2)
plot(to_spillover_ols_Russia, type = "l", las = 2,  ylab = 'To Spillover', main = 'Russia')
axis(1, at = c(2007:2022), las = 2)
plot(to_spillover_ols_UK, type = "l", las = 2,  ylab = 'To Spillover', main = 'UK')
axis(1, at = c(2007:2022), las = 2)
dev.off()
## FROM
png("desktop/typ/from_ols.png", width = 2824, height = 3010, res = 300)
par(mfrow=c(4,2))
plot(from_spillover_ols_usa, type = "l", las = 2,  ylab = 'From Spillover', main = 'USA')
axis(1, at = c(2007:2022), las = 2)
plot(from_spillover_ols_China, type = "l", las = 2,  ylab = 'From Spillover', main = 'China')
axis(1, at = c(2007:2022), las = 2)
plot(from_spillover_ols_Germany, type = "l", las = 2,  ylab = 'From Spillover', main = 'Germany')
axis(1, at = c(2007:2022), las = 2)
plot(from_spillover_ols_Greece, type = "l", las = 2,  ylab = 'From Spillover', main = 'Greece')
axis(1, at = c(2007:2022), las = 2)
plot(from_spillover_ols_SK, type = "l", las = 2,  ylab = 'From Spillover', main = 'South Korea')
axis(1, at = c(2007:2022), las = 2)
plot(from_spillover_ols_mx, type = "l", las = 2, ylab = 'From Spillover', main = 'Mexico')
axis(1, at = c(2007:2022), las = 2)
plot(from_spillover_ols_Russia, type = "l", las = 2,  ylab = 'From Spillover', main = 'Russia')
axis(1, at = c(2007:2022), las = 2)
plot(from_spillover_ols_UK, type = "l", las = 2,  ylab = 'From Spillover', main = 'UK')
axis(1, at = c(2007:2022), las = 2)
dev.off()
##NET
png("desktop/typ/net_ols.png", width = 2824, height = 3010, res = 300)
par(mfrow=c(4,2))
plot(net_spillover_ols_usa, type = "l", las = 2,  ylab = 'Net Spillover', main = 'USA')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
plot(net_spillover_ols_China, type = "l", las = 2,  ylab = 'Net Spillover', main = 'China')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
plot(net_spillover_ols_Germany, type = "l", las = 2,  ylab = 'Net Spillover', main = 'Germany')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
plot(net_spillover_ols_Greece, type = "l", las = 2,  ylab = 'Net Spillover', main = 'Greece')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
plot(net_spillover_ols_SK, type = "l", las = 2,  ylab = 'Net Spillover', main = 'South Korea')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
plot(net_spillover_ols_mx, type = "l", las = 2, ylab = 'Net Spillover', main = 'Mexico')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
plot(net_spillover_ols_Russia, type = "l", las = 2,  ylab = 'Net Spillover', main = 'Russia')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
plot(net_spillover_ols_UK, type = "l", las = 2,  ylab = 'Net Spillover', main = 'UK')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
dev.off()
##bilateral##
png("desktop/typ/pairwise1-ols.png", width = 2824, height = 3010, res = 300)
par(mfrow=c(2,1))
plot(cnct_ols_US_CN, type = "l", main = 'US-China', las = 2, ylab = '')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
plot(cnct_ols_US_MX, type = "l", main = 'US-Mexico', las = 2, ylab = '')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
dev.off()

png("desktop/typ/pairwise2-ols.png", width = 2824, height = 3010, res = 300)
par(mfrow=c(2,1))
plot(cnct_ols_US_RU, type = "l", main = 'US-Russia', las = 2, ylab = '')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
plot(cnct_ols_CN_HK, type = "l", main = 'China-Hong Kong', las = 2, ylab = '')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
dev.off()

png("desktop/typ/pairwise3-ols.png", width = 2824, height = 3010, res = 300)
par(mfrow=c(2,1))
plot(cnct_ols_GR_GE, type = "l", main = 'Germany-Greece', las = 2, ylab = '')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
plot(cnct_ols_CN_SK, type = "l", main = 'China-South Korea', las = 2, ylab = '')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
dev.off()


