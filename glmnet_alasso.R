library(glmnet)
library(vars)
library(readxl)
library(ggplot2)
library(ggrepel)                   

epu <- read_excel("~/Desktop/typ/data/data.xlsx")
epu <- epu[, 2:25]
epu <- as.matrix(epu)
# epu <- scale(epu)
View(epu)

### generate the ols coefficients ###
model_ols <- VAR(epu, ic = 'AIC')
coef_ols <- Bcoef(model_ols)[,-25] # the last column is the const

### generate the alasso coefficients ###
data=embed(epu,2)
one <- matrix(1, ncol = 24, nrow = 24)
alasso_weight <- one / coef_ols

coef_alasso <- matrix(NA, 24, 25)
res_alasso <- NULL
for (i in 1:24)
{
  model_alasso <- glmnet(x = data[,25:48], y = data[,i], lambda = 2.871421, alpha = 1, penalty.factor = alasso_weight[i,]) # in penalty.factor, the intercept is asked to be removed
  coef_list <- coef(model_alasso)
  coef_alasso[i, ] = as.matrix(coef_list)
  pred <- predict(model_alasso, newx= data[,25:48])
  res <- data[ ,i] - pred
  res_alasso <- cbind(res_alasso, res)
  Q_alasso <- (t(res_alasso)%*%res_alasso)/nrow(res_alasso)
}
# mean(res_alasso ^ 2)
# View(coef_alasso[,-1]) # the first vector is the intercept, remove it
# write.csv(coef_alasso[,-1],"desktop/typ/coef_alasso.csv")
# sum(coef_alasso[,-1] == 0) # 167 coefficients are shrunk to zeros
# diag(coef_alasso[,-1]) # see if the auto-regression term disappears
# dim(Q_alasso)

### connectedness ###
cnct_alasso <- GFEVD(coef_alasso[,-1], Q_alasso, n.ahead=6, standardize = T)$GFEVD * 100
# View(cnct_alasso)
rowSums(cnct_alasso)
colSums(cnct_alasso)

## static analysis ##
cnct_alasso_st <- cnct_analysis(cnct_alasso) # use cnct_analysis from the ols.R
cnct_alasso <- cbind(cnct_alasso, cnct_alasso_st$From)
cnct_alasso <- rbind(cnct_alasso, append(cnct_alasso_st$To, 0), append(cnct_alasso_st$Net, 0))

cnct_alasso_st$Sum
cnct_alasso[nrow(cnct_alasso),ncol(cnct_alasso)] <- cnct_alasso_st$Sum
# write.csv(cnct_alasso,"desktop/typ/cnct_alasso.csv")

## dynamic analysis ##
### generate the rolling window cnectedness tables
numofzero_dy <- NULL
cnct_alasso_rw <- array(NA, dim = c(24, 24, 172))
res_alasso_dy_collection <- array(NA, dim = c(50, 24, 172))
cnct_alasso_fromrw <- array(NA, dim =c(24,1,172))
cnct_alasso_torw <- array(NA, dim =c(24,1,172))
for (t in 1:172)
{
  coef_alasso_dy <- matrix(NA, 24, 25)
  res_alasso_dy <- NULL
  for (i in 1:24)
  {
    model_ols_dy <- VAR(epu[(0:49)+t, ], ic = 'AIC')
    coef_ols_dy <- Bcoef(model_ols_dy)[,-25]
    one <- matrix(1, ncol = 24, nrow = 24)
    alasso_weight_dy <- one / coef_ols_dy
    model_alasso_dy <- glmnet(x = data[(0:49) + t, 25:48], y = data[(0:49) + t, i], lambda = 3.320117, alpha = 1, penalty.factor = alasso_weight_dy[i,]) 
    # an optional value of hyperparameter is 26.94545, but it's coming from cv assuming non-time series process.
    coef_list_dy <- coef(model_alasso_dy) 
    coef_alasso_dy[i,] = as.matrix(coef_list_dy)
    pred_dy <- predict(model_alasso_dy, newx= data[(0:49) + t, 25:48])
    res_dy <- data[(0:49) + t,i] - pred_dy
    res_alasso_dy <- cbind(res_alasso_dy, res_dy)
  }
  numofzero <- sum(coef_alasso_dy[,-1] == 0)
  numofzero_dy <- rbind(numofzero_dy, numofzero/576) # 576 is the number of the slope coefficient
  Q_alasso_dy <- (t(res_alasso_dy)%*%res_alasso_dy) / nrow(res_alasso_dy)
  cnct_alasso_mid <- GFEVD(coef_alasso_dy[,-1], Q_alasso_dy, n.ahead=6, standardize = F)$GFEVD * 100
  res_alasso_dy_collection[,,t] <- as.matrix(res_alasso_dy)
  cnct_alasso_rw[,,t] <- as.matrix(cnct_alasso_mid)
  
  cnct_alasso_mid_dy <- cnct_analysis(cnct_alasso_mid)
  cnct_alasso_fromrw[,,t] <- as.matrix(cnct_alasso_mid_dy$From)
  cnct_alasso_torw[,,t] <- as.matrix(cnct_alasso_mid_dy$To)
}
cnct_alasso_rw[,,172]
mean(apply(res_alasso_dy_collection[,,1:172], 3 , function(x) mean(x ^ 2))) #mse of the rolling window

var(cnct_alasso_fromrw)
var(cnct_alasso_torw)

#check if there is a pattern in US residual
US_res <- NULL
for(i in 1:172){
  US_res <- rbind(US_res, res_alasso_dy_collection[1,1,i]^2)
}
US_res <- ts(US_res, frequency = 12, start = c(2003, 2))
plot(US_res)

### number of zeros
numofzero_dy <- ts(numofzero_dy, frequency = 12, start = c(2007, 3))
numofzero_mean <- NULL
for(i in 1:166){ 
  M <- mean(numofzero_dy[i:(i+5)])
  numofzero_mean[i] <- M
}
numofzero_mean <- ts(numofzero_mean,frequency = 12, start = c(2007, 9))
min(numofzero_dy)

### total, from, to, net spillover
zero <- ts(matrix(1,200) - 1, frequency = 12, start = c(2006, 1))
cnct_alasso_dy <- apply(cnct_alasso_rw, 3, cnct_analysis) #get the TO, FROM, TOTAL
##total
total_spillover_alasso <- matrix(NA, 172, 1)
for (j in 1:172){
  total_spillover1 <- cnct_alasso_dy[[j]]$Sum
  total_spillover_alasso[j,1] = as.numeric(total_spillover1)
} 
total_spillover_alasso <- ts(total_spillover_alasso, frequency = 12, start = c(2007, 3))

##TO
##USA to
to_spillover_alasso_usa <- matrix(NA, 172, 1)
for (j in 1:172){
  to_spillover1 <- cnct_alasso_dy[[j]]$To[1]
  to_spillover_alasso_usa[j,1] = as.numeric(to_spillover1)
} 
to_spillover_alasso_usa <- ts(to_spillover_alasso_usa, frequency = 12, start = c(2007, 3))
##China to
to_spillover_alasso_China <- matrix(NA, 172, 1)
for (j in 1:172){
  to_spillover1 <- cnct_alasso_dy[[j]]$To[7]
  to_spillover_alasso_China[j,1] = as.numeric(to_spillover1)
} 
to_spillover_alasso_China <- ts(to_spillover_alasso_China, frequency = 12, start = c(2007, 3))
##Germany to
to_spillover_alasso_Germany <- matrix(NA, 172, 1)
for (j in 1:172){
  to_spillover1 <- cnct_alasso_dy[[j]]$To[12]
  to_spillover_alasso_Germany[j,1] = as.numeric(to_spillover1)
} 
to_spillover_alasso_Germany <- ts(to_spillover_alasso_Germany, frequency = 12, start = c(2007, 3))
##Greece to 
to_spillover_alasso_Greece <- matrix(NA, 172, 1)
for (j in 1:172){
  to_spillover1 <- cnct_alasso_dy[[j]]$To[13]
  to_spillover_alasso_Greece[j,1] = as.numeric(to_spillover1)
} 
to_spillover_alasso_Greece <- ts(to_spillover_alasso_Greece, frequency = 12, start = c(2007, 3))
##SK to 
to_spillover_alasso_SK <- matrix(NA, 172, 1)
for (j in 1:172){
  to_spillover1 <- cnct_alasso_dy[[j]]$To[19]
  to_spillover_alasso_SK[j,1] = as.numeric(to_spillover1)
} 
to_spillover_alasso_SK <- ts(to_spillover_alasso_SK, frequency = 12, start = c(2007, 3))
##Mexico to
to_spillover_alasso_mx <- matrix(NA, 172, 1)
for (j in 1:172){
  to_spillover1 <- cnct_alasso_dy[[j]]$To[20]
  to_spillover_alasso_mx[j,1] = as.numeric(to_spillover1)
} 
to_spillover_alasso_mx <- ts(to_spillover_alasso_mx, frequency = 12, start = c(2007, 3))
##Russia to
to_spillover_alasso_Russia <- matrix(NA, 172, 1)
for (j in 1:172){
  to_spillover1 <- cnct_alasso_dy[[j]]$To[21]
  to_spillover_alasso_Russia[j,1] = as.numeric(to_spillover1)
} 
to_spillover_alasso_Russia <- ts(to_spillover_alasso_Russia, frequency = 12, start = c(2007, 3))
##UK to
to_spillover_alasso_UK <- matrix(NA, 172, 1)
for (j in 1:172){
  to_spillover1 <- cnct_alasso_dy[[j]]$To[24]
  to_spillover_alasso_UK[j,1] = as.numeric(to_spillover1)
} 
to_spillover_alasso_UK <- ts(to_spillover_alasso_UK, frequency = 12, start = c(2007, 3))


##From
##USA from
from_spillover_alasso_usa <- matrix(NA, 172, 1)
for (j in 1:172){
  from_spillover1 <- cnct_alasso_dy[[j]]$From[1]
  from_spillover_alasso_usa[j,1] = as.numeric(from_spillover1)
} 
from_spillover_alasso_usa <- ts(from_spillover_alasso_usa, frequency = 12, start = c(2007, 3))
##China from
from_spillover_alasso_China <- matrix(NA, 172, 1)
for (j in 1:172){
  from_spillover1 <- cnct_alasso_dy[[j]]$From[7]
  from_spillover_alasso_China[j,1] = as.numeric(from_spillover1)
} 
from_spillover_alasso_China <- ts(from_spillover_alasso_China, frequency = 12, start = c(2007, 3))
##Germany from
from_spillover_alasso_Germany <- matrix(NA, 172, 1)
for (j in 1:172){
  from_spillover1 <- cnct_alasso_dy[[j]]$From[12]
  from_spillover_alasso_Germany[j,1] = as.numeric(from_spillover1)
} 
from_spillover_alasso_Germany <- ts(from_spillover_alasso_Germany, frequency = 12, start = c(2007, 3))
##Greece from
from_spillover_alasso_Greece <- matrix(NA, 172, 1)
for (j in 1:172){
  from_spillover1 <- cnct_alasso_dy[[j]]$From[13]
  from_spillover_alasso_Greece[j,1] = as.numeric(from_spillover1)
} 
from_spillover_alasso_Greece <- ts(from_spillover_alasso_Greece, frequency = 12, start = c(2007, 3))
##SK from
from_spillover_alasso_SK <- matrix(NA, 172, 1)
for (j in 1:172){
  from_spillover1 <- cnct_alasso_dy[[j]]$From[19]
  from_spillover_alasso_SK[j,1] = as.numeric(from_spillover1)
} 
from_spillover_alasso_SK <- ts(from_spillover_alasso_SK, frequency = 12, start = c(2007, 3))
##Mexico from
from_spillover_alasso_mx <- matrix(NA, 172, 1)
for (j in 1:172){
  from_spillover1 <- cnct_alasso_dy[[j]]$From[20]
  from_spillover_alasso_mx[j,1] = as.numeric(from_spillover1)
} 
from_spillover_alasso_mx <- ts(from_spillover_alasso_mx, frequency = 12, start = c(2007, 3))
##Russia from
from_spillover_alasso_Russia <- matrix(NA, 172, 1)
for (j in 1:172){
  from_spillover1 <- cnct_alasso_dy[[j]]$From[21]
  from_spillover_alasso_Russia[j,1] = as.numeric(from_spillover1)
} 
from_spillover_alasso_Russia <- ts(from_spillover_alasso_Russia, frequency = 12, start = c(2007, 3))
##UK from
from_spillover_alasso_UK <- matrix(NA, 172, 1)
for (j in 1:172){
  from_spillover1 <- cnct_alasso_dy[[j]]$From[24]
  from_spillover_alasso_UK[j,1] = as.numeric(from_spillover1)
} 
from_spillover_alasso_UK <- ts(from_spillover_alasso_UK, frequency = 12, start = c(2007, 3))

##Net
zero <- ts(matrix(1,200) - 1, frequency = 12, start = c(2006, 1))
##USA net
net_spillover_alasso_usa <- matrix(NA, 172, 1)
for (j in 1:172){
  net_spillover1 <- cnct_alasso_dy[[j]]$Net[1]
  net_spillover_alasso_usa[j,1] = as.numeric(net_spillover1)
} 
net_spillover_alasso_usa <- ts(net_spillover_alasso_usa, frequency = 12, start = c(2007, 3))
##China net
net_spillover_alasso_China <- matrix(NA, 172, 1)
for (j in 1:172){
  net_spillover1 <- cnct_alasso_dy[[j]]$Net[7]
  net_spillover_alasso_China[j,1] = as.numeric(net_spillover1)
} 
net_spillover_alasso_China <- ts(net_spillover_alasso_China, frequency = 12, start = c(2007, 3))
##Germany net
net_spillover_alasso_Germany <- matrix(NA, 172, 1)
for (j in 1:172){
  net_spillover1 <- cnct_alasso_dy[[j]]$Net[12]
  net_spillover_alasso_Germany[j,1] = as.numeric(net_spillover1)
} 
net_spillover_alasso_Germany <- ts(net_spillover_alasso_Germany, frequency = 12, start = c(2007, 3))
##Greece net
net_spillover_alasso_Greece <- matrix(NA, 172, 1)
for (j in 1:172){
  net_spillover1 <- cnct_alasso_dy[[j]]$Net[13]
  net_spillover_alasso_Greece[j,1] = as.numeric(net_spillover1)
} 
net_spillover_alasso_Greece <- ts(net_spillover_alasso_Greece, frequency = 12, start = c(2007, 3))
##SK net
net_spillover_alasso_SK <- matrix(NA, 172, 1)
for (j in 1:172){
  net_spillover1 <- cnct_alasso_dy[[j]]$Net[19]
  net_spillover_alasso_SK[j,1] = as.numeric(net_spillover1)
} 
net_spillover_alasso_SK <- ts(net_spillover_alasso_SK, frequency = 12, start = c(2007, 3))
##Mexico net
net_spillover_alasso_mx <- matrix(NA, 172, 1)
for (j in 1:172){
  net_spillover1 <- cnct_alasso_dy[[j]]$Net[20]
  net_spillover_alasso_mx[j,1] = as.numeric(net_spillover1)
} 
net_spillover_alasso_mx <- ts(net_spillover_alasso_mx, frequency = 12, start = c(2007, 3))
##Russia net
net_spillover_alasso_Russia <- matrix(NA, 172, 1)
for (j in 1:172){
  net_spillover1 <- cnct_alasso_dy[[j]]$Net[21]
  net_spillover_alasso_Russia[j,1] = as.numeric(net_spillover1)
} 
net_spillover_alasso_Russia <- ts(net_spillover_alasso_Russia, frequency = 12, start = c(2007, 3))
##UK net
net_spillover_alasso_UK <- matrix(NA, 172, 1)
for (j in 1:172){
  net_spillover1 <- cnct_alasso_dy[[j]]$Net[24]
  net_spillover_alasso_UK[j,1] = as.numeric(net_spillover1)
} 
net_spillover_alasso_UK <- ts(net_spillover_alasso_UK, frequency = 12, start = c(2007, 3))


### pairwise spillover
## this part only plots one line, which is the difference between pairwise spillovers of each economies
##us-cn, us-mexico
cnct_alasso_US_CN <- apply(cnct_alasso_rw, 3, cnct_pairwise_analysis, i = 1, j = 7)
cnct_alasso_US_CN <- ts(cnct_alasso_US_CN, frequency = 12, start = c(2007, 3))
cnct_alasso_US_MX <- apply(cnct_alasso_rw, 3, cnct_pairwise_analysis, i = 1, j = 20)
cnct_alasso_US_MX <- ts(cnct_alasso_US_MX, frequency = 12, start = c(2007, 3))

##us-russia, cn-hk
cnct_alasso_US_RU <- apply(cnct_alasso_rw, 3, cnct_pairwise_analysis, i = 1, j = 21)
cnct_alasso_US_RU <- ts(cnct_alasso_US_RU, frequency = 12, start = c(2007, 3))
cnct_alasso_CN_HK <- apply(cnct_alasso_rw, 3, cnct_pairwise_analysis, i = 7, j = 14)
cnct_alasso_CN_HK <- ts(cnct_alasso_CN_HK, frequency = 12, start = c(2007, 3))

##gr-ge, cn-sk
cnct_alasso_GR_GE <- apply(cnct_alasso_rw, 3, cnct_pairwise_analysis, i = 12, j = 13)
cnct_alasso_GR_GE <- ts(cnct_alasso_GR_GE, frequency = 12, start = c(2007, 3))
cnct_alasso_CN_SK <- apply(cnct_alasso_rw, 3, cnct_pairwise_analysis, i = 7, j = 19)
cnct_alasso_CN_SK <- ts(cnct_alasso_CN_SK, frequency = 12, start = c(2007, 3))

##us-japan, cn-sk
cnct_alasso_US_JP <- apply(cnct_alasso_rw, 3, cnct_pairwise_analysis, i = 1, j = 18)
cnct_alasso_US_JP <- ts(cnct_alasso_US_JP, frequency = 12, start = c(2007, 3))
cnct_alasso_US_SK <- apply(cnct_alasso_rw, 3, cnct_pairwise_analysis, i = 1, j = 19)
cnct_alasso_US_SK <- ts(cnct_alasso_US_SK, frequency = 12, start = c(2007, 3))

## this part draws two lines, which are pairwise spillover of each country
cnct_alasso_UStoCN <- ts(cnct_alasso_rw[7,1,], frequency = 12, start = c(2007, 3))
cnct_alasso_CNtoUS <- ts(cnct_alasso_rw[1,7,], frequency = 12, start = c(2007, 3))
cnct_alasso_UStoMX <- ts(cnct_alasso_rw[20,1,], frequency = 12, start = c(2007, 3))
cnct_alasso_MXtoUS <- ts(cnct_alasso_rw[1,20,], frequency = 12, start = c(2007, 3))
cnct_alasso_UStoRU <- ts(cnct_alasso_rw[21,1,], frequency = 12, start = c(2007, 3))
cnct_alasso_RUtoUS <- ts(cnct_alasso_rw[1,21,], frequency = 12, start = c(2007, 3))
cnct_alasso_CNtoHK <- ts(cnct_alasso_rw[14,7,], frequency = 12, start = c(2007, 3))
cnct_alasso_HKtoCN <- ts(cnct_alasso_rw[7,14,], frequency = 12, start = c(2007, 3))
cnct_alasso_CNtoSK <- ts(cnct_alasso_rw[19,7,], frequency = 12, start = c(2007, 3))
cnct_alasso_SKtoCN <- ts(cnct_alasso_rw[7,19,], frequency = 12, start = c(2007, 3))
cnct_alasso_GRtoGE <- ts(cnct_alasso_rw[13,12,], frequency = 12, start = c(2007, 3))
cnct_alasso_GEtoGR <- ts(cnct_alasso_rw[12,13,], frequency = 12, start = c(2007, 3))
cnct_alasso_GRtoSP <- ts(cnct_alasso_rw[22,12,], frequency = 12, start = c(2007, 3))
cnct_alasso_SPtoGR <- ts(cnct_alasso_rw[12,22,], frequency = 12, start = c(2007, 3))
cnct_alasso_GRtoFR <- ts(cnct_alasso_rw[11,12,], frequency = 12, start = c(2007, 3))
cnct_alasso_FRtoGR <- ts(cnct_alasso_rw[12,11,], frequency = 12, start = c(2007, 3))





#test pairwise connectedness
png("desktop/typ/test.png", width = 4854, height = 3000, res = 350)
cnct_alasso_1_2 <- apply(cnct_alasso_rw, 3, cnct_pairwise_analysis, i = 12, j =22)
cnct_alasso_1_2 <- ts(cnct_alasso_1_2, frequency = 12, start = c(2007, 3))
plot(cnct_alasso_1_2, type = "l", main = '1-2', las = 2, ylab = '')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
dev.off()


### decide the best lambda ###
##static
bestlambda <- matrix(NA, 24, 1)
for (i in 1:24){
  cvfit <- cv.glmnet(x = data[,25:48], y = data[,i], penalty.factor = alasso_weight[i,], alpha = 1)
  bestlambda[i,1] <- cvfit$lambda.min
}
bestlambda
plot(bestlambda)
mean(bestlambda)
sd(bestlambda)

onebyone <- cv.glmnet(x = data[,25:48], y = data[,24], penalty.factor = alasso_weight[24,], alpha = 1)
plot(onebyone)

fit_usa <- glmnet(x = data[,25:48], y = data[,1], penalty.factor = alasso_weight[1,], alpha = 1)
plot(fit_usa, xvar = "lambda", label = TRUE)

##dynamic
min_lambda <- NULL
for (t in 1:172)
{
  lambda_min <- matrix(NA, 1, 24)
  for (i in 1:24)
  {
    model_ols_dy_cv <- VAR(epu[(0:49)+t, ], ic = 'AIC')
    coef_ols_dy_cv <- Bcoef(model_ols_dy_cv)[,-25]
    one <- matrix(1, ncol = 24, nrow = 24)
    alasso_weight_dy_cv <- one / coef_ols_dy_cv
    model_alasso_test_cv <- cv.glmnet(x = data[(0:49) + t, 25:48], y = data[(0:49) + t, i], alpha = 1, penalty.factor = alasso_weight_dy_cv[i,])
    lambda_min[, i] <- model_alasso_test_cv$lambda.min
  }
   min_lambda <- rbind(min_lambda, lambda_min)
}
apply(min_lambda, 2, mean)
mean(min_lambda)
sd(min_lambda)

## cross validation for time series
model_ols_cv <- VAR(epu, ic = 'AIC')
coef_ols_cv <- Bcoef(model_ols_cv)[,-25]
one <- matrix(1, ncol = 24, nrow = 24)
alasso_weight_cv <- one / coef_ols_cv
K = 10 # 10 fold
size_train = floor(nrow(epu)/K) - 1
mse_total <- NULL
for (t in 1:K)
{
  mse_by_train_group <- NULL
  for (j in seq(from = 0, to = 5, by = 0.1))
  {
    mse_by_hyper <- NULL
    coef_alasso_cv <- matrix(NA,24,25)
    for (i in 1:ncol(epu))
    {
      # print("Before glmnet call")  # Add this line
      model_alasso_cv <- glmnet(x = data[(1:(t*size_train)), 25:48], y = data[(1:(t*size_train)), i], lambda = exp(j), alpha = 1, penalty.factor = alasso_weight_cv[i,])
      # print("After glmnet call")  # Add this line
      coef_list_cv <- coef(model_alasso_cv) 
      coef_alasso_cv[i,] = as.matrix(coef_list_cv)
      pred_cv <- predict(model_alasso_cv, newx= data[(t*size_train+1):(t*size_train+10), 25:48])
      res_cv <- data[(t*size_train+1):(t*size_train+10),i] - pred_cv
      mse_individual_cv <- sum(res_cv^2) / 10
      mse_by_hyper <- rbind(mse_by_hyper, mse_individual_cv)
    }
    mse_cv <- mean(mse_by_hyper) #derive the mse of the whole regression of a certain j and t
    mse_by_train_group <- rbind(mse_by_train_group,mse_cv)#stack elements into a column according to j
  }
  mse_total <- cbind(mse_total,mse_by_train_group)#stack columns into a matrix according to t
}
means_across_group <- rowMeans(mse_total)#average within t
# print(means_across_group)
smallest_location <- which(means_across_group == min(means_across_group))
print(smallest_location)#the location marks the j, with j you can derive the optimal lambda

###symposium###
symposium_epu <- ts(epu, frequency = 12, start = c(2003, 1))

png("desktop/us&canada.png", width = 3000, height = 2000, res = 350)
plot(symposium_epu[,5], las = 2, type = "l", ylab = "EPU", main = "EPU of the US and Canada", col = "blue")
lines(symposium_epu[,1], col = "red")
legend("topleft", c("USA","Canada"), fill=c("red","blue"))
axis(1, at = c(2003:2022), las = 2)
dev.off()

png("desktop/us&mexico.png", width = 3000, height = 2000, res = 350)
plot(symposium_epu[,1], las = 2, type = "l", ylab = "EPU", main = "EPU of the US and Mexico", col = "red", ylim =c(0, 600))
lines(symposium_epu[,20], col = "green")
legend("topleft", c("USA","Mexico"), fill=c("red","green"))
axis(1, at = c(2003:2022), las = 2)
dev.off()
png("desktop/uscnada_sym.png", width = 3000, height = 2000, res = 350)
plot(cnct_alasso_1_2, type = "l", main = 'US-Canada', las = 2, ylab = 'aLasso')
lines(zero, lty=2)
axis(1, at = c(2007:2022), las = 2)
dev.off()
  ### function ###
# LassoVAR = function(y,p) {
#   y = as.matrix(y)
#   p = p
#   res1 = coef1 = NULL
#   k = ncol(y)
#   for (i in 1:k) {
#     yx = embed(y,p+1)
#     x = yx[,-c(1:k)]
#     z = y[-c(1:p),i]
#     mod = cv.glmnet(x, z, alpha=1, type.measure="mae")
#     pred = predict(mod, s = mod$lambda.min, newx=yx[,-c(1:k)])
#     coef2 = predict(mod, type = "coefficients", s = mod$lambda.min)[-1] # remove the intercept term
#     res = z - pred
#     coef1 = rbind(coef1,coef2) # each row is the coefficients
#     res1 = cbind(res1,res)
#   }
#   Q = (t(res1)%*%res1)/nrow(res1)
#   results = list(B=coef1,Q=Q)
# }
# abc <- LassoVAR(epu,1)
# dim(abc$Q)
