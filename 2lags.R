####use 2 lags in AdaLasso####

##generate x and y##
epu_2lags = embed(epu,3)
dependent_epu = epu_2lags[,1:24]
independent_epu = epu_2lags[,25:72]


### STATIC ANALYSIS
bestlambda <- matrix(NA, 24, 1)
for (i in 1:24){
  cvfit <- cv.glmnet(x = independent_epu, y = dependent_epu[,i], penalty.factor = vol_alasso_weight_cv[i,], alpha = 1)
  bestlambda[i,1] <- cvfit$lambda.min
}
# bestlambda
# plot(bestlambda)
mean(bestlambda) # 4.761924
# sd(bestlambda)

## generate parameters for GFEVD function
coef_alasso <- matrix(NA, 24, 49)
res_alasso <- NULL
for (i in 1:24)
{
  model_alasso <- glmnet(x = independent_epu, y = dependent_epu[,i], lambda = 4.761924, alpha = 1, penalty.factor = vol_alasso_weight_cv[i,]) # in penalty.factor, the intercept is asked to be removed
  coef_list <- coef(model_alasso)
  coef_alasso[i, ] = as.matrix(coef_list)
  pred <- predict(model_alasso, newx= epu_2lags[,25:72])
  res <- epu_2lags[ ,i] - pred
  res_alasso <- cbind(res_alasso, res)
  Q_alasso <- (t(res_alasso)%*%res_alasso)/nrow(res_alasso)
}
GVD_nor_preweighted_static <- GFEVD(coef_alasso[,-1], Q_alasso, n.ahead=6, standardize = T)$GFEVD * 100
# write.csv(GVD_nor_preweighted_static,"desktop/var(2)/GVD_nor_preweighted_static.csv")

# TO
static_TO <- Calculate_To(GVD_nor_preweighted_static, matrix(static_weight_nor[,1])) #成功
View(static_TO$weight)
# write.csv(static_TO,"desktop/var(2)/static_TO.csv")

# FROM
static_FROM <- cnct_analysis(GVD_nor_preweighted_static) # use cnct_analysis from the ols.R
static_FROM <- matrix(static_FROM$From)
# write.csv(static_FROM,"desktop/var(2)/static_FROM.csv")

# TOTAL
static_TOTAL1 <- give_weight_updated(GVD_nor_preweighted_static, static_weight_nor, 1)
# write.csv(static_TOTAL1,"desktop/var(2)/weighted_GVD.csv")
static_TOTAL <- cnct_analysis_updated(static_TOTAL1)$Sum * 24

# quadrants
lags2_static_quadrants <- read_excel("~/Desktop/var(2)/2lags_static_quadrant.xlsx")
png("desktop/var(2)/2lags_static_quadrants.png", width = 3570, height = 1865, res = 350)
ggplot(lags2_static_quadrants, aes(y=From, x=To, label=country)) +
  geom_text_repel(aes(y=From, x=To, label=country)) +
  geom_point() +
  scale_y_continuous(breaks=seq(0,100,10),limits=c(0,100)) + 
  scale_x_continuous(breaks=seq(0,200,10), limits=c(0,200)) +
  # #lims(y=c(70,110),x=c(-50,250)) +
  #coord_fixed() +  # 
  geom_vline(xintercept = mean(lags2_static_quadrants$To)) + geom_hline(yintercept = mean(lags2_static_quadrants$From)) + 
  geom_abline(intercept = 0)
dev.off()

### DYNAMIC ANALYSIS

## cross validation for time series
epu_ols_cv <- VAR(epu, p = 2, ic = 'AIC')
epu_coef_ols_cv <- Bcoef(epu_ols_cv)[,-49]
one <- matrix(1, ncol = 48, nrow = 24)
epu_alasso_weight_cv <- one / epu_coef_ols_cv
#if I want the length of the rolling window in cv similar with that in the estimates, the size is 222, size of rw in estimate is 60, the K should be 3.7, rounded up to 4
K = 4 # 4 fold
size_train = floor(nrow(epu)/K) - 1 #size_train = 54
mse_total <- NULL
for (t in 1:K)
{
  mse_by_train_group <- NULL
  for (j in seq(from = 0, to = 5, by = 0.1)) # lambda move by 0.1 from 0 to 5
  {
    mse_by_hyper <- NULL
    vol_coef_alasso_cv <- matrix(NA,24,49)
    for (i in 1:ncol(epu))
    {
      # print("Before glmnet call")  # Add this line
      model_alasso_cv <- glmnet(x = independent_epu[(((t-1)*size_train)+1):(t*size_train),], y = dependent_epu[(((t-1)*size_train)+1):(t*size_train), i], lambda = exp(j), alpha = 1, penalty.factor = epu_alasso_weight_cv[i,])
      # print("After glmnet call")  # Add this line
      coef_list_cv <- coef(model_alasso_cv) 
      vol_coef_alasso_cv[i,] = as.matrix(coef_list_cv)
      pred_cv <- predict(model_alasso_cv, newx= independent_epu[(t*size_train+1):(t*size_train+4),]) # making predictions, but cannot exceed the size of independent_epu, so max is 4, 220 - (4*54)
      res_cv <- dependent_epu[(t*size_train+1):(t*size_train+4),i] - pred_cv
      mse_individual_cv <- sum(res_cv^2) / 4
      mse_by_hyper <- rbind(mse_by_hyper, mse_individual_cv)
    }
    mse_cv <- mean(mse_by_hyper) # calculate the mean of residuals from 4 predicitions
    mse_by_train_group <- rbind(mse_by_train_group,mse_cv) # stack the mean of residuals into a column according to the ascending order of hyperparameter
  }
  mse_total <- cbind(mse_total,mse_by_train_group) # each row represents each hyperparameter, each column represents each fold
}
means_across_group <- rowMeans(mse_total) #calculate row means of every column, meaning calculate the average of each hyperparameter cross folds
# print(means_across_group)
smallest_location <- which(means_across_group == min(means_across_group))
print(smallest_location) # here it is 36, so optimal lambda is exp(3.6)

## dynamic analysis ##
### generate the rolling window cnectedness tables
num_nonzero_coef_1st_lag <- NULL
num_nonzero_coef_2nd_lag <- NULL
GVD_nor_alasso_rw <- array(NA, dim = c(24, 24, 161)) # normalized GVD
res_alasso_dy_collection <- array(NA, dim = c(60, 24, 161)) #residual matrix
# cnct_alasso_fromrw <- array(NA, dim =c(24,1,161))
# cnct_alasso_torw <- array(NA, dim =c(24,1,161))
for (t in 1:161)
{
  coef_alasso_dy <- matrix(NA, 24, 49)
  res_alasso_dy <- NULL
  for (i in 1:24)
  {
    model_ols_dy <- VAR(epu[(0:59)+t, ], p = 2, ic = 'AIC')
    coef_ols_dy <- Bcoef(model_ols_dy)[,-49]
    one <- matrix(1, ncol = 48, nrow = 24)
    alasso_weight_dy <- one / coef_ols_dy
    model_alasso_dy <- glmnet(x = epu_2lags[(0:59) + t, 25:72], y = epu_2lags[(0:59) + t, i], lambda = 36.5982344437, alpha = 1, penalty.factor = alasso_weight_dy[i,]) 
    # an optional value of hyperparameter is 26.94545, but it's coming from cv assuming non-time series process.
    coef_list_dy <- coef(model_alasso_dy) 
    coef_alasso_dy[i,] = as.matrix(coef_list_dy)
    pred_dy <- predict(model_alasso_dy, newx= epu_2lags[(0:59) + t, 25:72])
    res_dy <- epu_2lags[(0:59) + t,i] - pred_dy
    res_alasso_dy <- cbind(res_alasso_dy, res_dy)
  }
  numofzero_1st <- sum(coef_alasso_dy[,1:24] == 0)
  numofzero_2nd <- sum(coef_alasso_dy[,24:48] == 0)
  num_nonzero_coef_1st_lag <- rbind(num_nonzero_coef_1st_lag, 1 - numofzero_1st/576) # 1152 is the number of the slope coefficient
  num_nonzero_coef_2nd_lag <- rbind(num_nonzero_coef_2nd_lag, 1 - numofzero_2nd/576) # 1152 is the number of the slope coefficient
  Q_alasso_dy <- (t(res_alasso_dy)%*%res_alasso_dy) / nrow(res_alasso_dy)
  cnct_alasso_mid <- GFEVD(coef_alasso_dy[,-1], Q_alasso_dy, n.ahead=9, standardize = F)$GFEVD * 100 # horizon= 3, 6, 9
  # res_alasso_dy_collection[,,t] <- as.matrix(res_alasso_dy)
  GVD_nor_alasso_rw[,,t] <- as.matrix(cnct_alasso_mid)
  
  # cnct_alasso_mid_dy <- cnct_analysis(cnct_alasso_mid)
  # cnct_alasso_fromrw[,,t] <- as.matrix(cnct_alasso_mid_dy$From)
  # cnct_alasso_torw[,,t] <- as.matrix(cnct_alasso_mid_dy$To)
}
GVD_nor_alasso_rw[,,161]

### number of zeros
num_nonzero_coef_1st_lag <- ts(num_nonzero_coef_1st_lag, frequency = 12, start = c(2008, 2))
num_nonzero_coef_2nd_lag <- ts(num_nonzero_coef_2nd_lag, frequency = 12, start = c(2008, 2))


# calculate weights
rolling_weight_GDP2 <- read_excel("~/Desktop/typ/data/GDP_interpolate1.xlsx") # GDP from July 02 to June 21, sample size 228
rolling_weight_GDP3 <- rolling_weight_GDP2[7:228, 1:24] # weight started from Jan 03
Moving_Average_GDP <- weight_moving_average(rolling_weight_GDP3[3:222,],60)
Moving_Average_GDP_nor <- Moving_Average_GDP / apply(Moving_Average_GDP, 1, sum)
arr_Moving_Average_GDP_nor <- array(t(Moving_Average_GDP_nor), dim=c(nrow(t(Moving_Average_GDP_nor)), 1, ncol(t(Moving_Average_GDP_nor))))

#unweighted cnct
GVD_nor_alasso_dy <- apply(GVD_nor_alasso_rw, 3, cnct_analysis) 

# give weight
# TO
GVD_nor_weighted_To <- Rolling_Calculate_To(GVD_nor_alasso_rw, arr_Moving_Average_GDP_nor) #成功

GVD_nor_unweighted_To <- array(NA, dim = c(24, 1, 161))
to_spillover_test <- matrix(NA, 24, 1)
for (j in 1:161){
  to_spillover_test <- matrix(GVD_nor_alasso_dy[[j]]$To)
  GVD_nor_unweighted_To[,,j] <- to_spillover_test
} 


# FROM
GVD_nor_unweighted_From <- array(NA, dim = c(24, 1, 161))
from_spillover_test <- matrix(NA, 24, 1)
for (j in 1:161){
  from_spillover_test <- matrix(GVD_nor_alasso_dy[[j]]$From)
  GVD_nor_unweighted_From[,,j] <- from_spillover_test
} 

# Total
GVD_normalized_weighted_rw <- give_weight_rw(GVD_nor_alasso_rw, Moving_Average_GDP_nor) 
GVD_nor_weighted_dy <- apply(GVD_normalized_weighted_rw, 3, cnct_analysis)

Weighted_2lags_Total <- matrix(NA, 161, 1)
for (j in 1:161){
  Total_spillover1 <- GVD_nor_weighted_dy[[j]]$Sum
  Weighted_2lags_Total[j,1] = as.numeric(Total_spillover1)
} 
Weighted_2lags_Total <- ts(Weighted_2lags_Total, frequency = 12, start = c(2008, 2))

#OLS with sampling size of 60
cnct_ols_rw <- array(NA, dim = c(24, 24, 163))
res_ols_dy_collection <- array(NA, dim = c(59, 24, 163))
cnct_ols_fromrw <- array(NA, dim =c(24,1,163))
cnct_ols_torw <- array(NA, dim =c(24,1,163))
for (t in 1:163)
{
  model_ols_dy <- VAR(epu[(0:59)+t,], ic = 'AIC')
  coef_ols_dy <- Bcoef(model_ols_dy)[,-25]
  res_ols_dy <- resid(model_ols_dy) 
  Q_ols_dy <- (t(res_ols_dy)%*%res_ols_dy)/nrow(res_ols_dy)
  cnct_ols_mid <- GFEVD(coef_ols_dy, Q_ols_dy, n.ahead=9, standardize = T)$GFEVD * 100
  res_ols_dy_collection[,,t] <- as.matrix(res_ols_dy)
  cnct_ols_rw[,,t] <- as.matrix(cnct_ols_mid)
  
  cnct_ols_mid_dy <- cnct_analysis(cnct_ols_mid)
  cnct_ols_fromrw[,,t] <- as.matrix(cnct_ols_mid_dy$From)
  cnct_ols_torw[,,t] <- as.matrix(cnct_ols_mid_dy$To)
}

# cnct_ols_rw[,,161]
cnct_ols_dy <- apply(cnct_ols_rw, 3, cnct_analysis)
#re-run 88-93
##total
sum_spillover_ols <- matrix(NA, 163, 1)
for (j in 1:163){
  sum_spillover1 <- cnct_ols_dy[[j]]$Sum
  sum_spillover_ols[j,1] = as.numeric(sum_spillover1)
} 
sum_spillover_ols <- ts(sum_spillover_ols, frequency = 12, start = c(2007, 12))

#re-run 95-210
##TO
##USA to
to_spillover_ols_usa <- matrix(NA, 163, 1)
for (j in 1:163){
  to_spillover1 <- cnct_ols_dy[[j]]$To[1]
  to_spillover_ols_usa[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_usa <- ts(to_spillover_ols_usa, frequency = 12, start = c(2007, 12))
##China to
to_spillover_ols_China <- matrix(NA, 163, 1)
for (j in 1:163){
  to_spillover1 <- cnct_ols_dy[[j]]$To[7]
  to_spillover_ols_China[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_China <- ts(to_spillover_ols_China, frequency = 12, start = c(2007, 12))
##Germany to
to_spillover_ols_Germany <- matrix(NA, 163, 1)
for (j in 1:163){
  to_spillover1 <- cnct_ols_dy[[j]]$To[12]
  to_spillover_ols_Germany[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_Germany <- ts(to_spillover_ols_Germany, frequency = 12, start = c(2007, 12))
##Greece to 
to_spillover_ols_Greece <- matrix(NA, 163, 1)
for (j in 1:163){
  to_spillover1 <- cnct_ols_dy[[j]]$To[13]
  to_spillover_ols_Greece[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_Greece <- ts(to_spillover_ols_Greece, frequency = 12, start = c(2007, 12))
##SK to 
to_spillover_ols_SK <- matrix(NA, 163, 1)
for (j in 1:163){
  to_spillover1 <- cnct_ols_dy[[j]]$To[19]
  to_spillover_ols_SK[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_SK <- ts(to_spillover_ols_SK, frequency = 12, start = c(2007, 12))
##Mexico to
to_spillover_ols_mx <- matrix(NA, 163, 1)
for (j in 1:163){
  to_spillover1 <- cnct_ols_dy[[j]]$To[20]
  to_spillover_ols_mx[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_mx <- ts(to_spillover_ols_mx, frequency = 12, start = c(2007, 12))
##Russia to
to_spillover_ols_Russia <- matrix(NA, 163, 1)
for (j in 1:163){
  to_spillover1 <- cnct_ols_dy[[j]]$To[21]
  to_spillover_ols_Russia[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_Russia <- ts(to_spillover_ols_Russia, frequency = 12, start = c(2007, 12))
##UK to
to_spillover_ols_UK <- matrix(NA, 163, 1)
for (j in 1:163){
  to_spillover1 <- cnct_ols_dy[[j]]$To[24]
  to_spillover_ols_UK[j,1] = as.numeric(to_spillover1)
} 
to_spillover_ols_UK <- ts(to_spillover_ols_UK, frequency = 12, start = c(2007, 12))


##From
##USA from
from_spillover_ols_usa <- matrix(NA, 163, 1)
for (j in 1:163){
  from_spillover1 <- cnct_ols_dy[[j]]$From[1]
  from_spillover_ols_usa[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_usa <- ts(from_spillover_ols_usa, frequency = 12, start = c(2007, 12))
##China from
from_spillover_ols_China <- matrix(NA, 163, 1)
for (j in 1:163){
  from_spillover1 <- cnct_ols_dy[[j]]$From[7]
  from_spillover_ols_China[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_China <- ts(from_spillover_ols_China, frequency = 12, start = c(2007, 12))
##Germany from
from_spillover_ols_Germany <- matrix(NA, 163, 1)
for (j in 1:163){
  from_spillover1 <- cnct_ols_dy[[j]]$From[12]
  from_spillover_ols_Germany[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_Germany <- ts(from_spillover_ols_Germany, frequency = 12, start = c(2007, 12))
##Greece from
from_spillover_ols_Greece <- matrix(NA, 163, 1)
for (j in 1:163){
  from_spillover1 <- cnct_ols_dy[[j]]$From[13]
  from_spillover_ols_Greece[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_Greece <- ts(from_spillover_ols_Greece, frequency = 12, start = c(2007, 12))
##SK from
from_spillover_ols_SK <- matrix(NA, 163, 1)
for (j in 1:163){
  from_spillover1 <- cnct_ols_dy[[j]]$From[19]
  from_spillover_ols_SK[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_SK <- ts(from_spillover_ols_SK, frequency = 12, start = c(2007, 12))
##Mexico from
from_spillover_ols_mx <- matrix(NA, 163, 1)
for (j in 1:163){
  from_spillover1 <- cnct_ols_dy[[j]]$From[20]
  from_spillover_ols_mx[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_mx <- ts(from_spillover_ols_mx, frequency = 12, start = c(2007, 12))
##Russia from
from_spillover_ols_Russia <- matrix(NA, 163, 1)
for (j in 1:163){
  from_spillover1 <- cnct_ols_dy[[j]]$From[21]
  from_spillover_ols_Russia[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_Russia <- ts(from_spillover_ols_Russia, frequency = 12, start = c(2007, 12))
##UK from
from_spillover_ols_UK <- matrix(NA, 163, 1)
for (j in 1:163){
  from_spillover1 <- cnct_ols_dy[[j]]$From[24]
  from_spillover_ols_UK[j,1] = as.numeric(from_spillover1)
} 
from_spillover_ols_UK <- ts(from_spillover_ols_UK, frequency = 12, start = c(2007, 12))

## LINES
# Australia, weighted
Australia_TO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_weighted_To[2,,j]
  Australia_TO[j,1] = as.numeric(To_spillover1)
} 
Australia_TO <- ts(Australia_TO, frequency = 12, start = c(2008, 2))

Australia_FROM <- matrix(NA, 161, 1)
for (j in 1:161){
  From_spillover1 <- GVD_nor_unweighted_From[2,,j]
  Australia_FROM[j,1] = as.numeric(From_spillover1)
} 
Australia_FROM <- ts(Australia_FROM, frequency = 12, start = c(2008, 2))



# US, weighted
US_TO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_weighted_To[1,,j]
  US_TO[j,1] = as.numeric(To_spillover1)
} 
US_TO <- ts(US_TO, frequency = 12, start = c(2008, 2))

US_FROM <- matrix(NA, 161, 1)
for (j in 1:161){
  From_spillover1 <- GVD_nor_unweighted_From[1,,j]
  US_FROM[j,1] = as.numeric(From_spillover1)
} 
US_FROM <- ts(US_FROM, frequency = 12, start = c(2008, 2))

US_unTO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_unweighted_To[1,,j]
  US_unTO[j,1] = as.numeric(To_spillover1)
} 
US_unTO <- ts(US_unTO, frequency = 12, start = c(2008, 2))

# China, weighted
China_TO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_weighted_To[7,,j]
  China_TO[j,1] = as.numeric(To_spillover1)
} 
China_TO <- ts(China_TO, frequency = 12, start = c(2008, 2))

China_FROM <- matrix(NA, 161, 1)
for (j in 1:161){
  From_spillover1 <- GVD_nor_unweighted_From[7,,j]
  China_FROM[j,1] = as.numeric(From_spillover1)
} 
China_FROM <- ts(China_FROM, frequency = 12, start = c(2008, 2))

China_unTO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_unweighted_To[7,,j]
  China_unTO[j,1] = as.numeric(To_spillover1)
} 
China_unTO <- ts(China_unTO, frequency = 12, start = c(2008, 2))

# Germany, weighted
Germany_TO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_weighted_To[12,,j]
  Germany_TO[j,1] = as.numeric(To_spillover1)
} 
Germany_TO <- ts(Germany_TO, frequency = 12, start = c(2008, 2))

Germany_FROM <- matrix(NA, 161, 1)
for (j in 1:161){
  From_spillover1 <- GVD_nor_unweighted_From[12,,j]
  Germany_FROM[j,1] = as.numeric(From_spillover1)
} 
Germany_FROM <- ts(Germany_FROM, frequency = 12, start = c(2008, 2))

Germany_unTO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_unweighted_To[12,,j]
  Germany_unTO[j,1] = as.numeric(To_spillover1)
} 
Germany_unTO <- ts(Germany_unTO, frequency = 12, start = c(2008, 2))


# Greece, weighted
Greece_TO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_weighted_To[13,,j]
  Greece_TO[j,1] = as.numeric(To_spillover1)
} 
Greece_TO <- ts(Greece_TO, frequency = 12, start = c(2008, 2))

Greece_FROM <- matrix(NA, 161, 1)
for (j in 1:161){
  From_spillover1 <- GVD_nor_unweighted_From[13,,j]
  Greece_FROM[j,1] = as.numeric(From_spillover1)
} 
Greece_FROM <- ts(Greece_FROM, frequency = 12, start = c(2008, 2))

Greece_unTO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_unweighted_To[13,,j]
  Greece_unTO[j,1] = as.numeric(To_spillover1)
} 
Greece_unTO <- ts(Greece_unTO, frequency = 12, start = c(2008, 2))


# UK, weighted
UK_TO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_weighted_To[24,,j]
  UK_TO[j,1] = as.numeric(To_spillover1)
} 
UK_TO <- ts(UK_TO, frequency = 12, start = c(2008, 2))

UK_FROM <- matrix(NA, 161, 1)
for (j in 1:161){
  From_spillover1 <- GVD_nor_unweighted_From[24,,j]
  UK_FROM[j,1] = as.numeric(From_spillover1)
} 
UK_FROM <- ts(UK_FROM, frequency = 12, start = c(2008, 2))

UK_unTO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_unweighted_To[24,,j]
  UK_unTO[j,1] = as.numeric(To_spillover1)
} 
UK_unTO <- ts(UK_unTO, frequency = 12, start = c(2008, 2))


# Russia, weighted
Russia_TO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_weighted_To[21,,j]
  Russia_TO[j,1] = as.numeric(To_spillover1)
} 
Russia_TO <- ts(Russia_TO, frequency = 12, start = c(2008, 2))

Russia_FROM <- matrix(NA, 161, 1)
for (j in 1:161){
  From_spillover1 <- GVD_nor_unweighted_From[21,,j]
  Russia_FROM[j,1] = as.numeric(From_spillover1)
} 
Russia_FROM <- ts(Russia_FROM, frequency = 12, start = c(2008, 2))

Russia_unTO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_unweighted_To[21,,j]
  Russia_unTO[j,1] = as.numeric(To_spillover1)
} 
Russia_unTO <- ts(Russia_unTO, frequency = 12, start = c(2008, 2))


# SK, weighted
SK_TO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_weighted_To[19,,j]
  SK_TO[j,1] = as.numeric(To_spillover1)
} 
SK_TO <- ts(SK_TO, frequency = 12, start = c(2008, 2))

SK_FROM <- matrix(NA, 161, 1)
for (j in 1:161){
  From_spillover1 <- GVD_nor_unweighted_From[19,,j]
  SK_FROM[j,1] = as.numeric(From_spillover1)
} 
SK_FROM <- ts(SK_FROM, frequency = 12, start = c(2008, 2))

SK_unTO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_unweighted_To[19,,j]
  SK_unTO[j,1] = as.numeric(To_spillover1)
} 
SK_unTO <- ts(SK_unTO, frequency = 12, start = c(2008, 2))


# Japan, weighted
Japan_TO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_weighted_To[18,,j]
  Japan_TO[j,1] = as.numeric(To_spillover1)
} 
Japan_TO <- ts(Japan_TO, frequency = 12, start = c(2008, 2))

Japan_FROM <- matrix(NA, 161, 1)
for (j in 1:161){
  From_spillover1 <- GVD_nor_unweighted_From[18,,j]
  Japan_FROM[j,1] = as.numeric(From_spillover1)
} 
Japan_FROM <- ts(Japan_FROM, frequency = 12, start = c(2008, 2))

png("desktop/var(2)/dynamic/Japan.png", width = 4854, height = 3000, res = 350)
plot(Japan_TO, type = "l", las = 2, lty=1, ylab = '', main = 'Japan Spillover')
lines(Japan_FROM, lty=3)
legend("topright", legend=c("TO","FROM"), lty=c(1,3))
abline(v=2008.75, col="red")
text(x=2009,y=240,"Sep,2008")
text(x=2009,y=230,"Lehman's collapse")
abline(v=2019.583333, col="red")
text(x=2019,y=270,"Jul,2019")
text(x=2019,y=260,"Japan–South Korea trade dispute")
abline(v=2020.166667, col="red")
text(x=2021,y=250,"Feb,2020")
text(x=2021,y=240,"Covid-19")
axis(1, at = c(2007:2022), las = 2)
dev.off()

# Mexico, weighted
Mexico_TO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_weighted_To[20,,j]
  Mexico_TO[j,1] = as.numeric(To_spillover1)
} 
Mexico_TO <- ts(Mexico_TO, frequency = 12, start = c(2008, 2))

Mexico_FROM <- matrix(NA, 161, 1)
for (j in 1:161){
  From_spillover1 <- GVD_nor_unweighted_From[20,,j]
  Mexico_FROM[j,1] = as.numeric(From_spillover1)
} 
Mexico_FROM <- ts(Mexico_FROM, frequency = 12, start = c(2008, 2))

Mexico_unTO <- matrix(NA, 161, 1)
for (j in 1:161){
  To_spillover1 <- GVD_nor_unweighted_To[20,,j]
  Mexico_unTO[j,1] = as.numeric(To_spillover1)
} 
Mexico_unTO <- ts(Mexico_unTO, frequency = 12, start = c(2008, 2))

#pairwise
GVD_nor_alasso_UStoCN <- ts(GVD_nor_alasso_rw[7,1,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_CNtoUS <- ts(GVD_nor_alasso_rw[1,7,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_UStoMX <- ts(GVD_nor_alasso_rw[20,1,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_MXtoUS <- ts(GVD_nor_alasso_rw[1,20,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_UStoRU <- ts(GVD_nor_alasso_rw[21,1,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_RUtoUS <- ts(GVD_nor_alasso_rw[1,21,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_CNtoHK <- ts(GVD_nor_alasso_rw[14,7,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_HKtoCN <- ts(GVD_nor_alasso_rw[7,14,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_CNtoSK <- ts(GVD_nor_alasso_rw[19,7,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_SKtoCN <- ts(GVD_nor_alasso_rw[7,19,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_GRtoGE <- ts(GVD_nor_alasso_rw[13,12,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_GEtoGR <- ts(GVD_nor_alasso_rw[12,13,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_GRtoSP <- ts(GVD_nor_alasso_rw[22,12,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_SPtoGR <- ts(GVD_nor_alasso_rw[12,22,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_GRtoFR <- ts(GVD_nor_alasso_rw[11,12,], frequency = 12, start = c(2008, 2))
GVD_nor_alasso_FRtoGR <- ts(GVD_nor_alasso_rw[12,11,], frequency = 12, start = c(2008, 2))

#find white noise
GVD_nor_alasso_CroatiatoUS <- ts(GVD_nor_alasso_rw[7,13,], frequency = 12, start = c(2008, 2))

png("desktop/var(2)/dynamic/pairwise/pairwise_gr_cn_lag2.png", width = 4854, height = 3000, res = 350)
plot(GVD_nor_alasso_CroatiatoUS, type = "l",  main = 'Greece-China', las = 2, ylab = '')
axis(1, at = c(2007:2022), las = 2)
dev.off()
