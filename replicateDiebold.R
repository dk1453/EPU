volatility <- read_excel("~/Desktop/typ/data/replication/StockVolatility.xlsx")
volatility <- volatility[, 2:14]
volatility <- as.matrix(volatility)
######DIEBOLD'S#########
### generate the ols coefficients ###
vol_ols <- VAR(volatility, p = 3)
vol_coef_ols <- Bcoef(vol_ols)[,-40]

### connectedness ###
vol_res_ols <- resid(vol_ols) #each country a column
mean(vol_res_ols ^ 2) #mse of the static analysis
Q_ols <- (t(vol_res_ols)%*%vol_res_ols)/nrow(vol_res_ols)
dim(Q_ols) # dim of Q has to be consistant with the dim of coef(withou const)

vol_cnct_ols <- GFEVD(vol_coef_ols, Q_ols, n.ahead=12, standardize = T)$GFEVD * 100 # std or not, no difference, * 100 to let row sum to be 100, it's normalized because it's default, this is why std doesn't have effect.
vol_cnct_ols_nonnormalize <- GFEVD(vol_coef_ols, Q_ols, n.ahead=12, normalize=F, standardize = F)$GFEVD * 100
total_sum1 <- c(sum(vol_cnct_ols_nonnormalize), apply(vol_cnct_ols_nonnormalize, 1, sum))
rowSums(vol_cnct_ols_nonnormalize)
colSums(vol_cnct_ols_nonnormalize)
write.csv(vol_cnct_ols_nonnormalize,"desktop/vol_cnct_ols_nonnormalize1.csv")

View(vol_cnct_ols)
rowSums(vol_cnct_ols)
colSums(vol_cnct_ols)
## static analysis ##
vol_cnct_ols_st <- cnct_analysis(vol_cnct_ols)
vol_cnct_ols <- cbind(vol_cnct_ols, vol_cnct_ols_st$From)
vol_cnct_ols <- rbind(vol_cnct_ols, append(vol_cnct_ols_st$To, 0), append(vol_cnct_ols_st$Net, 0))
vol_cnct_ols[nrow(vol_cnct_ols),ncol(vol_cnct_ols)] <- vol_cnct_ols_st$Sum
# write.csv(vol_cnct_ols,"desktop/typ/rep/vol_cnct_ols.csv")

######MINE#########
##generate x and y##
volatility_lag3 = embed(volatility,4)
dependent = volatility_lag3[,1:13]
independent = volatility_lag3[,14:52]

## cross validation for time series
vol_ols_cv <- VAR(volatility, p = 3)
vol_coef_ols_cv <- Bcoef(vol_ols_cv)[,-40]
one <- matrix(1, ncol = 39, nrow = 13)
vol_alasso_weight_cv <- one / vol_coef_ols_cv
K = 10 # 10 fold
size_train = floor(nrow(volatility)/K) - 1
mse_total <- NULL
for (t in 1:K)
{
  mse_by_train_group <- NULL
  for (j in seq(from = 0, to = 5, by = 0.1))
  {
    mse_by_hyper <- NULL
    vol_coef_alasso_cv <- matrix(NA,13,40)
    for (i in 1:ncol(volatility))
    {
      # print("Before glmnet call")  # Add this line
      model_alasso_cv <- glmnet(x = independent[(1:(t*size_train)),], y = dependent[(1:(t*size_train)), i], lambda = exp(j), alpha = 1, penalty.factor = vol_alasso_weight_cv[i,])
      # print("After glmnet call")  # Add this line
      coef_list_cv <- coef(model_alasso_cv) 
      vol_coef_alasso_cv[i,] = as.matrix(coef_list_cv)
      pred_cv <- predict(model_alasso_cv, newx= independent[(t*size_train+1):(t*size_train+10),])
      res_cv <- dependent[(t*size_train+1):(t*size_train+10),i] - pred_cv
      mse_individual_cv <- sum(res_cv^2) / 10
      mse_by_hyper <- rbind(mse_by_hyper, mse_individual_cv)
    }
    mse_cv <- mean(mse_by_hyper)
    mse_by_train_group <- rbind(mse_by_train_group,mse_cv)
  }
  mse_total <- cbind(mse_total,mse_by_train_group)
}
means_across_group <- rowMeans(mse_total)
# print(means_across_group)
smallest_location <- which(means_across_group == min(means_across_group))
print(smallest_location) # here it is 9, so optimal lambda is exp(0.9)
##here I got the optimal lambda wrong, this is the lambda for rolling analysis, the static analysis should be like this:

bestlambda <- matrix(NA, 13, 1)
for (i in 1:13){
  cvfit <- cv.glmnet(x = independent, y = dependent[,i], penalty.factor = vol_alasso_weight_cv[i,], alpha = 1)
  bestlambda[i,1] <- cvfit$lambda.min
}
# bestlambda
# plot(bestlambda)
mean(bestlambda) #4.813555 以下部分更新后没run

### generate the alasso coefficients ###
vol_coef_ols <- vol_coef_ols_cv
one <- matrix(1, ncol = 39, nrow = 13)
vol_alasso_weight <- one / vol_coef_ols

vol_coef_alasso <- matrix(NA, 13, 40)
res_alasso <- NULL
for (i in 1:13)
{
  model_alasso <- glmnet(x = independent, y = dependent[,i], lambda = 4.823555, alpha = 1, penalty.factor = vol_alasso_weight[i,]) # in penalty.factor, the intercept is asked to be removed
  coef_list <- coef(model_alasso)
  vol_coef_alasso[i, ] = as.matrix(coef_list)
  pred <- predict(model_alasso, newx= independent)
  res <- dependent[ ,i] - pred
  res_alasso <- cbind(res_alasso, res)
  Q_alasso <- (t(res_alasso)%*%res_alasso)/nrow(res_alasso)
}

### connectedness ###
vol_cnct_alasso <- GFEVD(vol_coef_alasso[,-1], Q_alasso, n.ahead=12, standardize = T)$GFEVD * 100
# View(cnct_alasso)
rowSums(vol_cnct_alasso)
colSums(vol_cnct_alasso)

## static analysis ##
vol_cnct_alasso_st <- cnct_analysis(vol_cnct_alasso) # use cnct_analysis from the ols.R
vol_cnct_alasso <- cbind(vol_cnct_alasso, vol_cnct_alasso_st$From)
vol_cnct_alasso <- rbind(vol_cnct_alasso, append(vol_cnct_alasso_st$To, 0), append(vol_cnct_alasso_st$Net, 0))

vol_cnct_alasso_st$Sum
vol_cnct_alasso[nrow(vol_cnct_alasso),ncol(vol_cnct_alasso)] <- vol_cnct_alasso_st$Sum
View(vol_cnct_alasso)

market_cap <- rbind(98.5, 186,32,175,86,170,62.5,23,53.5,129,95.5,30.15,23.95)
vol_cnct_alasso_preweighted <- vol_cnct_alasso[1:13,1:13]
vol_cnct_a_weighted <- give_weight(vol_cnct_alasso_preweighted, market_cap, 1)

# write.csv(vol_cnct_a_weighted,"desktop/typ/rep/vol_cnct_a_weighted.csv")

