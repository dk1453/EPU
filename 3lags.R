# generate the weight vector
GDP <- read_excel("/Users/kangdi/Desktop/typ/data/GDP.xlsx")
GDP <- GDP[-1,-1]
static_weight <- colMeans(GDP)
static_weight_nor <- (t(static_weight) / apply(t(static_weight), 1, sum))
static_weight_nor <- t(static_weight_nor)

epu_3lags = embed(epu,4)
dependent3_epu = epu_3lags[,1:24]
independent3_epu = epu_3lags[,25:96]

ols3_cv <- VAR(epu, p = 3, ic = 'AIC')
coef_ols3_cv <- Bcoef(ols3_cv)[,-73]
one3 <- matrix(1, ncol = 72, nrow = 24)
alasso_weight3_cv <- one3 / coef_ols3_cv

bestlambda3 <- matrix(NA, 24, 1)
for (i in 1:24){
  cvfit3 <- cv.glmnet(x = independent3_epu, y = dependent3_epu[,i], penalty.factor = alasso_weight3_cv[i,], alpha = 1)
  bestlambda3[i,1] <- cvfit3$lambda.min
}
# bestlambda
# plot(bestlambda)
mean(bestlambda3) # 5.135885


coef_alasso <- matrix(NA, 24, 73)
res_alasso <- NULL
for (i in 1:24)
{
  model_alasso <- glmnet(x = independent3_epu, y = dependent3_epu[,i], lambda = 5.135885, alpha = 1, penalty.factor = alasso_weight3_cv[i,]) # in penalty.factor, the intercept is asked to be removed
  coef_list <- coef(model_alasso)
  coef_alasso[i, ] = as.matrix(coef_list)
  pred <- predict(model_alasso, newx= epu_3lags[,25:96])
  res <- epu_3lags[ ,i] - pred
  res_alasso <- cbind(res_alasso, res)
  Q_alasso <- (t(res_alasso)%*%res_alasso)/nrow(res_alasso)
}
GVD_nor_preweighted3_static <- GFEVD(coef_alasso[,-1], Q_alasso, n.ahead=6, standardize = T)$GFEVD * 100
# write.csv(GVD_nor_preweighted3_static,"desktop/var(2)/GVD_nor_preweighted3_static.csv")

# TO
static3_TO <- Calculate_To(GVD_nor_preweighted3_static, matrix(static_weight_nor[,1])) #成功
# write.csv(static3_TO,"desktop/var(2)/static3_TO.csv")

# FROM
static3_ <- cnct_analysis(GVD_nor_preweighted3_static) # use cnct_analysis from the ols.R
static3_TO_unweighted <- matrix(static3_$To)
static3_FROM <- matrix(static3_$From)
# write.csv(static3_FROM,"desktop/var(2)/static3_FROM.csv")
# write.csv(static3_TO_unweighted,"desktop/var(2)/static3_TO_unweighted.csv")

# TOTAL
static3_matrix_weighted <- give_weight_updated(GVD_nor_preweighted3_static, static_weight_nor, 1)
# write.csv(static3_matrix_weighted,"desktop/var(2)/weighted3_GVD.csv")
static3_TOTAL <- cnct_analysis_updated(static3_matrix_weighted)$Sum * 24

# quadrants
lags3_static_quadrants <- read_excel("~/Desktop/var(2)/3lags_static_quadrant.xlsx")
png("desktop/var(2)/3lags_static_quadrants.png", width = 3570, height = 1865, res = 350)
ggplot(lags3_static_quadrants, aes(y=From, x=To, label=country)) +
  geom_text_repel(aes(y=From, x=To, label=country)) +
  geom_point() +
  scale_y_continuous(breaks=seq(0,100,10),limits=c(0,100)) + 
  scale_x_continuous(breaks=seq(0,200,10), limits=c(0,200)) +
  # #lims(y=c(70,110),x=c(-50,250)) +
  #coord_fixed() +  # 
  geom_vline(xintercept = mean(lags3_static_quadrants$To)) + geom_hline(yintercept = mean(lags3_static_quadrants$From)) + 
  geom_abline(intercept = -2.881531)
dev.off()

lags3_static_quadrants_unweighted <- read_excel("~/Desktop/var(2)/3lags_static_quadrant_unweighted.xlsx")
png("desktop/var(2)/lags3_static_quadrants_unweighted.png", width = 3570, height = 1865, res = 350)
ggplot(lags3_static_quadrants_unweighted, aes(y=From, x=To, label=country)) +
  geom_text_repel(aes(y=From, x=To, label=country)) +
  geom_point() +
  scale_y_continuous(breaks=seq(0,100,10),limits=c(0,100)) + 
  scale_x_continuous(breaks=seq(0,200,10), limits=c(0,200)) +
  # #lims(y=c(70,110),x=c(-50,250)) +
  #coord_fixed() +  # 
  geom_vline(xintercept = mean(lags3_static_quadrants_unweighted$To)) + geom_hline(yintercept = mean(lags3_static_quadrants_unweighted$From)) + 
  geom_abline(intercept = 0)
dev.off()