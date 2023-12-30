#static_weight <- read_excel("~/Desktop/typ/data/static_weight.xlsx")# cannot use it as it is derived from EWMA
# View(static_weight)
# View(cnct_alasso_preweight)
# View(epu)



test_to <- Calculate_To(GVD_alasso_nor_preweight, matrix(static_weight_nor[,1])) #成功
test_from <- cnct_alasso_st$From

test_matrix <- give_weight_updated(GVD_alasso_unnor_preweighted, static_weight_nor, 1)
test_total <- cnct_analysis_updated(test_matrix)$Sum

test_quadrants <- read_excel("~/Desktop/test.xlsx")
png("desktop/test_quadrants.png", width = 3570, height = 1865, res = 350)
ggplot(test_quadrants, aes(y=From, x=To, label=country)) +
  geom_text_repel(aes(y=From, x=To, label=country)) +
  geom_point() +
  scale_y_continuous(breaks=seq(0,100,10),limits=c(0,100)) + 
  scale_x_continuous(breaks=seq(0,200,10), limits=c(0,200)) +
  # #lims(y=c(70,110),x=c(-50,250)) +
  #coord_fixed() +  # 
  geom_vline(xintercept = mean(test_quadrants$To)) + geom_hline(yintercept = mean(test_quadrants$From)) + 
  geom_abline(intercept = 0)
dev.off()


# write.csv(GVD_alasso_nor_preweight,"desktop/GVD_alasso_nor_preweight.csv")



### PREPARE THE DATA 
GVD_alasso_nor_preweight <- cnct_alasso[1:24,1:24] #normalized GVD
GVD_alasso_unnor_preweighted <- GFEVD(coef_alasso[,-1], Q_alasso, n.ahead=6, normalize=F)$GFEVD * 100 #unnormalized GVD

### weight the GVD
unnormalized_weighted <- give_weight(GVD_alasso_unnor_preweighted, static_weight, 1)
normalized_weighted <- give_weight(GVD_alasso_nor_preweight, static_weight, 1)

#### PLOT QUADRANTS
# write.csv(unnormalized_weighted,"desktop/unnormalized_weighted.csv")
# write.csv(normalized_weighted, "desktop/normalized_weighted.csv")

Nor_Wei_quadrants <- read_excel("~/Desktop/Nor_Wei_quadrant.xlsx")
png("desktop/Nor_Wei_quadrants.png", width = 3570, height = 1865, res = 350)
ggplot(Nor_Wei_quadrants, aes(y=From, x=To, label=country)) +
  geom_text_repel(aes(y=From, x=To, label=country)) +
  geom_point() +
  scale_y_continuous(breaks=seq(0,0.5,0.1),limits=c(0,0.4)) + 
  scale_x_continuous(breaks=seq(0,0.5,0.1), limits=c(0,0.15)) +
  # #lims(y=c(70,110),x=c(-50,250)) +
  #coord_fixed() +  # 
  geom_vline(xintercept = mean(Nor_Wei_quadrants$To)) + geom_hline(yintercept = mean(Nor_Wei_quadrants$From)) + 
  geom_abline(intercept = 0)
dev.off()

Unnor_wei_quadrants <- read_excel("~/Desktop/Unnor_wei_quadrant.xlsx")
png("desktop/Unnor_wei_quadrants.png", width = 3570, height = 1865, res = 350)
ggplot(Unnor_wei_quadrants, aes(y=From, x=To, label=country)) +
  geom_text_repel(aes(y=From, x=To, label=country)) +
  geom_point() +
  scale_y_continuous(breaks=seq(0,0.5,0.1),limits=c(0,0.4)) + 
  scale_x_continuous(breaks=seq(0,0.5,0.1), limits=c(0,0.15)) +
  # #lims(y=c(70,110),x=c(-50,250)) +
  #coord_fixed() +  # 
  geom_vline(xintercept = mean(Unnor_wei_quadrants$To)) + geom_hline(yintercept = mean(Unnor_wei_quadrants$From)) + 
  geom_abline(intercept = 0)
dev.off()

nor_unwei_quadrants <- read_excel("~/Desktop/nor_unwei_quadrant.xlsx")
png("desktop/nor_unwei_quadrants.png", width = 3570, height = 1865, res = 350)
ggplot(nor_unwei_quadrants, aes(y=From, x=To, label=country)) +
  geom_text_repel(aes(y=From, x=To, label=country)) +
  geom_point() +
  scale_y_continuous(breaks=seq(0,0.5,0.1),limits=c(0,0.1)) + 
  scale_x_continuous(breaks=seq(0,0.5,0.1), limits=c(0,0.1)) +
  # #lims(y=c(70,110),x=c(-50,250)) +
  #coord_fixed() +  # 
  geom_vline(xintercept = mean(nor_unwei_quadrants$To)) + geom_hline(yintercept = mean(nor_unwei_quadrants$From)) + 
  geom_abline(intercept = 0)
dev.off()

### Rolling Analysis
#this part is for normalized GVD, to do unnormalized GVD, rerun line 57 to 91 in glmnet_alasso, set line 84 to normalize=F.
rolling_weight_GDP2 <- read_excel("~/Desktop/typ/data/GDP_interpolate1.xlsx") # GDP from July 02 to June 21, sample size 228
rolling_weight_GDP3 <- rolling_weight_GDP2[7:228, 1:24] # weight started from Jan 03
Moving_Average_GDP <- weight_moving_average(rolling_weight_GDP3,60)


rolling_weight_GDP_nor <- rolling_weight_GDP / apply(rolling_weight_GDP, 1, sum)

normalized_weighted_rw <- give_weight_rw(cnct_alasso_rw, rolling_weight_GDP)
normalized_weighted_dy <- apply(normalized_weighted_rw, 3, cnct_analysis_updated)

# 计算rolling window的To
arr_rolling_weight_GDP <- array(t(rolling_weight_GDP_nor), dim=c(nrow(t(rolling_weight_GDP_nor)), 1, ncol(t(rolling_weight_GDP_nor))))

arraytest_To <- Rolling_Calculate_To(cnct_alasso_rw, arr_rolling_weight_GDP) #成功
### rolling window 的 from
arraytest_From <- array(NA, dim = c(24, 1, 172))
from_spillover_test <- matrix(NA, 24, 1)
for (j in 1:172){
  from_spillover_test <- matrix(cnct_alasso_dy[[j]]$From)
  arraytest_From[,,j] <- from_spillover_test
} 

### ploting quadrants for the first and last window
normalized_weighted_firstwindow <- cnct_analysis_updated(normalized_weighted_rw[,,1])
normalized_weighted_lastwindow <- cnct_analysis_updated(normalized_weighted_rw[,,172])
normalized_weighted_FirstWindow <- cbind(normalized_weighted_firstwindow$From, normalized_weighted_firstwindow$To)
normalized_weighted_LastWindow <- cbind(normalized_weighted_lastwindow$From, normalized_weighted_lastwindow$To)
# write.csv(normalized_weighted_FirstWindow,"desktop/normalized_weighted_FirstWindow.csv")
# write.csv(normalized_weighted_LastWindow, "desktop/normalized_weighted_LastWindow.csv")

# plot
normalized_weighted_FirstWindow_quadrants <- read_excel("~/Desktop/normalized_weighted_FirstWindow.xlsx")
png("desktop/normalized_weighted_FirstWindow_quadrants.png", width = 3579, height = 1865, res = 350)
ggplot(normalized_weighted_FirstWindow_quadrants, aes(y=From, x=To, label=country)) +
  geom_text_repel(aes(y=From, x=To, label=country)) +
  geom_point() +
  scale_y_continuous(breaks=seq(0,0.5,0.1),limits=c(0,0.4)) + 
  scale_x_continuous(breaks=seq(0,0.5,0.1), limits=c(0,0.15)) +
  geom_vline(xintercept = mean(normalized_weighted_FirstWindow_quadrants$To)) + geom_hline(yintercept = mean(normalized_weighted_FirstWindow_quadrants$From)) + 
  geom_abline(intercept = 0)
dev.off()

normalized_weighted_LastWindow_quadrants <- read_excel("~/Desktop/normalized_weighted_LastWindow.xlsx")
png("desktop/normalized_weighted_LastWindow_quadrants.png", width = 3570, height = 1865, res = 350)
ggplot(normalized_weighted_LastWindow_quadrants, aes(y=From, x=To, label=country)) +
  geom_text_repel(aes(y=From, x=To, label=country)) +
  geom_point() +
  scale_y_continuous(breaks=seq(0,0.5,0.1),limits=c(0,0.4)) + 
  scale_x_continuous(breaks=seq(0,0.5,0.1), limits=c(0,0.15)) +
  geom_vline(xintercept = mean(normalized_weighted_LastWindow_quadrants$To)) + geom_hline(yintercept = mean(normalized_weighted_LastWindow_quadrants$From)) + 
  geom_abline(intercept = 0)
dev.off()

### ploting quadrants for the windows after GFC, first window used as before GFC
normalized_weighted_after <- cnct_analysis_updated(normalized_weighted_rw[,,23])
normalized_weighted_After <- cbind(normalized_weighted_after$From, normalized_weighted_after$To)
# write.csv(normalized_weighted_After,"desktop/normalized_weighted_After.csv")

# plot
normalized_weighted_After_quadrants <- read_excel("~/Desktop/normalized_weighted_After.xlsx")
png("desktop/normalized_weighted_After_quadrants.png", width = 3579, height = 1865, res = 350)
ggplot(normalized_weighted_After_quadrants, aes(y=From, x=To, label=country)) +
  geom_text_repel(aes(y=From, x=To, label=country)) +
  geom_point() +
  scale_y_continuous(breaks=seq(0,0.5,0.1),limits=c(0,0.4)) + 
  scale_x_continuous(breaks=seq(0,0.5,0.1), limits=c(0,0.15)) +
  geom_vline(xintercept = mean(normalized_weighted_After_quadrants$To)) + geom_hline(yintercept = mean(normalized_weighted_After_quadrants$From)) + 
  geom_abline(intercept = 0)
dev.off()

