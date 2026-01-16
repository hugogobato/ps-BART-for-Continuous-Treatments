# Calculate SE and AE for 'empirical_ATE_cover' and 'theoretical_CATE_cover'
ps_bart$SE_empirical_ATE_cover <- (0.95 - ps_bart$empirical_ATE_cover)^2
ps_bart$AE_empirical_ATE_cover <- abs(0.95 - ps_bart$empirical_ATE_cover)

ps_bart$SE_theoretical_CATE_cover <- (0.95 - ps_bart$theoretical_CATE_cover)^2
ps_bart$AE_theoretical_CATE_cover <- abs(0.95 - ps_bart$theoretical_CATE_cover)

# Calculate the mean and standard deviation for each column
mean_values <- colMeans(ps_bart, na.rm = TRUE)
std_values <- apply(ps_bart, 2, sd, na.rm = TRUE)

print(ps_bart)
