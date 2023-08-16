#Script for real-time forecasting

# To cite:  
# Cathal Mills (2023)
# URL: https://github.com/cathalmills/peru_dengue

#Real-time forecasting involves for months in 2018 to 2021,
  # omission of a month's current and future data,
  # both surveillance and climate data,
  # and forecasting the incidence using the learned model
  # and climate data with a lag of between 1-4 months.

#References: 
## Rachel Lowe (2021) - https://github.com/drrachellowe/hydromet_dengue


#Extract last month of 2017 (last time point for climatic data used to forecast)
first_time_2018 <- head(ptl_dt[which(YEAR == 8), ]$TIME, 1)-1
head(ptl_dt)
first_time_2018
for(i in first_time_2018:(nrow(ptl_dt)/3 - 1)){
  print(i)
  time_ids <- seq(i*3+1, i*3+3)
  print(time_ids)
  s <- 5000
  rt_forecast_dt <- data.table(ptl_dt)
  rt_forecast_dt[, IND:= seq(1, nrow(rt_forecast_dt))]
  true_cases <- c(rt_forecast_dt[which(TIME == i + 1), ]$CASES)
  true_dirs <- c(rt_forecast_dt[which(TIME == i + 1), ]$DIR_POP_INTERP)
  pop_offsets <- c(rt_forecast_dt[which(TIME == i + 1), ]$POP_OFFSET_INTERP)
  
  rt_forecast_dt[which(TIME == i + 1), CASES:= NA]
  rt_forecast_dt[which(TIME == i + 1), DIR_POP_INTERP:= NA]
  rt_forecast_dt[which(TIME == i + 1), DIR:= NA]
  max_ind <- max(rt_forecast_dt[which(TIME == i + 1),]$IND)
  rt_forecast_dt <- rt_forecast_dt[c(1:max_ind),]
  # print(time_ids)
  
  
  maximum_climate_lag <- 4
  rt_forecast_lag_tmax <- tsModel::Lag(climate_dt$tmax, 
                                               group = climate_dt$REGION, 
                                               k = 1:maximum_climate_lag)
  rt_forecast_lag_tmax <- rt_forecast_lag_tmax[13:nrow(rt_forecast_lag_tmax),]
  rt_forecast_lag_tmax <- rt_forecast_lag_tmax[1:nrow(rt_forecast_dt),]
  
  # Minimum temperature (Tmin)
  rt_forecast_lag_tmin <- tsModel::Lag(climate_dt$tmin, 
                                               group = climate_dt$REGION, 
                                               k = 1:maximum_climate_lag)
  rt_forecast_lag_tmin <- rt_forecast_lag_tmin[13:nrow(rt_forecast_lag_tmin),]
  rt_forecast_lag_tmin <- rt_forecast_lag_tmin[1:nrow(rt_forecast_dt),]
  
  # rt_forecast_lag_tmin <- rt_forecast_lag_tmin[climate_dt$TIME > 4,]
  # rt_forecast_lag_tmin
  
  # Precipitation (prec)
  rt_forecast_lag_prec <- tsModel::Lag(climate_dt$prec, 
                                               group = climate_dt$REGION, 
                                               k = 1:maximum_climate_lag)
  rt_forecast_lag_prec <- rt_forecast_lag_prec[13:nrow(rt_forecast_lag_prec),]
  rt_forecast_lag_prec <- rt_forecast_lag_prec[1:nrow(rt_forecast_dt),]
  
  # rt_forecast_lag_prec <- rt_forecast_lag_prec[climate_dt$TIME > 4,]
  # dim(rt_forecast_lag_prec)
  # SPI-6
  rt_forecast_lag_spi <- tsModel::Lag(climate_dt$SPI_6, 
                                              group = climate_dt$REGION, 
                                              k = 1:maximum_climate_lag)
  rt_forecast_lag_spi <- rt_forecast_lag_spi[13:nrow(rt_forecast_lag_spi),]
  rt_forecast_lag_spi <- rt_forecast_lag_spi[1:nrow(rt_forecast_dt),]
  
  # rt_forecast_lag_spi <- rt_forecast_lag_spi[climate_dt$TIME > 4, ]
  
  
  
  # rt_forecast_lag_oni <- rt_forecast_lag_oni[climate_dt$TIME > 4,]
  
  rt_forecast_lag_icen <- tsModel::Lag(climate_dt$E_INDEX, 
                                               k = 1:maximum_climate_lag)
  rt_forecast_lag_icen <- rt_forecast_lag_icen[13:nrow(rt_forecast_lag_icen),]
  rt_forecast_lag_icen <- rt_forecast_lag_icen[1:nrow(rt_forecast_dt),]
  dim(rt_forecast_lag_icen)
  lagknot <- 2
  #Cross basis matrices----
  #Via the cross basis (a bi-dimensional functional space) we are specifying simultaneously the 
  # relationships in the dimensions of the predictor and lags, respectively.
  tmax_basis <- crossbasis(rt_forecast_lag_tmax, 
                                   argvar=list(fun = "bs", knots = equalknots(climate_dt$tmax, 2)),
                                   arglag=list(fun="bs", knots =lagknot))
  print(dim(tmax_basis))
  tmin_basis <- crossbasis(rt_forecast_lag_tmin, 
                                   argvar=list(fun = "bs", knots = equalknots(climate_dt$tmin, 2)),
                                   arglag=list(fun="bs", knots = lagknot))
  prec_basis <- crossbasis(rt_forecast_lag_prec, 
                                   argvar=list(fun = "bs", knots = equalknots(climate_dt$prec, 2)),
                                   arglag=list(fun="bs", knots = lagknot))
  spi_basis <- crossbasis(rt_forecast_lag_spi, 
                                  argvar=list(fun = "bs", knots = equalknots(climate_dt$SPI_6, 2)),
                                  arglag=list(fun="bs", knots = lagknot))
  icen_basis <- crossbasis(rt_forecast_lag_icen,  
                                   argvar=list(fun = "bs", knots = equalknots(climate_dt$E_INDEX, 2)),
                                   arglag=list(fun="bs", knots = lagknot))
  
  colnames(tmax_basis) = paste0("tmax_basis.", colnames(tmax_basis))
  colnames(prec_basis) = paste0("prec_basis.", colnames(prec_basis))
  colnames(icen_basis) = paste0("icen_basis.", colnames(icen_basis))
  colnames(spi_basis) = paste0("spi_basis.", colnames(spi_basis))
  
  baseline_rt_forecast_formula <- CASES ~ 1 + 
    f(MONTH, replicate = RGN_IND,
      model = "rw1", cyclic = TRUE, constr = TRUE,
      scale.model = TRUE,  hyper = prior.prec)+
    f(RGN_IND, model = "bym2",
      hyper = prior.prec, scale.model = TRUE,
      graph = file.path(peru.inla.data.in.dir, "nbr_piura_tumbes_lambayeque.graph")
    )+
    f(YEAR, replicate = RGN_IND, model = "iid")+
    SEASON+RSI_DIR_POP_INTERP_LAG
  rt_forecast_climate_formula_tmax_spi_prec_icen <- update.formula(baseline_rt_forecast_formula, ~. + tmax_basis + prec_basis+
                                                                         icen_basis + spi_basis)
  
  tmp_climate_cv_fit <- run_model_func(data = rt_forecast_dt, formula = rt_forecast_climate_formula_tmax_spi_prec_icen)
  xx <- inla.posterior.sample(s, tmp_climate_cv_fit)
  xx.s <- inla.posterior.sample.eval(function(...) c(theta[1], Predictor[time_ids]), xx)

  y.pred <- matrix(NA, 3, s)
  error.pred <- matrix(NA, 3, s)
  abs_error.pred <- matrix(NA, 3, s)
  dir.pred <- matrix(NA, 3, s)
  dir_error.pred <- matrix(NA, 3, s)
  dir_abs_error.pred <- matrix(NA, 3, s)

  for(s.idx in 1:s) {
    xx.sample <- xx.s[, s.idx]
    y.pred[, s.idx] <- rnbinom(3, mu = exp(xx.sample[-1]), size = xx.sample[1])
    error.pred[, s.idx] <- true_cases - y.pred[, s.idx]
    abs_error.pred[, s.idx] <- abs(true_cases - y.pred[, s.idx])
    
    dir.pred[, s.idx] <- y.pred[, s.idx]/pop_offsets
    dir_error.pred[, s.idx] <- true_dirs - dir.pred[, s.idx]
    dir_abs_error.pred[, s.idx] <- abs(true_dirs - dir.pred[, s.idx])
    
  }
  #Probabilistic probability of surpassing outbreak thresholds:
  outbreak_prob <- rep(0, 3)
  for(j in 1:nrow(dir.pred)){
    outbreak_prob[j] <-  length(which(dir.pred[j,] >= 100))/s
  }
  
  outbreak_prob_50 <- rep(0, 3)
  for(j in 1:nrow(dir.pred)){
    outbreak_prob_50[j] <-  length(which(dir.pred[j,] >= 50))/s
  }
  
  outbreak_prob_150 <- rep(0, 3)
  for(j in 1:nrow(dir.pred)){
    outbreak_prob_150[j] <-  length(which(dir.pred[j,] >= 150))/s
  }
  
  outbreak_prob_200 <- rep(0, 3)
  for(j in 1:nrow(dir.pred)){
    outbreak_prob_200[j] <-  length(which(dir.pred[j,] >= 200))/s
  }
  outbreak_prob_250 <- rep(0, 3)
  for(j in 1:nrow(dir.pred)){
    outbreak_prob_250[j] <-  length(which(dir.pred[j,] >= 250))/s
  }
  
  tmp <- ptl_dt[, list(UPPER_DIR = quantile(DIR_POP_INTERP, 0.95)), by = "REGION"]
  upper_vals <- tmp$UPPER_DIR
  
  upper_quantile_threshold_prob <- rep(0, 3)
  for(j in 1:nrow(dir.pred)){
    upper_quantile_threshold_prob[j] <-  length(which(dir.pred[j,] >= upper_vals[j]))/s
  }
  
  morans_i <- c(moran.test(apply(y.pred, 1, quantile ,probs = c(0.5)),
                           nb2listw(nbr_piura_tumbes), randomisation =  FALSE)$p.value)
  rt_forecast_preds_dt <- data.table(TIME = rep(i + 1, each = 3),
                                     MORANS = rep(morans_i, each = 3),
                                     CI_L = apply(y.pred, 1, quantile, probs = c(0.025)),
                                     CI_U = apply(y.pred, 1, quantile, probs = c(0.975)),
                                     MEDIAN = apply(y.pred, 1, quantile ,probs = c(0.5)),
                                     MEAN = apply(y.pred, 1, mean),
                                     ERROR_CI_L = apply(error.pred, 1, quantile, probs = c(0.025)),
                                     ERROR_CI_U = apply(error.pred, 1, quantile, probs = c(0.975)),
                                     ERROR_MEDIAN = apply(error.pred, 1, quantile ,probs = c(0.5)),
                                     ERROR_MEAN = apply(error.pred, 1, mean),
                                     ABS_ERROR_CI_L = apply(abs_error.pred, 1, quantile, probs = c(0.025)),
                                     ABS_ERROR_CI_U = apply(abs_error.pred, 1, quantile, probs = c(0.975)),
                                     ABS_ERROR_MEDIAN = apply(abs_error.pred, 1, quantile ,probs = c(0.5)),
                                     ABS_ERROR_MEAN = apply(abs_error.pred, 1, mean),
                                     
                                     OUTBREAK_PROB = outbreak_prob,
                                     OUTBREAK_PROB_50 = outbreak_prob_50,
                                     OUTBREAK_PROB_150 = outbreak_prob_150,
                                     OUTBREAK_PROB_200 = outbreak_prob_200,
                                     OUTBREAK_PROB_250 = outbreak_prob_250,
                                     UPPER_QUANTILE_THRESHOLD_PROB = upper_quantile_threshold_prob,
                                     
                                     DIR_CI_L = apply(dir.pred, 1, quantile, probs = c(0.025)),
                                     DIR_CI_U = apply(dir.pred, 1, quantile, probs = c(0.975)),
                                     DIR_MEDIAN = apply(dir.pred, 1, quantile ,probs = c(0.5)),
                                     DIR_MEAN = apply(dir.pred, 1, mean),
                                     
                                     DIR_ERROR_CI_L = apply(dir_error.pred, 1, quantile, probs = c(0.025)),
                                     DIR_ERROR_CI_U = apply(dir_error.pred, 1, quantile, probs = c(0.975)),
                                     DIR_ERROR_MEDIAN = apply(dir_error.pred, 1, quantile ,probs = c(0.5)),
                                     DIR_ERROR_MEAN = apply(dir_error.pred, 1, mean),
                                     DIR_ABS_ERROR_CI_L = apply(dir_abs_error.pred, 1, quantile, probs = c(0.025)),
                                     DIR_ABS_ERROR_CI_U = apply(dir_abs_error.pred, 1, quantile, probs = c(0.975)),
                                     DIR_ABS_ERROR_MEDIAN = apply(dir_abs_error.pred, 1, quantile ,probs = c(0.5)),
                                     DIR_ABS_ERROR_MEAN = apply(dir_abs_error.pred, 1, mean))
  print("Is the error here?")
  # ls[[fold]] <- rt_forecast_preds_dt
  saveRDS(rt_forecast_preds_dt, file = file.path(peru.inla.data.out.dir, paste0("rt_forecast_2018_2021_piura_tumbes_lambayeque_preds_dt", i, ".RDS")))
}



#Bind forecasted month's DIRs into single data.table
ovr_rt_forecast_preds_dt <- NULL
for(i in (first_time_2018):(nrow(ptl_dt)/3 - 1)){
  rt_forecast_preds_dt <- readRDS(file = file.path(peru.inla.data.out.dir, paste0("rt_forecast_2018_2021_piura_tumbes_lambayeque_preds_dt", i, ".RDS")))
  
  ovr_rt_forecast_preds_dt <- rbind(ovr_rt_forecast_preds_dt, rt_forecast_preds_dt)
}

rt_forecast_piura_tumbes_posterior_pred_dt <- plot_cv_posterior_preds(ovr_rt_forecast_preds_dt, ptl_dt[which(TIME > first_time_2018)])[[1]]
rt_forecast_piura_tumbes_posterior_pred_dt[, length(which((CASES) >= q2.5 & (CASES) <= q97.5))/length(CASES)]
rt_forecast_piura_tumbes_posterior_pred_dt[, length(which((CASES) > q97.5))/length(CASES)]
rt_forecast_piura_tumbes_posterior_pred_dt[, mean(q97.5 - q2.5)]
cor(rt_forecast_piura_tumbes_posterior_pred_dt$MEDIAN,
    rt_forecast_piura_tumbes_posterior_pred_dt$CASES)
R2(rt_forecast_piura_tumbes_posterior_pred_dt$MEDIAN,
   rt_forecast_piura_tumbes_posterior_pred_dt$CASES)
MAE(rt_forecast_piura_tumbes_posterior_pred_dt$MEDIAN,
    rt_forecast_piura_tumbes_posterior_pred_dt$CASES)
rt_forecast_piura_tumbes_posterior_pred_dt[, MAE(MEDIAN, CASES), by = "REGION"]


dir_rt_forecast_piura_tumbes_posterior_pred_dt <- plot_dir_cv_posterior_preds(ovr_rt_forecast_preds_dt, ptl_dt[which(TIME > first_time_2018)])[[1]]
dir_rt_forecast_piura_tumbes_posterior_pred_dt[, length(which((DIR) >= q2.5 & (DIR) <= q97.5))/length(DIR)]
dir_rt_forecast_piura_tumbes_posterior_pred_dt[, length(which((DIR) > q97.5))/length(DIR)]
head(dir_rt_forecast_piura_tumbes_posterior_pred_dt[order(DIR, decreasing = TRUE)], 10)
cor(dir_rt_forecast_piura_tumbes_posterior_pred_dt$MEDIAN,
    dir_rt_forecast_piura_tumbes_posterior_pred_dt$DIR)
R2(dir_rt_forecast_piura_tumbes_posterior_pred_dt$MEDIAN,
   dir_rt_forecast_piura_tumbes_posterior_pred_dt$DIR)
MAE(dir_rt_forecast_piura_tumbes_posterior_pred_dt$MEDIAN,
    dir_rt_forecast_piura_tumbes_posterior_pred_dt$DIR)
dir_rt_forecast_piura_tumbes_posterior_pred_dt[, MAE(MEDIAN, DIR), by = "REGION"]
dir_rt_forecast_piura_tumbes_posterior_pred_dt[, length(which(ABS_ERROR_MEDIAN <= 50))/
                                                 nrow(dir_rt_forecast_piura_tumbes_posterior_pred_dt)]
hist_median_rt_forecast_errors <- ggplot(dir_rt_forecast_piura_tumbes_posterior_pred_dt, aes(x = ABS_ERROR_MEDIAN))+
  geom_histogram(fill = "steelblue", color = "black", boundary = 0, closed = "left",
                 bins = 15)+theme_bw()+
  labs(x = "Posterior Median Absolute Error", y = "Frequency")+
  theme(text = element_text(size = 28),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        # panel.grid.minor.y = element_blank(),
        # panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 30, hjust = 0.5, face = "bold"),
        plot.subtitle = element_blank(),
        axis.title=element_text(size=25), 
        legend.text=element_text(size=25),
        legend.position = "none", legend.key.size = unit(1.5, 'cm'))
ggsave(hist_median_rt_forecast_errors, file = file.path(peru.inla.data.out.dir, "hist_median_rt_forecast_errors.pdf"), h = 12, w = 18)


#Confusion Matrices for Outbreak Detection ----
dir_rt_forecast_piura_tumbes_posterior_pred_dt[, OUTBREAK_150_PRED:= ifelse(OUTBREAK_PROB_150 >= 0.5, 1, 0)]
dir_rt_forecast_piura_tumbes_posterior_pred_dt
dir_rt_forecast_piura_tumbes_posterior_pred_dt[, OUTBREAK_150_TRUE:= ifelse(DIR >= 150, 1, 0)]

confusion_mat_rt_forecast_150 <- confusionMatrix(data=as.factor(dir_rt_forecast_piura_tumbes_posterior_pred_dt$OUTBREAK_150_PRED), 
                                                 reference = as.factor(dir_rt_forecast_piura_tumbes_posterior_pred_dt$OUTBREAK_150_TRUE),
                                                 positive = "1")
confusion_mat_rt_forecast_150$byClass
accuracy_rt_forecast_150 <- unname(confusion_mat_rt_forecast_150$overall[1])
precision_rt_forecast_150 <- unname(confusion_mat_rt_forecast_150$byClass[3])

true_pos_rt_forecast_150 <- unname(confusion_mat_rt_forecast_150$byClass[1]) #Sensitivity/Hit Rate
fpr_rt_forecast_150 <- 1-unname(confusion_mat_rt_forecast_150$byClass[2]) #1-TNR = 1-Specificity
true_neg_rt_forecast_150 <- unname(confusion_mat_rt_forecast_150$byClass[2]) #Sensitivity/Hit Rate
fnr_rt_forecast_150 <- 1 - unname(confusion_mat_rt_forecast_150$byClass[1]) #Miss Rate
roc_rt_forecast_150 <- roc(dir_rt_forecast_piura_tumbes_posterior_pred_dt$OUTBREAK_150_TRUE, 
                           dir_rt_forecast_piura_tumbes_posterior_pred_dt$OUTBREAK_150_PRED)
auc_rt_forecast_150 <- auc(roc_rt_forecast_150)
auc_rt_forecast_150
sens_spec_rt_forecast_150_table <-data.table("Accuracy" = accuracy_rt_forecast_150,
                                             "Precision" = precision_rt_forecast_150,
                                             "True Positive" = true_pos_rt_forecast_150,
                                             "False Positive" = fpr_rt_forecast_150,
                                             "True Negative Rate" = true_neg_rt_forecast_150,
                                             "False Negative Rate" = fnr_rt_forecast_150,
                                             "AUC" = auc_rt_forecast_150)
sens_spec_rt_forecast_150_table



dir_rt_forecast_piura_tumbes_posterior_pred_dt[, OUTBREAK_50_PRED:= ifelse(OUTBREAK_PROB_50 >= 0.5, 1, 0)]
dir_rt_forecast_piura_tumbes_posterior_pred_dt
dir_rt_forecast_piura_tumbes_posterior_pred_dt[, OUTBREAK_50_TRUE:= ifelse(DIR >= 50, 1, 0)]

confusion_mat_rt_forecast_50 <- confusionMatrix(data=as.factor(dir_rt_forecast_piura_tumbes_posterior_pred_dt$OUTBREAK_50_PRED), 
                                                reference = as.factor(dir_rt_forecast_piura_tumbes_posterior_pred_dt$OUTBREAK_50_TRUE),
                                                positive = "1")
confusion_mat_rt_forecast_50$byClass
accuracy_rt_forecast_50 <- unname(confusion_mat_rt_forecast_50$overall[1])
precision_rt_forecast_50 <- unname(confusion_mat_rt_forecast_50$byClass[3])
true_pos_rt_forecast_50 <- unname(confusion_mat_rt_forecast_50$byClass[1]) #Sensitivity/Hit Rate
fpr_rt_forecast_50 <- 1-unname(confusion_mat_rt_forecast_50$byClass[2]) #1-TNR = 1-Specificity
true_neg_rt_forecast_50 <- unname(confusion_mat_rt_forecast_50$byClass[2]) #Sensitivity/Hit Rate
fnr_rt_forecast_50 <- 1 - unname(confusion_mat_rt_forecast_50$byClass[1]) #Miss Rate
roc_rt_forecast_50 <- roc(dir_rt_forecast_piura_tumbes_posterior_pred_dt$OUTBREAK_50_TRUE, 
                          dir_rt_forecast_piura_tumbes_posterior_pred_dt$OUTBREAK_50_PRED)
auc_rt_forecast_50 <- auc(roc_rt_forecast_50)
auc_rt_forecast_50
rt_forecast_dt
sens_spec_rt_forecast_50_table <-data.table("Accuracy" = accuracy_rt_forecast_50,
                                            "Precision" = precision_rt_forecast_50,
                                            "True Positive" = true_pos_rt_forecast_50,
                                            "False Positive" = fpr_rt_forecast_50,
                                            "True Negative Rate" = true_neg_rt_forecast_50,
                                            "False Negative Rate" = fnr_rt_forecast_50,
                                            "AUC" = auc_rt_forecast_50)
xtable(rbind(sens_spec_rt_forecast_50_table, sens_spec_rt_forecast_150_table))



#Greater than 50
outbreak_50_rt_forecast_dt <- dir_rt_forecast_piura_tumbes_posterior_pred_dt[which(DIR >= 50), ]
outbreak_50_rt_forecast_dt[, IND:= seq(1, nrow(outbreak_50_rt_forecast_dt))]
max_ind <- max(outbreak_50_rt_forecast_dt$IND)
outbreak_50_rt_forecast_dt[which(OUTBREAK_PROB_50 <= 0.5)]
outbreak_50_rt_forecast_plot <- ggplot(outbreak_50_rt_forecast_dt, aes(x = IND, y = DIR))+
  geom_col(aes(fill = OUTBREAK_PROB_50), width = 0.7)+
  geom_hline(aes(yintercept = 50),
             linetype = "longdash", linewidth = 2.0,
             colour = "forest green")+
  coord_flip()+
  scale_x_continuous(breaks = seq(0, max_ind, by = 2))+
  scale_y_continuous(breaks = seq(0, 1000, by = 100))+
  geom_text(aes(label = label_percent(0.1)(OUTBREAK_PROB_50)),
            hjust = 0.55, fontface = "bold",
            size = 8)+
  scale_fill_gradient(
    name = expression("Prob(DIR" >= "50) "),
    low = "#FEED99",
    high = "#AF3301",
    guide = "colourbar",
    breaks = c(0.0, 0.25, 0.5, 0.75),
    labels = percent(0.25*0:3))+theme_bw()+
  labs(x = "Outbreak Number")+
  theme(text = element_text(size = 28),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 30, hjust = 0.5, face = "bold"),
        plot.subtitle = element_blank(),
        axis.title=element_text(size=25), 
        legend.text=element_text(size=25),
        legend.position = "bottom", legend.key.size = unit(1.5, 'cm'))

outbreak_50_rt_forecast_plot
ggsave(outbreak_50_rt_forecast_plot, 
       file = file.path(peru.inla.data.out.dir, "outbreak_50_rt_forecast_plot.pdf"), h = 14, w = 20)
ggsave(outbreak_50_rt_forecast_plot, 
       file = file.path(peru.inla.data.out.dir, "outbreak_50_rt_forecast_plot.png"), h = 14, w = 20)



#Greater than 150
outbreak_150_rt_forecast_dt <- dir_rt_forecast_piura_tumbes_posterior_pred_dt[which(DIR >= 150), ]
outbreak_150_rt_forecast_dt[, IND:= seq(1, nrow(outbreak_150_rt_forecast_dt))]
outbreak_150_rt_forecast_dt[which(OUTBREAK_PROB_150 <= 0.5)]
max_ind <- max(outbreak_150_rt_forecast_dt$IND)
outbreak_150_rt_forecast_plot
outbreak_150_rt_forecast_plot <- ggplot(outbreak_150_rt_forecast_dt, aes(x = IND, y = DIR))+
  geom_col(aes(fill = OUTBREAK_PROB_150), width = 0.7)+
  geom_hline(aes(yintercept = 150),
             linetype = "longdash", linewidth = 2.0,
             colour = "forest green")+
  coord_flip()+
  scale_x_continuous(breaks = seq(0, max_ind, by = 2))+
  scale_y_continuous(breaks = seq(0, 1000, by = 100))+
  geom_text(aes(label = label_percent(0.1)(OUTBREAK_PROB_150)),
            hjust = 0.85, fontface = "bold",
            size = 8)+
  scale_fill_gradient(
    name = expression("Prob(DIR" >= "150) "),
    low = "#FEED99",
    high = "#AF3301",
    guide = "colourbar",
    breaks = c(0.0, 0.25, 0.5, 0.75),
    labels = percent(0.25*0:3))+theme_bw()+
  labs(x = "Outbreak Number")+
  theme(text = element_text(size = 28),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 30, hjust = 0.5, face = "bold"),
        plot.subtitle = element_blank(),
        axis.title=element_text(size=25), 
        legend.text=element_text(size=25),
        legend.position = "bottom", legend.key.size = unit(1.5, 'cm'))
outbreak_150_rt_forecast_plot
grid.arrange(outbreak_50_rt_forecast_plot,
             outbreak_150_rt_forecast_plot, nrow = 1)
ggsave(outbreak_150_rt_forecast_plot, 
       file = file.path(peru.inla.data.out.dir, "outbreak_150_rt_forecast_plot.pdf"), h = 14, w = 20)
ggsave(outbreak_150_rt_forecast_plot, 
       file = file.path(peru.inla.data.out.dir, "outbreak_150_rt_forecast_plot.png"), h = 14, w = 20)

outbreak_50_150_rt_forecast_grid <- ggarrange(plotlist = list(outbreak_50_rt_forecast_plot, outbreak_150_rt_forecast_plot), nrow = 1, ncol = 2,common.legend = FALSE,
                                              legend = "bottom")
ggsave(outbreak_50_150_rt_forecast_grid, 
       file = file.path(peru.inla.data.out.dir, "outbreak_50_150_rt_forecast_grid.pdf"), h = 14, w = 24)
ggsave(outbreak_50_150_rt_forecast_grid, 
       file = file.path(peru.inla.data.out.dir, "outbreak_50_150_rt_forecast_grid.png"), h = 14, w = 24)

false_outbreak_50_rt_forecast_dt
false_outbreak_50_rt_forecast_dt <- dir_rt_forecast_piura_tumbes_posterior_pred_dt[which(DIR < 50), ]
mean(false_outbreak_50_rt_forecast_dt$OUTBREAK_PROB_50)
max_false_outbreak_50_rt_forecast_prob <- 
  false_outbreak_50_rt_forecast_dt[, max(OUTBREAK_PROB_50)]
max_false_outbreak_50_rt_forecast_prob
false_outbreak_50_rt_forecast_dt[which.max(OUTBREAK_PROB_50),]

max_false_outbreak_50_rt_forecast_prob
false_outbreak_50_rt_forecast_dt[, IND:= seq(1, nrow(false_outbreak_50_rt_forecast_dt))]
false_outbreak_50_rt_forecast_plot <- ggplot(false_outbreak_50_rt_forecast_dt, aes(x = IND, y = DIR))+
  geom_col(aes(fill = OUTBREAK_PROB_50), width = 0.7)+
  coord_flip()+
  scale_x_continuous(breaks = seq(0, 170, by = 10))+
  scale_y_continuous(breaks = seq(0, 60, by = 10),
                     limits = c(0, 50))+
  # geom_text(aes(label = label_percent()(OUTBREAK_PROB)),
  #           hjust = 1.75, fontface = "bold",
  #           size = 7)+
  scale_fill_gradient(
    name = expression("Prob(DIR" >= "50) "),
    low = "#FEED99",
    high = "#AF3301",
    guide = "colourbar",
    breaks = c(0.0, 0.25, 0.5, 0.75),
    labels = percent(0.25*0:3))+theme_bw()+
  labs(x = "Non-Outbreak Number")+
  theme(text = element_text(size = 28),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 30, hjust = 0.5, face = "bold"),
        plot.subtitle = element_blank(),
        axis.title=element_text(size=25), 
        legend.text=element_text(size=25),
        legend.position = "right", legend.key.size = unit(1.5, 'cm'))
false_outbreak_50_rt_forecast_plot
ggsave(false_outbreak_50_rt_forecast_plot, 
       file = file.path(peru.inla.data.out.dir, "false_outbreak_50_rt_forecast_plot.pdf"), h = 14, w = 20)
ggsave(false_outbreak_50_rt_forecast_plot, 
       file = file.path(peru.inla.data.out.dir, "false_outbreak_50_rt_forecast_plot.png"), h = 14, w = 20)


false_outbreak_150_rt_forecast_dt
false_outbreak_150_rt_forecast_dt <- dir_rt_forecast_piura_tumbes_posterior_pred_dt[which(DIR < 150), ]
mean(false_outbreak_150_rt_forecast_dt$OUTBREAK_PROB_150)

max_false_outbreak_150_rt_forecast_prob <- 
  false_outbreak_150_rt_forecast_dt[, max(OUTBREAK_PROB_150)]
max_false_outbreak_150_rt_forecast_prob
false_outbreak_150_rt_forecast_dt[, IND:= seq(1, nrow(false_outbreak_150_rt_forecast_dt))]
false_outbreak_150_rt_forecast_dt$IND
false_outbreak_150_rt_forecast_plot <- ggplot(false_outbreak_150_rt_forecast_dt, aes(x = IND, y = DIR))+
  geom_col(aes(fill = OUTBREAK_PROB_150), width = 0.7)+
  coord_flip()+
  scale_x_continuous(breaks = seq(0, 170, by = 10))+
  scale_y_continuous(breaks = seq(0, 150, by = 15),
                     limits = c(0, 150))+
  scale_fill_gradient(
    name = expression("Prob(DIR" >= "150) "),
    low = "#FEED99",
    high = "#AF3301",
    guide = "colourbar",
    breaks = c(0.0, 0.25, 0.5, 0.75),
    labels = percent(0.25*0:3))+theme_bw()+
  labs(x = "Non-Outbreak Number")+
  theme(text = element_text(size = 28),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(size = 30, hjust = 0.5, face = "bold"),
        plot.subtitle = element_blank(),
        axis.title=element_text(size=25), 
        legend.text=element_text(size=25),
        legend.position = "right", legend.key.size = unit(1.5, 'cm'))
false_outbreak_150_rt_forecast_plot
ggsave(false_outbreak_150_rt_forecast_plot, 
       file = file.path(peru.inla.data.out.dir, "false_outbreak_150_rt_forecast_plot.pdf"), h = 14, w = 20)
ggsave(false_outbreak_150_rt_forecast_plot, 
       file = file.path(peru.inla.data.out.dir, "false_outbreak_150_rt_forecast_plot.png"), h = 14, w = 20)

false_outbreak_50_150_rt_forecast_grid <- ggarrange(false_outbreak_50_rt_forecast_plot, false_outbreak_150_rt_forecast_plot, nrow = 1)
ggsave(false_outbreak_50_150_rt_forecast_grid, 
       file = file.path(peru.inla.data.out.dir, "false_outbreak_50_150_rt_forecast_grid.pdf"), h = 14, w = 24)
ggsave(false_outbreak_50_150_rt_forecast_grid, 
       file = file.path(peru.inla.data.out.dir, "false_outbreak_50_150_rt_forecast_grid.png"), h = 14, w = 24)




#Various visualisations of forecasted values vs observed values
rt_forecast_piura_tumbes_posterior_pred_plot <- plot_cv_posterior_preds(ovr_rt_forecast_preds_dt, ptl_dt[which(TIME > first_time_2018)])[[2]]
rt_forecast_piura_tumbes_posterior_pred_plot

dir_rt_forecast_piura_tumbes_posterior_pred_plot <- plot_dir_cv_posterior_preds(ovr_rt_forecast_preds_dt, ptl_dt[which(TIME > first_time_2018)])[[2]]
dir_rt_forecast_piura_tumbes_posterior_pred_plot
dir_rt_forecast_piura_tumbes_posterior_pred_plot
ggsave(dir_rt_forecast_piura_tumbes_posterior_pred_plot, 
       file = file.path(peru.inla.data.out.dir, "dir_rt_forecast_piura_tumbes_posterior_pred_plot.pdf"), h = 12, w = 20)
ggsave(dir_rt_forecast_piura_tumbes_posterior_pred_plot, 
       file = file.path(peru.inla.data.out.dir, "dir_rt_forecast_piura_tumbes_posterior_pred_plot.png"), h = 12, w = 20)

facet_rt_forecast_piura_tumbes_posterior_pred_plot <- plot_cv_posterior_preds(ovr_rt_forecast_preds_dt, ptl_dt[which(TIME > first_time_2018)])[[3]]
facet_rt_forecast_piura_tumbes_posterior_pred_plot
dir_facet_rt_forecast_piura_tumbes_posterior_pred_plot <- plot_dir_cv_posterior_preds(ovr_rt_forecast_preds_dt, ptl_dt[which(TIME > first_time_2018)])[[3]]
dir_facet_rt_forecast_piura_tumbes_posterior_pred_plot

errorbar_facet_rt_forecast_piura_tumbes_posterior_pred_plot <- plot_cv_posterior_preds(ovr_rt_forecast_preds_dt, ptl_dt[which(TIME > first_time_2018)])[[7]]
errorbar_facet_rt_forecast_piura_tumbes_posterior_pred_plot

dir_errorbar_facet_rt_forecast_piura_tumbes_posterior_pred_plot <- plot_dir_cv_posterior_preds(ovr_rt_forecast_preds_dt, ptl_dt[which(TIME > first_time_2018)])[[7]]
dir_errorbar_facet_rt_forecast_piura_tumbes_posterior_pred_plot

facet_no_ci_rt_forecast_piura_tumbes_posterior_pred_plot <- plot_cv_posterior_preds(ovr_rt_forecast_preds_dt, ptl_dt[which(TIME > first_time_2018)])[[4]]
facet_no_ci_rt_forecast_piura_tumbes_posterior_pred_plot
