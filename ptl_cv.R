#Script for leave-one-time-point-out cross-validation

# To cite:  
# Cathal Mills (2023)
# URL: https://github.com/cathalmills/peru_dengue

#Leave-one-time-point-out cross-validation involves omission of
  # single month's (three) observations and deriving posterior predictive
  # estimates for the month

#References: 
## Rachel Lowe (2021) - https://github.com/drrachellowe/hydromet_dengue

baseline_formula <- CASES ~ 1 + 
  f(MONTH, replicate = RGN_IND,
    model = "rw1", cyclic = TRUE, constr = TRUE,
    scale.model = TRUE,  hyper = prior.prec)+
  f(YEAR, replicate = RGN_IND, model = "iid")+
  f(RGN_IND, model = "bym2",
    hyper = prior.prec, scale.model = TRUE,
    graph = file.path(peru.inla.data.in.dir, "nbr_piura_tumbes_lambayeque.graph")
  )+
  SEASON+RSI_DIR_POP_INTERP_LAG
tmax_basis <- crossbasis(lag_natural_tmax, 
                                 argvar=list(fun = "bs", knots = equalknots(climate_dt$tmax, 2)),
                                 arglag=list(fun="bs", knots =lagknot))
prec_basis <- crossbasis(lag_natural_prec, 
                                 argvar=list(fun = "bs", knots = equalknots(climate_dt$prec, 2)),
                                 arglag=list(fun="bs", knots = lagknot))
icen_basis <- crossbasis(lag_natural_icen, 
                                 argvar=list(fun = "bs", knots = equalknots(climate_dt$E_INDEX, 2)),
                                 arglag=list(fun="bs", knots = lagknot))
spi_basis <- crossbasis(lag_natural_spi, 
                                argvar=list(fun = "bs", knots = equalknots(climate_dt$SPI_6, 2)),
                                arglag=list(fun="bs", knots = lagknot))
colnames(tmax_basis)
colnames(tmax_basis) = paste0("tmax_basis.", colnames(tmax_basis))
colnames(prec_basis) = paste0("prec_basis.", colnames(prec_basis))
colnames(icen_basis) = paste0("icen_basis.", colnames(icen_basis))
colnames(spi_basis) = paste0("spi_basis.", colnames(spi_basis))

climate_formula_tmax_spi_prec_icen <- update.formula(baseline_formula, ~. + tmax_basis + prec_basis+
                                                           icen_basis + spi_basis)


run_model_func <- function(data, formula = climate_formula_tmax_spi_prec_icen){
  model <- inla(formula = formula, 
                data = data, family = "nbinomial", offset = log(POP_OFFSET_INTERP),
                verbose = FALSE,
                control.inla = list(strategy = 'adaptive'), 
                control.compute = list(waic = TRUE, dic = TRUE, 
                                       cpo = TRUE, config = TRUE,
                                       return.marginals = TRUE),
                control.fixed = list(correlation.matrix = TRUE, 
                                     prec.intercept = 1, prec = 1),
                control.predictor = list(link = 1, compute = TRUE)       
  )
  model <- inla.rerun(model)
  return(model)
}



for(i in 1:(nrow(ptl_dt)/3 - 1)){
  time_ids <- seq(i*3+1, i*3+3)
  
  s <- 5000 #Number of posterior samples
  
  #Set up data.table for cross-validation
  cv_dt <- data.table(ptl_dt)
  setkeyv(cv_dt, c("TIME", "REGION"))
  cv_dt[, IND:= seq(1, nrow(cv_dt))]
  
  #Extract true cases and DIRs for comparison to posterior predictions
  true_cases <- c(cv_dt[which(TIME == i + 1), ]$CASES)
  true_dirs <- c(cv_dt[which(TIME == i + 1), ]$DIR_POP_INTERP)
  pop_offsets <- c(cv_dt[which(TIME == i + 1), ]$POP_OFFSET_INTERP)
  
  #Set Cases to NA to ensure we "leave" them out 
  cv_dt[which(TIME == i + 1), CASES:= NA] 
  setkeyv(cv_dt, c("TIME", "REGION"))
  tmp_climate_cv_fit <- run_model_func(data = cv_dt, formula = climate_nat_rgn_formula_tmax_spi_prec_icen)
  xx <- inla.posterior.sample(s, tmp_climate_cv_fit)
  xx.s <- inla.posterior.sample.eval(function(...) c(theta[1], Predictor[time_ids]), xx)
  
  #Set up matrices for posterior predictions
  #3 rows correspond to departments and s=5000 columns
  y.pred <- matrix(NA, 3, s)
  error.pred <- matrix(NA, 3, s)
  abs_error.pred <- matrix(NA, 3, s)
  dir.pred <- matrix(NA, 3, s)
  dir_error.pred <- matrix(NA, 3, s)
  dir_abs_error.pred <- matrix(NA, 3, s)
  # morans_error.pred <- matrix(NA, 1, s)
  for(s.idx in 1:s) {
    xx.sample <- xx.s[, s.idx]
    y.pred[, s.idx] <- rnbinom(3, mu = exp(xx.sample[-1]), size = xx.sample[1])
    error.pred[, s.idx] <- true_cases - y.pred[, s.idx]
    abs_error.pred[, s.idx] <- abs(true_cases - y.pred[, s.idx])
    
    dir.pred[, s.idx] <- y.pred[, s.idx]/pop_offsets
    dir_error.pred[, s.idx] <- true_dirs - dir.pred[, s.idx]
    dir_abs_error.pred[, s.idx] <- abs(true_dirs - dir.pred[, s.idx])
    
  }
  
  #Probabilistic estimates of breaching DIR thresholds
  
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
  
  #Setting up data.table to save results
  preds_dt <- data.table(TIME = rep(i + 1, each = 3),
                         MORANS = rep(morans_i, each = 3),
                         # MORANS = rep(apply(morans_error.pred, 1, quantile, probs = c(0.025)), each = 3),
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
                         
                         DIR_CI_L = apply(dir.pred, 1, quantile, probs = c(0.025)),
                         DIR_CI_U = apply(dir.pred, 1, quantile, probs = c(0.975)),
                         DIR_MEDIAN = apply(dir.pred, 1, quantile ,probs = c(0.5)),
                         DIR_MEAN = apply(dir.pred, 1, mean),
                         
                         OUTBREAK_PROB = outbreak_prob,
                         OUTBREAK_PROB_50 = outbreak_prob_50,
                         OUTBREAK_PROB_150 = outbreak_prob_150,
                         OUTBREAK_PROB_200 = outbreak_prob_200,
                         OUTBREAK_PROB_250 = outbreak_prob_250,
                         UPPER_QUANTILE_THRESHOLD_PROB = upper_quantile_threshold_prob,
                         
                         
                         DIR_ERROR_CI_L = apply(dir_error.pred, 1, quantile, probs = c(0.025)),
                         DIR_ERROR_CI_U = apply(dir_error.pred, 1, quantile, probs = c(0.975)),
                         DIR_ERROR_MEDIAN = apply(dir_error.pred, 1, quantile ,probs = c(0.5)),
                         DIR_ERROR_MEAN = apply(dir_error.pred, 1, mean),
                         DIR_ABS_ERROR_CI_L = apply(dir_abs_error.pred, 1, quantile, probs = c(0.025)),
                         DIR_ABS_ERROR_CI_U = apply(dir_abs_error.pred, 1, quantile, probs = c(0.975)),
                         DIR_ABS_ERROR_MEDIAN = apply(dir_abs_error.pred, 1, quantile ,probs = c(0.5)),
                         DIR_ABS_ERROR_MEAN = apply(dir_abs_error.pred, 1, mean))
  saveRDS(preds_dt, file = file.path(peru.inla.data.out.dir, paste0("2010_2021_loocv_piura_tumbes_lambayeque_preds_dt", i, ".RDS")))
}


#Bind all LOOCV posterior prediction results into single dt
ovr_preds_dt <- NULL
for(i in 1:(nrow(nat_rgn_inla_df)/3 - 1)){
  preds_dt <- readRDS(file = file.path(peru.inla.data.out.dir, paste0("2010_2021_loocv_piura_tumbes_lambayeque_preds_dt", i, ".RDS")))
  ovr_preds_dt <- rbind(ovr_preds_dt, preds_dt)
}





#Plotting Functions----
#Functions for generating a range of plots for posterior predictions
#1) On scale of cases
plot_cv_posterior_preds <- function(cv_data, data){
  setkeyv(data, c("TIME", "REGION"))
  cv_post_preds <- data.table(MEDIAN = (cv_data$MEDIAN),
                              MEAN = (cv_data$MEAN),
                              q2.5 = (cv_data$CI_L),
                              q97.5 = (cv_data$CI_U),
                              ERROR_MEDIAN = (cv_data$ERROR_MEDIAN),
                              ERROR_MEAN = (cv_data$ERROR_MEAN),
                              ERROR_q2.5 = (cv_data$ERROR_CI_L),
                              ERROR_q97.5 = (cv_data$ERROR_CI_U),
                              ABS_ERROR_MEDIAN = (cv_data$ABS_ERROR_MEDIAN),
                              ABS_ERROR_MEAN = (cv_data$ABS_ERROR_MEAN),
                              ABS_ERROR_q2.5 = (cv_data$ABS_ERROR_CI_L),
                              ABS_ERROR_q97.5 = (cv_data$ABS_ERROR_CI_U),
                              OUTBREAK_PROB = cv_data$OUTBREAK_PROB,
                              DIR = data$DIR,
                              CASES = data$CASES,
                              REGION = data$REGION,
                              TIME = cv_data$TIME,
                              MONTH = data$MONTH,
                              YEAR = data$YEAR,
                              POP = data$POP,
                              IND = seq(1, nrow(data))
  )
  posterior_pred_plot <- ggplot(cv_post_preds)+
    geom_line(aes(x = IND, y = CASES, col = "Observed"))+ 
    geom_ribbon(aes(x = IND, ymin = q2.5, ymax = q97.5), alpha = 0.2)+
    geom_line(aes(x = IND, y = MEDIAN, col = "Estimated"))+
    theme_bw()+
    theme(legend.position = "bottom")+
    ylim(c(0,5000))

  posterior_pred_plot2 <- ggplot(cv_post_preds)+
    geom_line(aes(x = MONTH, y = CASES, col = "Observed"))+ 
    geom_ribbon(aes(x = MONTH, ymin = q2.5, ymax = q97.5), alpha = 0.2)+
    geom_line(aes(x = MONTH, y = MEDIAN, col = "Estimated"))+
    theme_bw()+
    facet_wrap(REGION ~ YEAR, scales = "free_y",
               ncol= 12)+
    theme(legend.position = "bottom")+
    coord_cartesian(xlim = c(1, 12))+
    scale_x_continuous(breaks=seq(1,12))
  
  posterior_pred_plot3 <- ggplot(cv_post_preds)+
    geom_line(aes(x = MONTH, y = CASES, col = "Observed"))+ 
    # geom_ribbon(aes(x = MONTH, ymin = q2.5, ymax = q97.5), alpha = 0.2)+
    geom_line(aes(x = MONTH, y = MEDIAN, col = "Estimated"))+
    theme_bw()+
    facet_wrap(REGION ~ YEAR, scales = "free_y", ncol= 12 )+
    theme(legend.position = "bottom")+
    coord_cartesian(xlim = c(1, 12))+
    scale_x_continuous(breaks=seq(1,12))
  posterior_pred_plot_errorbars <- ggplot(cv_post_preds)+
    geom_point(aes(x = MONTH, y = CASES, col = "Observed"))+ 
    geom_point(aes(x = MONTH, y = MEDIAN, col = "Estimated"))+ 
    geom_errorbar(aes(x = MONTH, ymin = q2.5, ymax = q97.5, col = "Estimated"), alpha = 0.6)+
    theme_bw()+
    facet_wrap(REGION ~ YEAR, scales = "free_y", ncol= 12 )+
    theme(legend.position = "bottom")+
    coord_cartesian(xlim = c(1, 12))+
    scale_x_continuous(breaks=seq(1,12))+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=25),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=25)+geom_text(size = 25),
          legend.position = "bottom",
          legend.title = element_blank())
  tmp_xy_plot <- ggplot(cv_post_preds, aes(x = MEDIAN, y = CASES))+ 
    geom_point(alpha = 0.3, size = 3)+
    geom_abline(intercept = 0, slope = 1)+
    theme_bw()+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=25),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=25)+geom_text(size = 25),
          legend.position = "bottom",
          legend.title = element_blank())+
    coord_cartesian(xlim=c(0, 1000), ylim=c(0, 1000))
  coloured_xy_plot <- ggplot(cv_post_preds, aes(x = MEDIAN, y = CASES, colour = REGION))+ 
    geom_point(alpha = 0.3, size = 3)+
    geom_abline(intercept = 0, slope = 1)+
    theme_bw()+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=25),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=25)+geom_text(size = 25),
          legend.position = "bottom",
          legend.title = element_blank())+
    coord_cartesian(xlim=c(0, 1000), ylim=c(0, 1000))
  
  
  table_errors_by_region <- cv_post_preds[, median(ERROR_MEDIAN), by = "REGION"]
  
  posterior_pred_plot_of_errors <- ggplot(cv_post_preds)+
    # geom_ribbon(aes(x = MONTH, ymin = q2.5, ymax = q97.5), alpha = 0.2)+
    geom_line(aes(x = MONTH, y = MEDIAN, col = REGION))+
    theme_bw()+
    facet_wrap(YEAR ~. , scales = "free_y", nrow= 5 )+
    theme(legend.position = "bottom")+
    coord_cartesian(xlim = c(1, 12))+
    scale_x_continuous(breaks=seq(1,12))
  
  
  two_facet_posterior_pred_plot_errorbars <- ggplot(cv_post_preds)+
    geom_point(aes(x = TIME, y = CASES, col = "Observed"))+ 
    geom_point(aes(x = TIME, y = MEDIAN, col = "Estimated"))+ 
    geom_errorbar(aes(x = TIME, ymin = q2.5, ymax = q97.5, col = "Estimated"), alpha = 0.6)+
    theme_bw()+
    facet_wrap(REGION ~ ., scales = "free_y", ncol= 1)+
    theme(legend.position = "bottom")+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=25),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=25)+geom_text(size = 25),
          legend.position = "bottom",
          legend.title = element_blank())
  
  two_facet_posterior_pred_plot_cis <- ggplot(cv_post_preds)+
    geom_line(aes(x = TIME, y = CASES, col = "Observed"))+ 
    geom_point(aes(x = TIME, y = MEDIAN, col = "Estimated"))+ 
    geom_ribbon(aes(x = TIME, ymin = q2.5, ymax = q97.5), alpha = 0.3)+
    theme_bw()+
    facet_wrap(REGION ~ ., scales = "free_y", ncol= 1)+
    theme(legend.position = "bottom")
  two_facet_posterior_pred_plot_no_cis <- ggplot(cv_post_preds)+
    geom_line(aes(x = TIME, y = CASES, col = "Observed"))+ 
    geom_point(aes(x = TIME, y = MEDIAN, col = "Estimated"))+ 
    theme_bw()+
    facet_wrap(REGION ~ ., scales = "free_y", ncol= 1)+
    theme(legend.position = "bottom")
  
  return(list(cv_post_preds, posterior_pred_plot, posterior_pred_plot2, posterior_pred_plot3,
              tmp_xy_plot, coloured_xy_plot, posterior_pred_plot_errorbars, table_errors_by_region,
              posterior_pred_plot_of_errors, two_facet_posterior_pred_plot_errorbars, two_facet_posterior_pred_plot_cis,
              two_facet_posterior_pred_plot_no_cis))
}

#2) On scale of Dengue Incidence Rates (DIRs)
plot_dir_cv_posterior_preds <- function(cv_data, data){
  setkeyv(data, c("TIME", "REGION"))
  cv_post_preds <- data.table(MEDIAN = (cv_data$DIR_MEDIAN),
                              MEAN = (cv_data$DIR_MEAN),
                              q2.5 = (cv_data$DIR_CI_L),
                              q97.5 = (cv_data$DIR_CI_U),
                              OUTBREAK_PROB = cv_data$OUTBREAK_PROB,
                              OUTBREAK_PROB_50 = cv_data$OUTBREAK_PROB_50,
                              OUTBREAK_PROB_150 = cv_data$OUTBREAK_PROB_150,
                              OUTBREAK_PROB_200 = cv_data$OUTBREAK_PROB_200,
                              OUTBREAK_PROB_250 = cv_data$OUTBREAK_PROB_250,
                              UPPER_QUANTILE_THRESHOLD_PROB  = cv_data$UPPER_QUANTILE_THRESHOLD_PROB,
                              ERROR_MEDIAN = (cv_data$DIR_ERROR_MEDIAN),
                              ERROR_MEAN = (cv_data$DIR_ERROR_MEAN),
                              ERROR_q2.5 = (cv_data$DIR_ERROR_CI_L),
                              ERROR_q97.5 = (cv_data$DIR_ERROR_CI_U),
                              ABS_ERROR_MEDIAN = (cv_data$DIR_ABS_ERROR_MEDIAN),
                              ABS_ERROR_MEAN = (cv_data$DIR_ABS_ERROR_MEAN),
                              ABS_ERROR_q2.5 = (cv_data$DIR_ABS_ERROR_CI_L),
                              ABS_ERROR_q97.5 = (cv_data$DIR_ABS_ERROR_CI_U),
                              # CASES = (nat_rgn_inla_df$CASES/nat_rgn_inla_df$POP*1e5),
                              DIR = data$DIR_POP_INTERP,
                              CASES = data$CASES,
                              REGION = data$REGION,
                              TIME = cv_data$TIME,
                              MONTH = data$MONTH,
                              YEAR = data$YEAR,
                              POP = data$POP,
                              IND = seq(1, nrow(data))
  )
  posterior_pred_plot <- ggplot(cv_post_preds)+
    geom_line(aes(x = IND, y = DIR, col = "Observed"))+ 
    geom_ribbon(aes(x = IND, ymin = q2.5, ymax = q97.5), alpha = 0.2)+
    geom_line(aes(x = IND, y = MEDIAN, col = "Estimated"))+
    theme_bw()+
    theme(legend.position = "bottom")+
    ylim(c(0,1000))
  # facet_wrap(REGION ~ YEAR, scales = "free_y")
  # 
  # posterior_pred_plot <- ggplot(in_sample_post_preds)+
  #   geom_ribbon(aes(x = MONTH, ymin = q2.5, ymax = q97.5, fill = "Estimated"), alpha = 0.4)+
  #   geom_line(aes(x = MONTH, y = MEDIAN, col = "Estimated"))+
  #   geom_line(aes(x = MONTH, y = DIR, col = "Observed"))+
  #   facet_wrap(REGION ~ YEAR, scales = "free_y")
  posterior_pred_plot2 <- ggplot(cv_post_preds)+
    geom_line(aes(x = MONTH, y = DIR, col = "Observed"))+ 
    geom_ribbon(aes(x = MONTH, ymin = q2.5, ymax = q97.5), alpha = 0.2)+
    geom_line(aes(x = MONTH, y = MEDIAN, col = "Estimated"))+
    theme_bw()+
    facet_wrap(REGION ~ YEAR, scales = "free_y",
               ncol= 12)+
    theme(legend.position = "bottom")+
    coord_cartesian(xlim = c(1, 12))+
    scale_x_continuous(breaks=seq(1,12, by = 2))+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 30, hjust = 0.5, face = "bold"),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=20)+geom_text(size = 20),
          legend.position = "bottom",
          legend.title = element_blank())
  
  posterior_pred_plot3 <- ggplot(cv_post_preds)+
    geom_line(aes(x = MONTH, y = DIR, col = "Observed"))+ 
    # geom_ribbon(aes(x = MONTH, ymin = q2.5, ymax = q97.5), alpha = 0.2)+
    geom_line(aes(x = MONTH, y = MEDIAN, col = "Estimated"))+
    theme_bw()+
    facet_wrap(REGION ~ YEAR, scales = "free_y", ncol= 12 )+
    theme(legend.position = "bottom")+
    coord_cartesian(xlim = c(1, 12))+
    scale_x_continuous(breaks=seq(1,12))+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(size = 30, hjust = 0.5, face = "bold"),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=20)+geom_text(size = 20),
          legend.position = "bottom",
          legend.title = element_blank())
  posterior_pred_plot_errorbars <- ggplot(cv_post_preds)+
    geom_point(aes(x = MONTH, y = DIR, col = "Observed"))+ 
    geom_point(aes(x = MONTH, y = MEDIAN, col = "Estimated"))+ 
    geom_errorbar(aes(x = MONTH, ymin = q2.5, ymax = q97.5, col = "Estimated"), alpha = 0.6)+
    theme_bw()+
    facet_wrap(REGION ~ YEAR, scales = "free_y", ncol= 12 )+
    theme(legend.position = "bottom")+
    coord_cartesian(xlim = c(1, 12))+
    scale_x_continuous(breaks=seq(1,12))+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=25),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=25)+geom_text(size = 25),
          legend.position = "bottom",
          legend.title = element_blank())
  tmp_xy_plot <- ggplot(cv_post_preds, aes(x = MEDIAN, y = DIR))+ 
    geom_point(alpha = 0.3, size = 3)+
    geom_abline(intercept = 0, slope = 1)+
    theme_bw()+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=25),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=25)+geom_text(size = 25),
          legend.position = "bottom",
          legend.title = element_blank())+
    coord_cartesian(xlim=c(0, 1000), ylim=c(0, 1000))
  coloured_xy_plot <- ggplot(cv_post_preds, aes(x = MEDIAN, y = DIR, colour = REGION))+ 
    geom_point(alpha = 0.4, size = 8)+
    geom_abline(intercept = 0, slope = 1)+
    theme_bw()+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=25),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=25)+geom_text(size = 25),
          legend.position = "bottom",
          legend.title = element_blank())+
    coord_cartesian(xlim=c(0, 1000), ylim=c(0, 1000))
  
  
  table_errors_by_region <- cv_post_preds[, median(ERROR_MEDIAN), by = "REGION"]
  
  
  posterior_pred_plot_of_errors <- ggplot(cv_post_preds)+
    # geom_ribbon(aes(x = MONTH, ymin = q2.5, ymax = q97.5), alpha = 0.2)+
    geom_line(aes(x = TIME, y = ABS_ERROR_MEDIAN))+
    theme_bw()+
    facet_wrap(REGION ~ . , scales = "free_y", nrow = 3 )+
    theme(legend.position = "bottom")+
    xlab("Month")+ylab("Median Absolute Error")+
    scale_x_continuous(breaks = seq(80, 140, by = 10))+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=25),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=25)+geom_text(size = 25),
          legend.position = "bottom",
          legend.title = element_blank())
  
  
  two_facet_posterior_pred_plot_errorbars <- ggplot(cv_post_preds)+
    geom_point(aes(x = TIME, y = DIR, col = "Observed"))+ 
    geom_point(aes(x = TIME, y = MEDIAN, col = "Estimated"))+ 
    geom_errorbar(aes(x = TIME, ymin = q2.5, ymax = q97.5, col = "Estimated"), alpha = 0.6)+
    theme_bw()+
    facet_wrap(REGION ~ ., scales = "free_y", ncol= 1)+
    theme(legend.position = "bottom")+
    xlab("Month")+
    scale_x_continuous(breaks = seq(0, 140, by = 20))+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=25),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=25)+geom_text(size = 25),
          legend.position = "bottom",
          legend.title = element_blank())
  
  two_facet_posterior_pred_plot_cis <- ggplot(cv_post_preds)+
    geom_point(aes(x = TIME, y = DIR, col = "Observed"))+ 
    geom_line(aes(x = TIME, y = MEDIAN, col = "Estimated"))+ 
    geom_point(aes(x = TIME, y = MEDIAN, col = "Estimated"))+ 
    geom_ribbon(aes(x = TIME, ymin = q2.5, ymax = q97.5), alpha = 0.18)+
    theme_bw()+
    xlab("Month")+
    scale_x_continuous(breaks = seq(0, 140, by = 20))+
    facet_wrap(REGION ~ ., scales = "free_y", ncol= 1)+
    theme(legend.position = "bottom")+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=25),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=25)+geom_text(size = 25),
          legend.position = "bottom",
          legend.title = element_blank())
  
  two_facet_posterior_pred_plot_no_cis <- ggplot(cv_post_preds)+
    geom_line(aes(x = TIME, y = DIR, col = "Observed"))+ 
    geom_point(aes(x = TIME, y = MEDIAN, col = "Estimated"))+ 
    theme_bw()+
    facet_wrap(REGION ~ ., scales = "free_y", ncol= 1)+
    theme(legend.position = "bottom")
  
  
  
  tmp_xy_plot_limited <- ggplot(cv_post_preds, aes(x = MEDIAN, y = DIR))+ 
    geom_point(alpha = 0.3, size = 3)+
    geom_abline(intercept = 0, slope = 1)+
    theme_bw()+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=25),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=25)+geom_text(size = 25),
          legend.position = "bottom",
          legend.title = element_blank())+
    coord_cartesian(xlim=c(0, 300), ylim=c(0, 300))
  coloured_xy_plot_limited <- ggplot(cv_post_preds, aes(x = MEDIAN, y = DIR, colour = REGION))+ 
    geom_point(alpha = 0.3, size = 3)+
    geom_abline(intercept = 0, slope = 1)+
    theme_bw()+
    labs(x = "Estimated DIR", y = "Observed DIR")+
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size=25),
          axis.text.y = element_text(size=25),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_blank(),
          plot.subtitle = element_blank(),
          axis.title=element_text(size=25), 
          legend.text=element_text(size=25)+geom_text(size = 25),
          legend.position = "bottom",
          legend.title = element_blank())+
    coord_cartesian(xlim=c(0, 300), ylim=c(0, 300))
  
  return(list(cv_post_preds, posterior_pred_plot, posterior_pred_plot2, posterior_pred_plot3,
              tmp_xy_plot, coloured_xy_plot, posterior_pred_plot_errorbars, table_errors_by_region,
              posterior_pred_plot_of_errors, two_facet_posterior_pred_plot_errorbars, two_facet_posterior_pred_plot_cis,
              two_facet_posterior_pred_plot_no_cis, tmp_xy_plot_limited, coloured_xy_plot_limited))
}



setkeyv(ptl_dt, c("TIME", "REGION"))
cv_piura_tumbes_posterior_pred_dt <- plot_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[1]]
cv_piura_tumbes_posterior_pred_dt[, length(which((CASES) >= q2.5 & (CASES) <= q97.5))/length(CASES)]
cv_piura_tumbes_posterior_pred_dt[, length(which((CASES) > q97.5))/length(CASES)]
cv_piura_tumbes_posterior_pred_dt[, length(which((CASES) > q97.5))]


cor(cv_piura_tumbes_posterior_pred_dt$MEDIAN,
    cv_piura_tumbes_posterior_pred_dt$CASES)
R2(cv_piura_tumbes_posterior_pred_dt$MEDIAN,
   cv_piura_tumbes_posterior_pred_dt$CASES)
MAE(cv_piura_tumbes_posterior_pred_dt$MEDIAN,
    cv_piura_tumbes_posterior_pred_dt$CASES)
cv_piura_tumbes_posterior_pred_dt[, MAE(MEDIAN, CASES), by = "REGION"]
cv_piura_tumbes_posterior_pred_plot <- plot_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[2]]
cv_piura_tumbes_posterior_pred_plot
ptl_dt[which(TIME > 1)]

dir_cv_piura_tumbes_posterior_pred_dt <- plot_dir_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[1]]
dir_cv_piura_tumbes_posterior_pred_dt[, length(which((DIR) >= q2.5 & (DIR) <= q97.5))/length(DIR)]
cor(dir_cv_piura_tumbes_posterior_pred_dt$MEDIAN,
    dir_cv_piura_tumbes_posterior_pred_dt$DIR)
R2(dir_cv_piura_tumbes_posterior_pred_dt$MEDIAN,
   dir_cv_piura_tumbes_posterior_pred_dt$DIR)
MAE(dir_cv_piura_tumbes_posterior_pred_dt$MEDIAN,
    dir_cv_piura_tumbes_posterior_pred_dt$DIR)


dir_cv_piura_tumbes_posterior_pred_dt[, length(which(ABS_ERROR_MEDIAN <= 50))/
                                        nrow(dir_cv_piura_tumbes_posterior_pred_dt)]
hist_median_cv_errors <- ggplot(dir_cv_piura_tumbes_posterior_pred_dt, aes(x = ABS_ERROR_MEDIAN))+
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
hist_median_cv_errors
ggsave(hist_median_cv_errors, file = file.path(peru.inla.data.out.dir, "hist_median_cv_errors.pdf"), h = 12, w = 18)


facet_cv_piura_tumbes_posterior_pred_plot <- plot_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[3]]
facet_cv_piura_tumbes_posterior_pred_plot
dev.off()
dir_facet_cv_piura_tumbes_posterior_pred_plot <- plot_dir_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[3]]
dir_facet_cv_piura_tumbes_posterior_pred_plot

facet_no_ci_cv_piura_tumbes_posterior_pred_plot <- plot_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[4]]
facet_no_ci_cv_piura_tumbes_posterior_pred_plot

errorbar_facet_cv_piura_tumbes_posterior_pred_plot <- plot_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[7]]
errorbar_facet_cv_piura_tumbes_posterior_pred_plot
gc()
dir_errorbar_facet_cv_piura_tumbes_posterior_pred_plot <- plot_dir_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[7]]
dir_errorbar_facet_cv_piura_tumbes_posterior_pred_plot


two_facet_cv_piura_tumbes_posterior_pred_plot <- plot_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[10]]
two_facet_cv_piura_tumbes_posterior_pred_plot


dir_two_facet_cv_piura_tumbes_posterior_pred_plot <- plot_dir_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[10]]
dir_two_facet_cv_piura_tumbes_posterior_pred_plot <- set_palette(dir_two_facet_cv_piura_tumbes_posterior_pred_plot, "Set2")
dir_two_facet_cv_piura_tumbes_posterior_pred_plot

two_facet_cis_cv_piura_tumbes_posterior_pred_plot <- plot_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[11]]
two_facet_cis_cv_piura_tumbes_posterior_pred_plot

dir_two_facet_cis_cv_piura_tumbes_posterior_pred_plot <- plot_dir_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[11]]
dir_two_facet_cis_cv_piura_tumbes_posterior_pred_plot <- set_palette(dir_two_facet_cis_cv_piura_tumbes_posterior_pred_plot, "Set2")
dir_two_facet_cis_cv_piura_tumbes_posterior_pred_plot
ggsave(dir_two_facet_cis_cv_piura_tumbes_posterior_pred_plot, 
       file = file.path(peru.inla.data.out.dir, "dir_two_facet_cis_cv_piura_tumbes_posterior_pred_plot.pdf"), h = 12, w = 20)
ggsave(dir_two_facet_cis_cv_piura_tumbes_posterior_pred_plot, 
       file = file.path(peru.inla.data.out.dir, "dir_two_facet_cis_cv_piura_tumbes_posterior_pred_plot.png"), h = 12, w = 20)

#Above Line = Under-predict
#Under Line = Over-predict
xy_ci_cv_piura_tumbes_posterior_pred_plot <- plot_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[5]]
xy_ci_cv_piura_tumbes_posterior_pred_plot

dir_xy_ci_cv_piura_tumbes_posterior_pred_plot <- plot_dir_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[5]]
dir_xy_ci_cv_piura_tumbes_posterior_pred_plot

coloured_xy_cv_piura_tumbes_posterior_pred_plot <- plot_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[6]]
coloured_xy_cv_piura_tumbes_posterior_pred_plot

dir_coloured_xy_cv_piura_tumbes_posterior_pred_plot <- plot_dir_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[6]]
dir_coloured_xy_cv_piura_tumbes_posterior_pred_plot
ggsave(dir_coloured_xy_cv_piura_tumbes_posterior_pred_plot, 
       file = file.path(peru.inla.data.out.dir, "dir_coloured_xy_cv_piura_tumbes_posterior_pred_plot.pdf"), h = 12, w = 20)
ggsave(dir_coloured_xy_cv_piura_tumbes_posterior_pred_plot, 
       file = file.path(peru.inla.data.out.dir, "dir_coloured_xy_cv_piura_tumbes_posterior_pred_plot.png"), h = 12, w = 20)
dir_cv_plot_grid <- ggarrange(dir_two_facet_cis_cv_piura_tumbes_posterior_pred_plot,
                              dir_coloured_xy_cv_piura_tumbes_posterior_pred_plot, nrow = 1,
                              widths = c(1.5, 1), heights = c(2,1 ))
dir_cv_plot_grid
ggsave(dir_cv_plot_grid, 
       file = file.path(peru.inla.data.out.dir, "dir_cv_plot_grid.pdf"), h = 10, w = 24)
ggsave(dir_cv_plot_grid, 
       file = file.path(peru.inla.data.out.dir, "dir_cv_plot_grid.png"), h = 10, w = 24)

grid.arrange(dir_two_facet_cis_cv_piura_tumbes_posterior_pred_plot,
             dir_coloured_xy_cv_piura_tumbes_posterior_pred_plot, nrow = 1)
forecasted_errors_table_region <- plot_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[8]]
forecasted_errors_table_region

errors_cv_piura_tumbes_posterior_pred_plot <- plot_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[9]]
errors_cv_piura_tumbes_posterior_pred_plot

dir_errors_cv_piura_tumbes_posterior_pred_plot <- plot_dir_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[9]]
dir_errors_cv_piura_tumbes_posterior_pred_plot

table_errors_cv_piura_tumbes_posterior_pred <- plot_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[8]]
table_errors_cv_piura_tumbes_posterior_pred

dir_table_errors_cv_piura_tumbes_posterior_pred <- plot_dir_cv_posterior_preds(ovr_preds_dt, ptl_dt[which(TIME > 1)])[[8]]
dir_table_errors_cv_piura_tumbes_posterior_pred






#Confusion Matrices for Outbreak Detection ----

#Precision = TP/(TP+FP)
#False Positive Rate = 1 - TNR = 1 - Specificity

#Sensitivity = TPR = % Positives detected
#Specificity = True Negative Rate = % Negatives Detected
dir_cv_piura_tumbes_posterior_pred_dt[, OUTBREAK_100_TRUE:= ifelse(DIR >= 100, 1, 0)]
dir_cv_piura_tumbes_posterior_pred_dt[, OUTBREAK_100_PRED:= ifelse(OUTBREAK_PROB >= 0.5, 1, 0)]
confusion_mat_cv_100 <- confusionMatrix(data=as.factor(dir_cv_piura_tumbes_posterior_pred_dt$OUTBREAK_100_PRED), 
                                        reference = as.factor(dir_cv_piura_tumbes_posterior_pred_dt$OUTBREAK_100_TRUE),
                                        positive = "1")

accuracy_cv_100 <- unname(confusion_mat_cv_100$overall[1])
true_pos_cv_100 <- unname(confusion_mat_cv_100$byClass[1]) #Sensitivity/Hit Rate
fpr_cv_100 <- 1-unname(confusion_mat_cv_100$byClass[2]) #1-TNR = 1-Specificity
true_neg_cv_100 <- unname(confusion_mat_cv_100$byClass[2]) #Sensitivity/Hit Rate
fnr_cv_100 <- 1 - unname(confusion_mat_cv_100$byClass[1]) #Miss Rate
roc_cv_100 <- roc(dir_cv_piura_tumbes_posterior_pred_dt$OUTBREAK_100_TRUE, 
                  dir_cv_piura_tumbes_posterior_pred_dt$OUTBREAK_100_PRED)
auc_cv_100 <- auc(roc_cv_100)
auc_cv_100
sens_spec_cv_100_table <-data.table("Accuracy" = accuracy_cv_100,
                                    "True Positive" = true_pos_cv_100,
                                    "False Positive" = fpr_cv_100,
                                    "True Negative Rate" = true_neg_cv_100,
                                    "False Negative Rate" = fnr_cv_100,
                                    "AUC" = auc_cv_100)
sens_spec_cv_100_table
dir_cv_piura_tumbes_posterior_pred_dt[, OUTBREAK_150_TRUE:= ifelse(DIR >= 150, 1, 0)]
dir_cv_piura_tumbes_posterior_pred_dt[, OUTBREAK_150_PRED:= ifelse(OUTBREAK_PROB_150 >= 0.5, 1, 0)]
confusion_mat_cv_150 <- confusionMatrix(data=as.factor(dir_cv_piura_tumbes_posterior_pred_dt$OUTBREAK_150_PRED), 
                                        reference = as.factor(dir_cv_piura_tumbes_posterior_pred_dt$OUTBREAK_150_TRUE),
                                        positive = "1")
accuracy_cv_150 <- unname(confusion_mat_cv_150$overall[1])
precision_cv_150 <- unname(confusion_mat_cv_150$byClass[3])

true_pos_cv_150 <- unname(confusion_mat_cv_150$byClass[1]) #Sensitivity/Hit Rate
fpr_cv_150 <- 1-unname(confusion_mat_cv_150$byClass[2]) #1-TNR = 1-Specificity
true_neg_cv_150 <- unname(confusion_mat_cv_150$byClass[2]) #Sensitivity/Hit Rate
fnr_cv_150 <- 1 - unname(confusion_mat_cv_150$byClass[1]) #Miss Rate
roc_cv_150 <- roc(dir_cv_piura_tumbes_posterior_pred_dt$OUTBREAK_150_TRUE, 
                  dir_cv_piura_tumbes_posterior_pred_dt$OUTBREAK_150_PRED)
auc_cv_150 <- auc(roc_cv_150)
auc_cv_150
sens_spec_cv_150_table <- data.table("Accuracy" = accuracy_cv_150,
                                     "Precision" = precision_cv_150,
                                     "True Positive" = true_pos_cv_150,
                                     "False Positive" = fpr_cv_150,
                                     "True Negative Rate" = true_neg_cv_150,
                                     "False Negative Rate" = fnr_cv_150,
                                     "AUC" = auc_cv_150)



dir_cv_piura_tumbes_posterior_pred_dt[, OUTBREAK_50_TRUE:= ifelse(DIR >= 50, 1, 0)]
dir_cv_piura_tumbes_posterior_pred_dt[, OUTBREAK_50_PRED:= ifelse(OUTBREAK_PROB_50 >= 0.5, 1, 0)]
confusion_mat_cv_50 <- confusionMatrix(data=as.factor(dir_cv_piura_tumbes_posterior_pred_dt$OUTBREAK_50_PRED), 
                                       reference = as.factor(dir_cv_piura_tumbes_posterior_pred_dt$OUTBREAK_50_TRUE),
                                       positive = "1")
confusion_mat_cv_50$byClass
accuracy_cv_50 <- unname(confusion_mat_cv_50$overall[1])
precision_cv_50 <- unname(confusion_mat_cv_50$byClass[3])

true_pos_cv_50 <- unname(confusion_mat_cv_50$byClass[1]) #Sensitivity/Hit Rate
fpr_cv_50 <- 1-unname(confusion_mat_cv_50$byClass[2]) #1-TNR = 1-Specificity
true_neg_cv_50 <- unname(confusion_mat_cv_50$byClass[2]) #Sensitivity/Hit Rate
fnr_cv_50 <- 1 - unname(confusion_mat_cv_50$byClass[1]) #Miss Rate
roc_cv_50 <- roc(dir_cv_piura_tumbes_posterior_pred_dt$OUTBREAK_50_TRUE, 
                 dir_cv_piura_tumbes_posterior_pred_dt$OUTBREAK_50_PRED)
auc_cv_50 <- auc(roc_cv_50)
auc_cv_50
sens_spec_cv_50_table <-data.table("Accuracy" = accuracy_cv_50,
                                   "Precision" = precision_cv_50, 
                                   "True Positive" = true_pos_cv_50,
                                   "False Positive" = fpr_cv_50,
                                   "True Negative Rate" = true_neg_cv_50,
                                   "False Negative Rate" = fnr_cv_50,
                                   "AUC" = auc_cv_50)
sens_spec_cv_50_table
xtable(rbind(sens_spec_cv_50_table, sens_spec_cv_150_table))

#Outbreaks truly greater than 50
outbreak_50_cv_dt <- dir_cv_piura_tumbes_posterior_pred_dt[which(DIR >= 50), ]
outbreak_50_cv_dt[, IND:= seq(1, nrow(outbreak_50_cv_dt))]
median_outbreak_50_cv <- median(outbreak_50_cv_dt$OUTBREAK_PROB_50)
median_outbreak_50_cv
outbreak_50_cv_plot <- ggplot(outbreak_50_cv_dt, aes(x = IND, y = DIR))+
  geom_col(aes(fill = OUTBREAK_PROB_50), width = 0.7)+
  coord_flip()+
  geom_hline(aes(yintercept = 50),
             linetype = "longdash", linewidth = 2.0,
             colour = "forest green")+
  scale_y_continuous(breaks = seq(0, 1000, by = 100))+
  geom_text(aes(label = label_percent()(OUTBREAK_PROB)),
            hjust = 0, fontface = "bold",
            size = 3.5)+
  scale_fill_gradient(
    name = expression("Prob(DIR" >= "50) "),
    low = "#FEED99",
    high = "#AF3301",
    guide = "colourbar",
    breaks = c(0.0, 0.25, 0.5, 0.75),
    labels = percent(0.25*0:3)
  )+theme_bw()+
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
        legend.position = "right", legend.key.size = unit(1.5, 'cm'))
outbreak_50_cv_plot
ggsave(outbreak_50_cv_plot, 
       file = file.path(peru.inla.data.out.dir, "outbreak_50_cv_plot.pdf"), h = 14, w = 20)
ggsave(outbreak_50_cv_plot, 
       file = file.path(peru.inla.data.out.dir, "outbreak_50_cv_plot.png"), h = 14, w = 20)



#Outbreaks truly greater than 150
outbreak_150_cv_dt <- dir_cv_piura_tumbes_posterior_pred_dt[which(DIR >= 150), ]
outbreak_150_cv_dt[, IND:= seq(1, nrow(outbreak_150_cv_dt))]
outbreak_150_cv_dt[(which(OUTBREAK_PROB_150 >= 0.5))]
mean_outbreak_150_cv <- mean(outbreak_150_cv_dt$OUTBREAK_PROB_150)
mean_outbreak_150_cv
outbreak_150_cv_dt
outbreak_150_cv_plot <- ggplot(outbreak_150_cv_dt, aes(x = IND, y = DIR))+
  geom_col(aes(fill = OUTBREAK_PROB_150), width = 0.7)+
  coord_flip()+  geom_hline(aes(yintercept = 150),
                            linetype = "longdash", linewidth = 2.0,
                            colour = "forest green")+
  scale_y_continuous(breaks = seq(0, 1000, by = 100))+
  geom_text(aes(label = label_percent()(OUTBREAK_PROB)),
            hjust = 1.25, fontface = "bold",
            size = 7)+
  scale_fill_gradient(
    name = expression("Prob(DIR" >= "150) "),
    low = "#FEED99",
    high = "#AF3301",
    guide = "colourbar",
    breaks = c(0.0, 0.25, 0.5, 0.75),
    labels = percent(0.25*0:3)
  )+theme_bw()+
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
outbreak_150_cv_plot
ggsave(outbreak_150_cv_plot, 
       file = file.path(peru.inla.data.out.dir, "outbreak_150_cv_plot.pdf"), h = 14, w = 20)
ggsave(outbreak_150_cv_plot, 
       file = file.path(peru.inla.data.out.dir, "outbreak_150_cv_plot.png"), h = 14, w = 20)

#Plot for side-by-side outbreak probabilities
outbreak_50_150_cv_grid <- ggarrange(plotlist = list(outbreak_50_cv_plot, outbreak_150_cv_plot), nrow = 1, ncol = 2,common.legend = FALSE,
                                     legend = "bottom")
outbreak_50_150_cv_grid
ggsave(outbreak_50_150_cv_grid, 
       file = file.path(peru.inla.data.out.dir, "outbreak_50_150_cv_grid.pdf"), h = 16, w = 24)
ggsave(outbreak_50_150_cv_grid, 
       file = file.path(peru.inla.data.out.dir, "outbreak_50_150_cv_grid.png"), h = 16, w = 24)

#False Outbreaks----
false_outbreak_50_cv_dt <- dir_cv_piura_tumbes_posterior_pred_dt[which(DIR < 50), ]
mean_false_outbreak_50_cv_prob <- mean(false_outbreak_50_cv_dt$OUTBREAK_PROB_50)
mean_false_outbreak_50_cv_prob
max_false_outbreak_50_cv_prob <- 
  false_outbreak_50_cv_dt[, max(OUTBREAK_PROB_50)]
max_false_outbreak_50_cv_prob
false_outbreak_50_cv_dt[which.max(OUTBREAK_PROB_50),]

false_outbreak_50_cv_dt[, IND:= seq(1, nrow(false_outbreak_50_cv_dt))]
false_outbreak_50_cv_dt
false_outbreak_50_cv_plot <- ggplot(false_outbreak_50_cv_dt, aes(x = IND, y = DIR))+
  geom_col(aes(fill = OUTBREAK_PROB_50), width = 0.7)+
  coord_flip()+
  scale_x_continuous(breaks = seq(0, 360, by = 25))+
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
    labels = percent(0.25*0:3)
  )+theme_bw()+
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
false_outbreak_50_cv_plot
ggsave(false_outbreak_50_cv_plot, 
       file = file.path(peru.inla.data.out.dir, "false_outbreak_50_cv_plot.pdf"), h = 14, w = 20)
ggsave(false_outbreak_50_cv_plot, 
       file = file.path(peru.inla.data.out.dir, "false_outbreak_50_cv_plot.png"), h = 14, w = 20)


false_outbreak_150_cv_dt
false_outbreak_150_cv_dt <- dir_cv_piura_tumbes_posterior_pred_dt[which(DIR < 150), ]
mean(false_outbreak_150_cv_dt$OUTBREAK_PROB_50)

max_false_outbreak_150_cv_prob <- 
  false_outbreak_150_cv_dt[, max(OUTBREAK_PROB_150)]
max_false_outbreak_150_cv_prob
false_outbreak_150_cv_dt[, IND:= seq(1, nrow(false_outbreak_150_cv_dt))]
false_outbreak_150_cv_dt$IND
false_outbreak_150_cv_plot <- ggplot(false_outbreak_150_cv_dt, aes(x = IND, y = DIR))+
  geom_col(aes(fill = OUTBREAK_PROB_150), width = 0.7)+
  coord_flip()+
  scale_x_continuous(breaks = seq(0, 400, by = 25))+
  scale_y_continuous(breaks = seq(0, 150, by = 15),
                     limits = c(0, 150))+
  scale_fill_gradient(
    name = expression("Prob(DIR" >= "150) "),
    low = "#FEED99",
    high = "#AF3301",
    guide = "colourbar",
    breaks = c(0.0, 0.25, 0.5, 0.75),
    labels = percent(0.25*0:3)
  )+theme_bw()+
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
false_outbreak_150_cv_plot
ggsave(false_outbreak_150_cv_plot, 
       file = file.path(peru.inla.data.out.dir, "false_outbreak_150_cv_plot.pdf"), h = 14, w = 20)
ggsave(false_outbreak_150_cv_plot, 
       file = file.path(peru.inla.data.out.dir, "false_outbreak_150_cv_plot.png"), h = 14, w = 20)

false_outbreak_50_150_cv_grid <- ggarrange(false_outbreak_50_cv_plot, false_outbreak_150_cv_plot, nrow = 1)
false_outbreak_50_150_cv_grid
ggsave(false_outbreak_50_150_cv_grid, 
       file = file.path(peru.inla.data.out.dir, "false_outbreak_50_150_cv_grid.pdf"), h = 14, w = 24)
ggsave(false_outbreak_50_150_cv_grid, 
       file = file.path(peru.inla.data.out.dir, "false_outbreak_50_150_cv_grid.png"), h = 14, w = 24)

