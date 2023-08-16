g = inla.read.graph(file.path(peru.inla.data.in.dir, "nbr_piura_tumbes_lambayeque.graph"))

#Specifying PC Priors (Simpson, 2017) ----
prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(0.5, 0.01) )) 


#Setting up lagged variables ----
maximum_climate_lag <- 4

# Maximum temperature (Tmax)
lag_tmax <- tsModel::Lag(ptl_climate_dt$tmax, 
                                 group = ptl_climate_dt$REGION, 
                                 k = 0:maximum_climate_lag)
#Only keep from Time Point 5 (May 2010) onwards -
  # Reason = We use a 3-month period to evaluate momentum (RSI) and the
        # RSI predictor is lagged by 1 month -> May RSI predictor is the value from April
lag_tmax <- lag_tmax[13:nrow(lag_tmax),]

# Minimum temperature (Tmin)
lag_tmin <- tsModel::Lag(ptl_climate_dt$tmin, 
                                 group = ptl_climate_dt$REGION, 
                                 k = 0:maximum_climate_lag)
lag_tmin <- lag_tmin[13:nrow(lag_tmin),]

#TMAX-PREC
lag_tmax_prec <- tsModel::Lag(ptl_climate_dt$tmax_prec, 
                                      group = ptl_climate_dt$REGION, 
                                      k = 0:maximum_climate_lag)
lag_tmax_prec <- lag_tmax_prec[13:nrow(lag_tmax_prec),]


#SPI-PREC Interaction
lag_spi_prec <- tsModel::Lag(ptl_climate_dt$spi_prec, 
                                     group = ptl_climate_dt$REGION, 
                                     k = 0:maximum_climate_lag)
lag_spi_prec <- lag_spi_prec[13:nrow(lag_spi_prec),]

# Precipitation (prec)
lag_prec <- tsModel::Lag(ptl_climate_dt$prec, 
                                 group = ptl_climate_dt$REGION, 
                                 k = 0:maximum_climate_lag)
lag_prec <- lag_prec[13:nrow(lag_prec),]

# SPI-6
lag_spi <- tsModel::Lag(ptl_climate_dt$SPI_6, 
                                group = ptl_climate_dt$REGION, 
                                k = 0:maximum_climate_lag)
lag_spi <- lag_spi[13:nrow(lag_spi),]


# ONI
lag_oni <- tsModel::Lag(ptl_climate_dt$ANOM, 
                                k = 0:maximum_climate_lag)
lag_oni <- lag_oni[13:nrow(lag_oni),]


#ICEN
lag_icen <- tsModel::Lag(ptl_climate_dt$E_INDEX, 
                                 k = 0:maximum_climate_lag)
lag_icen <- lag_icen[13:nrow(lag_icen),]





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

lagknot <- 2
setkeyv(ptl_dt, c("TIME", "REGION"))

#Cross basis matrices----
#Via the cross basis (a bi-dimensional functional space) we are specifying simultaneously the 
# relationships in the dimensions of the predictor and lags, respectively.
tmax_basis <- crossbasis(lag_natural_tmax, 
                                 argvar=list(fun = "bs", knots = equalknots(ptl_climate_dt$tmax, 2)),
                                 arglag=list(fun="bs", knots =lagknot))
tmin_basis <- crossbasis(lag_natural_tmin, 
                                 argvar=list(fun = "bs", knots = equalknots(ptl_climate_dt$tmin, 2)),
                                 arglag=list(fun="bs", knots =lagknot))

tmax_prec_basis <- crossbasis(lag_natural_tmax_prec, 
                                      argvar=list(fun = "bs", knots = equalknots(ptl_climate_dt$tmax_prec, 2)),
                                      arglag=list(fun="bs", knots =lagknot))

spi_prec_basis <- crossbasis(lag_natural_spi_prec, 
                                     argvar=list(fun = "bs", knots = equalknots(ptl_climate_dt$spi_prec, 2)),
                                     arglag=list(fun="bs", knots =lagknot))

prec_basis <- crossbasis(lag_natural_prec, 
                                 argvar=list(fun = "bs", knots = equalknots(ptl_climate_dt$prec, 2)),
                                 arglag=list(fun="bs", knots = lagknot))

spi_basis <- crossbasis(lag_natural_spi, 
                                argvar=list(fun = "bs", knots = equalknots(ptl_climate_dt$SPI_6, 2)),
                                arglag=list(fun="bs", knots = lagknot))

oni_basis <- crossbasis(lag_natural_oni,  
                                argvar=list(fun = "bs", knots = equalknots(ptl_climate_dt$ANOM, 2)),
                                arglag=list(fun="bs", knots = lagknot))
icen_basis <- crossbasis(lag_natural_icen,  
                                 argvar=list(fun = "bs", knots = equalknots(ptl_climate_dt$E_INDEX, 2)),
                                 arglag=list(fun="bs", knots = lagknot))

#Column Names of cross basis matrices
colnames(tmax_basis) = paste0("tmax_basis.", colnames(tmax_basis))
colnames(tmax_prec_basis) = paste0("tmax_prec_basis.", colnames(tmax_prec_basis))
colnames(spi_prec_basis) = paste0("spi_prec_basis.", colnames(spi_prec_basis))
colnames(oni_basis) = paste0("oni_basis.", colnames(oni_basis))
colnames(tmin_basis) = paste0("tmin_basis.", colnames(tmin_basis))
colnames(prec_basis) = paste0("prec_basis.", colnames(prec_basis))
colnames(spi_basis) = paste0("spi_basis.", colnames(spi_basis))
colnames(icen_basis) = paste0("icen_basis.", colnames(icen_basis))


ptl_climate_formula_tmax <- update.formula(baseline_formula, ~. + tmax_basis)
ptl_climate_formula_tmin <- update.formula(baseline_formula, ~. + tmin_basis)
ptl_climate_formula_prec <- update.formula(baseline_formula, ~. + prec_basis)
ptl_climate_formula_oni <- update.formula(baseline_formula, ~. + oni_basis)
ptl_climate_formula_spi <- update.formula(baseline_formula, ~. + spi_basis)
ptl_climate_formula_icen <- update.formula(baseline_formula, ~. + icen_basis)
ptl_climate_formula_tmax_prec <- update.formula(baseline_formula, ~. + tmax_basis + prec_basis)
ptl_climate_formula_icen_prec <- update.formula(baseline_formula, ~. + icen_basis+ prec_basis)
ptl_climate_formula_tmax_spi <- update.formula(baseline_formula, ~. + tmax_basis + spi_basis)
ptl_climate_formula_icen_spi <- update.formula(baseline_formula, ~. + icen_basis+ spi_basis)
ptl_climate_formula_tmax_icen <- update.formula(baseline_formula, ~. + tmax_basis + icen_basis)
ptl_climate_formula_icen_icen <- update.formula(baseline_formula, ~. + icen_basis + icen_basis)
ptl_climate_formula_tmax_prec_icen <- update.formula(baseline_formula, ~. + tmax_basis + prec_basis+icen_basis)
ptl_climate_formula_tmax_prec_icen_spi <- update.formula(baseline_formula, ~. + tmax_basis + prec_basis+icen_basis+
                                                               spi_basis)
ptl_climate_formula_tmax_icen_spi <- update.formula(baseline_formula, ~. + tmax_basis +icen_basis+
                                                          spi_basis)
ptl_climate_formula_tmax_prec_icen_spi_spi_prec <- update.formula(baseline_formula, ~. + tmax_basis + prec_basis+icen_basis+spi_basis+
                                                                spi_prec)


ptl_climate_formulae <- c(ptl_climate_formula_tmax, ptl_climate_formula_icen, ptl_climate_formula_prec,
                              ptl_climate_formula_oni, ptl_climate_formula_spi, ptl_climate_formula_icen, 
                              ptl_climate_formula_tmax_prec, ptl_climate_formula_icen_prec,
                              ptl_climate_formula_tmax_spi, ptl_climate_formula_icen_spi,
                              ptl_climate_formula_tmax_icen, ptl_climate_formula_icen_icen,
                              ptl_climate_formula_tmax_prec_icen, ptl_climate_formula_tmax_prec_icen_spi,
                              ptl_climate_formula_tmax_icen_spi, ptl_climate_formula_tmax_prec_icen_spi_spi_prec)

run_inla_mods <- function(formulae, family, data)
{
  inla_base_mod_summary <- data.table(MOD = paste0("Mod_", seq(1, length(formulae))),
                                      DIC = rep(0, length(formulae)),
                                      WAIC = rep(0, length(formulae)),
                                      p_eff = rep(0, length(formulae)),
                                      LOG_SCORE = rep(0, length(formulae)))
  for(i in 1:length(formulae)){
    print(paste0("We're now on Model ", i))
    formula_in_q <- formulae[[i]]
    if(family == "nbinomial")
    {inla_mod_fit <- inla(formula = formula_in_q, 
                          data = data, family = family, offset = log(POP_OFFSET_INTERP),
                          verbose = FALSE,
                          control.inla = list(strategy = 'adaptive'), 
                          control.compute = list(waic = TRUE, dic = TRUE, 
                                                 cpo = TRUE, config = TRUE,
                                                 return.marginals = TRUE),
                          control.fixed = list(correlation.matrix = TRUE, 
                                               prec.intercept = 1, prec = 1),
                          control.predictor = list(link = 1, compute = TRUE), 
    )
    }
    inla_mod_fit <- inla.rerun(inla_mod_fit)
    #Save Model Fit
    saveRDS(inla_mod_fit, file = file.path(peru.inla.data.out.dir, paste0("Mod_", i,".RDS")))
    #Update model summary table
    inla_base_mod_summary[which(MOD == paste0("Mod_", i)),
                          DIC:= inla_mod_fit$dic$dic]
    inla_base_mod_summary[which(MOD == paste0("Mod_", i)),
                          WAIC:= inla_mod_fit$waic$waic]
    inla_base_mod_summary[which(MOD == paste0("Mod_", i)),
                          p_eff:= inla_mod_fit$waic$p.eff]
    inla_base_mod_summary[which(MOD == paste0("Mod_", i)),
                          LOG_SCORE:= - mean(log(inla_mod_fit$cpo$cpo))]
    inla_base_mod_summary[which(MOD == paste0("Mod_", i)),
                          MAE:= 
                            MAE(inla_mod_fit$summary.fitted.values$`0.5quant`,
                                data$DIR_POP_INTERP)]
  }
  return(inla_base_mod_summary)
}

climate_mods_summary <- run_inla_mods(climate_formulae,
                                              family = "nbinomial",
                                              data = ptl_dt)

climate_mods_summary

process_posterior_preds <- function(inla_mod_fit){
  in_sample_post_preds <- data.table(MEDIAN = (inla_mod_fit$summary.fitted.values$`0.5quant`),
                                     q2.5 = (inla_mod_fit$summary.fitted.values$`0.025quant`),
                                     q97.5 = (inla_mod_fit$summary.fitted.values$`0.975quant`),
                                     # CASES = (inla_df$CASES),
                                     DIR_POP_INTERP = inla_df$DIR_POP_INTERP,
                                     REGION = inla_df$REGION,
                                     MONTH = inla_df$MONTH,
                                     YEAR = inla_df$YEAR,
                                     POP = inla_df$POP,
                                     TIME = inla_df$TIME,
                                     IND = seq(1, nrow(inla_df))
  )
  # # print(in_sample_post_preds)
  # posterior_pred_plot <- ggplot(in_sample_post_preds)+
  #   geom_line(aes(x = TIME, y = DIR_POP_INTERP, col = "Observed"),
  #             linewidth = 1.2)+ 
  #   # geom_point(aes(x = TIME, y = DIR_POP_INTERP, col = "Observed"))+ 
  #   geom_ribbon(aes(x = TIME, ymin = q2.5, ymax = q97.5), alpha = 0.2)+
  #   geom_line(aes(x = TIME, y = MEDIAN, col = "Estimated"),
  #             linewidth = 1.2)+
  #   # geom_point(aes(x = TIME, y = MEDIAN, col = "Estimated"))+
  #   scale_x_continuous(breaks = seq(0, 140, by = 20))+
  #   facet_wrap( .~ REGION, nrow = 3,
  #               scales = "free_y")+theme_bw()+
  #   xlab("Month")+ylab("DIR")+
  #   theme(text = element_text(size = 25),
  #         axis.text.x = element_text(size=25),
  #         axis.text.y = element_text(size=25),
  #         panel.grid.minor.y = element_blank(),
  #         panel.grid.minor.x = element_blank(),
  #         panel.grid.major.y = element_blank(),
  #         panel.grid.major.x = element_blank(),
  #         plot.title = element_blank(),
  #         plot.subtitle = element_blank(),
  #         axis.title=element_text(size=25), 
  #         legend.text=element_text(size=25)+geom_text(size = 25),
  #         legend.position = "bottom",
  #         legend.title = element_blank())
  return(list(in_sample_post_preds))
}


ptl_climate_mod_fit <- inla(formula = ptl_climate_formulae[[14]], 
                                data = ptl_dt, family = "nbinomial", offset = log(POP_OFFSET_INTERP),
                                verbose = FALSE,
                                control.inla = list(strategy = 'adaptive'), 
                                control.compute = list(waic = TRUE, dic = TRUE, 
                                                       cpo = TRUE, config = TRUE,
                                                       return.marginals = TRUE),
                                control.fixed = list(correlation.matrix = TRUE, 
                                                     prec.intercept = 1, prec = 1),
                                control.predictor = list(link = 1, compute = TRUE)
)

climate_posterior_pred_dt <- process_posterior_preds(ptl_climate_mod_fit)[[1]]
climate_mod_fit$dic$dic
climate_mod_fit$waic$waic





#Supplementary Material Figures For Model Fit 
#Marginal Effects ----
marginal <- inla.smarginal(ptl_climate_mod_fit$marginals.fixed$SEASON)
marginal <- data.frame(marginal)
ggplot(marginal, aes(x = x, y = y)) + geom_line() +
  labs(x = expression(beta[1]), y = "Density") +
  geom_vline(xintercept = 0, col = "black") + theme_bw()


marginal <- inla.smarginal(ptl_climate_mod_fit$marginals.fixed$RSI_DIR_POP_INTERP_LAG)
marginal <- data.frame(marginal)
ggplot(marginal, aes(x = x, y = y)) + geom_line() +
  labs(x = expression(beta[1]), y = "Density") +
  geom_vline(xintercept = 0, col = "black") + theme_bw()

#Plotting Random Effects ----
#Yearly IID Rand Effects are per region
yearly_random_effects <-data.table(ptl_climate_mod_fit$summary.random$YEAR)
yearly_random_effects[, YEAR:= ID+2009]
yearly_random_effects[, REGION:= rep(unique(ptl_dt$REGION), 
                                     each = length(unique(YEAR)))]
yearly_random_effects <- subset(yearly_random_effects, select = c("mean", "0.025quant", "0.5quant", "0.975quant", "YEAR", "REGION"))
setnames(yearly_random_effects, colnames(yearly_random_effects), toupper(c("mean", "CI_L", "MEDIAN", "CI_U", "YEAR", "REGION")))
yearly_rand_effect_plot <- ggplot(yearly_random_effects) + geom_point(aes(x = YEAR, y = MEDIAN, col = REGION),
                                                                      size = 1.4)+
  geom_errorbar(aes(x = YEAR, ymin = CI_L, ymax = CI_U, col = REGION), alpha = 0.6,
                linewidth = 1.2)+
  theme_bw()+
  geom_hline(aes(yintercept = 0),
             linewidth = 0.5, linetype = "longdash")+
  facet_wrap(REGION ~. )+
  labs(x = "Year", y = "Effect Size (Log)",
       title = "Out-of-Sample Prevalence: REACT vs Estimated",
       subtitle = "Rounds 15-19")+
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
        legend.position = "none",
        legend.title = element_blank())+
  scale_x_continuous(breaks = seq(2010, 2021, by = 2))
yearly_rand_effect_plot
ggsave(yearly_rand_effect_plot, 
       file = file.path(peru.inla.data.out.dir, "yearly_rand_effect_plot.pdf"), h = 12, w = 22)
ggsave(yearly_rand_effect_plot, 
       file = file.path(peru.inla.data.out.dir, "yearly_rand_effect_plot.png"), h = 12, w = 22)


#Cyclic RW(1) Effects
monthly_random_effects <- data.table(ptl_climate_mod_fit$summary.random$MONTH)
monthly_random_effects[, MONTH:= ID]
monthly_random_effects[, REGION:= rep(unique(ptl_dt$REGION), 
                                      each = length(unique(MONTH)))]
monthly_random_effects <- subset(monthly_random_effects, select = c("mean", "0.025quant", "0.5quant", "0.975quant", "MONTH", "REGION"))
setnames(monthly_random_effects, colnames(monthly_random_effects), toupper(c("mean", "CI_L", "MEDIAN", "CI_U", "MONTH", "REGION")))
monthly_rand_effect_plot <- ggplot(monthly_random_effects) + geom_point(aes(x = MONTH, y = MEDIAN, col = REGION),
                                                                        size = 1.4)+
  geom_errorbar(aes(x = MONTH, ymin = CI_L, ymax = CI_U, col = REGION), alpha = 0.6,
                linewidth = 1.2)+
  theme_bw()+
  facet_wrap(REGION ~. )+
  geom_hline(aes(yintercept = 0),
             linewidth = 0.5, linetype = "longdash")+
  labs(x = "Month", y = "Effect Size (Log)",
       title = "Out-of-Sample Prevalence: REACT vs Estimated",
       subtitle = "Rounds 15-19")+
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
        legend.position = "none",
        legend.title = element_blank())+
  scale_x_continuous(breaks = seq(1, 12, by = 2))
monthly_rand_effect_plot
ggsave(monthly_rand_effect_plot, 
       file = file.path(peru.inla.data.out.dir, "monthly_rand_effect_plot.pdf"), h = 12, w = 22)
ggsave(monthly_rand_effect_plot, 
       file = file.path(peru.inla.data.out.dir, "monthly_rand_effect_plot.png"), h = 12, w = 22)


#Have 2x Number of areas
#First half correspond to combined effects 
#BYM2 rand effects include structured and unstructured effects

regional_bym2_random_effects <- data.table(ptl_climate_mod_fit$summary.random$RGN_IND)
regional_bym2_random_effects <- subset(regional_bym2_random_effects, select = c("mean", "0.025quant", "0.5quant", "0.975quant"))
setnames(regional_bym2_random_effects, colnames(regional_bym2_random_effects), toupper(c("mean", "CI_L", "MEDIAN", "CI_U")))
regional_bym2_random_effects <- regional_bym2_random_effects[1:3, ]
regional_bym2_random_effects[, REGION:= unique(ptl_dt$REGION)]
regional_bym2_random_effects[, IND:= seq(1,3)]

bym2_rand_effect_plot <- ggplot(regional_bym2_random_effects) + geom_point(aes(x = IND, y = MEDIAN, col = REGION),
                                                                           size = 1.4)+
  geom_errorbar(aes(x = IND, ymin = CI_L, ymax = CI_U, col = REGION), alpha = 0.6,
                linewidth = 1.2)+
  theme_bw()+
  geom_hline(aes(yintercept = 0),
             linewidth = 0.5, linetype = "longdash")+
  coord_flip()+
  labs(x = "", y = "Effect Size (Log)",
       title = "Out-of-Sample Prevalence: REACT vs Estimated",
       subtitle = "Rounds 15-19")+
  theme(text = element_text(size = 25),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=25),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.title=element_text(size=25), 
        legend.text=element_text(size=25),
        legend.position = "bottom",
        legend.title = element_blank())
bym2_rand_effect_plot
ggsave(bym2_rand_effect_plot, 
       file = file.path(peru.inla.data.out.dir, "bym2_rand_effect_plot.pdf"), h = 14, w = 20)
ggsave(bym2_rand_effect_plot, 
       file = file.path(peru.inla.data.out.dir, "bym2_rand_effect_plot.png"), h = 14, w = 20)



#Hyperparamaters
ggplot() + geom_line(data = as.data.frame(climate_mod_fit$marginals.hyperpar$`Phi for RGN_IND`), 
                     aes(x = x, y = y)) + theme_bw() + 
  ggtitle("Posterior of sd of the mixing parameter") 
                                





#DLNM Plotting ----
#Function for RR on Log Scale

dlnm_func <- function(){
  
  #Citation: 
  ## Rachel Lowe (2021) - https://github.com/drrachellowe/hydromet_dengue
  
  coef <- climate_mod_fit$summary.fixed$mean
  vcov <- climate_mod_fit$misc$lincomb.derived.covariance.matrix
  indt <- grep("tmax_basis", climate_mod_fit$names.fixed)
  predt <- crosspred(tmax_basis, coef = coef[indt], vcov=vcov[indt,indt],
                     model.link = "log", bylag = 0.5, cen = mean(ptl_dt$tmax)) 
  y <- predt$predvar
  x <- seq(0, 4, 0.1)
  z <- t(log(predt$matRRfit))
  pal <- rev(brewer.pal(11, "PRGn"))
  levels <- seq(-1.1, 1.1, by = 0.05)
  col1 <- colorRampPalette(pal[1:6])
  
  col2 <- colorRampPalette(pal[6:11])
  cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))
  
  pdf(file.path(peru.inla.data.out.dir, "log_rr_new_plots_four_pages.pdf"), h = 6, w = 12)
  
  filled.contour(x,y,z,
                 plot.title = title(main = "Maximum Temperature", cex.lab = 1.25,
                                    xlab = "Lag (months)", 
                                    ylab = expression(paste("Maximum Temperature (",degree,"C)"))),
                 col = cols,levels = levels,
                 plot.axes = { axis(1, at = 0:4, c(0:4),
                                    cex.axis = 1.25) 
                   axis(2, cex.axis = 1.25)},
                 key.title = title(main = "log(RR)", cex = 1.1),
                 cex.lab = 1.5,
                 cex.axis = 1.5)
  
  indt <- grep("prec_basis", climate_mod_fit$names.fixed)
  # extract predictions from the tmax DLNM centred on overall mean tmax (19 deg C)
  predt <- crosspred(prec_basis, coef = coef[indt], vcov=vcov[indt,indt],
                     model.link = "log", bylag = 0.2, cen = (mean(ptl_dt$prec)),
                     by = 10) 

  y <- predt$predvar
  x <- seq(0, 4, 0.2)
  z <- t(log(predt$matRRfit))
  col_val <- max(abs(z))
  pal <- rev(brewer.pal(11, "PRGn"))
  levels <- seq(-1.1, 1.1, by = 0.05)
  col1 <- colorRampPalette(pal[1:6])
  
  col2 <- colorRampPalette(pal[6:11])
  cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))
  
  filled.contour(x,y,z,
                 col = cols,levels = levels,
                 plot.axes = { axis(1, at = 0:6, c(0:6), cex.axis = 1.25) 
                   axis(2, cex.axis = 1.25)},
                 key.title = title(main = "log(RR)", cex = 1.1),
                 cex.lab = 1.1, plot.title = title(main = "Monthly Precipitation", cex.lab = 1.25,
                                                   xlab = "Lag (months)", ylab = expression(paste("Precipitation (mm)"))))

  indt <- grep("spi_basis", climate_mod_fit$names.fixed)
  # extract predictions from the oni DLNM centred on overall mean SPI (19 deg C)
  predt <- crosspred(spi_basis, coef = coef[indt], vcov=vcov[indt,indt],
                     model.link = "log", bylag = 0.1, cen = (mean(ptl_dt$SPI_6))) 
  y <- predt$predvar
  x <- seq(0, 4, 0.1)
  z <- t(log(predt$matRRfit))
  pal <- rev(brewer.pal(11, "PRGn"))
  levels <- seq(-1.1, 1.1, by = 0.05)
  col1 <- colorRampPalette(pal[1:6])
  
  col2 <- colorRampPalette(pal[6:11])
  cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))
  
  # pdf(file.path(peru.inla.data.out.dir, "spi.pdf"), h = 10, w = 14)
  filled.contour(x,y,z,
                 plot.title = title(main = "Standardized Precipitation Index", cex.lab = 1.25,
                                    xlab = "Lag (months)", 
                                    ylab = "SPI-6"),
                 col = cols,levels = levels,
                 plot.axes = { axis(1, at = 0:4, c(0:4),
                                    cex.axis = 1.25) 
                   axis(2, cex.axis = 1.25)},
                 key.title = title(main = "log(RR)", cex = 1.1),
                 cex.lab = 1.5,
                 cex.axis = 1.5)
  
  indt <- grep("icen_basis", climate_mod_fit$names.fixed)
  predt <- crosspred(icen_basis, coef = coef[indt], vcov=vcov[indt,indt],
                     model.link = "log", bylag = 0.1, cen = round(mean(ptl_dt$E_INDEX), 0) )
  y <- predt$predvar
  x <- seq(0, 4, 0.1)
  z <- t(log(predt$matRRfit))
  pal <- rev(brewer.pal(11, "PRGn"))
  levels <- seq(-1.1, 1.1, by = 0.05)
  col1 <- colorRampPalette(pal[1:6])
  col2 <- colorRampPalette(pal[6:11])
  cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))
  

  filled.contour(x,y,z,
                 plot.title = title(main = "ICEN E-Index", cex.lab = 1.25,
                                    xlab = "Lag (months)", 
                                    ylab = "ICEN E-Index"),
                 col = cols,levels = levels,
                 plot.axes = { axis(1, at = 0:4, c(0:4),
                                    cex.axis = 1.25) 
                   axis(2, cex.axis = 1.25)},
                 key.title = title(main = "log(RR)", cex = 1.1),
                 cex.lab = 1.5,
                 cex.axis = 1.5)
  
  dev.off()
  
}
dlnm_func()




#Cumulative Risk Plotting Function
dlnm_slice_rr_func <- function(){
  coef <- climate_mod_fit$summary.fixed$mean
  vcov <- climate_mod_fit$misc$lincomb.derived.covariance.matrix
  
  #tmax
  indt <- grep("tmax_basis", climate_mod_fit$names.fixed)
  pred2 <- crosspred(tmax_basis, coef = coef[indt], vcov=vcov[indt,indt],
                     model.link = "log", bylag = 0.25, cen = round(mean(ptl_dt$tmax), 0), by = 1,
                     cumul = TRUE)
  
  pdf(file.path(peru.inla.data.out.dir, "tmax_dlnm_cumulative.pdf"), h = 6, w = 12)
  plot(pred2, "slices", var= "28",col=1, pch=19, ylab="Cumulative Risk", xlab="Lag (months)",
       ylim = c(0, 2.1), main = "Maximum Temperature")
  lines(pred2, "slices", var= "25" ,col= "purple", pch=19, ylab="Cumulative Risk", xlab="Lag (months)")
  lines(pred2, "slices", var= "31" ,col= "orange", pch=19, ylab="Cumulative Risk", xlab="Lag (months)")
  
  legend("topright", legend=c(expression(paste("28",degree, "C")), 
                              expression(paste("25",degree, "C")),
                              expression(paste("31",degree, "C"))),
         col=c("black","purple", "orange"), lty=1,
         title = "Maximum Temperature")
  dev.off()
  
  #prec
  indt <- grep("prec_basis", climate_mod_fit$names.fixed)
  pred2 <- crosspred(prec_basis, coef = coef[indt], vcov=vcov[indt,indt],
                     model.link = "log", bylag = 0.25, cen = round(mean(ptl_dt$prec), 0), by = 1,
                     cumul = TRUE)
  
  pdf(file.path(peru.inla.data.out.dir, "prec_dlnm_cumulative.pdf"), h = 6, w = 12)
  

  plot(pred2, "slices", var= "22",col= "black", pch=19, ylab="Cumulative Risk", xlab="Lag (months)",
       ylim = c(0, 2), main = "Monthly Precipitation")
  lines(pred2, "slices", var= "1" ,col= "purple", pch=19, ylab="Cumulative Risk", xlab="Lag (months)")
  lines(pred2, "slices", var= "125", col="orange", pch=19, ylab="Cumulative Risk", xlab="Lag (months)")
  
  
  legend("topright", legend=c(paste("22mm"), 
                              paste("1mm"),
                              paste("125mm")),
         col=c("black","purple", "orange"), lty=1,
         title = "Precipitation")
  
  dev.off()
  
  indt <- grep("spi_basis", climate_mod_fit$names.fixed)
  pred2 <- crosspred(spi_basis, coef = coef[indt], vcov=vcov[indt,indt],
                     model.link = "log", bylag = 0.25, cen = round(mean(ptl_dt$SPI_6), 0), by = 0.1,
                     cumul = TRUE)
  #Cumulative Risk for different lags
  pdf(file.path(peru.inla.data.out.dir, "spi_dlnm_cumulative.pdf"), h = 6, w = 12)
  
  plot(pred2, "slices", var= "0",col= "black", pch=19, ylab="Cumulative Risk", xlab="Lag (months)", xlim = c(0, 4), main = "SPI-6")
  lines(pred2, "slices", var= "-1" ,col= "purple", pch=19, ylab="Cumulative Risk", xlab="Lag (months)")
  lines(pred2, "slices", var= "1", col="orange", pch=19, ylab="Cumulative Risk", xlab="Lag (months)")
  
  
  legend("topright", legend=c(paste("0"), 
                              paste("-1"),
                              paste("1")),
         col=c("black","purple", "orange"), lty=1,
         title = "SPI-6")
  dev.off()
  
  #ICEN
  indt <- grep("icen_basis", climate_mod_fit$names.fixed)
  pred2 <- crosspred(icen_basis, coef = coef[indt], vcov=vcov[indt,indt],
                     model.link = "log", bylag = 0.25, cen = round(mean(ptl_dt$E_INDEX), 0), by = 0.1,
                     cumul = TRUE)
  #Cumulative Risk for different lags
  pdf(file.path(peru.inla.data.out.dir, "icen_dlnm_cumulative.pdf"), h = 6, w = 12)
  plot(pred2, "slices", var= "0",col= "black", pch=19, ylab="Cumulative Risk", xlab="Lag (months)",
       ylim = c(0, 2.5), main = "ICEN")
  lines(pred2, "slices", var= "-1.4" ,col= "purple", pch=19, ylab="Cumulative Risk", xlab="Lag (months)")
  # lines(pred2, "slices", var= "1" , col= "green", pch=19, ylab="RR", xlab="Lag (months)")
  lines(pred2, "slices", var= "1.7", col="orange", pch=19, ylab="Cumulative Risk", xlab="Lag (months)")
  
  legend("topright", legend=c(paste("0"), 
                              paste("-1.4"),
                              paste("1.7")),
         col=c("black","purple", "orange"), lty=1,
         title = "E-Index")
  dev.off()
}
