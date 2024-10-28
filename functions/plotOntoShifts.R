inverse_logit = function(TL, muAlpha, muBeta, sigmaAlpha, sigmaBeta) {
  
  # performs inverse-logit transformation using parameter values
  # from logistic model
  
  # parameter values -----------------------------------------------------------
  lowerAlpha = muAlpha - sigmaAlpha
  upperAlpha = muAlpha + sigmaAlpha
  lowerBeta  = muBeta - sigmaBeta
  upperBeta  = muBeta + sigmaBeta
  
  # inverse-logit transform ----------------------------------------------------
  mean  = 1/(1 + exp(-(muAlpha + muBeta*TL)))
  upper = 1/(1 + exp(-(upperAlpha + upperBeta*TL)))
  lower = 1/(1 + exp(-(lowerAlpha + lowerBeta*TL)))
  
  # create output
  output = data.frame(TL, mean, upper, lower)
  
  return(output)
  
}

plotOntoShifts = function(outputs, prey, minTL, maxTL, plotNcols) {
  
  # Plots relationships between total lenght (TL) and the probabilities 
  # of observing certain prey types
  
  # load required packages, install if not detected ----------------------------
  library(ggmcmc)
  library(tidyr)
  library(ggplot2)
  
  # return error message if lists are NOT of same length -----------------------
  if (length(outputs) != length(prey)) {
    stop("Please provide lists with same dimensions")
  }
  
  nModels = length(outputs) # number of models
  
  # create summarized dataframes with parameter values -------------------------
  modelDFs = list()
  for(i in 1:nModels) {
    
    modelDFs[[i]] = ggs(outputs[[i]]) %>%
      group_by(Parameter) %>%
      filter(Parameter %in% c("alpha", "beta")) %>%
      summarise(
        mean = median(value),
        sd   = sd(value)
        ) %>%
      mutate(Prey = paste0(prey[[i]]))
      
  }
  
  # create dataframes to be used for plotting ----------------------------------
  
  TL = seq(minTL, maxTL, .1)
  
  plotInputs = list()
  for(i in 1:nModels) {
    
    plotInputs[[i]] = inverse_logit(
                                    TL         = TL,
                                    muAlpha    = as.numeric(modelDFs[[i]][1,2]),
                                    muBeta     = as.numeric(modelDFs[[i]][2,2]),
                                    sigmaAlpha = as.numeric(modelDFs[[i]][1,3]/2),
                                    sigmaBeta  = as.numeric(modelDFs[[i]][2,3]/2)
                                    ) %>%
      mutate(Prey = paste0(prey[[i]]))
    
  }

  TL_dfs = do.call("rbind", plotInputs) # joining dataframes
  
  # ploting --------------------------------------------------------------------

   TL_dfs %>%
    ggplot() +
    geom_line(aes(x = TL, y = mean)) +
    geom_ribbon(aes(x = TL, ymin = lower, ymax = upper),
                alpha = 0.2, fill = "green") +
    facet_wrap(~Prey, scales = "free_y", ncol = plotNcols) + 
    theme_custom() +
    xlab("Total length (mm)") +
    ylab("PoE (%)")

}