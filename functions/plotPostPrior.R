plotPostPrior = function (res, priorList, postList, parNames, plotTitle, phiLimit) {
  
  # extract priors from model
  priors         = rstan::extract(res, pars = priorList)
  priorDF        = data.frame(priors[[1]], priors[[2]], priors[[3]])
  names(priorDF) = parNames
  priorDF$phi    = ifelse(priorDF$phi > phiLimit, NA, priorDF$phi)
  priorDF$phi    = ifelse(priorDF$phi < 0, NA, priorDF$phi)
  
  # extract posteriors from model
  posts         = rstan::extract(res, pars = postList)
  postDF        = data.frame(posts[[1]], posts[[2]], posts[[3]])  
  names(postDF) = parNames
  postDF$phi    = ifelse(postDF$phi > phiLimit, NA, postDF$phi)
  postDF$phi    = ifelse(postDF$phi < 0, NA, postDF$phi)
  
  # convert to long format
  require(dplyr);require(tidyr)
  priorDF = priorDF %>% pivot_longer(1:ncol(priorDF), names_to = "Parameter", values_to = "Values") %>%
    mutate(Source = "Prior")
  postDF  = postDF  %>% pivot_longer(1:ncol(postDF),  names_to = "Parameter", values_to = "Values") %>%
    mutate(Source = "Posterior")
  
  mergedDF = rbind(priorDF, postDF)
  
  # plotting
  mergedDF$Parameter = factor(mergedDF$Parameter,
                              labels = c(bquote(alpha),
                                         bquote("beta[TL]"),
                                         bquote(phi)))
  
  require(ggplot2)
  mergedDF %>% 
    ggplot() +
    stat_density(aes(x = Values,
                     y = ..scaled..,
                     color = Source,
                     linetype = Source),
                 position = "dodge",
                 geom = "density",
                 alpha = 0.3,
                 linewidth = 0.7) +
    facet_wrap(~Parameter, scales = "free", labeller = label_parsed) + 
    theme_custom() +
    xlab("Parameter value") + ylab("Scaled density") +
    theme(legend.title = element_blank()) +
    scale_color_manual(values = c("lightgreen", "lightseagreen")) +
    ggtitle(plotTitle)

}
