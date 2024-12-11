plotPostPrior_logit = function (resList, priorList, postList,
                                parNames, preyNames, plotTitle) {

  priors = list()
  nPrey  = length(resList)
  nPar   = length(parNames)
  
  # extract priors from models
  priors = list()
  for (i in 1:nPrey) priors[[i]] = rstan::extract(resList[[i]], pars = priorList)
  
  # mutate prey and parameter names
  for (i in 1:nPrey) {
    for (j in 1:nPar) {
      priors[[i]][[j]] = as.data.frame(priors[[i]][[j]]) %>%
        mutate(Prey      = preyNames[i]) %>%
        mutate(Parameter = parNames[j])
    }
  }
  
  # bind dataframes inside nested list
  priors = do.call(Map, c(f = rbind, priors))
  priors = do.call(rbind, priors)
  priors = priors %>% mutate(Source = "Prior")
  names(priors) = c("Value", "Prey", "Parameter", "Source")
  
  # extract posteriors from models
  posts = list()
  for (i in 1:nPrey) posts[[i]] = rstan::extract(resList[[i]], pars = postList)
  
  # mutate prey and parameter names
  for (i in 1:nPrey) {
    for (j in 1:nPar) {
      posts[[i]][[j]] = as.data.frame(posts[[i]][[j]]) %>%
        mutate(Prey      = preyNames[i]) %>%
        mutate(Parameter = parNames[j])
    }
  }
  
  # bind dataframes inside nested list
  posts = do.call(Map, c(f = rbind, posts))
  posts = do.call(rbind, posts)
  posts = posts %>% mutate(Source = "Posterior")
  names(posts) = c("Value", "Prey", "Parameter", "Source")
  
  # plotting
  mergedDF = rbind(priors, posts)
  mergedDF$Parameter = factor(mergedDF$Parameter,
                              labels = c(bquote(alpha),
                                         bquote("beta[TL]")))
  mergedDF %>%
    ggplot() +
    stat_density(aes(x = Value,
                     y = ..scaled..,
                     color = Source),
                 position = "dodge",
                 geom = "density",
                 alpha = 0.3,
                 linewidth = 0.7)+
    facet_grid(Prey~Parameter, scales = "free", labeller = label_parsed) +
    theme_custom() +
    xlab("Parameter value") + ylab("Scaled density") +
    theme(legend.title = element_blank()) +
    scale_color_manual(values = c("lightgreen", "lightseagreen")) +
    ggtitle(plotTitle)
  
} # end of function