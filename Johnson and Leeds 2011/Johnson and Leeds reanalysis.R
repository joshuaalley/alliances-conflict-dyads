# Joshua Alley
# Texas A&M University
# Replication and check of Johnson and Leeds 2011 
# model of alliance participation and conflict onset


# Load packages
library(here)
library(haven)
library(MASS)
library(dplyr)
library(ggplot2)
library(separationplot)
library(sandwich)
library(lmtest)



# Set working directory to current folder 
setwd(here::here())
getwd()


# Load data
jl.2011 <- read.csv("Johnson & Leeds 2011/johnson-leeds-2011.csv")

# Replicate their probit results
jl.2011.probit <- glm(dispute ~ ptargdef + pchalally + pchaloff + pchalneu +
                        ln_distance + capprop + jdem + s_un_glo + 
                        peaceyrs + peaceyrs2 + peaceyrs3,
                      family = binomial(link = "probit"),
                      data = jl.2011)


### Robustness Check: cluster robust variance estimator 

# Start by creating a dataframe with complete cases
jl.2011.comp <- select(jl.2011, dispute, ptargdef, pchalally, pchaloff, pchalneu,
                       ln_distance, capprop, jdem, s_un_glo, 
                       peaceyrs, peaceyrs2, peaceyrs3, challenger, target, ddyad)
jl.2011.comp <- jl.2011.comp[complete.cases(jl.2011.comp), ]


# Run the model: use a logit in place of probit. 
# Similar inferences: marginal increase in p-value for challenger and target defense pacts
jl.2011.model <- glm(dispute ~ ptargdef + pchalally + pchaloff + pchalneu +
                       ln_distance + capprop + jdem + s_un_glo + 
                       peaceyrs + peaceyrs2 + peaceyrs3,
                     family = binomial(link = "logit"),
                     data = jl.2011.comp, x = T)
summary(jl.2011.model)


# Implement Aronow et al clustered standard errors 
# accounts for participation of states in multiple dyads 
# Cluster robust variance estimation function
robust.se.nodfc <- function(model, cluster){
  require(sandwich)
  require(lmtest)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- model$rank
  dfc <- 1
  uj <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum))
  rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
  rcse.se <- coeftest(model, rcse.cov)
  return(list(rcse.cov, rcse.se))
}


# Create an index for the dyad members
jl.index.sc <- unique(c(jl.2011.comp$challenger, jl.2011.comp$target))
jl.index <- index.sc[order(index.sc)]
jl.dyad.mat <- jl.2011.comp[,c("challenger","target")]


# Correct the standard errors
# Dyadic cluster robust via multiway decomposition
for(i in 1:length(jl.index)){
  iUp <- jl.index[i]
  jl.clusUp <- apply(jl.dyad.mat, 
                     1,
                     function(x)as.numeric(iUp %in% x))
  jl.clusIndexUp <- jl.clusUp*-99 + (1-jl.clusUp)*1:nrow(jl.dyad.mat)
  if(i==1){jl.dcrUp <- robust.se.nodfc(jl.2011.model, jl.clusIndexUp)[[1]]}
  if(i>1){jl.dcrUp <- jl.dcrUp + robust.se.nodfc(jl.2011.model, jl.clusIndexUp)[[1]]}	
  cat(paste(iUp, " ")); flush.console()
}

# substract naive CR:
jl.dcrUp2 <- jl.dcrUp - robust.se.nodfc(jl.2011.model, jl.2011.comp$ddyad)[[1]]
# substract HR:
jl.Vhat <- jl.dcrUp2 - (length(index)-2)*vcovHC(jl.2011.model,type="HC0")
jl.Vhat <- as.matrix(jl.Vhat)

# Standard errors and test statistics
jl.dcr_se <- sqrt(diag(jl.Vhat))
jl.stat.sig <- abs(coef(jl.2011.model)/sqrt(diag(jl.Vhat))) > 1.64

jl.V1 <- vcov(jl.2011.model)
jl.se_nc <- sqrt(diag(jl.V1))


# Export the results
results.jl.corr <- cbind(jl.2011.model$coefficients, jl.dcr_se, jl.stat.sig)
results.jl.corr

library(xtable)
xtable(jl.2011.model)
xtable(results.jl.corr)

# Correcting the standard errors means that Target and challenger defensive pacts 
# are no longer associated with MID initiation at conventional levels of statistical significance
