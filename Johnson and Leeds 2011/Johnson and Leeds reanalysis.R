# Joshua Alley
# Texas A&M University
# Replication and check of Johnson and Leeds 2011 
# model of alliance participation and conflict onset


# Load packages
library(here)
library(arm)
library(MASS)
library(dplyr)
library(ggplot2)
library(separationplot)
library(sandwich)
library(lmtest)
library(xtable)
library(rstan)
library(shinystan)



# Set working directory to current folder 
setwd(here::here())
getwd()


# Set up STAN guidelines
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Set seed
set.seed(12)

# Load data
jl.2011 <- read.csv("Johnson and Leeds 2011/johnson-leeds-2011.csv")


### Robustness Check: cluster robust variance estimator 

# Start by creating a dataframe with complete cases
jl.2011.comp <- select(jl.2011, dispute, ptargdef, pchalally, pchaloff, pchalneu,
                       ln_distance, capprop, jdem, s_un_glo, 
                       peaceyrs, peaceyrs2, peaceyrs3, challenger, target, ddyad)
jl.2011.comp <- jl.2011.comp[complete.cases(jl.2011.comp), ]


# Run the model: use a logit link instead of probit. 
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
print(results.jl.corr)


xtable(jl.2011.model)
xtable(results.jl.corr)

# Correcting the standard errors means that Target and challenger defensive pacts 
# are no longer associated with MID initiation at conventional levels of statistical significance






### Fit varying intercept model 

# STAN Model

# Pull IVs and rescale so posterior variances are not wildly different
ivs.jl <- jl.2011.comp[, 2: 12]
ivs.jl[1: ncol(ivs.jl)] <- lapply(ivs.jl[1: ncol(ivs.jl)], 
                                          function(x) rescale(x, binary.inputs = "0/1"))
ivs.jl <- as.matrix(ivs.jl)

# create a challenger state index variable
jl.2011.comp$challenger.id <- jl.2011.comp %>% group_indices(challenger)
# Create a target index variable 
jl.2011.comp$target.id <- jl.2011.comp %>% group_indices(target)

# Create a dyad index variable
jl.2011.comp$dyad.id <- jl.2011.comp %>% group_indices(ddyad)


# Assign data and place in list
stan.data.jl <- list(N = nrow(jl.2011.comp), y = jl.2011.comp$dispute, 
                         X = ivs.jl, K = ncol(ivs.jl),
                         chall = jl.2011.comp$challenger.id, 
                         C = length(unique(jl.2011.comp$challenger.id)), 
                         targ = jl.2011.comp$target.id,
                         T = length(unique(jl.2011.comp$target.id)),
                         dyad = jl.2011.comp$dyad.id,
                         D = length(unique(jl.2011.comp$dyad.id))
)

# compile STAN code
model.1 <- stan_model(file = "VI Stan Model.stan")

# Run variational Bayes approximation: this did not converge
vb.jl <- vb(model.1, data = stan.data.jl, seed = 12)
launch_shinystan(vb.jl)


# Regular STAN: this will take a long time. 
# system.time(
#  vi.model.jl <- sampling(model.1,
#                              data = stan.data.jl, 
#                              iter = 2000, warmup = 1000, chains = 4
#  )
#)

# Check convergence diagnostics
#check_hmc_diagnostics(vi.model.jl)

vi.summary.jl <- summary(vb.jl, pars = c("beta", "alpha"), probs = c(0.025, 0.975))$summary
rownames(vi.summary.jl) <- c("ptargdef", "pchalally", "pchaloff", "pchalneu",
                               "ln_distance", "capprop", "jdem", "s_un_glo", 
                               "peaceyrs", "peaceyrs2", "peaceyrs3", 
                            "constant")
print(vi.summary.jl)


