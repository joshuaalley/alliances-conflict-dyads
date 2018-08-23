# Joshua Alley
# Texas A&M University
# Replication and check of Johnson and Leeds 2017 
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


# Load data
kv.2017 <- read.csv("Kenwick and Vasquez 2017/kenwick-vasquez-2017.csv")


### Robustness Check: cluster robust variance estimator 

# Start by creating a dataframe with complete cases
kv.2017.comp <- select(kv.2017, dispute, fbdef, nfbdef, pchaloff, pchalneu,
                       ln_distance, capprop, jdem, s_un_glo, 
                       peaceyrs, peaceyrs2, peaceyrs3, challenger, target, ddyad)
kv.2017.comp <- kv.2017.comp[complete.cases(kv.2017.comp), ]


# Run the model on complete data
kv.2017.model <- glm(dispute ~ fbdef + nfbdef + pchaloff +
                       pchalneu + ln_distance + capprop + jdem + s_un_glo +
                       peaceyrs + peaceyrs2 + peaceyrs3,
                     family = binomial(link = "logit"),
                     data = kv.2017.comp, x = T)
summary(kv.2017.model) # This replicates their results w/o corrected standard errors


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
kv.index.sc <- unique(c(kv.2017.comp$challenger, kv.2017.comp$target))
kv.index <- index.sc[order(index.sc)]
kv.dyad.mat <- kv.2017.comp[,c("challenger","target")]


# Correct the standard errors
# Dyadic cluster robust via multiway decomposition
for(i in 1:length(kv.index)){
  iUp <- kv.index[i]
  kv.clusUp <- apply(kv.dyad.mat, 
                     1,
                     function(x)as.numeric(iUp %in% x))
  kv.clusIndexUp <- kv.clusUp*-99 + (1-kv.clusUp)*1:nrow(kv.dyad.mat)
  if(i==1){kv.dcrUp <- robust.se.nodfc(kv.2017.model, kv.clusIndexUp)[[1]]}
  if(i>1){kv.dcrUp <- kv.dcrUp + robust.se.nodfc(kv.2017.model, kv.clusIndexUp)[[1]]}	
  cat(paste(iUp, " ")); flush.console()
}

# substract naive CR:
kv.dcrUp2 <- kv.dcrUp - robust.se.nodfc(kv.2017.model, kv.2017.comp$ddyad)[[1]]
# substract HR:
kv.Vhat <- kv.dcrUp2 - (length(index)-2)*vcovHC(kv.2017.model,type="HC0")
kv.Vhat <- as.matrix(kv.Vhat)

# Standard errors and test statistics
kv.dcr_se <- sqrt(diag(kv.Vhat))
kv.stat.sig <- abs(coef(kv.2017.model)/sqrt(diag(kv.Vhat))) > 1.64

kv.V1 <- vcov(kv.2017.model)
kv.se_nc <- sqrt(diag(kv.V1))


# Export the results
results.kv.corr <- cbind(kv.2017.model$coefficients, kv.dcr_se, kv.stat.sig)
print(results.kv.corr)


xtable(kv.2017.model)
xtable(results.kv.corr)





### Estimate model with varying intercepts for states, and dyads
# Pull IVs and rescale so posterior variances are not wildly different
ivs.kv <- kv.2017.comp[, 2: 12]
ivs.kv[1: ncol(ivs.kv)] <- lapply(ivs.kv[1: ncol(ivs.kv)], 
                                  function(x) rescale(x, binary.inputs = "0/1"))
ivs.kv <- as.matrix(ivs.kv)

# create a challenger state index variable
kv.2017.comp$challenger.id <- kv.2017.comp %>% group_indices(challenger)
# Create a target index variable 
kv.2017.comp$target.id <- kv.2017.comp %>% group_indices(target)

# Create a dyad index variable
kv.2017.comp$dyad.id <- kv.2017.comp %>% group_indices(ddyad)


# Assign data and place in list
stan.data.kv <- list(N = nrow(kv.2017.comp), y = kv.2017.comp$dispute, 
                     X = ivs.kv, K = ncol(ivs.kv),
                     chall = kv.2017.comp$challenger.id, 
                     C = length(unique(kv.2017.comp$challenger.id)), 
                     targ = kv.2017.comp$target.id,
                     T = length(unique(kv.2017.comp$target.id)),
                     dyad = kv.2017.comp$dyad.id,
                     D = length(unique(kv.2017.comp$dyad.id))
)

# compile STAN code
model.1 <- stan_model(file = "VI Stan Model.stan")

# Run variational Bayes approximation: this did not converge
vb.kv <- vb(model.1, data = stan.data.kv, seed = 12)
launch_shinystan(vb.kv)


# Regular STAN: this will take a long time. 
# system.time(
#  vi.model.kv <- sampling(model.1,
#                              data = stan.data.kv, 
#                              iter = 2000, warmup = 1000, chains = 4
#  )
#)

# Check convergence diagnostics
#check_hmc_diagnostics(vi.model.kv)

vi.summary.kv <- summary(vb.kv, pars = c("beta", "alpha"), probs = c(0.025, 0.975))$summary
rownames(vi.summary.kv) <- c("fbdef", "nfbdef", "pchaloff", "pchalneu",
                             "ln_distance", "capprop", "jdem", "s_un_glo", 
                             "peaceyrs", "peaceyrs2", "peaceyrs3", 
                             "constant")
print(vi.summary.kv)

