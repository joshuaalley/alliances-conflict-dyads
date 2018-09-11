# Joshua Alley
# Texas A&M University
# Replication and check of Benson's 2011 model of conflict onset


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
benson.2011 <- read.csv("Benson 2011/Benson 2011 data.csv")


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


# Create data frame of complete observations 
# Necessary for to correct for dyadic clustering
benson.2011.comp <- select(benson.2011, ddyad, year, ccode1, ccode2, mzinit, jointdem, cont,
                           caprat2, sglo, uncondcompel, condcompel, uncondeter,
                             condeter, probabilistic, peaceyears, peace2, peace3)

benson.2011.comp <- benson.2011.comp[complete.cases(benson.2011.comp), ]

# Fit the model
fit.benson <- glm(mzinit ~ jointdem + cont + caprat2 + sglo + uncondcompel + condcompel + uncondeter + 
             condeter + probabilistic + peaceyears + peace2 + peace3,
           family = binomial(link = "logit"), x = T, data = benson.2011.comp)
# Create an index for the dyad members
index.sc <- unique(c(benson.2011.comp$ccode1, benson.2011.comp$ccode2))
index <- index.sc[order(index.sc)]
dyad.mat <- benson.2011.comp[,c("ccode1","ccode2")]


# Correct the standard errors
# Dyadic cluster robust via multiway decomposition
for(i in 1:length(index)){
  iUp <- index[i]
  clusUp <- apply(dyad.mat, 
                  1,
                  function(x)as.numeric(iUp %in% x))
  clusIndexUp <- clusUp*-99 + (1-clusUp)*1:nrow(dyad.mat)
  if(i==1){dcrUp <- robust.se.nodfc(fit.benson, clusIndexUp)[[1]]}
  if(i>1){dcrUp <- dcrUp + robust.se.nodfc(fit.benson, clusIndexUp)[[1]]}	
  cat(paste(iUp, " ")); flush.console()
}

# substract naive CR:
dcrUp2 <- dcrUp - robust.se.nodfc(fit.benson, benson.2011.comp$ddyad)[[1]]
# substract HR:
Vhat <- dcrUp2 - (length(index)-2)*vcovHC(fit.benson,type="HC0")
Vhat <- as.matrix(Vhat)

# Standard errors and test statistics
dcr_se <- sqrt(diag(Vhat))
stat.sig <- abs(coef(fit.benson)/sqrt(diag(Vhat))) > 1.64

V1 <- vcov(fit.benson)
se_nc <- sqrt(diag(V1))


# Export the results
results.benson.corr <- cbind(fit.benson$coefficients, dcr_se, stat.sig)
print(results.benson.corr)


xtable(fit.benson)
xtable(results.benson.corr)

# Corrrecting for state membership across dyads increases the standard errors,
# so conditional and unconditional deterrent treaties have no effect on conflict risk




### Fit varying intercept model 

# STAN Model

# Pull IVs and rescale so posterior variances are not wildly different
ivs.benson <- benson.2011.comp[, 6: ncol(benson.2011.comp)]
ivs.benson[1: ncol(ivs.benson)] <- lapply(ivs.benson[1: ncol(ivs.benson)], 
                     function(x) rescale(x, binary.inputs = "0/1"))
ivs.benson <- as.matrix(ivs.benson)

# create a challenger state index variable
benson.2011.comp$challenger.id <- benson.2011.comp %>% group_indices(ccode1)
# Create a target index variable 
benson.2011.comp$target.id <- benson.2011.comp %>% group_indices(ccode2)

# Create a dyad index variable
benson.2011.comp$dyad.id <- benson.2011.comp %>% group_indices(ddyad)


# Assign data and place in list
stan.data.benson <- list(N = nrow(benson.2011.comp), y = benson.2011.comp$mzinit, 
                         X = ivs.benson, K = ncol(ivs.benson),
                         chall = benson.2011.comp$challenger.id, 
                         C = length(unique(benson.2011.comp$challenger.id)), 
                         targ = benson.2011.comp$target.id,
                         T = length(unique(benson.2011.comp$target.id)),
                         dyad = benson.2011.comp$dyad.id,
                         D = length(unique(benson.2011.comp$dyad.id))
                         )

# compile STAN code
model.1 <- stan_model(file = "VI Stan Model.stan")

# Run variational Bayes approximation
vb.benson <- vb(model.1, data = stan.data.benson, seed = 12)
launch_shinystan(vb.benson)

# Regular STAN: 10 hour runtime on my laptop
system.time(
  vi.model.benson <- sampling(model.1,
                          data = stan.data.benson, 
                          iter = 2000, warmup = 1000, chains = 4
  )
)

# Check convergence diagnostics
check_hmc_diagnostics(vi.model.benson)

vi.summary.benson <- summary(vi.model.benson, pars = c("beta", "alpha"), probs = c(0.025, 0.975))$summary
rownames(vi.summary.benson) <- c("jointdem", "cont", "caprat2", "sglo", 
                            "uncondcompel", "condcompel", "uncondeter", 
                              "condeter", "probabilistic", 
                              "peaceyears", "peace2", "peace3", 
                              "constant")
print(vi.summary.benson)





