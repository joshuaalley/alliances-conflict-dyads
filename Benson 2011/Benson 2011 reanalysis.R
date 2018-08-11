# Joshua Alley
# Texas A&M University
# Replication and check of Benson 2011's model of conflict onset


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
benson.2011 <- read_dta("Benson 2011/ajps2011 replication.dta")


# Logit model 5 from first table (w/o cluster robust se)
benson.m5 <- glm(mzinit ~ jointdem + cont + caprat2 + sglo + uncondcompel + condcompel + uncondeter + 
                   condeter + probabilistic + peaceyears + peace2 + peace3,
                 family = binomial(link = "logit"), data = benson.2011, x = TRUE)
summary(benson.m5)

# separation plot to assess model fit
# Get the predicted probabilities
benson.pred.prob <- as.numeric(predict.glm(benson.m5, type = "response"))
benson.outcome <- as.numeric(benson.m5$y)
separationplot(benson.pred.prob, benson.outcome, type = "bands")


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
benson.2011.comp <- select(benson.2011, ddyad, ccode1, ccode2, mzinit, jointdem, cont,
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


library(xtable)
xtable(fit.benson)
xtable(results.benson.corr)

# Corrrecting for state membership across dyads increases the standard errors,
# so conditional and unconditional deterrent treaties have no effect on conflict risk










