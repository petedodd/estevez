## read in arguments
## R --slave --vanilla --args < inference.R ZWE
args <- commandArgs(trailingOnly = TRUE)
cniso <- as.character(args[1])

useage <- TRUE  #whether to use age-based likelihood or not
makePR <- FALSE #don't include unnecessary PR ODEs for inference

cat('cniso  = ',cniso,'...\n')
## stop('testing')

## libraries/dependencies
library(here)
library(reticulate)
zs <- import('zeus')

tgt <- 'c'
source(here("R/modelprep.R")) #above is actually in here NOTE est vn!
source(here("R/dataLL.R"))    #for dataLL & prior
load(here('data/irrAOdata.Rdata')) #data for bilinear interpolations

## set country here
cniso <- glue(cniso)
cz <- c('ETH','KEN','LSO','MOZ','MWI','NGA',
        'SWZ','TZA','UGA','ZAF','ZMB','ZWE')
## cniso <- cz[1]

## prep for this country
aimpars <- getCountryAimParmsInc(cniso,eyear = 2019)
## aimout <- AimDyns(aimpars,graph=FALSE,fullhist = TRUE) #
AO <- irrAOdata[[cniso]]                #bilinear interpolation for this country

##prepare ODE parms
estinputs <- make.estinputs(AO,aimpars,contacts = TRUE)

## -- cdrt with mortality
mortmod <- estebin(user=estinputs)
out <- runwithpars(list(),estinputs)

## ====================== inference ==========
## NOTE changes to controlled with flag
## modelprep/runwithpars
## modelprep/multirun
## dataLL/ps2l

## test priors
u0 <- rep(0.5,length(hyperparms)) #unit cube space
(psv <- uv2ps(u0))                #parameter space (vector)
psl <- ps2l(psv)                #parameter space (list)

## ------------- likelihoods ------------

## make the data LL function
dll <- makeDataLL(cniso,out[,'t'])
dll(out,age=FALSE) #test
dll(out,age=TRUE,show=TRUE) #test

## ss <- qfun(0.5,hyperparms[['sigma']])
## dll(out,list(sigma=1/3),show=TRUE) #test


## likelihood operating on parameter space vector
LLpsv <- function(x){
  pl <- ps2l(x)
  tmp <- runwithpars(pl,estinputs)
  ans <- dll(tmp,age=useage)
  if(!is.finite(ans)) ans <- -Inf #safety
  ans
}

## test
LLpsv(psv)


## likelihood operating on parameter space vector WITH PRIOR
LLpsvWithPrior <- function(x){
  y <- LLpsv(x) + prlogd(x)
  ifelse(is.finite(y),y,-1e+10)
}

LLpsvWithPrior(psv)


## stop()
## --------- inference work --------
ndim <- length(hyperparms)
print(ndim)
(NLIVE <- 2*max(50*ndim,ndim*(ndim+1)/2))
## also likes nlive > ndim*(ndim+1)/2

## set up zeus sampler
nwalkers <- 50

## sample from prior to start
S <- matrix(nrow=nwalkers,ncol=ndim)
for(i in 1:nwalkers) S[i,] <- uv2ps(0.5+runif(ndim)/1e2)

## ## test S
## SL <- rep(NA,nwalkers)
## for(i in 1:nwalkers) SL[i] <- LLpsvWithPrior(S[i,])
## summary(SL)

## make sampler
sampler <- zs$EnsembleSampler(as.integer(nwalkers),
                              as.integer(ndim),
                              LLpsvWithPrior)

## ---------------------------------------
## run sampler - shaping to be up to 20 hrs or so
nsteps <- 2000
tt <- system.time({
  sampler$run_mcmc(start=S,nsteps=as.integer(nsteps))
})

## save stats
catfn <- glue(here('plots/stats_{cniso}.txt'))
stats <- c(time=tt[3]/60,ESS=sampler$ess)
cat(stats,file=catfn,append = FALSE);cat("\n",file=catfn,append = TRUE)

## save sample
sampe <- sampler$get_chain(flat=TRUE,discard=as.integer(nsteps/2))
sampe <- as.data.table(sampe)
lbz <- names(hyperparms)
names(sampe)[1:length(lbz)] <- lbz
fno <- glue(here('data/sampe_{cniso}.Rdata'))
save(sampe,file=fno)   ## save out

## save arrays also
fno <- glue(here('data/chains/CH_{cniso}.Rdata'))
CH <- sampler$chain
if(!file.exists(here('data/chains'))) dir.create(here('data/chains'))
save(CH,file=fno)   ## save out

