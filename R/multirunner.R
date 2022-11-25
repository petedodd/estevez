## to do MR data afterwards
## use:
## R --slave --vanilla --args < multirunner.R ETH
library(here)
source(here('R/argumenthandler.R')) #parse arguments

## cniso <- 'ETH'
## stop('testing')

## neater version of running
useage <- TRUE #NOTE don't think matters here
makePR <- TRUE #whether to include PR code or not
tgt <- 'c'
source(here("R/modelprep.R"))   #above is actually in here
source(here("R/dataLL.R"))      #for dataLL & prior
source(here("R/plotters.R"))    #for utilities to grab PR
load(here('data/irrAOdata.Rdata')) #data for bilinear interpolations
gh <- function(x) glue(here(x))
setDTthreads(threads = 1) #because running script already multithreaded - otherwise slower


## set country here
## cniso <- 'SWZ'
## cniso <- 'LSO'
## cniso <- 'NGA'
cniso <- glue(cniso)
## cz <- c('ETH','KEN','LSO','MOZ','MWI','NGA','SWZ','TZA','UGA','ZAF','ZMB','ZWE')

## prep for this country
aimpars <- getCountryAimParmsInc(cniso,eyear = 2019)
AO <- irrAOdata[[cniso]]
estinputs <- make.estinputs(AO,aimpars,contacts=TRUE) #prepare ODE parms
mortmod <- estebin(user=estinputs)
out <- runwithpars(list(),estinputs)

OD <- out2df(out)
save(OD,file=gh('data/OD_{cniso}.Rdata'))

## read in and save out
fno <- glue(here('data/sampe_{cniso}{sensitivity.analysis}.Rdata'))
load(file=fno)

## correlations
scr <- cor(sampe)
fno <- glue(here('data/scr_{cniso}{sensitivity.analysis}.Rdata'))
save(scr,file=fno)

## multiple runs
test <- sampe[sample(nrow(sampe),300)]
MR <- multirun(test,estinputs,bothout=TRUE,YR = 2019)
fno <- glue(here('data/MR_{cniso}{sensitivity.analysis}.Rdata'))
save(MR,file=fno)
