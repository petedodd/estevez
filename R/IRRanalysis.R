## R --slave --vanilla --args < IRRanalysis.R iso3
args <- commandArgs(trailingOnly = TRUE)
cniso <- as.character(args[1])

## libraries
library(here)
library(data.table)
library(ggplot2)
library(lhs)
library(glue)

gh <- function(x) glue(here(x))


source(here('R/AIMdynInc.R'))              #AIM runner
source(here('R/getCountryParmsInc.R'))  #parameter constructer
load(here('indata/XD.Rdata'))           #data for AIM


## create any missing directories
fn <- gh('plots/IRRplots')
if(!file.exists(fn)){cat('creating ',fn,'...\n'); dir.create(fn);}
fn <- gh('data/IRRdata')
if(!file.exists(fn)) {cat('creating ',fn,'...\n'); dir.create(fn);}


## bilinear interpolation functions
bl <- function(F11,F12,F21,F22,x1,x2,y1,y2,x,y){
  A <- (x2-x1)*(y2-y1)
  (F11*(x2-x)*(y2-y) + F12*(x2-x)*(y-y1) + F21*(x-x1)*(y2-y) + F22*(x-x1)*(y-y1)) / A
}

## i  = alpha, j = HR
## list of lists

## biliniear interpolator, including finding the right square
BLI <- function(LL,                     #list of lists (to handle matrices/arrays)
                xz,yz,                  #knots
                x,y                     #new point
                ){
  i <- sum(x>xz) ; j <- sum(y>yz)
  i2 <- i+1 ; j2 <- j+1
  bl(LL[[i]][[j]], LL[[i]][[j2]], LL[[i2]][[j]], LL[[i2]][[j2]],
     xz[i],xz[i2],yz[j],yz[j2],
     x,y)
}


## fixing odd last-time NA
safety1 <- function(x){
  bad <- which(!is.finite(x))
  x[bad] <- x[bad-1] #LOCF
  x
}
safety <- function(X) apply(X,2,safety1)
checkfix <- function(X){
  nr <- length(X)
  nc <- length(X[[1]])
  for(i in 1:nr){
    for(j in 1:nc){
      tmp <- X[[i]][[j]]
      if(any(!is.finite(tmp))){
        X[[i]][[j]] <- safety(tmp)
      }
    }
  }
  X
}
checkfixes <- function(X,NMZ){
  for(nm in NMZ) X[[nm]] <- checkfix(X[[nm]])
  X
}

## snippet for near IRR testing
hrz <- seq(.1,.8,by=.1)
az <- seq(.1,.8,by=.05)

## ETH KEN LSO MOZ MWI NGA SWZ TZA UGA ZAF ZMB ZWE
cz <- c('ETH','KEN','LSO','MOZ','MWI','NGA','SWZ','TZA','UGA','ZAF','ZMB','ZWE')


testpart <- TRUE

if(cniso=="ALL"){
  ## loop over all and join
  irrAOdata <- list()
  for(i in 1:length(cz)){
    print(cz[i])
    fn <- glue(here('data/IRRdata/AO{cz[i]}.Rdata'))
    load(fn)
    AO <- checkfixes(AO,c('irrhf','irrhm','irraf','irram'))
    irrAOdata[[cz[i]]] <- AO
  }
  save(irrAOdata,file=here('data/irrAOdata.Rdata')) #NOTE large added to repo

} else {

  ## compute data
  cat("......",cniso,"\n")
  ## get country to run
  fn <- glue(here('data/IRRdata/AO{cniso}.Rdata'))
  aimpars <- getCountryAimParmsInc(cniso,eyear=2019)
  AO <- AimDyns(aimpars,graph=FALSE,alph=az,HR=hrz)
  save(AO,file=fn)
  ## load(fn)

  ## for one country, do test of approximation
  if(testpart){

    ## length(hrz)*length(az)                  #120

    set.seed(123)
    L <-randomLHS(50,2)
    tz <- seq(from=min(AO$yrz),to=max(AO$yrz),by=1)

    colnames(L) <- c('alph','HR')
    L[,'alph'] <- min(az) + diff(range(az)) * L[,'alph']
    L[,'HR'] <- min(hrz) + diff(range(hrz)) * L[,'HR']
    ## plot(L)

    ansaf <- ansam <- anshf <- anshm <- list()
    tstaf <- tstam <- tsthf <- tsthm <- list()
    for(i in 1:nrow(L)){
      print(i)
      ## run main model
      tst <- AimDyns(aimpars,graph=FALSE,fullhist = TRUE,alph=L[i,'alph'],HR=L[i,'HR'])
      ansaf[[i]] <- tst$irraf[time %in% tz]
      ansam[[i]] <- tst$irram[time %in% tz]
      anshf[[i]] <- tst$irrhf[time %in% tz]
      anshm[[i]] <- tst$irrhm[time %in% tz]
      ## approximations
      tst <- BLI(AO$irrhf,AO$alph,AO$HR,L[i,'alph'],L[i,'HR']) #HF
      tst <- as.data.table(exp(tst)); tst[,time:=tz]
      tsthf[[i]] <- tst
      tst <- BLI(AO$irrhm,AO$alph,AO$HR,L[i,'alph'],L[i,'HR']) #HM
      tst <- as.data.table(exp(tst)); tst[,time:=tz]
      tsthm[[i]] <- tst
      tst <- BLI(AO$irraf,AO$alph,AO$HR,L[i,'alph'],L[i,'HR']) #AF
      tst <- as.data.table(exp(tst)); tst[,time:=tz]
      tstaf[[i]] <- tst
      tst <- BLI(AO$irram,AO$alph,AO$HR,L[i,'alph'],L[i,'HR']) #AM
      tst <- as.data.table(exp(tst)); tst[,time:=tz]
      tstam[[i]] <- tst
    }

    ## add replicate label
    for( i in 1:nrow(L)){
      ansaf[[i]][,c('repn','calculation'):=.(i,'exact')]
      ansam[[i]][,c('repn','calculation'):=.(i,'exact')]
      anshf[[i]][,c('repn','calculation'):=.(i,'exact')]
      anshm[[i]][,c('repn','calculation'):=.(i,'exact')]
      ## approximations
      tsthf[[i]][,c('repn','calculation'):=.(i,'approximate')]
      tsthm[[i]][,c('repn','calculation'):=.(i,'approximate')]
      tstaf[[i]][,c('repn','calculation'):=.(i,'approximate')]
      tstam[[i]][,c('repn','calculation'):=.(i,'approximate')]
    }

    ## combining
    ansaf <- rbindlist(ansaf)
    ansam <- rbindlist(ansam)
    anshf <- rbindlist(anshf)
    anshm <- rbindlist(anshm)
    ## approximations
    tsthf <- rbindlist(tsthf)
    tsthm <- rbindlist(tsthm)
    tstaf <- rbindlist(tstaf)
    tstam <- rbindlist(tstam)

    ## annotating
    ansaf[,c('art','sex'):=.('yes','F')]
    ansam[,c('art','sex'):=.('yes','M')]
    anshf[,c('art','sex'):=.('no','F')]
    anshm[,c('art','sex'):=.('no','M')]
    tsthf[,c('art','sex'):=.('no','F')]
    tsthm[,c('art','sex'):=.('no','M')]
    tstaf[,c('art','sex'):=.('yes','F')]
    tstam[,c('art','sex'):=.('yes','M')]

    ## all together
    ALL <- rbindlist(list(ansaf,ansam,anshf,anshm,
                          tsthf,tsthm,tstaf,tstam
                          ))

    ALLM <- melt(ALL,id.vars = c('time','repn','calculation','art','sex'))
    ALLM <- ALLM[value>1e-6]
    ALLM <- dcast(ALLM, time + art + sex + repn + variable ~ calculation,value.var = 'value')

    ggplot(ALLM,aes(exact,approximate),) +
      geom_point(shape=1,alpha=0.3) +
      geom_abline(intercept = 0,slope=1,col=2) +
      coord_fixed() +
      theme_classic() + ggpubr::grids()

    nrow(ALLM) # 77K

    ggsave(file=gh('plots/IRRplots/IRRcompare{cniso}.pdf'),h=6,w=6)
    ggsave(file=gh('plots/IRRplots/IRRcompare{cniso}.png'),h=6,w=6)

    ALLM[,mean(abs(1-approximate/exact))]     #0.01
    save(ALLM,file=gh('data/IRRdata/IRRapprox{cniso}.Rdata'))

    ## ---- single parm example--------
    aimout <- AimDynInc(aimpars,graph=FALSE,fullhist = TRUE,
                     alph=0.36,HR=0.33)

    aprx <- BLI(AO$irrhf,AO$alph,AO$HR,0.36,0.33) #HF
    aprx <- as.data.table(exp(aprx)); aprx[,time:=tz]
    aprxm <- melt(aprx,'time')

    comp <- aimout$irrhf
    comps <- comp[time %in% tz]
    compsm <- melt(comps,'time')


    SV <- ggplot(compsm,aes(time,value,col=variable,group=variable)) +
      geom_line() +
      geom_point(data=aprxm) +
      ylab('Incidence rate ratio') + xlab('Year')+
      theme_classic() + ggpubr::grids()
    SV

    ggsave(SV,file=gh('plots/IRRplots/srIRR_{cniso}.pdf'),w=7,h=5)
    ggsave(SV,file=gh('plots/IRRplots/srIRR_{cniso}.png'),w=7,h=5)

  }

}

