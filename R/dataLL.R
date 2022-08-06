## making data likelihood & hyperparms

## === load data
load(here('indata/fitdata/NT.Rdata'))
load(here('indata/fitdata/NAS.Rdata'))
load(here('indata/fitdata/E.Rdata'))
load(here('indata/fitdata/P.Rdata'))
load(here('indata/fitdata/VR.Rdata'))
load(here('indata/fitdata/PA.Rdata'))
load(here('indata/fitdata/HF.Rdata'))

nagz <- c('0-4','5-14','15-24','25-34','35-44','45-54','55-64','65+')

## select and aggregate columns for incidence by age
getRateAge <- function(OM,term="IX|IH|IA"){ #NoN for notifications
  vrz <- grep(term,colnames(OM),value = TRUE)
  tmp <- OM[,vrz]
  nmz <- as.integer(gsub(".+\\[(\\d+),\\d+\\]","\\1",colnames(tmp)))
  ## key
  LK <- list()
  LK[['0-4']] <- (which(nmz==1))
  LK[['5-14']] <- (which(nmz %in% 2:3))
  LK[['15-24']] <- (which(nmz %in% 4:5))
  LK[['25-34']] <- (which(nmz %in% 6:7))
  LK[['35-44']] <- (which(nmz %in% 8:9))
  LK[['45-54']] <- (which(nmz %in% 10:11))
  LK[['55-64']] <- (which(nmz %in% 12:13))
  LK[['65+']] <- (which(nmz %in% 14:length(agz)))
  ## apply aggregations
  tmpo <- lapply(LK,function(x) rowSums(tmp[,x]))
  ## answer structure
  MO <- matrix(nrow=length(tmpo[[1]]),ncol=length(tmpo))
  colnames(MO) <- names(tmpo)
  for(nm in names(tmpo)) MO[,nm] <- tmpo[[nm]]
  MO <- cbind(t=OM[,'t'],MO)
  MO
}

## ## ## test
## test <- getRateAge(out)
## test <- getRateAge(out,term="NoN")
## head(test)

## function factory for data LL
## NOTE using midpoints of years
makeDataLL <- function(cniso,tz){
  ## utility:
  F <- function(x) which.min(abs(tz-x))
  ## notifications by time
  cat('...making notification data...\n')
  ntmp <- NT[cniso][!is.na(c_newinc),.(year,c_newinc)]
  whont <- unlist(lapply(ntmp$year + 0.5,F)) #which entries to match
  notes <- ntmp$c_newinc
  notesfac <- max(notes) #NOTE heuristic scaling
  notelen <- length(notes)
  natmp <- NAS[cniso][,.(value=sum(value)),by=.(year,acat)]
  setkey(natmp,year)
  natmp <- natmp[ntmp]
  natmp <- dcast(natmp,year ~ acat,value.var = 'value')
  nmtmp <- as.matrix(natmp[,..nagz]) #matrix by ages
  nmtmp[is.na(nmtmp)] <- 0
  tota <- rowSums(nmtmp,na.rm=TRUE)
  kf <- tota/notes #age-known frac
  ## HIV start/end years
  Hlgt <- HF[iso3==cniso,logitH]
  HlgtV <- HF[iso3==cniso,logitH2]
  hstarty <- HF[iso3==cniso][1]$ny
  hendy <- HF[iso3==cniso][1]$xy;
  hiz <- which(out[,'t']>hstarty & out[,'t']<hendy)

  if(cniso %in% P$iso3){
    cat('...making prevalence data...\n')
    ## ptmp <- P[cniso][case.type=='a' & sex=='all' & acat=='all']
    ptmp <- P[cniso][acat=='15+'& sex=='all' & case.type=='b']
    pys <- ptmp[,year]
    pwhot <- unlist(lapply(pys + 0.5,F))
    pmn <- ptmp[,prev]
    psd <- ptmp[,(prev.hi-prev.lo)/3.92]
    parr <- PA[iso3==cniso,RR]
    parr.sd <- PA[iso3==cniso,RR.sd]
  } else {ptmp <- NULL}
  ## TODO test 1st line always OK, including multi years

  function(out,## P, #removed as now no parameters needed here
           show=FALSE,age=TRUE){
    ## function(out,P,show=FALSE,age=TRUE){
    LL <- c(notes=0,prev=0,hiv=0,preva=0)
    ## notifications over time
    ## print(notes)
    enotes <- out[whont,'Ntot']
    ## LL[1] <-  - sum((enotes - notes)^2/(2*(notesfac*P$sigma)^2)) -
    ##   (notelen) * log(P$sigma)
    E <- sum((enotes - notes)^2)/(2*notesfac^2)
    ## age dependence
    if(age){
      ## --- age unknown part
      ## ku * notes = ku * enotes + e
      ## var(e) = ku * v
      ## p^2(x-y)^2/(2*p*blah)
      Eu <- sum((1-kf)*(enotes - notes)^2)/(2*notesfac^2)
      ## --- age known part
      modelnotes <- getRateAge(out,term="NoN")[whont,-1]
      SE <- (kf*modelnotes - nmtmp)^2 #multiplied by known frac
      ## want var multiplied by 1/ncol
      SSE <- rowSums(SE) / (2*notesfac^2/ncol(SE))
      Ek <- sum(SSE[kf>0]/kf[kf>0]) #see above notes
      if(show) print(c(Ek,Eu))
      E <- Ek + Eu
    }
    LL[1] <-  - (ig$alpha + (notelen+1)/2)*log(E+ig$beta)
    ## HIV frac
    hmn <- mean(out[hiz,'NtotH']/out[hiz,'Ntot'])
    hmn <- logit(hmn)
    LL['hiv'] <- -(Hlgt-hmn)^2/2*HlgtV

    ## point prevalence
    if(!is.null(ptmp)){
      eprev <- 1e5 * out[pwhot,'prevtot15plus'] / out[pwhot,'pop15plus']
      LL[2] <- -sum((eprev-pmn)^2/(2*psd^2))
      if(age){
        llar <- out[pwhot,c('prev2534','prev3544',
                            'prev4554','prev5564',
                            'prev65pl')]/(out[pwhot,'prev1524']+1e-6)
        lla <- -sum(((llar-PA[iso3==cniso,RR])/
                     PA[iso3==cniso,RR.sd])^2)/2
        LL['preva'] <- lla
      }
    }
    if(show) print(LL)
    sum(LL)
  }

}

## === priors etc

## uses hyerparms in dataLL - name determines distribution etc
qfun <- function(u,L){
  if(names(L)[1]=='meanlog') x <- qlnorm(u,L[[1]],L[[2]])
  if(names(L)[1]=='shape1') x <- qbeta(u,L[[1]],L[[2]])
  if(names(L)[1]=='mean') x <- qnorm(u,L[[1]],L[[2]])
  if(names(L)[1]=='shape') x <- qgamma(u,L[[1]],scale=L[[2]])
  if(names(L)[1]=='omega') x <- qsn(u,omega=L[[1]],alpha=L[[2]])
  x
}

## ## durations
## tst <- rlnorm(1e4,meanlog=1.1,sdlog=0.2)
## c(mean(tst),sd(tst))
## tst <- rlnorm(1e4,meanlog=-1.386,sdlog=0.75)
## c(mean(tst),sd(tst))

checker <- function(nm){
  print(qfun(c(.025,.25,.5,.75,.975),hyperparms[[nm]]))
  qplot(qfun(runif(1e4),hyperparms[[nm]]))
}

betlist <- list(meanlog=log(20),sdlog=0.75)         #bet
if(useage) betlist <- list(meanlog=log(10),sdlog=0.75)         #bet

## NOTE new hyperparms
hyperparms <- list(
  bet=betlist,                             #beta
  ari0=list(meanlog=log(3e-2),sdlog=0.75),         #ari0
  ## ari0=list(shape1=1,shape2=25),             #ari0
  ## ari0=list(shape1=1.5,shape2=50),             #ari0
  ## bet=list(meanlog=1.68,sdlog=0.37),         #bet
  ## ari0=list(shape1=5,shape2=150),             #ari0
  arig=list(meanlog=0.62, sdlog=0.068),      #arig Ragonnet
  pp=list(meanlog=-2.837,sdlog=0.32),            #pp Ragonnet
  eps=list(meanlog=-6.89,sdlog=0.58),           #eps Ragonnet
  alph=list(meanlog=-1.02,sdlog=0.219),      #alph: log(.36),sqrt(.048)
  HR=list(meanlog=-1.05,sdlog=0.115),      #HR (ART)
  CDR=list(shape1=2,shape2=2),               #CDR 4 6
  drnX=list(meanlog=1.1,sdlog=0.2),        #durnX log(3)
  ## drnH=list(meanlog=-1.386,sdlog=0.75),  #durnH log(.25)
  drnH=list(shape=7.374824,scale_alpha_binned=0.06497331), #Ku see sigma.R NOTE changed
  v=list(shape1=20.7,shape2=77.9),          #protn Andrews
  txf=list(shape1=2.71,shape2= 87.55),      #CFRs
  cfrn=list(shape1=25.48, shape2= 33.78),
  cfrpn=list(shape1=23.68,shape2= 6.68),
  cfrpp=list(shape1=11.88, shape2= 12.37),
  rel=list(meanlog=-3.95,sdlog=0.27),      #relapse Crampin NOTE notin
  ## sigma=list(shape=3.626696,scale=0.01040186) #gamm * notes v 2010
  ## sigma=list(shape=2.599076,scale=0.02578648), #gamm * notes v 2000
  ## sigma=list(shape=2.599076,scale=0.02578648*2), #gamm * notes v 2000
  ## cdrdt=list(mean=0,sd=5e-2)       #cdrdt trend
  ## ecdrdt=list(mean=0,sd=1e-1)       #exp cdr trend
  ## cdrdt=list(omega=5e-2,alpha=2)       #cdrdt trend skew normal
  ## cdrdt=list(shape=1,scale=2e-2)       #cdrdt trend exp
  cdrdt=list(shape=0.5,scale=4e-2)       #cdrdt trend exp NOTE was above
)
kidhps <- list(
  OR04=list(meanlog=-0.9995206,sdlog=0.6258456),        #OR CDR u5
  OR514=list(meanlog=-0.5668002,sdlog=0.4577016),       #OR CDR 514
  pp04=list(shape1=5.152793,shape2=21.96717)   #u5 progn
)
hyperparms <- c(hyperparms,kidhps)
## ig <- list(alpha=2,beta=0.05^2)
ig <- list(alpha=5,beta=(5-1)*0.05^2) #NOTE changed

## sigma=list(meanlog=10,sdlog=1),
## list(meanlog=log(300),sdlog=1.5),     #sigma
## cdrdt=list(meanlog=log(5e-2),sdlog=1)       #cdrdt trend

## checker('bet')
## checker('ari0')
## checker('CDR')
## hyperparms[['sigma']] <- list(shape=1,scale=1)
## checker('sigma')
## checker('rel')
## checker('txf')
## checker('pp04')
## checker('OR04')
## checker('OR514')
## checker('cdrdt')
## checker('drnH')

## tz <- seq(from=0,to=40,by=.1)
## y <- -5e-2*tz*exp(-4e-2*tz)
## ## plot(tz,y)
## plot(tz,1/(exp(-y)+1))



## map unit cube to parameter space using quantile functions
uv2ps <- function(u){
  for(i in 1:length(hyperparms)){
    u[i] <- qfun(u[i],hyperparms[[i]])
  }
  u[10] <- u[10]+0.1 #shift for this (durn H)
  u
}


logdfun <- function(x,L){
  if(names(L)[1]=='meanlog') y <- dlnorm(x,L[[1]],L[[2]],log = TRUE)
  if(names(L)[1]=='shape1') y <- dbeta(x,L[[1]],L[[2]],log = TRUE)
  if(names(L)[1]=='mean') y <- dnorm(x,L[[1]],L[[2]],log = TRUE)
  if(names(L)[1]=='shape') y <- dgamma(x,L[[1]],scale=L[[2]],log = TRUE)
  if(names(L)[1]=='omega') y <- dsn(x,omega=L[[1]],
                                    alpha=L[[2]],
                                    log = TRUE)
  y
}

## logdfun(10,hyperparms[['bet']])
## logdfun(0,hyperparms[['cdrdt']])


prlogd <- function(x){
  for(i in 1:length(hyperparms)){
    x[i] <- logdfun(x[i],hyperparms[[i]])
  }
  sum(x)
}

## parameter space vector to list for running
ps2l <- function(x){
  ## NOTE safety
  if(x[6] > max(AO$alph)) x[6] <- .999* max(AO$alph)
  if(x[6] < min(AO$alph)) x[6] <- 1.001* min(AO$alph)
  if(x[7] > max(AO$HR)) x[7] <- .999* max(AO$HR)
  if(x[7] < min(AO$HR)) x[7] <- 1.001* min(AO$HR)
  ## make list
  y <- as.list(x[1:(length(hyperparms))])
  names(y) <- names(hyperparms)
  if(is.null(y[['drnA']]) &&
     !is.null(y[['drnH']]) && !is.null(y[['drnX']])){
    y[['drnA']] <- sqrt(abs(y[['drnH']]*y[['drnX']])) #geom mean
  }

  ##logit-linear trend
  if(is.null(y[['ecdrdt']])) y[['ecdrdt']] <- 0
  if(is.null(y[['cdrdt']])) y[['cdrdt']] <- 0
  DT <- (estinputs$ttq - min(estinputs$ttq))
  lgtcdr.traj <- logit(y[['CDR']]) +
    DT * y[['cdrdt']] * exp(-DT * y[['ecdrdt']])
  y[['ecdrdt']] <- y[['cdrdt']] <- NULL
  y[['CDR']] <- invlogit(lgtcdr.traj)

  y
}
