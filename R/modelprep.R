## trying to implement as simplified odin model
##
## library(here)
library(odin)
library(glue)
library(data.table)

## source files
source(here('R/plotters.R'))                    #plots and data utils
source(here('R/getCountryParmsInc.R'))          #parameter constructer
source(here('R/AIMdynInc.R'))                      #dynamics

## WPP demographic estimates
load(here('indata/ND80.Rdata'))       #pre-calculated outrates
load(here('indata/N80_simple.Rdata')) #population by age/sex/year
load(here('indata/BL.Rdata')) #births data

## other external data
load(here('indata/XD.Rdata')) #data from SPECTRUM/AIM model
load(here('indata/CML.Rdata')) #contact data from Prem et al


## load/compile model


fn <- here('R/estebinxkcdneat.R') #vn recovery + mixing
## estebin <- odin::odin(fn,target=tgt)
## NOTE to allow PR code to be compliled in only when needed
OF <- paste(readLines(fn), collapse="\n")
if(makePR){
  cat('***Adding in PR code!***\n')
  fn <- here('R/estebinxkcdPR19.R')
  PR19 <- paste(readLines(fn), collapse="\n")
  OF <- paste(OF,PR19,collapse="\n")
}
estebin <- odin::odin(OF,target=tgt)

cat('--- Mix version ready -----\n')


## ====== param functions ========

## logit <- function(x) log(x/(1-x))
logit <- function(x) ifelse(x<1e-8,-19,
                     ifelse(x>1-1e-8,19,log(x/(1-x)))) #safer
invlogit <- function(x) 1/(1+exp(-x))

odemopars <- function(cn,eyear=2010){
  tc <- 1970:eyear                         #coarse times NOTE
  ## ND[iso3==cntry & Year %in% tc]
  ## men
  omegaL <- ND80[iso3==cn & Year %in% tc,.(Year,AgeGrp,omegaM)]
  omegaW <- dcast(omegaL,Year~AgeGrp,value.var = 'omegaM')
  omega <- as.matrix(omegaW)
  omega <- omega[,-1]                     #
  row.names(omega) <- tc
  omegaM <- omega                             #NOTE
  bzm <- BL[iso3==cn & Year %in% tc,BirthsMale] #NOTE
  pinitM <- N80[iso3==cn & Year==min(tc),PopMale] #NOTE
  ## women
  omegaL <- ND80[iso3==cn & Year %in% tc,.(Year,AgeGrp,omegaF)]
  omegaW <- dcast(omegaL,Year~AgeGrp,value.var = 'omegaF')
  omega <- as.matrix(omegaW)
  omega <- omega[,-1]                     #
  row.names(omega) <- tc
  omegaF <- omega                             #NOTE
  bzf <- BL[iso3==cn & Year %in% tc,BirthsFemale] #NOTE
  pinitF <- N80[iso3==cn & Year==min(tc),PopFemale] #NOTE
  ## make parameter object
  list(ttp=tc,
       popdatF=omegaF,
       popdatM=omegaM,
       BF=bzf,BM=bzm,
       popinitM=pinitM,
       popinitF=pinitF
       )
}


irrpars <- function(aimd,n=17){
  nend <- n+1
  tc <- seq(from=min(aimd$yrz),to=max(aimd$yrz),by=1)
  if(is.null(aimd$alph) | (length(aimd$alph)==1 & length(aimd$HR)==1)){
      AF <- as.matrix(aimd$irraf[time %in% tc])[,-nend]
      AM <- as.matrix(aimd$irram[time %in% tc])[,-nend]
      HF <- as.matrix(aimd$irrhf[time %in% tc])[,-nend]
      HM <- as.matrix(aimd$irrhm[time %in% tc])[,-nend]
  } else {
      AF <- exp(aimd$irraf[[1]][[1]])
      AM <- exp(aimd$irram[[1]][[1]])
      HF <- exp(aimd$irrhf[[1]][[1]])
      HM <- exp(aimd$irrhm[[1]][[1]])
  }
  list(irrhfD=HF,irrhmD=HM,irrafD=AF,irramD=AM)
}


## parameter maker
make.estinputs <- function(aimout,aimpars,contacts=FALSE){
  ## fixed
  nage <- 17                              #number of age categories
  CA <- matrix(0,ncol=81,nrow=17)
  for(i in 1:80) CA[((i-1) %/% 5)+1,i] <- 1
  CA[nrow(CA),81] <- 1
  colnames(CA) <- XD$nmage
  xz <- seq(from=0,to=80,by=5)
  rownames(CA) <- paste0('[',xz,',',c(xz[-1]-1,Inf),')')
  CAV <- colSums((1:nrow(CA)) * CA)
  CAV
  agz <<- rownames(CA)
  agz[4:10]
  ## mapping from year band ages to five years
  nmmap <- list()
  k <- 0
  for(i in 1:nrow(CA)){
    for(j in which(CA[i,]>0)){
      k <- k+1
      nmmap[[k]] <- data.table(age=colnames(CA)[j],ages=rownames(CA)[i])
    }
  }
  nmmap <- rbindlist(nmmap)

  ## making data for odin
  agen <- names(aimout$moUm)[1:nage]
  tc <- aimout$moUm$time                  #times
  tcs <- min(aimpars$yrz)-1 + 1:length(aimpars$nbmv)
  summary(aimout$moUm[,..agen])

  odepz <- odemopars(aimpars$iso3,eyear = aimpars$eyear)

  ## tc <- 1970:aimpars$eyear

  ## mortality 
  moUmD <- as.matrix(aimout$moUm[,..agen])
  moHmD <- as.matrix(aimout$moHm[,..agen])
  moAmD <- as.matrix(aimout$moAm[,..agen])
  moUfD <- as.matrix(aimout$moUf[,..agen])
  moHfD <- as.matrix(aimout$moHf[,..agen])
  moAfD <- as.matrix(aimout$moAf[,..agen])
  moTmD <- as.matrix(odepz$popdatM)
  moTfD <- as.matrix(odepz$popdatF)
  moUmD <- moTmD
  moUfD <- moTfD

  popinitF <- odepz$popinitF*1e3
  popinitM <- odepz$popinitM*1e3


  ## HIV
  hitD <- aimpars$hivtarget
  hitD2 <- c(0,diff(hitD))                #deriv correction
  hsrD <- aimpars$srv
  hivbyage <- as.data.table(aimpars$AR)                              #aggregate
  hivbyage[,cav:=CAV]
  hivbyage <- hivbyage[,.(M=mean(M),F=mean(F)),by=cav]
  hivbyage <- as.matrix(hivbyage[,.(M,F)])
  ## ART by sex/age
  hatD <- aimpars$arttarget
  hatD2 <- c(0,diff(hatD))                #deriv correction
  ## make finer grained ART target
  hatfun <- splinefun(tcs,sqrt(hatD),method='natural')
  hatDL <- hatfun(tc)^2
  artsyear <- max(tcs[hatD==0])
  hatDL[tc<=artsyear] <- 0
  ##plot(tcs,hatD);lines(tc,hatDL,col=2); points(tcs[hatD==0],hatD[hatD==0],pch='x',col=4)


  ## population data
  WW <- as.data.table(aimout$W)
  names(WW) <- c('time','sex','age','CD4','ART','value')
  WW[,year:=aimpars$syear + (as.numeric(as.character(time))-1)*aimpars$tstep,]
  WW[,time:=NULL]
  WW <- merge(WW,nmmap,by='age',all.x=TRUE)
  WW[,ages:=factor(ages,levels=unique(ages),ordered=TRUE)]
  WW[,art:=ART]
  WW[art!='none',art:='ART+ve']
  WW[art=='none',art:='ART-ve']
  WW[,hiv:=CD4]
  WW[hiv!='hiv-ve',hiv:='hiv+ve']
  WWP <- WW[,.(pop=sum(value)),by=.(year,sex,ages,hiv,art)]

  ## ART starts
  tmp <- melt(aimout$ram,id='time')
  tmp[,sex:='M']
  tmp2 <- melt(aimout$raf,id='time')
  tmp2[,sex:='F']
  tmp <- rbind(tmp,tmp2)
  names(tmp)[2] <- 'ages'
  names(tmp)[1] <- 'year'
  tmp[,ages:=factor(ages,levels=unique(ages),ordered=TRUE)]

  ## rates
  tmp <- merge(tmp,
               WWP[hiv=='hiv+ve' & art=='ART-ve',.(year,sex,ages,pop)],
               by=c('year','ages','sex'))
  tmp[,vr:=1e3*value/(pop+1e-10)]
  tmp[,vrn:=vr/(max(vr)),by=.(year)]#,sex)] # normalize
  tmp2 <- tmp[is.finite(vrn),.(rra=mean(vrn,na.rm=TRUE)),by=.(sex,ages)]
  tmp2c <- dcast(tmp2,ages~sex,value.var='rra')
  rra <- as.matrix(tmp2c[,.(M,F)])         #relative rate of ART starts
  dim(moUfD)
  length(tcs)

  ## put into parameter object
  parms <- list(ttp=tc,                    #times for interpolation
                ttq=tcs,
                tscale=5e-2,              #for pusuit
                popinitM=popinitM,
                popinitF=popinitF,
                ## BM=aimpars$nbmv,#births
                ## BF=aimpars$nbfv,
                BF=odepz$BF*1e3,
                BM=odepz$BM*1e3,
                muxfD=moUfD,   #mortality
                muxmD=moUmD,
                muhfD=moHfD,
                muhmD=moHmD,
                muafD=moAfD,
                muamD=moAmD,
                hitD=hitD,               #HIV prevalence
                ## hitD2=hitD2,               #HIV prevalence
                hsrD=hsrD,               #HIV SR over time
                hivbyage=hivbyage,       #distribution
                hatD=hatD,               #ART targets
                rra=rra                  #ART IRR by age
                )
  parms$hatD <- parms$hatD + hatD2 * parms$tscale 
  parms$hitD <- parms$hitD + hitD2 * parms$tscale 
  irrp <- irrpars(aimout)

  parms <- c(parms,irrp) #join IRR
  agzm <- 2.5 + seq(from=0,by=5,len=length(agz))

  ## sensitivity analysis for relative infectiousness in 10-14 year olds
  RI <- 0.0
  if(sensitivity.analysis=='INF1014'){
    cat('---!! setting infectiousness in 10-14 yos  to 0.5 !!---\n')
    RI <- 0.5
  }

  ## CDR over time (ttq)
  parms$CDRdata <- rep(0.6,length(parms$ttq))
  ## NOTE could do via qfun and hyerparms
  ## TB parameters
  tbparz2L <- list(drnX = 3,
               drnH = 0.25,
               drnA = 0.25,             #careful
               pp = 2.7e-4 * 365,       #raggonet
               arig = 5.4e-3 * 365,     #raggonet
               ari0 = 5e-2,
               eps = 3.3e-6 * 365,
               v=0.3,
               bet=10,
               ## CDR=0.6,
               txf = 0.1,   #ATT + (all HIV)
               cfrn = 0.5,  #ATT-, HIV-
               cfrpn = 0.9, #ATT-, HIV+/ART-
               cfrpp = 0.5, ##ATT-, HIV+/ART+
               rel=2e-2,    #relapse
               OR04=0.4,    #OR CDR u5
               OR514=0.6,   #OR CDR 514
               pp04=0.2,    #u5 progn
               rinf1014=RI, #10-14 relative infectionsness
               am=agzm)
  parms <- c(parms,tbparz2L)
  if(contacts){
    MM <- CML[[aimpars$iso3]]
    MM <- MM / mean(rowSums(MM)) #normalize
    parms[['MM']] <- MM
  }

  return(parms)
}


## bilinear interpolation
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
## NOTE edge cases? just use truncated distribution - alwasy strictly inside square


## function to correct alpha/HR
amendAHR <- function(ao,     #aimout object with grid data for alpha/HR
                     op,     #ode parameter object
                     alph,HR #new values of alpha and HR
                     ){
    ## errors if out of range
    if(HR > max(ao$HR)) cat('HR out of range! ',HR,' >', max(ao$HR),'\n')
    if(HR < min(ao$HR)) cat('HR out of range! ',HR,' <', min(ao$HR),'\n')
    if(alph > max(ao$alph)) cat('alph out of range! ',alph,' >', max(ao$alph),'\n')
    if(alph < min(ao$alph)) cat('alph out of range! ',alph,' <', min(ao$alph),'\n')
    ## approximations
    tst <- BLI(ao$irrhf,ao$alph,ao$HR,alph,HR) #HF
    op$irrhfD <- exp(tst)
    tst <- BLI(ao$irrhm,ao$alph,ao$HR,alph,HR) #HM
    op$irrhmD <- exp(tst)
    tst <- BLI(ao$irraf,ao$alph,ao$HR,alph,HR) #AF
    op$irrafD <- exp(tst)
    tst <- BLI(ao$irram,ao$alph,ao$HR,alph,HR) #AM
    op$irramD <- exp(tst)
    return(op)
}


## other parameter changes (simple wrapper)
## for changing in parameters
changepars <- function(bpar,L){
    for(nm in names(L))
        bpar[[nm]] <- L[[nm]]
    bpar
}


## use a parameter list and run model
runwithpars <- function(parmlist,baseparms){
  if(length(parmlist)>0){
    aph <- 0.34; hr <- 0.3;
    if('HR' %in% names(parmlist)){
      hr <- parmlist$HR
      parmlist$HR <- NULL
    }
    if('alph' %in% names(parmlist)){
      aph <- parmlist$alph
      parmlist$alph <- NULL
    }
    if('sigma' %in% names(parmlist)) parmlist$sigma <- NULL
    if('cdrdt' %in% names(parmlist)) parmlist$cdrdt <- NULL
    if('CDR' %in% names(parmlist)){
      ##needs handling differently diff models
      if(length(parmlist$CDR)==1){ #single CDR extended
        baseparms$CDRdata <- rep(parmlist$CDR,
                                 length(baseparms$ttq))
      } else {                            #full length CDR data
        baseparms$CDRdata <- parmlist$CDR #NOTE required to get right in ps2l
      }
      baseparms$CDRdata <- baseparms$CDRdata * 0.9 #NOTE cap
      parmlist$CDR <- NULL
    }
    baseparms <- amendAHR(AO,baseparms,aph,hr)
    baseparms <- changepars(baseparms,parmlist) #tweak parms
  }
  ltz <- seq(from=min(baseparms$ttp),
             to=max(baseparms$ttp),
             length.out = 1*500)
  mortmod$set_user(user=baseparms)
  system.time({out <- mortmod$run(ltz,method='rk4')})
  out
}



## term="IX|IH|IA", NoN for notifications
## getRateAge() <- in dataLL.R
## ## multirunner
nmvec2dt <- function(x) data.table(value=x,variable=names(x)) #utility

multirun <- function(DF,baseparms,bothout=FALSE,YR=2015){
  ansts <- ansss <- PR <- list()
  ARIS <- ARIA <- ARIF <- list()
  for(i in 1:nrow(DF)){
    if(!i%%10) print(i)
    tmp <- DF[i]
    out <- runwithpars(tmp,baseparms)
    ## timeseries
    outr <- out[,c('t','Itot','Ntot','NtotH',
                   'poptot','poph','popart',
                   'poph1549',"pop1549",
                   'Ipn','Ipp',
                   'prevtot15plus','pop15plus',
                   'mort','mortA','mortH')]
    outd <- as.data.table(outr)
    outd[,id:=i]
    ansts[[i]] <- outd
    ## -- snapshots
    ## by age
    ii <- which.min(abs(out[,'t']-YR))[1]
    allinc <- nmvec2dt(getRateAge(out,term="IX|IH|IA")[ii,-1])
    hinc <- nmvec2dt(getRateAge(out,term="IH|IA")[ii,-1])
    allnote <- nmvec2dt(getRateAge(out,term="NoN")[ii,-1])
    allinc[,c('qty','type'):=.('incidence','total')];
    hinc[,c('qty','type'):=.('incidence','HIV')];
    allnote[,c('qty','type'):=.('notifications','total')]
    outdat <- rbindlist(list(allinc,hinc,allnote))
    outdat[,id:=i]
    ansss[[i]] <- outdat
    ## --- ARI stuff
    if(bothout){
      OD <- out2df(out)
      ## scalar ari by time
      tmp <- OD[variable=='aris']
      tmp[,id:=i]
      ARIS[[i]] <- tmp
      ## ari by age
      tmp <- OD[year==2019][grep('ariv',variable),
                            .(year,value=mean(value)),by=age]
      tmp$age <- factor(tmp$age,levels=tmp$age,ordered = TRUE)
      tmp[,id:=i]
      ARIA[[i]] <- tmp
      ## ari from
      tmp <- OD[year==2019][grep('arif',variable),
                            .(year,value=mean(value)),by=age]
      tmp$age <- factor(tmp$age,levels=tmp$age,ordered = TRUE)
      tmp[,value:=value/sum(value)]
      tmp[,id:=i]
      ARIF[[i]] <- tmp
    }
    ## PR stuff
    if(makePR){
      tmp <- getPRdata(out)
      tmp[,id:=i]
      PR[[i]] <- tmp
    }
  }
  ansts <- rbindlist(ansts)
  ansss <- rbindlist(ansss)
  if(bothout){
    ARIS <- rbindlist(ARIS)
    ARIA <- rbindlist(ARIA)
    ARIF <- rbindlist(ARIF)
  }
  if(makePR)
    PR <- rbindlist(PR)
  ## tmp <- ansss[,.(value=sum(value)),by=.(id,age,hstate)]
  ## tmp1 <- tmp[,.(value=sum(value)),by=.(id,age)]
  ## tmp2 <- tmp[hstate!='X',.(value=sum(value)),by=.(id,age)]
  ## tmp1[,type:='total']; tmp2[,type:='HIV'];
  ## ansss <- rbind(tmp1,tmp2)
  setnames(ansss,'variable','age')
  ansss$age <- factor(ansss$age,levels=nagz,ordered=TRUE)
  ro <- ansts
  if(bothout==TRUE){
    ro <- list(ts=ansts,ss=ansss)
    ro$ARIS <- ARIS
    ro$ARIA <- ARIA
    ro$ARIF <- ARIF
    if(makePR) ro$pr <- PR
  }
  ro #return object
}

## Kish ESS calculator, from log weights
calcess <- function(lwt){
  wt <- lwt - max(lwt)
  wt <- exp(wt)
  wt <- wt/sum(wt)
  1/sum(wt^2)
}
