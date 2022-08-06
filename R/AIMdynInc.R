##' A function for running the AIM dynamics (R version)
##'
##' Content to be written
##' @title AIM dynamics
##' @param P 
##' @param graph 
##' @return X an array describing the end population state
##' @author Pete Dodd
##' @export

AimDynInc <- function(P,graph=FALSE,fullhist=FALSE,
                      alph=0.36,HR=0.3,gyr=Inf #which year to graph snapshots from: the end if == Inf
                      ){
  list2env(P,envir=environment())

  ## records
  popt <- poph <- pop1549 <- poph1549 <- popa <- hit <- hat <- yrz  #total pop, hiv
  W <- NULL
  if(fullhist) W <- array(0,dim = c(ntime,XD[['dz']]),
                          dimnames=list(paste0(1:ntime),
                                        XD[['nmsex']],XD[['nmage']],
                                        XD[['nmhiv']],XD[['nmart']]))

  ## initialize X
  X <- dI <- dA <- mu <- hivin <-
    mr <- ap <- array(0,dim=XD[['dz']],
                      dimnames=list(XD[['nmsex']],
                                    XD[['nmage']],
                                    XD[['nmhiv']],
                                    XD[['nmart']]))
  man <- hmredfac <- X
  K <- array(1,dim=XD[['dz']][1:3],## dimnames=list(nmsex,nmage,nmhiv)
             dimnames=list(XD[['nmsex']],
                           XD[['nmage']],
                           XD[['nmhiv']]))
  K[,,1] <- 0 #a mask to zero any HIV-ves

  ap[,,-1,2] <- 2; ap[,,-1,3] <- 4;        #ART aging rates (out)
  X[1,,1,1] <- UPD$bp[year==syear & sex==1 & age<81,value] #men
  X[2,,1,1] <- UPD$bp[year==syear & sex==2 & age<81,value] #women
  if( is.na(sum(X)) ) cat('popn pbm! 1 \n') #
  ## X[,,1,1] <- 1                           #some HIV-ve folk for testing..
  moUm <- moHm <- moAm <- moUf <- moHf <- moAf <- list() #records
  rim <- rif <- ram <- raf <- list() #records
  irrhm <- irrhf <- irram <- irraf <- list() #IRR records
  ## alph <- 0.36; HR <- 0.3                    #for IRRs

  if(fullhist){
    CA <- matrix(0,ncol=81,nrow=17)
    for(i in 1:80) CA[((i-1) %/% 5)+1,i] <- 1
    CA[nrow(CA),81] <- 1
    colnames(CA) <- XD[['nmage']]
    xz <- seq(from=0,to=80,by=5)
    rownames(CA) <- paste0('[',xz,',',c(xz[-1]-1,Inf),')')
    ## IRR matrices
    IRR <- array(1,dim=XD[['dz']],
                 dimnames=list(XD[['nmsex']],XD[['nmage']],
                               XD[['nmhiv']],XD[['nmart']]))
    IRR[,,">=500",] <- exp((1e3-750)*alph/100)
    IRR[,,"350-499",] <- exp((1e3-(350+500)/2)*alph/100)
    IRR[,,"250-349",] <- exp((1e3-(250+350)/2)*alph/100)
    IRR[,,"200-249",] <- exp((1e3-(200+250)/2)*alph/100)
    IRR[,,"100-199",] <- exp((1e3-(100+200)/2)*alph/100)
    IRR[,,"50-99",] <- exp((1e3-(50+100)/2)*alph/100)
    IRR[,,"<50",] <- exp((1e3-(50)/2)*alph/100)
    IRR[,,,"[0,6)m"] <- IRR[,,,"[0,6)m"]*(1-(1-HR)*.25)
    IRR[,,,"[7,12)m"] <- IRR[,,,"[7,12)m"]*(1-(1-HR)*.75)
    IRR[,,,"[12,Inf)m"] <- IRR[,,,"[12,Inf)m"]*HR
  }

  j <- 0
  ## ## time loop
  yrzs <- unique(floor(yrz))
  man <- ma #mortality rate on ART now
  for(i in 1:ntime){

    ## cat('timestep ',i,'\n')
    ## introduce a lookup for new year data so not done by 10
    if(!(i-1)%%(1/tstep)){           #NB assumes 1/tstep is integer
      k <- 1
      j <- j+1
      yr <- floor(yrz[i])

      ## ## new external inputs from getcountryparms
      ## HIV & ART targets
      htarget <- hivtarget[j]
      htarget2 <- hivtarget2[j]
      atarget <- arttarget[j]
      atarget2 <- arttarget2[j]
      ## births
      nbm <- nbmv[j]
      nbf <- nbfv[j]
      ## HIV sex disaggregation
      HSR <- srv[j]
      ## migration rates
      mr['M',,1,1] <- migrates[j,1,]
      mr['F',,1,1] <- migrates[j,2,]
      ## mortality rates (check fortran syntax when porting)
      mu['M',,,] <- mortrates[j,1,]
      mu['F',,,] <- mortrates[j,2,]

      ## next 4 state dependent  = not to shift
      noma <- apply(X=X[1,,,],MARGIN=1,FUN = sum) #number of men at each age
      nofa <- apply(X=X[2,,,],MARGIN=1,FUN = sum) #number of women at each agex

      ## recalcaulate
      fracm <- (AR[,1]*noma/(sum(AR[,1]*noma)+1e-10))
      fracm <- (1-HSR) * fracm / sum(fracm)
      fracf <- (AR[,2]*nofa/(sum(AR[,2]*nofa)+1e-10))
      fracf <- HSR * fracf / sum(fracf)
      hivin['M',,,1] <- fracm * cdinit['M',,,1] #into
      hivin['M',,1,1] <- - fracm                #outof
      hivin['F',,,1] <- fracf * cdinit['F',,,1] #into
      hivin['F',,1,1] <- - fracf                #outof

      ## update ART mortality according to year [RG]
      if(yr %in% y1){
        man <- ma #mortality rate on ART now
        whoy <- which(y1==yr)
        man[,,,2:3] <- man[,,,2:3] * rr1[whoy] #ART < 12 mo
        man[,,,4] <- man[,,,4] * rr2[whoy] #ART >= 12 mo
      }

      if(any(X<0)) X[X<0] <- 0 #safety
    }

    ## recording
    if(fullhist){
      tmp1 <- (X-X*mu) - mr*tstep #background mortality/migration:NOTE +/-
      tmp1 <- tmp1/tstep          #rate
    }

    ## ------- background mortality ---------
    X <- X * mu                     #todo: adjust to include migration
    ## ------- migration ---------
    X <- X + mr * tstep             #todo: figure out if per capita or what
    ## -------- aging ---------
    dX <- X * tstep
    X <- X * (1-tstep)
    X[,2:dim(X)[2],,] <- X[,2:dim(X)[2],,] + dX[,1:(dim(X)[2]-1),,]
    ## ------ new births --------
    X['M',1,1,1] <- X['M',1,1,1] + tstep * nbm #new boys
    X['F',1,1,1] <- X['F',1,1,1] + tstep * nbf  #new girls
    ## ------- new infections ------
    ## interpolate target & relax towards - works well
    htgt <- htarget2*k/npy + htarget*(npy-k)/npy #interpolate
    ## newI <- max(htgt*sum(X[,XD[['i1549']],,]) - sum(X[,XD[['i1549']],-1,]),0) * (1-exp(-20*tstep))
    speed <- 20
    newI <- max(htgt*sum(X[,XD[['io15']],,]) - sum(X[,XD[['io15']],-1,]),0) * (1-exp(-speed*tstep)) #1549 vs wholejj
    dI <- newI * hivin           #disaggregate
    X <- X + dI                             #population shifts

    ## --------- starting ART.--------
    ## [HIV+ not on ART (K nixes H-ves)] / [all HIV+]
    hmredfac[,,,1] <- K*X[,,,1] /
      (K*X[,,,1]+X[,,,2]+X[,,,3]+X[,,,4] + 1e-10) #check
    ## data to target
    atgt <- atarget2*k/npy + atarget*(npy-k)/npy #interpolate
    newA <- 0
    if(atgt>0 & TRUE){
      newA <- max(atgt*sum(X[,XD[['io15']],-1,]) -
                  sum(X[,XD[['io15']],,2:4]),
                  0) * (1-exp(-20*tstep))
      ## start pattern
      P1 <- X
      P1[,,1:2,] <- 0;            #not eligible if HIV- or >500
      P1[,,,2:4] <- 0;            #not eligible if on ART
      P1 <- P1 / (sum(P1)+1e-10)
      P2 <- X * mh * hmredfac   #proportional to HIV mortality
      P2 <- P2 / (sum(P2)+1e-10)
      artin <- -.5 * (P1*(1-alc) + P2*alc) # [RG6]
      artin[,,,2] <- -artin[,,,1] #going into 2=ART
      ## NB updating this here avoids having negative folk in low CD4
      dA <- newA * artin           #disaggregate
      X <- X + dA                             #population shifts
    }
    ## interpolator within year
    k <- k+1                     #iterate interpolator

    ## ------- progression ----------
    ## HIV
    dX <- X * tstep * pp
    X <- X - dX                       #out from the old
    X[,,1+(2:7),] <- X[,,1+(2:7),] + dX[,,2:7,] #in with the new
    ## ART
    dX <- X * tstep * ap
    X <- X - dX
    X[,,,3:4] <- X[,,,3:4] +  dX[,,,2:3]
    ## -------- HIV mortality ----------
    X <- X * (1 - tstep * mh * hmredfac) #[RG1]
    ## -------- ART mortality ----------
    X <- X * (1 - tstep * man) #[RG2]
    ## ---------record for checks -----------
    popt[i] <- sum(X)                 #population total
    poph[i] <- sum(X[,,-1,])          #hiv total
    ## poph1549[i] <- sum(X[,XD[['i1549']],-1,]) #hiv in 15-49
    ## pop1549[i] <- sum(X[,XD[['i1549']],,])  #15-49 population
    ## NOTE despite names, now o15
    poph1549[i] <- sum(X[,XD[['io15']],-1,]) #hiv in >=15
    pop1549[i] <- sum(X[,XD[['io15']],,])  #>=15
    popa[i] <- sum(X[,,,2:4])       #ART

    ## recording
    if(fullhist){
      RI <- IRR * X                   #relative incidences, to be averaged by H/A/age
      ## background mortality
      tmpm <- CA %*% apply(X=X[1,,1,],MARGIN=1,FUN = sum)
      tmpf <- CA %*% apply(X=X[2,,1,],MARGIN=1,FUN = sum)
      moUm[[i]] <- data.table(t((CA %*% apply(X=tmp1[1,,1,],
                                              MARGIN=1,FUN = sum))/(tmpm + 1e-15)))
      moUf[[i]] <- data.table(t((CA %*% apply(X=tmp1[2,,1,],
                                              MARGIN=1,FUN = sum))/(tmpf + 1e-15)))
      ## HIV
      tmp2 <- X * mh + X * ma
      tmpm <- CA %*% apply(X=X[1,,2:8,1],MARGIN=1,FUN = sum)
      tmpf <- CA %*% apply(X=X[2,,2:8,1],MARGIN=1,FUN = sum)
      moHm[[i]] <- data.table(t((CA %*% apply(X=tmp2[1,,2:8,1],
                                              MARGIN=1,FUN = sum))/(tmpm + 1e-15)))
      moHf[[i]] <- data.table(t((CA %*% apply(X=tmp2[2,,2:8,1],
                                              MARGIN=1,FUN = sum))/(tmpf + 1e-15)))
      irrhm[[i]] <- data.table(t((CA %*% apply(X=RI[1,,2:8,1],
                                               MARGIN=1,FUN = sum))/(tmpm + 1e-15)))
      irrhf[[i]] <- data.table(t((CA %*% apply(X=RI[2,,2:8,1],
                                               MARGIN=1,FUN = sum))/(tmpf + 1e-15)))
      ## ART
      tmpm <- CA %*% apply(X=X[1,,2:8,2:4],MARGIN=1,FUN = sum)
      tmpf <- CA %*% apply(X=X[2,,2:8,2:4],MARGIN=1,FUN = sum)
      moAm[[i]] <- data.table(t((CA %*% apply(X=tmp2[1,,2:8,2:4],
                                              MARGIN=1,FUN = sum))/(tmpm + 1e-15)))
      moAf[[i]] <- data.table(t((CA %*% apply(X=tmp2[2,,2:8,2:4],
                                              MARGIN=1,FUN = sum))/(tmpf + 1e-15)))
      irram[[i]] <- data.table(t((CA %*% apply(X=RI[1,,2:8,2:4],
                                               MARGIN=1,FUN = sum))/(tmpm + 1e-15)))
      irraf[[i]] <- data.table(t((CA %*% apply(X=RI[2,,2:8,2:4],
                                               MARGIN=1,FUN = sum))/(tmpf + 1e-15)))
      ## HIV/ART incidence NOTE not per capita
      rim[[i]] <- data.table(t((CA %*% apply(X=dI[1,,2:8,],MARGIN=1,FUN = sum))))
      rif[[i]] <- data.table(t((CA %*% apply(X=dI[2,,2:8,],MARGIN=1,FUN = sum))))
      ram[[i]] <- data.table(t((CA %*% apply(X=dA[1,,2:8,2:4],MARGIN=1,FUN = sum))))
      raf[[i]] <- data.table(t((CA %*% apply(X=dA[2,,2:8,2:4],MARGIN=1,FUN = sum))))
      hit[i] <- newI; hat[i] <- newA;
    }

    if(fullhist) W[i,,,,] <- X
    ## error:
    if( is.na(sum(X)) ) {print(yr);print(1e2*mean(is.na(X))); stop(paste0('popn pbm! t=',i,'\n'))}
    if(yr<gyr) XG <- X #capture from a particular year

  }

  if(fullhist){
    ## reformat
    moUm <- rbindlist(moUm); moHm <- rbindlist(moHm); moAm <- rbindlist(moAm)
    moUf <- rbindlist(moUf); moHf <- rbindlist(moHf); moAf <- rbindlist(moAf)
    rim <- rbindlist(rim); rif <- rbindlist(rif); ram <- rbindlist(ram); raf <- rbindlist(raf);
    irraf <- rbindlist(irraf); irram <- rbindlist(irram);
    irrhf <- rbindlist(irrhf); irrhm <- rbindlist(irrhm);
    moUm[,time:=syear+(1:nrow(moUm))*tstep-tstep]; moHm[,time:=syear+(1:nrow(moUm))*tstep-tstep];
    moAm[,time:=syear+(1:nrow(moUm))*tstep-tstep]; moUf[,time:=syear+(1:nrow(moUm))*tstep-tstep];
    moHf[,time:=syear+(1:nrow(moUm))*tstep-tstep]; moAf[,time:=syear+(1:nrow(moUm))*tstep-tstep];
    rim[,time:=syear+(1:nrow(moUm))*tstep-tstep]; rif[,time:=syear+(1:nrow(moUm))*tstep-tstep];
    ram[,time:=syear+(1:nrow(moUm))*tstep-tstep]; raf[,time:=syear+(1:nrow(moUm))*tstep-tstep];
    irraf[,time:=syear+(1:nrow(moUm))*tstep-tstep]; irram[,time:=syear+(1:nrow(moUm))*tstep-tstep];
    irrhf[,time:=syear+(1:nrow(moUm))*tstep-tstep]; irrhm[,time:=syear+(1:nrow(moUm))*tstep-tstep]
  }

  ## graphing
  if(graph!=FALSE){
    dfn <- graph + 'tmp.Rdata'
    save(P,XG,yrz,popt,poph,poph1549,pop1549,popa,graph,file=dfn)
    AIMplotterInc(P,XG,yrz,popt,poph,poph1549,pop1549,popa,graph,gyr)
  }

  if(fullhist)
    return(list(W=W,
                moUm=moUm,moHm=moHm,moAm=moAm,
                moUf=moUf,moHf=moHf,moAf=moAf,
                rim=rim,rif=rif,ram=ram,raf=raf,
                irraf=irraf, irram=irram,
                irrhf=irrhf, irrhm=irrhm,
                hit=hit,hat=hat,yrz=yrz))
  else
    return(X)
}

AIMplotterInc <- function(P,X,yrz,
                          popt,poph,poph1549,pop1549,popa,
                          graph,gyr){

  ## population data
  tmpP <- N80[iso3==P$iso3 & Year %in% yrz,
             .(value=1e3*sum(PopTotal)),by=.(year=Year)]
  modpop <- data.table(yrz,popt,poph,poph1549,pop1549,popa)
  tmp <- P$Hdata[year<max(modpop$yrz)]

  ## --------- pop growth
  ggplot(data=modpop,aes(x=yrz,y=popt)) +
    geom_line(color=2) +
    scale_y_continuous(label=absspace)+
    xlab('year')  + ylab('Total population') +
    geom_point(data=tmpP,aes(x=year,y=value))+ ggtitle(P$iso3)+
    theme_light()
  fn <- graph + P$iso3 + '_chk_totpop.pdf'
  ggsave(fn,width=7,height=7)
  cat('totpop saved...\n--- see:',fn,'\n')

  ## -----------------
  ## 1549 vs whole - now >15 despite names
  ggplot(data=modpop,aes(x=yrz,y=poph1549/pop1549)) +
    geom_line(color=2)+
    ## geom_point(data=P$data,aes(x=year,y=hivp*1e-2)) + #hp80pc
    geom_point(data=tmp,aes(x=year,y=hp80pc*1e-2)) + #
    scale_y_continuous(label=percent) +
    xlab('year') + ylab('HIV prevalence in 15-80 year olds') +
    ggtitle(P$iso3) + theme_light()
  ##+ geom_point(data=hatmp,aes(x=year,y=hp1549*1e-2),col=2)
  fn <- graph + P$iso3 + '_chk_hivp.pdf'
  ggsave(fn,width=7,height=7)
  cat('HIV check saved...\n--- see:',fn,'\n')
  ## -----------------

  aatmp <- P$Adata[year<max(modpop$yrz)]
  ggplot(data=modpop,aes(x=yrz,y=popa/(poph+1e-6))) +
    geom_line(color=2) + scale_y_continuous(label=percent) +
    xlim(c(2000,max(modpop$yrz)))+
    xlab('year') + ylab('ART coverage in HIV') +
    geom_point(data=aatmp,aes(x=year,y=art*1e-2))+
    ggtitle(P$iso3) + theme_light()

  fn <- graph + P$iso3 + '_chk_artp.pdf'
  ggsave(fn,width=7,height=7)
  cat('ART check saved...\n--- see:',fn,'\n')

  ## ## -----------------
  tmp <- reshape2::melt(X)
  names(tmp) <- c('sex','age','hiv','art','value')
  tmp <- as.data.table(tmp)
  tmpart <- tmp[,.(value=sum(value)),by=.(age,hiv,art)]

  tmp <- tmp[,.(value=sum(value)),by=.(sex,age,hiv)]
  tmp$value <- pmax(tmp$value,0)
  tmp[,HIV:=!(hiv=='hiv-ve')]
  tmpa <- tmp[,.(value=sum(value)),by=.(sex,age)]
  tmpah <- tmp[,.(value=sum(value)),by=.(sex,age,HIV)]
  tmpah <- tmpah[(order(HIV))]        #check pltfm indpt! todo!!
  tmpah$HIV <- factor(tmpah$HIV,levels=c('FALSE','TRUE'),ordered=TRUE)

  ## -----------------
  ggplot(data=tmpah,aes(x=age,fill=HIV,order=(HIV)))+
    geom_bar(data=tmpah[sex=='M'],aes(y=value),stat='identity')+
    geom_bar(data=tmpah[sex=='F'],aes(y=-value),stat='identity') +
    coord_flip()+geom_abline(intercept=0,slope=0) +
    ylab('Number') + xlab('Age')+
    theme_light()+
    theme(axis.text.y = element_text(size = rel(0.5)))+
    scale_y_continuous(label=absspace) +
    ggtitle(paste0(P$iso3,': ',gyr))
  fn <- graph + P$iso3 + '_chk_agehiv.pdf'
  ggsave(fn,width=7*1,height=7)
  cat('HIV/age saved...\n--- see:',fn,'\n--- see:',fn,'\n')

  ## -----------------
  tmpart$hiv <- factor(tmpart$hiv,
                       levels=rev(c(levels(tmpart$hiv)[-1],
                                    levels(tmpart$hiv)[1])),
                       ordered=TRUE)
  tmpart$hiv <- factor(tmpart$hiv,
                       levels=rev(levels(tmpart$hiv)),ordered=TRUE)

  ggplot(data=tmpart[hiv!='hiv-ve',],aes(x=age,fill=hiv))+
    geom_bar(aes(y=value),stat='identity')+
    geom_abline(intercept=0,slope=0) + ylab('Number') +
    xlab('Age') + facet_wrap(~art,scales='free_y') +
    scale_y_continuous(label=absspace)+
    theme_light()+
    theme(axis.text.x = element_text(size = rel(0.5),angle = 90)) +
    labs(fill='CD4 category')+ ggtitle(paste0(P$iso3,': ',gyr))
  fn <- graph + P$iso3 + '_chk_agecd4.pdf'
  ggsave(fn,width=7*2,height=7)
  cat('ART check saved...\n--- see:',fn,'\n--- see:',fn,'\n')

}

absspace <- function(x,...) {             #works
  format(abs(x), ..., big.mark=" ",scientific = FALSE, trim = TRUE)
}


## wrapper that either runs singly or not
AimDyns <- function(P,graph=FALSE,fullhist=FALSE,
                    alph=0.36,HR=0.3){
  if(length(alph) ==1 & length(HR)==1 ){
    ans <- AimDynInc(P=P,graph=graph,fullhist = fullhist,alph=alph,HR=HR)
    if(fullhist){
      ans$alph <- alph
      ans$HR <- HR
    }
  } else {
    ## ans <- AimDyn(P=P,graph=FALSE, fullhist = TRUE,alph=alph[1],HR=HR[1])
    ans <- AimDynInc(P=P,graph=FALSE, fullhist = TRUE,
                  alph=mean(alph),HR=mean(HR))
    tz <- seq(from=min(ans$irraf$time),
              to=max(ans$irraf$time),by=1) #1 a year
    agz <- names(ans$irraf)[1:17]   #age names
    irraf <- irram <- irrhf <- irrhm <- list()
    for(i in 1:length(alph)){
      irraf1 <- irram1 <- irrhf1 <- irrhm1 <- list()
      for(j in 1:length(HR)){
        print(c(i,j))
        tmp  <- AimDynInc(P=P,graph=FALSE, fullhist = TRUE,
                       alph=alph[i],HR=HR[j])
        irraf1[[j]] <- log(as.matrix(tmp$irraf[time %in% tz, ..agz]) + 1e-10)
        irram1[[j]] <- log(as.matrix(tmp$irram[time %in% tz, ..agz]) + 1e-10)
        irrhf1[[j]] <- log(as.matrix(tmp$irrhf[time %in% tz, ..agz]) + 1e-10)
        irrhm1[[j]] <- log(as.matrix(tmp$irrhm[time %in% tz, ..agz]) + 1e-10)
      }
      irraf[[i]] <- irraf1
      irram[[i]] <- irram1
      irrhf[[i]] <- irrhf1
      irrhm[[i]] <- irrhm1
    }
    ans$alph <- alph; ans$HR <- HR
    ans$irraf <- irraf
    ans$irram <- irram
    ans$irrhf <- irrhf
    ans$irrhm <- irrhf
  }
  return(ans)
}
