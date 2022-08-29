## === libraries
library(data.table)
library(ggplot2)
library(ggthemes)
library(scales)
## library(ggdist)
library(ggpubr)
## library(GGally)
library(ggrepel)
library(grid)
library(ggthemes)
library(patchwork)


## === load data
load(here('indata/fitdata/NT.Rdata'))
load(here('indata/fitdata/NAS.Rdata'))
load(here('indata/fitdata/E.Rdata'))
load(here('indata/fitdata/P.Rdata'))
load(here('indata/fitdata/VR.Rdata'))


## === utilities
## graph utilities
absci_10 <- function(x) {
  x <- abs(x)
  parse(text=gsub("e", "%*%10^", scales::scientific_format()(x)))
}

absspace <- function(x,...) {
  format(abs(x), ..., big.mark=" ",scientific = FALSE, trim = TRUE)
}
rotx <- rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))


## ============== useful data ============
agz <- c("[0,4)", "[5,9)", "[10,14)", "[15,19)", "[20,24)", "[25,29)",
         "[30,34)", "[35,39)", "[40,44)", "[45,49)","[50,54)","[55,59)",
         "[60,64)",  "[65,69)",  "[70,74)",  "[75,79)",  "[80,Inf)")

nagz <- c("0-4", "5-14", "15-24","25-34","35-44","45-54","55-64","65+")
a1549 <- paste0(seq(from=15,by=5,to=45),'-',seq(from=19,by=5,to=49))

## cnz <- rev(c('ZWE','ZMB','ZAF','NGA','MWI','ETH'))
cz <- c('ETH','KEN','LSO','MOZ','MWI','NGA','SWZ',
        'TZA','UGA','ZAF','ZMB','ZWE')

## make Key for AIM HIV/ART data
cnm <-  c('Ethiopia','Kenya','Lesotho','Mozambique',
          'Malawi','Nigeria','Swaziland',
          'United Republic of Tanzania','Uganda',
          'South Africa','Zambia','Zimbabwe')

cknz <- data.table(iso3=cz,name=cnm) #TODO OR USE PREV


## cky <- unique(SD19::EAA[,.(Code,Country)]) #TODO
## for(cn in cname) if( !cn %in% cky$Country) print(cn)
## for(i in 1:length(cname)) cky[Country==cname[i],iso3:=cz[i]]
## cky <- cky[!is.na(iso3)]
## HA <- SD19::EAA[Code %in% cky$Code,.(year,Code,
## HA <- merge(HA,cky[,.(iso3,Code)],by='Code',all.x=TRUE)

## summary(XD)
## summary(XD[['ETH']])




NF <- N80[,.(PopTotal=sum(PopTotal)),by=.(year=Year,iso3)]
NF2 <- N80[AgeGrp %in% a1549,.(PopTotal1549=sum(PopTotal)),by=.(year=Year,iso3)]
NF <- merge(NF,NF2,by=c('iso3','year'))
NF[,F1549:=PopTotal1549/PopTotal]
NF[,c('PopTotal1549','PopTotal'):=NULL]

## HA <- merge(HA,NF,by=c('iso3','year'))
## HA[,unique(iso3)]


## parameter namings
useunicode <- TRUE #flag for use
modpmnmz <- c("bet","ari0","arig",
              "pp","eps","alph",
              "HR","CDR","drnX",
              "drnH","v","txf",
              "cfrn","cfrpn","cfrpp",
              "rel","cdrdt",
              "OR04","OR514","pp04")
## order to use
der <- c(1,11,2, #tx
         4,20,5,3,16, #nh
         6:7,       #hiv
         8,17:19,   #cdr
         9:10,      #dur
         12:15      #cfr
         )
parmlabs <- c("\u03B2","\u039B\u2080","\u03BA",
              "\u03B5","\u03BD","\u03B1",
              "\u03C1","K","D\u02E3",
              "D\u1D34","\u03A8","\u03B8",
              "\u03A6\u02E3","\u03A6\u1D34","\u03A6\u1D2C",
              "\u03C9","c",
              "OR\u2080\u2084","OR\u2085\u2081\u2084",
              "\u03B5\u2080\u2084")
## print(parmlabs)
parmkey <- data.frame(model=modpmnmz,
                      writeup=parmlabs,
                      stringsAsFactors = FALSE)

parmkey <- as.data.table(parmkey)
save(parmkey,file=here('data/parmkey.Rdata'))



## ====================== PLOTTING FUNCTIONS ===============
HML <- function(x){
  list(mid=quantile(x,0.5),lo=quantile(x,0.025),hi=quantile(x,1-0.025))
}


## ============== inference outputs ==============
MakeInfPlots <- function(cniso){
  ## load data
  sminf <- glue(here('data/sampe_{cniso}.Rdata'))
  load(sminf)
  if(useunicode){
    setnames(sampe,old=modpmnmz,new=parmlabs)
    setcolorder(sampe,der)
  }

  ## pairs plot
  PP <- ppairs(sampe)
  ppfn <- glue(here('plots/corner/npi_PP_{cniso}.pdf'))

  ## save
  ggsave(PP,file=ppfn,w=10,h=10,device = cairo_pdf)
}



prepOD <- function(OD,cniso){
  if(!'year' %in% names(OD))
    if('t' %in% names(OD))
      OD[,year:=t]
  PT <- OD[,.(year,value=poptot)]
  PH <- OD[,.(year,value=poph)]
  PA <- OD[,.(year,value=popart)]
  ## data to plot against
  PTD <- N80[iso3==cniso & Year<max(OD$year) & Year>1970]
  PTD <- PTD[,.(value=1e3*sum(PopTotal)),by=.(year=Year)]


  ## merge against pop
  popsr <- OD[,.(year,poptot,pop1549,pop15plus,poph,popart,
                 Itot,Ipn,Ipp,mort,mortA,mortH)]

  pops <- popsr[,.(poptot=median(poptot),
                  pop1549=median(pop1549),
                  pop15plus=median(pop15plus),
                  poph=mean(poph)## ,
                  ## Itot=mean(Itot),
                  ## Ipn=mean(Ipn),
                  ## Ipp=mean(Ipp),
                  ## mort=mean(mort),
                  ## mortA=mean(mortA),
                  ## mortH=mean(mortH)
                  ),
               by=.(year=round(year))]

  ## return
  list(pops=pops,popsr=popsr,PT=PT,PH=PH,PA=PA,PTD=PTD)
}

test <- prepOD(OD,cniso)

prepHdata <- function(cniso){
  ## HIV/ART data
  HAdata <- data.table(
    iso3=cniso,
    year=XD[[cniso]]$y2,
    mlhiv=XD[[cniso]]$mh,
    flhiv=XD[[cniso]]$fh,
    artm=XD[[cniso]]$ma,
    artf=XD[[cniso]]$fa
  )
  HAdata[,plhiv:=(mlhiv+flhiv)]
  HAdata[,art:=(artm*mlhiv+artf*flhiv)/(plhiv)]
  HAdata <- merge(HAdata,
                  data.table(iso3=cniso,
                             year=XD[[cniso]]$y3,
                             hp80pc=XD[[cniso]]$hp80pc),
                  by=c('iso3','year'),all=TRUE)

  HAdata
}

## === pop & HIV time-series
## MR or TBP type input
make.poph.plt <- function(cniso,  #country
                          OD,     #model run
                          SY=1980, #start year
                          na.rm=FALSE
                          ){
  ## prepare model input data
  odprep <- prepOD(OD,cniso)
  list2env(odprep,envir = environment())

  ## merge HIV data against pop
  HAdata <- prepHdata(cniso)
  HD <- merge(HAdata,pops,by='year',all.y=FALSE,all.x = TRUE)
  HDH <- HD[hp80pc>0,.(year,value=hp80pc*pop15plus)][!is.na(value)]
  HDA <- HD[art>0,.(year,value=art*poph/1e2)][!is.na(value)]

  ## plot
  Dt <- ggplot(PT,aes(year,value))+
    geom_line()+
    geom_line(data=PH,col=2)+
    geom_line(data=PA,col=3)+
    geom_point(data=HDA,col='green',shape=1)+
    geom_point(data=HDH,col='red',shape=1)+
    scale_y_sqrt(label=absspace)+
    scale_x_continuous(limits=c(SY,NA))+
    geom_point(data=PTD,pch=1)+
    ylab('Population')+xlab('Year')+
    theme_classic()+ggpubr::grids()
  Dt
}

## ## ## -----------------

## ## make.poph.plt(cniso,OD)

## === pop snapshot
make.popsnap.plt <- function(cniso,OD){
  ## TODO choose where tickmarks are!
  yr <- 2015
  myrs <- OD[,unique(year)] #model years
  myr <- myrs[which.min(abs(myrs-yr))] #model year closest yr
  ## model data
  DA <- OD[tbstate=='no TB' & year==myr & !is.na(age),
           .(value=sum(value)),by=.(age,sex)]
  DA[,Age:=(gsub(",","-",gsub("\\[|\\)","",age)))]
  DA[grepl('Inf',Age),Age:="80+"]
  setkey(DA,age)
  nags <- DA[agz][sex=='F']$Age
  DA$Age <- factor(DA$Age,levels=nags,ordered = TRUE)
  ## real data
  DAD <- N80[iso3==cniso & Year==yr]
  DAD <- melt(DAD[,.(Age=AgeGrp,M=PopMale,F=PopFemale)],id='Age')
  names(DAD)[names(DAD)=='variable'] <- 'sex'
  mx <- DA[,max(value)] #make breaks
  x <- floor(log10(mx))
  y <- ceiling(mx /10^x)
  y <- seq(from=1,to=4,length.out = 2) #not too many
  bks <- y*10^x;bks <- c(-bks,bks)

  ## graph
  Ds <- ggplot(data=DA,aes(x=Age,fill=sex))+
    geom_bar(data=DA,aes(y=ifelse(sex=='M',value,-value)),
             stat='identity') +
    coord_flip()+geom_abline(intercept=0,slope=0) +
    ylab('Population') + xlab('Age') +
    theme(axis.text.y = element_text(size = rel(0.6))) +
    scale_y_continuous(label=absci_10,breaks=bks)+
    geom_point(data=DAD,
               aes(y=ifelse(sex=='M',1e3*value,-1e3*value)),
               shape=1)+
    theme_minimal()+
    theme(legend.position = "none",
          axis.text.x=element_text(size=rel(0.75)))
  Ds
}


## ## make.popsnap.plt(cniso,OD)


## === prevalence
make.prev.plt <- function(cniso,MR,na.rm=FALSE){
  tmp <- P[cniso][sex=='all' & acat=='all' & case.type=='a']
  tmp[,type:='prevalence']
  clz <- c('prevalence'='purple')
  if('year' %in% names(MR) & !'t' %in% names(MR)) MR[,t:=year]
  MRL <- copy(MR)
  MRL[,type:="prevalence"]
  MRLM <- MRL[,.(value=median(1e5*prevtot15plus/pop15plus,na.rm=na.rm),
                 lo=quantile(1e5*prevtot15plus/pop15plus,0.025,na.rm=na.rm),
                 hi=quantile(1e5*prevtot15plus/pop15plus,0.975,na.rm=na.rm)),
              by=.(t,type)]
  Pplt <- ggplot(MRLM,aes(t,y=value,ymin=lo,ymax=hi,
                         col=type,fill=type))
  if('id' %in% names(MR)){
    Pplt <- Pplt + geom_ribbon(alpha=0.2,col=NA)+
      geom_line()
  } else {
    Pplt <- Pplt + geom_line()
  }
  Pplt <- Pplt +
    scale_fill_manual(values=clz,guide="none")+
    scale_color_manual(values=clz,guide="none")+
    scale_x_continuous(limits=c(1980,NA))+
    geom_pointrange(data=tmp,## size=1,
                    aes(x=year,y=prev,ymin=prev.lo,ymax=prev.hi),
                    col="purple")+
    expand_limits(y=0) +
    theme_classic()+
    ylab('Adult TB prevalence per 100,000') + xlab('Year')+
    ggpubr::grids()
  Pplt
}

## ## make.prev.plt(cniso,MR$ts)


## ## -- aggregate:
## ## ## new ages
## ## MR[age %in% agz[14:length(agz)],acat:='65+']
## ## MR[age %in% agz[12:13],acat:='55-64']
## ## MR[age %in% agz[10:11],acat:='45-54']
## ## MR[age %in% agz[8:9],acat:='35-44']
## ## MR[age %in% agz[6:7],acat:='25-34']
## ## MR[age %in% agz[4:5],acat:='15-24']
## ## MR[age %in% agz[2:3],acat:='5-14']
## ## MR[age %in% agz[1],acat:='0-4']
## ## MR <- MR[,.(value=sum(value)),by=.(id,acat,type)]
## ## MR$acat <- factor(MR$acat,levels=nagz,ordered=TRUE)

## ## ========  Age patterns

make.Age.inc.plt <- function(cniso,MR,narm=FALSE){
  clz <- c('total'='blue','HIV'='red')
  MR[,acat:=age]
  MRM <- MR[,.(value=median(value/epfac,na.rm=narm),
               lo=quantile(value/epfac,0.025,na.rm=narm),
               hi=quantile(value/epfac,0.975),na.rm=narm),
               by=.(acat,qty,type)]
  ## notification data
  tmpn <- NAS[cniso]
  tmpn <- tmpn[,.(value=sum(value)),by=.(year,acat)]
  tmpn <- tmpn[year==2015]
  tmpn[,c('type'):=NA]
  ## -- make plot
  Plt <- ggplot(MRM[qty!='notifications'],
                aes(acat,value,group=type,
                    col=type,fill=type))+
    expand_limits(y=0)+
    geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.2,color=NA)+
    geom_line()+geom_point()+
    scale_fill_manual(values = clz,guide="none")+
    scale_color_manual(values = clz,guide="none")+
    ylab('TB incidence')+ xlab('Age category')+
    scale_y_continuous(label=absspace)+
    geom_point(data=tmpn,aes(acat,value),shape=1,col=1,size=2)+
    theme_classic()+ggpubr::grids() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))
  Plt
}


## ## make.Age.inc.plt(cniso,MR$ss)

## ## MR$ss$age <- factor(MR$ss$age,levels=agz,ordered=TRUE)

## ## make.Age.inc.plt(cniso,MR$ss)

## ## nagz <- c('0-4','5-14','15-24','25-34','35-44','45-54','55-64','65+')

## ## tmp <- MR$ss
## ## tmp[,ag:='0-4']
## ## tmp[age %in% agz[2:3],ag:="5-14"]
## ## tmp[age %in% agz[4:5],ag:="15-24"]
## ## tmp[age %in% agz[6:7],ag:="25-34"]
## ## tmp[age %in% agz[8:9],ag:="35-44"]
## ## tmp[age %in% agz[10:11],ag:="45-54"]
## ## tmp[age %in% agz[12:13],ag:="55-64"]
## ## tmp[age %in% agz[14:length(agz)],ag:="65+"]
## ## ## tmp[,table(ag,age)]
## ## tmp <- tmp[,.(value=sum(value)),by=.(age=ag,type,id)]
## ## tmp$age <- factor(tmp$age,levels=nagz,ordered=TRUE)

## ## make.Age.inc.plt(cniso,tmp)

## ## === notifications & incidence
## ## WHO comparison TODO
## make.WHOcompare.data <- function(cniso,MR,na.rm=FALSE){
##   ## if('year' %in% names(MR) & !'t' %in% names(MR)) MR[,t:=year]
##   MRM <- melt(MR[,.(t,id,
##                     Itot=Itot/epfac,
##                     IHtot=(Ipn+Ipp)/epfac,
##                     Mtot=mort/epfac,
##                     MHtot=(mortH+mortA)/epfac)],id=c('id','t'))
##   yrs <- MR[,min(floor(t))]:MR[,max(floor(t))]
##   wyrs <- data.table(yrs=yrs) #nearest to mid pts
##   for(y in yrs) wyrs[yrs==y,
##                      tid:=which.min(abs(MR[id==1,t]-(y+0.5)))]
##   wyrs[,t:=MR[id==1][tid,t]]
##   MRM <- MRM[t %in% wyrs$t]
##   MRM <- merge(MRM,wyrs[,.(year=yrs,t)],by='t')
##   MRM <- MRM[,.(mid=median(value,na.rm=na.rm),
##                 hi=quantile(value,.975,na.rm=na.rm),
##                 lo=quantile(value,.025,na.rm=na.rm)),
##              by=.(year,variable)]
##   etmp <- melt(E[cniso][,.(year,
##                            Itot=e_inc_num,
##                            Itot.hi=e_inc_num_hi,
##                            Itot.lo=e_inc_num_lo,
##                            IHtot=e_inc_tbhiv_num,
##                            IHtot.hi=e_inc_tbhiv_num_hi,
##                            IHtot.lo=e_inc_tbhiv_num_lo,
##                            Mtot=e_mort_num,
##                            Mtot.lo=e_mort_num_lo,
##                            Mtot.hi=e_mort_num_hi,
##                            MHtot=e_mort_tbhiv_num,
##                            MHtot.lo=e_mort_tbhiv_num_lo,
##                            MHtot.hi=e_mort_tbhiv_num_hi
##                            )],
##                id='year')
##   etmp[,c('variable','type'):=tstrsplit(variable,"\\.")]
##   etmp[is.na(type),type:='mid']
##   etmp <- dcast(etmp,year+variable~type,value.var='value')
##   ## merge
##   etmp[,estimate:='WHO']; MRM[,estimate:='transmission model']
##   commonyrs <- intersect(etmp$year,MRM$year)
##   B <- rbind(etmp[year %in% commonyrs],
##              MRM[year %in% commonyrs],
##              use.names = TRUE)
##   B[,iso3:=cniso]
##   B
## }

## ## tt <- make.WHOcompare.data(cniso,MR$ts)

## time-series
make.NI.plt <- function(cniso,MR,WHO=FALSE,narm=FALSE){
  if('year' %in% names(MR) & !'t' %in% names(MR)) MR[,t:=year]
  MRM <- melt(MR[,.(t,
                    Itot=Itot/epfac,
                    Ntot=Ntot/epfac,
                    NtotH=NtotH/epfac)],id='t')
  MRMM <- MRM[,.(value=median(value,na.rm=narm),
                 hi=quantile(value,.975,na.rm=narm),
                 lo=quantile(value,.025,na.rm=narm)),
              by=.(t,variable)]
  ## print(MRMM)
  clz <- c('Itot'='blue','Ntot'='black',
           'NtotH'='red','NtotA'='green')
  tmp <- merge(NT[cniso],E[cniso],by='year',all=TRUE)
  etmp <- melt(tmp[,.(t=year,Itot=e_inc_num,NtotH=e_inc_tbhiv_num)],
               id='t')
  eih <- melt(tmp[,.(t=year,Itot=e_inc_num_hi)],id='t')
  eil <- melt(tmp[,.(t=year,Itot=e_inc_num_lo)],id='t')
  dtmp <- melt(tmp[,.(t=year,Ntot=c_newinc,
                      NtotH=newrel_hivpos,NtotA=newrel_art)],
               id='t')

  ## plot
  NIplt <- ggplot(MRMM,aes(t,y=value,
                           col=variable,fill=variable))
  if('id' %in% names(MR)){
    NIplt <- NIplt +
      geom_ribbon(aes(ymin=lo,ymax=hi),alpha=.2,col=NA)+geom_line()
  } else {
    NIplt <- NIplt + geom_line()
  }
  NIplt <- NIplt +
    scale_fill_manual(values = clz,guide="none")+
    scale_color_manual(values = clz,guide="none")+
    scale_y_continuous(label=absspace,limits=c(0,NA))+
    scale_x_continuous(limits=c(1980,NA))+
    ylab('TB cases')+xlab('Year')+
    theme_classic() + ggpubr::grids()+
    geom_point(data=dtmp[variable=='NtotA'],
               alpha=1,shape=1,size=2,col='green')+
    geom_point(data=dtmp[variable=='NtotH'],
               alpha=1,shape=1,size=2,col='red')+
    geom_point(data=dtmp[variable=='Ntot'],
               alpha=1,shape=1,size=2,col='black')
  if(WHO)
    NIplt <- geom_line(data=etmp[variable!='NtotH'],alpha=1,lty=2)+
      geom_line(data=eih,alpha=1,lty=3)+
      geom_line(data=eil,alpha=1,lty=3)
  NIplt
}

## ## make.NI.plt(cniso,MR$ts)

## ## make.NI.plt(cniso,MR$ts)

## ## MRreal <- copy(MR)
## ## MR <- MR$ts
## ## MR <- copy(MRreal)

## time-series
make.NIPC.plt <- function(cniso,MR,ttl='',
                          ylb='TB incidence per 100,000 per year',xlb='Year',
                          narm=FALSE){
  if('year' %in% names(MR) & !'t' %in% names(MR)){
    MR[,t:=year]
    MR[,year:=NULL]
  }
  MRM <- MR[,.(t,id,value=1e5*(Itot/epfac)/poptot)]
  MRMM <- MRM[,.(value=median(value,na.rm=narm),
                 hi=quantile(value,0.975,na.rm=narm),
                 lo=quantile(value,0.025,na.rm=narm)),
              by=.(t)] #,variable

  ## plot
  NIplt <- ggplot(MRMM,aes(t,y=value))
  if('id' %in% names(MR)){
    NIplt <- NIplt +
      geom_ribbon(aes(ymin=lo,ymax=hi),fill='blue',alpha=.2,col=NA)+
      geom_line(col='blue')
  } else {
    NIplt <- NIplt + geom_line(col='blue')
  }
  NIplt <- NIplt +
      scale_y_continuous(label=absspace,limits=c(0,NA))+
    scale_x_continuous(limits=c(1980,NA))+
    ylab(ylb)+xlab(xlb)+
    theme_classic() + ggpubr::grids()+ggtitle(ttl)
    NIplt
}

## load(here('data/MR_ETH.Rdata'))
## make.NIPC.plt(cniso,MR$ts)

## ## 1

## make.CDR.plt <- function(cniso,MR,sampe,narm=FALSE){
##   MRM <- melt(MR[,.(t,Itot=Itot/epfac,Ntot=Ntot/epfac,NtotH=NtotH/epfac)],id='t')
##   tmp <- merge(NT[cniso],E[cniso],by='year',all=TRUE)
##   etmp <- melt(tmp[,.(t=year,Itot=e_inc_num,NtotH=e_inc_tbhiv_num)],
##                id='t')
##   eih <- melt(tmp[,.(t=year,Itot=e_inc_num_hi)],id='t')
##   eil <- melt(tmp[,.(t=year,Itot=e_inc_num_lo)],id='t')
##   dtmp <- melt(tmp[,.(t=year,Ntot=c_newinc,
##                       NtotH=newrel_hivpos,NtotA=newrel_art)],
##                id='t')

##   ## also looking at CDR
##   cdrdat <- merge(MRM[variable=='Itot',
##                       .(inc=median(value,na.rm = narm),
##                         inc.lo=quantile(value,0.025,na.rm = narm),
##                         inc.hi=quantile(value,0.975,na.rm = narm)),
##                       by=.(t=round(t))],
##                   dtmp[variable=='Ntot',.(t,notes=value)],
##                   by='t',all=TRUE)
##   cdrdat[,c('CDR','CDR.lo','CDR.hi'):=.(notes/inc,
##                                         notes/inc.hi,
##                                         notes/inc.lo)]
##   ## who est
##   west <- merge(etmp[variable=='Itot',.(t,winc=value)],
##                 eil[variable=='Itot',.(t,wlo=value)],
##                 by='t')
##   west <- merge(west,
##                 eih[variable=='Itot',.(t,whi=value)],
##                 by='t')
##   west <- merge(west,
##                 dtmp[variable=='Ntot',.(t,notes=value)],
##                 by='t')
##   west[,c('wCDR','wCDR.lo','wCDR.hi'):=.(notes/winc,
##                                          notes/whi,
##                                          notes/wlo)]
##   ## CDR model - NOTE change wrt dataLL
##   ns <- 300
##   if('cdrdt' %in% names(sampe)){
##     if('ecdrdt' %in% names(sampe)){
##       cdrmod <- sampe[sample(nrow(sampe),ns),
##                       .(CDR,cdrdt,ecdrdt)]
##     } else {
##       cdrmod <- sampe[sample(nrow(sampe),ns),
##                       .(CDR,cdrdt)]
##       cdrmod[,ecdrdt:=0]
##     }
##   } else {
##     cdrmod <- sampe[sample(nrow(sampe),ns),
##                     .(CDR)]
##     cdrmod[,c('cdrdt','ecdrdt'):=0]
##   }
##   cdrmod[,id:=1:ns]
##   tz <- seq(from=1970,to=2020,by=1)
##   extra <- expand.grid(t=tz,id=1:ns)
##   cdrmod <- merge(cdrmod,extra,by='id',all=TRUE)
##   cdrmod[,dt:=t-min(t)]
##   cdrmod[,estcdr:=invlogit(logit(CDR) + dt*cdrdt*exp(-dt*ecdrdt))]

##   ## make plot
##   cdrplt <- ggplot(cdrdat,aes(t,CDR))+
##     stat_lineribbon(data=cdrmod,aes(t,estcdr),
##                     .width =.95,alpha=.2) +
##     geom_line()+geom_point()+ xlab('Year')+
##     geom_line(aes(t,CDR.lo),lty=2)+
##     geom_line(aes(t,CDR.hi),lty=2)+
##     geom_line(data=west,aes(t,wCDR.lo),lty=2,col=2)+
##     geom_line(data=west,aes(t,wCDR.hi),lty=2,col=2)+
##     geom_line(data=west,aes(t,wCDR),lty=1,col=2)+
##     scale_y_continuous(label=percent,limits=c(0,NA))+
##     theme(legend.position = 'none')
##   cdrplt

## }

## ## make.CDR.plt(cniso,MR,sampe)

## ## === mortality
## make.mort.plt <- function(cniso,MR,whoshow=FALSE){
##   ## data prep
##   MRM <- melt(MR[,.(t,mort=mort/epfac,morth=(mortH+mortA)/epfac)],id='t')
##   clz <- c('mort'='darkturquoise','morth'='red')
##   tmp <- E[cniso][,.(t=year,e_mort_num,e_mort_tbhiv_num,
##                      e_mort_num_lo,e_mort_num_hi,
##                      e_mort_tbhiv_num_lo,e_mort_tbhiv_num_hi)]
##   etmp <- melt(tmp,id='t')
##   names(etmp)[2] <- 'nm'
##   etmp[,variable:=ifelse(grepl('hiv',nm),'morth','mort')]
##   etmp[,tp:=ifelse(grepl('_hi|lo',nm),'hi','mid')]
##   etmp[grepl('lo',nm),tp:='lo']

##   ## plot
##   Mplt <- ggplot(MRM,aes(t,value,col=variable,fill=variable)) +
##     stat_lineribbon(.width =.95,alpha=.2) +
##     scale_fill_manual(values = clz,guide="none")+
##     scale_color_manual(values = clz,guide="none")+
##     scale_y_continuous(label=absspace,limits=c(0,NA))+
##     scale_x_continuous(limits=c(1980,NA))+
##     ylab('TB deaths')+xlab('Year')+
##     theme_classic() + ggpubr::grids()
##   if(whoshow)
##     Mplt <- Mplt+
##       geom_line(data=etmp[tp=='mid'],lty=2)+
##       geom_line(data=etmp[tp=='lo'],lty=3)+
##       geom_line(data=etmp[tp=='hi'],lty=3)
##   Mplt

## }

## ## make.mort.plt(cniso,MR)

