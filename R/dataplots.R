## === libraries
library(data.table)
library(ggplot2)
library(ggthemes)
library(scales)
library(ggpubr)
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

cz <- c('ETH','KEN','LSO','MOZ','MWI','NGA','SWZ',
        'TZA','UGA','ZAF','ZMB','ZWE')

## make Key for AIM HIV/ART data
cnm <-  c('Ethiopia','Kenya','Lesotho','Mozambique',
          'Malawi','Nigeria','Swaziland',
          'United Republic of Tanzania','Uganda',
          'South Africa','Zambia','Zimbabwe')

cknz <- data.table(iso3=cz,name=cnm) #TODO OR USE PREV


NF <- N80[,.(PopTotal=sum(PopTotal)),by=.(year=Year,iso3)]
NF2 <- N80[AgeGrp %in% a1549,.(PopTotal1549=sum(PopTotal)),by=.(year=Year,iso3)]
NF <- merge(NF,NF2,by=c('iso3','year'))
NF[,F1549:=PopTotal1549/PopTotal]
NF[,c('PopTotal1549','PopTotal'):=NULL]


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



prepOD <- function(OD,cniso){
  if(!'year' %in% names(OD))
    if('t' %in% names(OD))
      OD[,year:=t]

  ## totals
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
                  poph=mean(poph)
                  ),
               by=.(year=round(year))]

  ## return
  list(pops=pops,popsr=popsr,PT=PT,PH=PH,PA=PA,PTD=PTD)
}


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

## for legend
lls <- c('Total population','People living with HIV','People with HIV on ART')
AtmpD <- data.table(x=rep(c(0,1),3),y=rep(1:3,each=2),`A)`=rep(lls,each=2))
AtmpD$`A)` <- factor(AtmpD$`A)`,levels = lls,ordered = TRUE)
GL <- ggplot(AtmpD,aes(x,y,col=`A)`))+geom_point(shape=1)+geom_line()+
  scale_color_manual(values=c('black','red','green')) +
  theme(legend.title = element_text(face = "bold"))
lega <- as_ggplot(get_legend(GL))

lls <- c('Data','Model')
Atmp <- data.table(x=rep(c(0,1),2),y=rep(1:2,each=2),v1='Model',v2='Data')
GL <- ggplot(Atmp,aes(x,y))+geom_point(aes(col=v2),shape=1)+
  scale_color_manual(values=c('Data'='black')) +
  ggnewscale::new_scale_color() +
  geom_line(aes(col=v1))+
  scale_color_manual(values=c('Model'='black')) +
  theme(legend.title = element_blank())
lega2 <- as_ggplot(get_legend(GL))




## ## ## -----------------

## prepOD(MR$ts[id==1],cniso)$popsr

## === pop snapshot
make.popsnap.plt <- function(cniso,OD){
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

## make.poph.plt(cniso,MR$ts[id==1])
## ## make.popsnap.plt(cniso,OD)
## for legend
lls <- c('Male population','Female population')
BtmpD <- data.table(x=c(0,1),y=rep(1,each=2),`model:`=rep(lls,each=2),type=NA)
BA <- ggplot(BtmpD,aes(x,y,fill=`model:`))+geom_bar(stat='identity')+
  scale_fill_discrete(name='model:\n\n')

legb2 <- as_ggplot(get_legend(BA))
B2 <- data.table(x=c(0,1),y=rep(1,each=2),type='UN data')
BB <- ggplot(B2,aes(x,y)) +
  geom_point(aes(col=type),shape=1) +
  scale_color_manual(name='B)\n\n\n',values=c('black')) +
  theme(legend.title = element_text(face = "bold"))
legb1 <- as_ggplot(get_legend(BB))

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

## legend stuff
tmpd <- data.table(x=c(0,1),y=c(0,1))
tmpd[,c('ymin','ymax','f','v'):=.(c(-0.1,0.9),c(0.1,1.1),'Prevalence data',
                                  'TB prevalence')]
tmp <- ggplot(data=tmpd,aes(x=x,y=y,ymin=ymin,ymax=ymax,col=f))+
  geom_pointrange()+
  scale_color_manual(name='C)\n\n\n',values = 'purple') +
  theme(legend.title = element_text(face = "bold"))
legc1 <- as_ggplot(get_legend(tmp))
tmp <- ggplot(data=tmpd,aes(x=x,y=y,ymin=ymin,ymax=ymax,fill=v,col=v))+
  geom_line()+ geom_ribbon(alpha=0.2,col=NA)+
  scale_fill_manual(name='model:\n\n\n',values = 'purple')+
  scale_color_manual(name='model:\n\n\n',values = 'purple')
legc2 <- as_ggplot(get_legend(tmp))



## ## make.prev.plt(cniso,MR$ts)

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


## legend stuff
tmpd <- data.table(x=c(0,1,2),y=c(0,1,2))
tmpd[,v:=c('TB notifications','TB/HIV notifications','TB/HIV/ART notifications' )]
tmp <- ggplot(data=tmpd,aes(x=x,y=y,col=v))+
  geom_point(shape=1)+
  scale_color_manual(name='D) & E)',values = c('black','red','green')) +
  theme(legend.title = element_text(face = "bold"))
legd1 <- as_ggplot(get_legend(tmp))


tmpd <- data.table(x=rep(0:2,3),y=rep(0:2,3))
tmpd[,v:=rep(c('TB notifications','TB/HIV notifications','TB incidence'),each=3)]
tmpd[,c('ymin','ymax'):=.(0.9*y,1.1*y)]

tmp <- ggplot(data=tmpd,aes(x=x,y=y,
                            ymin=ymin,ymax=ymax,
                            fill=v,col=v,group=v))+
  geom_line()+geom_ribbon(alpha=0.2,color=NA)+
  scale_fill_manual(name='model:',
                    values = c('TB incidence'='blue',
                               'TB notifications'='black','TB/HIV notifications'='red'))+
  scale_color_manual(name='model:',
                    values = c('TB incidence'='blue',
                               'TB notifications'='black','TB/HIV notifications'='red'))
##   theme(legend.title = element_text(face = "bold"))

legd2 <- as_ggplot(get_legend(tmp))




## ## make.Age.inc.plt(cniso,MR$ss)

## ## === notifications & incidence
## ## WHO comparison 
make.WHOcompare.data <- function(cniso,MR,na.rm=FALSE){
  ## if('year' %in% names(MR) & !'t' %in% names(MR)) MR[,t:=year]
  MRM <- melt(MR[,.(t,id,
                    Itot=Itot/epfac,
                    IHtot=(Ipn+Ipp)/epfac,
                    Mtot=mort/epfac,
                    MHtot=(mortH+mortA)/epfac)],id=c('id','t'))
  yrs <- MR[,min(floor(t))]:MR[,max(floor(t))]
  wyrs <- data.table(yrs=yrs) #nearest to mid pts
  for(y in yrs) wyrs[yrs==y,
                     tid:=which.min(abs(MR[id==1,t]-(y+0.5)))]
  wyrs[,t:=MR[id==1][tid,t]]
  MRM <- MRM[t %in% wyrs$t]
  MRM <- merge(MRM,wyrs[,.(year=yrs,t)],by='t')
  MRM <- MRM[,.(mid=median(value,na.rm=na.rm),
                hi=quantile(value,.975,na.rm=na.rm),
                lo=quantile(value,.025,na.rm=na.rm)),
             by=.(year,variable)]
  etmp <- melt(E[cniso][,.(year,
                           Itot=e_inc_num,
                           Itot.hi=e_inc_num_hi,
                           Itot.lo=e_inc_num_lo,
                           IHtot=e_inc_tbhiv_num,
                           IHtot.hi=e_inc_tbhiv_num_hi,
                           IHtot.lo=e_inc_tbhiv_num_lo,
                           Mtot=e_mort_num,
                           Mtot.lo=e_mort_num_lo,
                           Mtot.hi=e_mort_num_hi,
                           MHtot=e_mort_tbhiv_num,
                           MHtot.lo=e_mort_tbhiv_num_lo,
                           MHtot.hi=e_mort_tbhiv_num_hi
                           )],
               id='year')
  etmp[,c('variable','type'):=tstrsplit(variable,"\\.")]
  etmp[is.na(type),type:='mid']
  etmp <- dcast(etmp,year+variable~type,value.var='value')
  ## merge
  etmp[,estimate:='WHO']; MRM[,estimate:='transmission model']
  commonyrs <- intersect(etmp$year,MRM$year)
  B <- rbind(etmp[year %in% commonyrs],
             MRM[year %in% commonyrs],
             use.names = TRUE)
  B[,iso3:=cniso]
  B
}

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

