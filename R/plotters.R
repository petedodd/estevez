## various plotting utilities & also utilities for processing output data


absspace <- function(x,...) {             #works
   format(abs(x), ..., big.mark=" ",scientific = FALSE, trim = TRUE)
}


##  ===================== DATA FORMATTING ============
## turn odin model output into sensible data table
out2df <- function(out){
  nmz <- grep('DA|DH|DX|XLL|HLL|ALL|XLR|HLR|ALR|XU|HU|AU|XR|HR|AR|aris|IreactX|IrelX|TA|TH|TX|Inn|Ipn|Ipp|Itot|prevtot15plus|prevTtot15plus|pop15plus|Ntot|NtotH|NtotA|poptot|popart|poph|poph1549|pop1549|mort|mortH',
              colnames(out),value=TRUE)
  nmz <- c('t',nmz)
  outdat <- out[,nmz]
  colnames(outdat)[1] <- 'year'
  outdat <- as.data.table(outdat)
  outdat <- melt(outdat,id.vars = 'year')
  outdat[,hstate:='X']                     #states
  outdat[grep('A',variable),hstate:='A']
  outdat[grep('H',variable),hstate:='H']
  outdat[,tbstate:='no TB']
  outdat[grepl('DA',variable) |
         grepl('DX',variable) |
         grepl('DH',variable),tbstate:='TB']
  outdat[,Lstate:='LTBI+']
  outdat[grepl('U',variable),Lstate:='LTBI-']
  outdat[grepl('LL',variable),Ltype:='distant']
  outdat[grepl('LR',variable),Ltype:='recent']
  outdat[grepl('XR',variable) |
         grepl('HR',variable) |
         grepl('AR',variable),Ltype:='recovered']
  outdat[Lstate=='LTBI-',Ltype:='LTBI-']
  outdat[,sex:='M']                       #sex
  outdat[grep(",2",variable),sex:='F']
  outdat[,acat:=as.integer(gsub(".+\\[(\\d+),\\d+\\]","\\1",variable))] #age
  outdat[,age:=agz[acat]]
  outdat
}

## make age/sex splits of TBI
makeTBA <- function(out){
    xn <- grep('IX',colnames(out),value=TRUE)
    hn <- grep('IH',colnames(out),value=TRUE)
    an <- grep('IA',colnames(out),value=TRUE)
    nmz <- c('t',xn,hn,an) #'Itot',
    outdat <- out[,nmz]
    colnames(outdat)[1] <- 'year'
    outdat <- as.data.table(outdat)
    outdat <- melt(outdat,id.vars = 'year')
    outdat[grep('X',variable),hstate:='X']                     #states
    outdat[grep('A',variable),hstate:='A']
    outdat[grep('H',variable),hstate:='H']
    outdat[,tbstate:='Incidence']
    outdat[grep(",1",variable),sex:='M']                       #sex
    outdat[grep(",2",variable),sex:='F']
    outdat[,acat:=as.integer(gsub(".+\\[(\\d+),\\d+\\]","\\1",variable))] #age
    outdat[,age:=agz[acat]]
    outdat
}

## odin output into timeseries relevant to HIV
makeHIVts <- function(out){
    HIV <- as.data.table(out[,c('t','Inn','Itot')])
    HIV[,TBhiv:=1-Inn/Itot]
    HIV
}


## odin output into TB timeseries
makeTBts <- function(outdat){
  TBP <- merge(outdat[tbstate=='TB',
                      .(TBnumber=sum(value)),by=.(year)],
               outdat[tbstate=='no TB',
                      .(pop=sum(value)),by=.(year)],
                 by='year'
               )
  TBP <- merge(TBP,
               outdat[tbstate=='no TB' & Lstate=='LTBI+',
                      .(L=sum(value)),by=.(year)],
               by='year'
               )
  tmp <- outdat[variable%in%c('IrelX','IreactX','aris',
                              'Itot','Inn','Ipn','Ipp',
                              "prevtot15plus","prevTtot15plus",
                              "pop15plus","poptot",'popart','poph',
                              "poph1549","pop1549",
                              'mort','morth',
                              'Ntot','NtotH','NtotA'),
                .(value,year,variable)]
  tmp <- dcast(tmp,year ~ variable,value.var = 'value')
  names(tmp)[names(tmp)=='aris'] <- 'ari'
  TBP <- merge(TBP,tmp,by='year')
  tmp <- outdat[grepl("^D",variable),.(prev=sum(value)),by=year]
  TBP <- merge(TBP,tmp,by='year')
  tmp <- outdat[grepl("^T",variable),.(tx=sum(value)),by=year]
  TBP <- merge(TBP,tmp,by='year')
  TBP
}


## ======================= PLOTS =====================
## HIV timeseries plots
HIVtsplots <- function(rootname,H){
    rfn <- glue(rootname)
    GP <- ggplot(data=H,
                 aes(t,TBhiv)) + geom_line() +
        xlab('Year') + ylab('Proportion of TB incidence HIV+ve') +
        scale_y_continuous(label=percent)
    fn <- rfn + 'HIVinTB.pdf'
    ggsave(GP,filename=fn); print(fn)
}



## TB timeseries plots
TBtsplots <- function(TBP, rootname){
    if(!is.null(rootname)){ rn <- glue(rootname); } else {rn <- glue('');}
    ## prevalence graph
    fn <- rn + 'TBprev.pdf'
    GP1 <- ggplot(TBP,aes(year,1e5*TBnumber/pop)) + geom_line() +
        expand_limits(y=0)+
        xlab('Year') + ylab('TB prevalence per 100,000')
    if(!is.null(rootname)) {ggsave(GP1,filename=fn,w=7,h=7); print(fn)}
    ## LTBI prevalence
    fn <- rn + 'LTBIprev.pdf'
    GP2 <- ggplot(TBP,aes(year,L/pop)) + geom_line() + xlab('Year') +
        ylab('LTBI prevalence') + scale_y_continuous(label=percent,limits = c(0,NA))
    if(!is.null(rootname)) {ggsave(GP2,filename=fn,w=7,h=7); print(fn)}
    ## TB incidence
    fn <- rn + 'TBI.pdf'
    GP3 <- ggplot(TBP,aes(year,1e5*Itot/pop)) + geom_line() +
        expand_limits(y=0)+
        xlab('Year') + ylab('TB incidence per 100,000')
    if(!is.null(rootname)) {ggsave(GP3,filename=fn,w=7,h=7); print(fn)}
    ## ARI
    fn <- rn + 'ARI.pdf'
    GP4 <- ggplot(TBP,aes(year,ari)) + geom_line() +
      scale_y_continuous(label=percent,limits = c(0,NA))+
      xlab('Year') + ylab('ARI')
    if(!is.null(rootname)) {ggsave(GP4,filename=fn,w=7,h=7); print(fn)}
    ## proportion
    fn <- rn + 'inctypes.pdf'
    tmp <- melt(TBP[,.(year,reactivation=IreactX,
                       relapse=IrelX,
                       recent=1-IrelX-IreactX)],
                id='year')
    GP5 <- ggplot(tmp,aes(year,value,col=variable)) +
      geom_line() +
      scale_y_continuous(label=percent,limits = c(0,NA))+
      xlab('Year') + ylab('Proportion of HIV-negative TB incidence')+
      theme(legend.position = 'top')
    if(!is.null(rootname)) {ggsave(GP5,filename=fn,w=7,h=7); print(fn)}
    if(is.null(rootname)) return(list(GP1,GP2,GP3,GP4,GP5))
}



ppairs <- function(D,alph=0.5){
  GP <- ggpairs(D,switch='both',
               lower=list(continuous=wrap('density',col='black',alpha=alph)),
               upper='blank')+
    theme_bw() + theme(panel.spacing.x=unit(0.1, "lines"),
                       panel.spacing.y=unit(0.1, "lines"),
                       axis.text.x = element_text(angle=45,hjust=1),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       strip.background = element_blank(),
                       strip.placement = "outside",
                       panel.border = element_rect(colour = "black")) +
    ggpubr::grids()
  GP
}

## demographic plots
demoplots <- function(rootname,outdat,yr){
    bnm <- glue(rootname)
    ## ages
    lagz <- seq(from=0,to=80,by=5)
    ragz <- c(lagz[2:length(lagz)]-1,Inf)
    agz <- paste(lagz,ragz,sep=',')
    agz <- paste0('[',agz,')')
    ## reduce data
    dyr <- outdat[which.min(abs(year-yr)),year] #closest year in data
    ## total
    tmp <- outdat[year==dyr,.(value=sum(value)),by=.(age,sex)]
    tmp$age <- factor(tmp$age,levels=agz,ordered = TRUE)
    ## LTBI
    tmpL <- outdat[year==yr,.(value=sum(value)),by=.(age,sex,Lstate)]
    tmpL$age <- factor(tmpL$age,levels=agz,ordered = TRUE)

    ## basic  demography
    fn <- bnm + 'demo.pdf'
    GP <- ggplot(data=tmp,aes(x=age,fill=sex))+
        geom_bar(data=tmp,aes(y=ifelse(sex=='M',value,-value)),stat='identity') +
        coord_flip()+geom_abline(intercept=0,slope=0) +
        ylab('Number') + xlab('Age') +
        theme(axis.text.y = element_text(size = rel(0.95))) +
        scale_y_continuous(label=absspace)
    ggsave(GP,file=fn); print(fn)
    
    ## LTBI1
    fn <- bnm + 'LTBItot.pdf'
    GP <- ggplot(data=tmpL,aes(x=age,fill=Lstate))+
        geom_bar(data=tmpL,aes(y=ifelse(sex=='M',value,-value)),stat='identity') +
        coord_flip()+geom_abline(intercept=0,slope=0) +
        ylab('Number') + xlab('Age') +
        theme(axis.text.y = element_text(size = rel(0.95))) +
        scale_y_continuous(label=absspace)
    ggsave(GP,file=fn); print(fn)
    ## LTBI1
    fn <- bnm + 'LTBIprop.pdf'
    GP <- ggplot(data=tmpL,aes(x=age,fill=Lstate))+
        geom_bar(data=tmpL,aes(y=ifelse(sex=='M',value,-value)),stat='identity',position='fill') +
        coord_flip()+geom_abline(intercept=0,slope=0) +
        ylab('Number') + xlab('Age') +
        theme(axis.text.y = element_text(size = rel(0.95))) +
        scale_y_continuous(label=absspace)
    ggsave(GP,file=fn); print(fn)
}

## TB demography
TBAplots <- function(rootname,outdat,yr){
    bnm <- glue(rootname)
    ## ages
    lagz <- seq(from=0,to=80,by=5)
    ragz <- c(lagz[2:length(lagz)]-1,Inf)
    agz <- paste(lagz,ragz,sep=',')
    agz <- paste0('[',agz,')')
    ## reduce data
    dyr <- outdat[which.min(abs(year-yr)),year] #closest year in data
    ## total
    tmp <- outdat[year==dyr,.(value=sum(value)),by=.(age,sex)]
    tmp$age <- factor(tmp$age,levels=agz,ordered = TRUE)
    tmp[,value:=value/sum(value)] #make proportions

    ## TBI  demography
    fn <- bnm + 'tbidemo.pdf'
    GP <- ggplot(data=tmp,aes(x=age,fill=sex))+
        geom_bar(data=tmp,aes(y=ifelse(sex=='M',value,-value)),
                 stat='identity') +
        coord_flip()+geom_abline(intercept=0,slope=0) +
        ylab('Proportion') + xlab('Age') +
        theme(axis.text.y = element_text(size = rel(0.95))) +
        scale_y_continuous(label=absspace)
    ggsave(GP,file=fn); print(fn)
}
