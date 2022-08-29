library(here)
library(glue)

useage <- TRUE #NOTE needs changing as beta prior is different
makePR <- TRUE
gh <- function(x) glue(here(x))
tgt <- 'r'
source(here('R/modelprep.R'))
## not really using yet!
source(here('R/plotters.R'))
source(here('R/dataplots.R'))
source(here('R/dataLL.R')) #contains epfac

## NOTE make directories if lacking
fn <- gh('plots/resfigs')
if(!file.exists(fn)) dir.create(fn)
fn <- gh('plots/figures')
if(!file.exists(fn)) dir.create(fn)


## cor stuff
SCL <- list()
for(cniso in cz){
  cat(cniso,'...\n')
  ## load data
  scrfn <- gh('data/scr_{cniso}.Rdata') #multirun data
  load(scrfn)
  nmz <- rownames(scr)
  scr <- as.data.table(scr)
  scr[,var:=nmz]
  scr <- melt(scr,id='var')
  scr[,iso3:=cniso]
  SCL[[cniso]] <- scr
}
SCL <- rbindlist(SCL)
save(SCL,file=gh('data/SCL.Rdata'))

SCL <- SCL[value<1]
mc <- SCL[,.(value=mean(value)),by=.(var,variable)]
mc <- mc[rev(order(abs(value)))]
mc <- mc[seq(from=1,to=nrow(mc),by=2)] #drop duplicates
mc[1:20]
fwrite(mc,file=gh('data/mc.csv'))

cat('========== HIV/PCI PLOTS ==========\n')

## loop to make data & also per cap inc plots
PCP <- HIVL <- HIVLF <- MRTS <- list()
for(cniso in cz){
  ## load data
  cat(cniso,'...\n')
  ## load data
  mrinf <- gh('data/MR_{cniso}.Rdata') #multirun data
  load(mrinf)
  tmp <- MR$ts
  tmp[,iso3:=cniso]
  MRTS[[cniso]] <- tmp
  tmp <- prepOD(MR$ts[id==1],cniso)$popsr
  tmp[,iso3:=cniso]
  HIVL[[cniso]] <- tmp
  tmp <- MR$ts
  tmp[,iso3:=cniso]
  HIVLF[[cniso]] <- tmp[,.(id,iso3,year=t,
                           hivintb=(Ipn+Ipp)/(Itot+1),
                           IRR=(pop15plus/(poph+1)-1) / (Itot/(Ipn+Ipp+1)-1))]
  PCP[[cniso]] <- make.NIPC.plt(cniso,MR$ts,cknz[iso3==cniso,name],
                                ylb='',xlb='',narm=TRUE) #TODO
}
HIVL <- rbindlist(HIVL)
HIVL[,IRR:=(pop15plus/(poph+1)-1) / (Itot/(Ipn+Ipp+1)-1)]
HIVL[,art:=popart/poph]
HIVL[,hivintb:=(Ipn+Ipp)/(Itot+1)]
HIVLF <- rbindlist(HIVLF)
MRTS <- rbindlist(MRTS)

GA <- ggarrange(plotlist = PCP,ncol=3,nrow=4)
GA <- annotate_figure(GA,
                left = textGrob("TB incidence per 100,000 per year",
                                rot = 90, vjust = 1,
                                gp = gpar(cex = 1.3)),
                bottom = textGrob("Year", gp = gpar(cex = 1.3))
                )


ggsave(GA,file=glue(here('plots/resfigs/TBPpercapita.pdf')),
       w=3*3,
       h=4*3)

ggsave(GA,file=glue(here('plots/figures/Figure4.pdf')),
       w=3*3,
       h=4*3)


## names(MRTS)
MRTS[,prev:=1e5*prevtot15plus/pop15plus]
MRTS <- MRTS[,.(iso3,id,t,Itot,mort,prev)]
MRTS <- melt(MRTS,id=c('t','id','iso3'))
MRTS <- MRTS[,HML(value),by=.(t,iso3,variable)]

save(MRTS,file=gh('data/MRTS.Rdata')) #edit?


MRTS[,year:=floor(t)]
tmp <- MRTS[year==2019]
tmp[,relprec:=1e2*(hi-lo)/(3.92*mid)]
tmp[variable=='Itot']
tmp[variable=='mort',range(relprec)]

fwrite(tmp,file=gh('data/RP19.csv'))

## uncertainty version of this
HIVLF <- HIVLF[,.(hivintb.mid=median(hivintb),
                  hivintb.lo=quantile(hivintb,0.025),
                  hivintb.hi=quantile(hivintb,0.975),
                  IRR.mid=median(IRR),
                  IRR.lo=quantile(IRR,0.025),
                  IRR.hi=quantile(IRR,0.975)),
               by=.(iso3,year)]
HIVLF <- merge(HIVLF,HIVL,by=c('iso3','year'))
HIVLF <- merge(HIVLF,cknz,by='iso3')

HPCu <- ggplot(data=HIVLF,aes(x=year,y=hivintb.mid,
                         ymin=hivintb.lo,ymax=hivintb.hi,
                         col=name,fill=name))+
  geom_ribbon(alpha=.2,col=NA)+
  geom_line()+
  scale_y_continuous(label=percent)+
  scale_x_continuous(limits=c(1980,2022))+
  scale_color_calc()+
  xlab('Year')+ylab('Proportion of incident TB living with HIV')+
  facet_wrap(~name)+
  theme_light()+
  theme(legend.position = 'none')
## HPCu

ggsave(HPCu,file=gh('plots/resfigs/HIVinTBu.pdf'),
       w=3*3,
       h=4*3)

## using medians
tmp <- HIVLF[year==max(year),.(name,hivintb.mid)]
tmp[,year:=2019]

HPC <- ggplot(HIVLF,aes(year,hivintb.mid,col=name))+
  geom_line()+
  geom_point(data=tmp)+
  geom_text_repel(data=tmp,aes(label=name),segment.linetype = 6)+
  scale_y_continuous(label=percent)+
  scale_x_continuous(limits=c(1980,2022))+
  scale_color_calc()+
  xlab('Year')+ylab('Proportion of incident TB living with HIV')+
  theme_classic()+ggpubr::grids()+
  theme(legend.position = 'none')

ggsave(HPC,file=gh('plots/resfigs/HIVinTB.pdf'),w=7,h=5)

save(HIVL,file=gh('data/HIVL.Rdata')) # NOTE
## save(HIVLF,file=gh('data/HIVLF.Rdata')) #
 
## loop to make ARI & also PR  plots
PR <- ARIS <- ARIA <- ARIF <- list()
narm <- FALSE
for(cniso in cz){
  ## load data cniso <- "ETH"
  cat(cniso,'...\n')
  ## load data
  mrinf <- gh('data/MR_{cniso}.Rdata') #multirun data
  load(mrinf)
  tmp <- MR$ARIS
  ARIS[[cniso]] <- tmp[year>=1980,.(iso3=cniso,
                          mid=median(value,na.rm=narm),
                          lo=quantile(value,0.025,na.rm=narm),
                          hi=quantile(value,0.975,na.rm=narm)),
                      by=.(year)]
  tmp <- MR$ARIA
  ARIA[[cniso]] <- tmp[,.(iso3=cniso,
                          mid=median(value,na.rm=narm),
                          lo=quantile(value,0.025,na.rm=narm),
                          hi=quantile(value,0.975,na.rm=narm)),
                       by=age]
  tmp <- MR$ARIF
  ARIF[[cniso]] <- tmp[,.(iso3=cniso,
                          mid=median(value,na.rm=narm),
                          lo=quantile(value,0.025,na.rm=narm),
                          hi=quantile(value,0.975,na.rm=narm)),
                       by=age]
  tmp <- MR$pr
  PR[[cniso]] <- tmp[,.(iso3=cniso,all=sum(all,na.rm=narm),
                        recent=sum(recent,na.rm=narm)),
                     by=.(id,hiv,age,year)]
}
ARIS <- rbindlist(ARIS)
ARIA <- rbindlist(ARIA)
ARIF <- rbindlist(ARIF)
PR <- rbindlist(PR)

## saving out
fn <- gh('data/ARIS.Rdata') #multirun data
save(ARIS,file=fn)
fn <- gh('data/ARIA.Rdata') #multirun data
save(ARIA,file=fn)
fn <- gh('data/ARIF.Rdata') #multirun data
save(ARIF,file=fn)
fn <- gh('data/PR.Rdata') #multirun data
save(PR,file=fn)

## load(fn)

## ARIS graph
ARIS <- merge(ARIS,cknz,by='iso3')
tmp <- ARIS[year==max(year),.(name,mid)]
tmp[,year:=2019]

ARISg <- ggplot(ARIS,aes(year,mid,col=name))+
  geom_line()+
  geom_point(data=tmp)+
  geom_text_repel(data=tmp,aes(label=name),segment.linetype = 6)+
  scale_y_sqrt(label=percent)+
  scale_x_continuous(limits=c(1980,2022))+
  scale_color_calc()+
  xlab('Year')+ylab('Annual risk of TB infection')+
  theme_classic()+ggpubr::grids()+
  theme(legend.position = 'none')

ggsave(ARISg,file=gh('plots/resfigs/ARISgByT.pdf'),w=7,h=5)


## ARIA graph
ARIA <- merge(ARIA,cknz,by='iso3')
tmp <- ARIA[age=="[20,24)",.(name,age,mid)]

ARIAg <- ggplot(ARIA,aes(age,mid,col=name,group=name))+
  geom_line()+geom_point()+
  geom_text_repel(data=tmp,aes(label=name),segment.linetype = 6)+
  scale_y_continuous(label=percent)+
  scale_color_calc()+
  xlab('Age')+ylab('Annual risk of TB infection in 2019')+
  theme_classic()+ggpubr::grids()+rotx+
  theme(legend.position = 'none')

ggsave(ARIAg,file=gh('plots/resfigs/ARIAgByAge.pdf'),w=7,h=5)

## ARIF graph
ARIF <- merge(ARIF,cknz,by='iso3')

tmp <- ARIF[age=="[15,19)",.(name,age,mid)]
ARIFg <- ggplot(ARIF,aes(age,mid,col=name,group=name))+
  geom_line()+geom_point()+
  geom_text_repel(data=tmp,aes(label=name),segment.linetype = 6)+
  scale_y_continuous(label=percent)+
  scale_color_calc()+
  xlab('Age')+ylab('Proportion of transmission in 2019')+
  theme_classic()+ggpubr::grids()+rotx+
  theme(legend.position = 'none')

ggsave(ARIFg,file=gh('plots/resfigs/ARIFgByAge.pdf'),w=7,h=5)

## uncertainty versions

ARIAgu <- ggplot(ARIA,aes(age,mid,ymin=lo,ymax=hi,
                          fill=name,col=name,group=name))+
  geom_ribbon(col=NA,alpha=0.2)+
  geom_line()+geom_point()+
  scale_y_continuous(label=percent)+
  scale_color_calc()+
  scale_fill_calc()+
  xlab('Age')+ylab('Annual risk of TB infection in 2019')+
  theme_light()+
  facet_wrap(~name,scales='free')+
  theme(legend.position = 'none')+ rotx

ggsave(ARIAgu,file=gh('plots/resfigs/ARIAgu.pdf'),
       w=4*3,
       h=4*3)


## ARIF graph

ARIFu <- ggplot(ARIF,aes(age,mid,
                         ymin=lo,ymax=hi,
                         fill=name,
                         col=name,group=name))+
  geom_ribbon(col=NA,alpha=0.2)+
  geom_line()+geom_point()+
  scale_y_continuous(label=percent)+
  scale_color_calc()+
  scale_fill_calc()+
  xlab('Age')+ylab('Proportion of transmission in 2019')+
  theme_light()+
  facet_wrap(~name,scales='free')+
  theme(legend.position = 'none')+rotx

ggsave(ARIFu,file=gh('plots/resfigs/ARIFu.pdf'),
       w=4*3,
       h=4*3)



PR[,hiv2:=ifelse(hiv=='X','HIV-ve','HIV+ve')]

## total
PRt <- PR[,.(all=sum(all,na.rm=narm),
             recent=sum(recent,na.rm=narm)),
          by=.(iso3,id,year)]
PRt <- PRt[,.(mid=median(recent/all,na.rm=narm),
              lo=quantile(recent/all,0.025,na.rm=narm),
              hi=quantile(recent/all,0.975,na.rm=narm)),
           by=.(iso3,year)]
PRt[,hiv:='all']

save(PRt,file=gh('data/PRt.Rdata'))

## HIV
PRh <- PR[,.(all=sum(all,na.rm=narm),
             recent=sum(recent,na.rm=narm)),
          by=.(iso3,id,hiv2,year)]
PRh <- PRh[,.(mid=median(recent/all,na.rm=narm),
              lo=quantile(recent/all,0.025,na.rm=narm),
              hi=quantile(recent/all,0.975,na.rm=narm)),
             by=.(iso3,hiv=hiv2,year)]

## age
PRa <- PR[,.(all=sum(all,na.rm=narm),
             recent=sum(recent,na.rm=narm)),
          by=.(iso3,id,age,year)]
PRa <- PRa[,.(mid=median(recent/all,na.rm=narm),
              lo=quantile(recent/all,0.025,na.rm=narm),
              hi=quantile(recent/all,0.975,na.rm=narm)),
           by=.(iso3,age,year)]
PRa$age <- factor(PRa$age,levels=agz,ordered = TRUE)

save(PRa,file=gh('data/PRa.Rdata'))


## age version
PRa <- PRa[year==2019]
PRa <- merge(PRa,cknz,by='iso3')
tmp <- PRa[age=="[10,14)",.(name,age,mid)]

PRag <- ggplot(PRa,aes(age,mid,col=name,group=name))+
  geom_line()+geom_point()+
  geom_text_repel(data=tmp,aes(label=name),segment.linetype = 6)+
  scale_y_continuous(label=percent,limits=c(0,1))+
  scale_color_calc()+
  xlab('Age')+
  ylab('Proportion recent incidence 2019')+
  ## Proportion of incidence in 2019 from transmision since 2017
  theme_classic()+ggpubr::grids()+rotx+
  theme(legend.position = 'none')

## PRag

ggsave(PRag,file=gh('plots/resfigs/PRag.pdf'),w=7,h=5)

## unc
PRagu <- ggplot(PRa,aes(age,mid,
                        ymin=lo,ymax=hi,fill=name,
                        col=name,group=name))+
  geom_ribbon(col=NA,alpha=0.2)+
  geom_line()+geom_point()+
  scale_y_continuous(label=percent,limits=c(0,1))+
  scale_color_calc()+
  scale_fill_calc()+
  xlab('Age')+
  ylab('Proportion recent incidence 2019')+
  theme_light()+
  facet_wrap(~name,scales='free')+
  theme(legend.position = 'none')+rotx

ggsave(PRagu,file=gh('plots/resfigs/PRagu.pdf'),
       w=4*3,
       h=4*3)


## total version
PRt <- PRt[year==2019]
PRt <- merge(PRt,cknz,by='iso3')

PRtg <- ggplot(PRt,aes(name,mid,ymin=lo,ymax=hi,col=name))+
  geom_errorbar(width=0)+geom_point()+
  scale_y_continuous(label=percent)+
  scale_color_calc()+
  xlab('Country')+
  ylab('Proportion recent incidence 2019')+
  coord_flip()+
  theme_classic()+ggpubr::grids()+
  theme(legend.position = 'none')

## PRtg

ggsave(PRtg,file=gh('plots/resfigs/PRtg.pdf'),w=6,h=6)


## TODO understand
## by HIV
PRh <- PRh[year==2019]
PRh <- merge(PRh,cknz,by='iso3')
PRta <- rbind(PRt,PRh)
psn <- position_dodge(0.5)

PRtga <- ggplot(PRta,aes(name,mid,ymin=lo,ymax=hi,col=name,shape=hiv))+
  geom_errorbar(width=0,position = psn)+geom_point(position = psn)+
  scale_y_continuous(label=percent)+
  scale_color_calc(guide='none')+
  xlab('Country')+
  ylab('Proportion recent incidence 2019')+
  coord_flip()+
  theme_classic()+ggpubr::grids()+
  theme(legend.position = 'top')

## PRtga

ggsave(PRtga,file=gh('plots/resfigs/PRtga.pdf'),w=6,h=6)



GA2 <- ggarrange(PRtg,PRag,labels=c('A','B'),widths = c(1,3),align='h')

## GA2

ggsave(GA2,file=gh('plots/resfigs/mg_PRg.pdf'),w=13,h=6)

## multigraph
GA3 <- ggarrange(ARIAg,ARIFg,labels=c('A','B'),widths = c(1,1),align='h')
GA2 <- ggarrange(PRtg,PRag,labels=c('C','D'),widths = c(1,3),align='h')
GA4 <- GA3/GA2

ggsave(GA4,file=gh('plots/resfigs/Fig3.pdf'),w=13,h=13)

nttl <- 'Proportion of incidence in 2019 from transmision since 2017'
PRtg <- PRtg + expand_limits(y=0.5)
PRtg <- PRtg+ylab(nttl)
PRag <- PRag+ylab(nttl)

GA1 <- ggarrange(HPC,labels='A')
GA3 <- ggarrange(ARIAg,ARIFg,labels=c('B','C'),widths = c(1,1),align='h')
GA2 <- ggarrange(PRtg,PRag,labels=c('D','E'),widths = c(1,1),align='h')
GA4a <- GA1/GA3/GA2

ggsave(GA4a,file=gh('plots/figures/Figure5.pdf'),w=13,h=18)

## ======================================
## --- all-doer
MakePlotsList <- function(CL,todo=1:6,na.rm=FALSE){
  k <- 1
  PltList <- list()
  for(cniso in CL){

    ## load data
    cat(cniso,'...\n')
    ## load data
    mrinf <- glue(here('data/MR_{cniso}.Rdata')) #multirun data
    load(mrinf)
    odfn <- glue(here('data/OD_{cniso}.Rdata'))
    load(file=odfn)
    sminf <- glue(here('data/sampe_{cniso}.Rdata')) #sample data
    load(sminf)

    ## ============== epi outputs ==============
    if(1 %in% todo){
      cat('...poph plot\n')
      Dt <- make.poph.plt(cniso,MR$ts[id==1])
      PltList[[k]] <- Dt; k <- k+1
    }

    if(2 %in% todo){
      Ds <- make.popsnap.plt(cniso,OD)
      PltList[[k]] <- Ds; k <- k+1
    }

    if(3 %in% todo){
      cat('...prev plot\n')
      Pplt <- make.prev.plt(cniso,MR$ts,na.rm = na.rm)
      PltList[[k]] <- Pplt; k <- k+1
    }

    if(4 %in% todo){
      cat('...NI plot\n')
      NIplt <- make.NI.plt(cniso,MR$ts,narm = na.rm)
      PltList[[k]] <- NIplt; k <- k+1
    }

    if(6 %in% todo){
      cat('...age plot\n')
      Aplt <- make.Age.inc.plt(cniso,MR$ss,narm = na.rm)
      PltList[[k]] <- Aplt; k <- k+1
    }

  }
  return(PltList)
}



## ----
## make data
ComparePlotData <- function(CL){
  CPD <- list()
  for(cn in CL){
    cat(cn,'...\n')
    ## load data
    sdfn <- glue(here('data/sampe_{cn}.Rdata'))
    load(sdfn)
    sampe <- melt(sampe)
    sampe[,iso3:=cn]
    CPD[[cn]] <- sampe
  }
  rbindlist(CPD)
}


## =========== prior vs posterior =========
cat('========== DOING BOXPLOTS ==========\n')
CPD <- ComparePlotData(cz) #make posterior data

  ## prior data
PR <- list()
for(nm in names(hyperparms)){
  cat(nm,'...\n')
  PR[[nm]] <- data.table(mid=qfun(0.5,hyperparms[[nm]]),
                         lo=qfun(0.25,hyperparms[[nm]]),
                         hi=qfun(0.75,hyperparms[[nm]]))
  PR[[nm]][,vrbl:=nm]
}
PR <- rbindlist(PR)
PR <- melt(PR,id='vrbl')
names(PR) <- c('variable','qty','value')
PR <- merge(PR,parmkey[,.(variable=model,varname=writeup)],
            by='variable')
PR$varname <- factor(PR$varname,
                     levels=parmkey[der,writeup],ordered=TRUE)

if(useunicode){
  CPD <- merge(CPD,parmkey[,.(variable=model,varname=writeup)],
               by='variable')
  CPD$varname <- factor(CPD$varname,
                        levels=parmkey[der,writeup],
                        ordered = TRUE)

  GPC <- ggplot(CPD,aes(iso3,value))+
    geom_boxplot(outlier.shape = NA)+
    geom_hline(data=PR,aes(yintercept=value,linetype=qty),col=2)+
    facet_wrap(~varname,scales = 'free')+
    theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
          legend.position = "none")
} else {
  GPC <- ggplot(CPD,aes(iso3,value))+
    geom_boxplot(outlier.shape = NA)+
    geom_hline(data=PR,aes(yintercept=value,linetype=qty),col=2)+
    facet_wrap(~variable,scales = 'free')+
    theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5),
          legend.position = "none")
}

## GPC
ggsave(GPC,
       file=glue(here('plots/resfigs/PostVprior.pdf')),
       w=12,h=8,device = cairo_pdf)## ,family='Japan1')

## ------ WHO comparisons ----------
cat('========== DOING WHO COMPARE ==========\n')

CompareWHOData <- function(CL){
  CPD <- list()
  for(cn in CL){
    cat(cn,'...\n')
    ## load data
    sdfn <- glue(here('data/MR_{cn}.Rdata'))
    load(sdfn)
    CPD[[cn]] <- make.WHOcompare.data(cn,MR$ts,na.rm = TRUE) #TODO
  }
  rbindlist(CPD)
}

WCD <- CompareWHOData(cz)
WCD <- WCD[year==2015]
WCD[,quantity:='Incidence']
WCD[grepl('M',variable),quantity:='Mortality']
WCD[,type:='Total']
WCD[grepl('H',variable),type:='HIV']
WCD <- WCD[,.(iso3,quantity,type,estimate,mid,lo,hi)]
WCDC <- dcast(WCD,
              iso3+quantity + type ~ estimate,
             value.var = c('mid','lo','hi'))

tp <- WCD[,max(hi)]*1.05
WCDC <- merge(WCDC,cknz,by='iso3')

GP <- ggplot(WCDC,aes(x=mid_WHO,y=`mid_transmission model`,
                      xmin=lo_WHO,ymin=`lo_transmission model`,
                      xmax=hi_WHO,ymax=`hi_transmission model`,
                      col=name,label=name))+
  geom_point()+
  geom_text_repel(show.legend = FALSE,max.overlaps=20)+
  geom_errorbar(width=0)+
  geom_errorbarh(height=0)+
  facet_grid(quantity ~ type)+
  scale_color_calc()+
  scale_y_sqrt(label=absspace,limits=c(0,tp))+
  scale_x_sqrt(label=absspace,limits=c(0,tp))+
  coord_fixed()+geom_abline(slope=1,intercept = 0,col=2)+
  xlab('WHO estimate (square root scale)') +
  ylab('Transmission model (square root scale)')+
  theme_light()+## ggpubr::grids()+
  theme(legend.position='none',
        legend.title = element_blank())

sz <- 10
ggsave(GP,file=here('plots/figures/Figure3.pdf'),w=sz,h=sz)
## h=1.1*sz) #if using legend

## ========== work ===============

cat('========== FIG 2 ==========\n')

PltList <- MakePlotsList(cz,na.rm=TRUE)


nc <- 5
lbz <- rep(NA,nc*length(cz))
lbz[seq(from=1,by=nc,to=length(lbz))] <- cz
GA <- ggarrange(plotlist=PltList,
                ncol=nc,nrow=length(cz),
                hjust=0, #-ve -> R
                vjust=1, #+ve -> down
                labels=lbz
                )

ggsave(GA,file=glue(here('plots/resfigs/Figure2.pdf')),
       w=nc*3,
       h=length(cz)*3)


## reorder
PltList2 <- list()
for(i in 1:6){
  ## top half
  PltList2[[i]] <- PltList[[(i-1)*5+1]]+ggtitle(cnm[i])
  PltList2[[6+i]] <- PltList[[(i-1)*5+2]]
  PltList2[[12+i]] <- PltList[[(i-1)*5+3]]
  PltList2[[18+i]] <- PltList[[(i-1)*5+4]]
  PltList2[[24+i]] <- PltList[[(i-1)*5+5]]
  ## bottom half
  PltList2[[30+i]] <- PltList[[(i-1)*5+31]]+ggtitle(cnm[i+6])
  PltList2[[36+i]] <- PltList[[(i-1)*5+32]]
  PltList2[[42+i]] <- PltList[[(i-1)*5+33]]
  PltList2[[48+i]] <- PltList[[(i-1)*5+34]]
  PltList2[[54+i]] <- PltList[[(i-1)*5+35]]
}


lbz2 <- rep(NA,60)
lbz2[c(1,31)] <- 'A)';
lbz2[c(1,31)+6] <- 'B)'
lbz2[c(1,31)+12] <- 'C)'
lbz2[c(1,31)+18] <- 'D)'
lbz2[c(1,31)+24] <- 'E)'

GA2 <- ggarrange(plotlist=PltList2,
                 ncol=6,nrow=10,
                 hjust=0, #-ve -> R
                 vjust=1, #+ve -> down
                 labels=lbz2
                )

ggsave(GA2,file=glue(here('plots/figures/Figure2.pdf')),
       w=6*3,
       h=10*3)






