## making the plots for SA results
library(here)
library(glue)


sensitivity.analysis <- ''
useage <- TRUE #NOTE needs changing as beta prior is different
makePR <- TRUE
gh <- function(x) glue(here(x))
tgt <- 'r'
source(here('R/modelprep.R'))
## not really using yet!
source(here('R/plotters.R'))
source(here('R/dataplots.R'))
source(here('R/dataLL.R')) #contains epfac
cknz[iso3=='SWZ',name:='Eswatini'] #correct



## SA: EPTB=hi/lo + INF1014
saz <- c('','EPTBlo','EPTBhi','INF1014')
fn <- gh('plots/resfigs/SA')
if(!file.exists(fn)) dir.create(fn)

## =============== data for PR, inc, mort, prop of transmission by age
## loop to make ARI & also PR  plots data & estimates
PRSA <- ARIFSA <- ESTA <- list()
for(sa in saz){ #loop over SAs
  cat('SA:',sa,'...\n')
  PR <- ARIF <- EST <- list()
  narm <- FALSE
  for(cniso in cz){
    ## load data cniso <- "ETH"
    cat('...',cniso,'...\n')
    ## load data
    mrinf <- gh('data/MR_{cniso}{sa}.Rdata') #multirun data
    load(mrinf)
    tmp <- MR$pr
    PR[[cniso]] <- tmp[,.(iso3=cniso,all=sum(all,na.rm=narm),
                          recent=sum(recent,na.rm=narm)),
                       by=.(id,hiv,age,year)]
    tmp <- MR$ARIF
    ARIF[[cniso]] <- tmp[,.(iso3=cniso,
                            mid=median(value,na.rm=narm),
                            lo=quantile(value,0.025,na.rm=narm),
                            hi=quantile(value,0.975,na.rm=narm)),
                         by=age]
    tmp <- MR$ts[abs(t-2019)<0.05,.(poptot,Itot,mort,Ipc=1e5*Itot/poptot,Mpc=1e5*mort/poptot)]
    EST[[cniso]] <- tmp[,.(iso3=cniso,
                           I.mid=median(Ipc,na.rm=narm),
                           I.lo=quantile(Ipc,0.025,na.rm=narm),
                           I.hi=quantile(Ipc,0.975,na.rm=narm),
                           M.mid=median(Mpc,na.rm=narm),
                           M.lo=quantile(Mpc,0.025,na.rm=narm),
                           M.hi=quantile(Mpc,0.975,na.rm=narm))]
  }
  PR <- rbindlist(PR)
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
  PRt <- PRt[year==2019]
  PRt <- merge(PRt,cknz,by='iso3')
  PRt[,satring:=sa]
  PRSA[[sa]] <- PRt #record this SA
  ## ARI
  ARIF <- rbindlist(ARIF)
  ARIF <- merge(ARIF,cknz,by='iso3')
  ARIF[,satring:=sa]
  ARIFSA[[sa]] <- ARIF
  ## ests
  EST <- rbindlist(EST)
  EST <- merge(EST,cknz,by='iso3')
  EST[,satring:=sa]
  ESTA[[sa]] <- EST
}
PRSA <- rbindlist(PRSA)
ARIFSA <- rbindlist(ARIFSA)
ESTA <- rbindlist(ESTA)
PRSA[,`sensitivity analysis`:='basecase']
PRSA[satring=='INF1014',`sensitivity analysis`:='infectious children 10-14 years']
PRSA[satring=='EPTBlo',`sensitivity analysis`:='70% PTB']
PRSA[satring=='EPTBhi',`sensitivity analysis`:='90% PTB']
PRSA$`sensitivity analysis` <- factor(PRSA$`sensitivity analysis`,
                                      levels=PRSA[,unique(`sensitivity analysis`)],ordered=TRUE)
ARIFSA[,`sensitivity analysis`:='basecase']
ARIFSA[satring=='INF1014',`sensitivity analysis`:='infectious children 10-14 years']
ARIFSA[satring=='EPTBlo',`sensitivity analysis`:='70% PTB']
ARIFSA[satring=='EPTBhi',`sensitivity analysis`:='90% PTB']
ARIFSA$`sensitivity analysis` <- factor(ARIFSA$`sensitivity analysis`,
                                      levels=ARIFSA[,unique(`sensitivity analysis`)],ordered=TRUE)
ESTA[,`sensitivity analysis`:='basecase']
ESTA[satring=='INF1014',`sensitivity analysis`:='infectious children 10-14 years']
ESTA[satring=='EPTBlo',`sensitivity analysis`:='70% PTB']
ESTA[satring=='EPTBhi',`sensitivity analysis`:='90% PTB']
ESTA$`sensitivity analysis` <- factor(ESTA$`sensitivity analysis`,
                                      levels=ESTA[,unique(`sensitivity analysis`)],ordered=TRUE)

## PRSA
## ARIFSA
## ESTSA

## ================= plots
## ----------------- PR

## plot
dg <- position_dodge(width = 0.5)
PRtgSA <- ggplot(PRSA,aes(name,mid,ymin=lo,ymax=hi,col=name,shape=`sensitivity analysis`))+
  ## geom_errorbar(width=0)+geom_point()+
  geom_pointrange(position = dg)+
  scale_y_continuous(label=percent)+
  scale_x_discrete(limits = rev)+
  scale_color_calc(guide='none')+
  xlab('Country')+
  ylab('Proportion recent incidence 2019')+
  coord_flip()+
  theme_classic()+ggpubr::grids()+
  theme(legend.position = 'top')
## PRtgSA

ggsave(PRtgSA,file=gh('plots/resfigs/SA/PRtgSA.pdf'),w=8,h=8)
ggsave(PRtgSA,file=gh('plots/resfigs/SA/PRtgSA.png'),w=8,h=8)

## stat:
tmp <- PRSA[,.(iso3,mid, satring)]
tmp[satring=='',satring:='basecase']
tmp <- dcast(tmp,iso3~satring,value.var='mid')
tmp <- tmp[,.(iso3,EPTBlo=EPTBlo/basecase,EPTBhi=EPTBhi/basecase,INF1014=INF1014/basecase)]
tmp <- tmp[,.(EPTBlo=median(EPTBlo),EPTBhi=median(EPTBhi),INF1014=median(INF1014))]
fwrite(tmp,file=gh('plots/resfigs/SA/PRtgSA.csv'))

## ----------------- tx by age

ARIFSA

## NOTE only do relevant age group
tmpa <- ARIFSA[age=='[10,14)' & satring =='INF1014']

ARIFsag <- ggplot(tmpa,aes(name,mid,ymin=lo,ymax=hi,col=name))+
  geom_pointrange()+
  scale_y_continuous(label=percent)+
  scale_x_discrete(limits = rev)+
  scale_color_calc(guide='none')+
  xlab('Country')+
  ylab('Proportion of transmission from children 10-14 years 2019')+
  coord_flip()+
  theme_classic()+ggpubr::grids()+
  theme(legend.position = 'top')

## ARIFsag

ggsave(ARIFsag,file=gh('plots/resfigs/SA/ARIFsa.pdf'),w=7,h=5)
ggsave(ARIFsag,file=gh('plots/resfigs/SA/ARIFsa.png'),w=7,h=5)

## stat:
tmp <- tmpa[,.(M=median(mid),H=max(mid),L=min(mid))]
fwrite(tmp,file=gh('plots/resfigs/SA/ARIFsa.csv'))


## ------------ incidence and mortality per capita?
ESTSA <- melt(ESTA[,.(country=name,`sensitivity analysis`,
                      I.mid,I.lo,I.hi,M.mid,M.lo,M.hi)],
              id=c('country','sensitivity analysis'))
ESTSA[,c('quantity','statistic'):=tstrsplit(variable,split='\\.')]
ESTSA[,quantity:=ifelse(quantity=='I','Incidence','Mortality')]
ESTSA <- dcast(data=ESTSA[,.(country,quantity,`sensitivity analysis`,value,statistic)],
               formula=country+`sensitivity analysis`+quantity~statistic,
               value.var = 'value')
ESTSA$`sensitivity analysis` <- factor(ESTSA$`sensitivity analysis`,
                                      levels=ESTA[,unique(`sensitivity analysis`)],ordered=TRUE)


dg <- position_dodge(width = 0.5)

ESTsp <- ggplot(ESTSA,aes(country,y=mid,ymin=lo,ymax=hi,
                 col=country,shape=`sensitivity analysis`))+
  geom_pointrange(position = dg)+
  scale_x_discrete(limits = rev)+
  scale_color_calc(guide='none')+
  xlab('Country')+
  ylab('Rate per 100,000 population in 2019')+
  coord_flip()+
  facet_wrap(~quantity,scales='free')+
  theme_light()+
  theme(legend.position = 'top')

ggsave(ESTsp,file=gh('plots/resfigs/SA/ESTsa.pdf'),w=12,h=8)
ggsave(ESTsp,file=gh('plots/resfigs/SA/ESTsa.png'),w=12,h=8)

## stat:
tmp <- ESTSA[,.(country,mid, sa=`sensitivity analysis`,quantity)]

tmp <- dcast(tmp,country + quantity~sa,value.var='mid')
tmp <- tmp[,.(country,EPTBlo=`70% PTB`/basecase,EPTBhi=`90% PTB`/basecase,
              INF1014=`infectious children 10-14 years`/basecase)]
tmp <- tmp[,.(EPTBlo=median(EPTBlo),EPTBhi=median(EPTBhi),INF1014=median(INF1014))]

fwrite(tmp,file=gh('plots/resfigs/SA/ESTsa.csv'))
