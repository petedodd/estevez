## data sources for fitting
library(data.table)
library(glue)
library(metafor)
library(here)

cz <- c('ETH','KEN','LSO','MOZ','MWI','NGA','SWZ','TZA','UGA','ZAF','ZMB','ZWE')


## --- where doing stuff
wdir <- glue(here('indata/fitdata/'))

## === notifications

## read in
ddir <- glue('~/Dropbox/Documents/WHO_TBreports/data2021/') #data dir
fn <- ddir + 'TB_notifications_2021-10-14.csv'
N <- fread(fn)

## names explore
nmz <- names(N)
grep('newinc',nmz,value=TRUE)
grep('ret',nmz,value=TRUE)
grep('hiv',nmz,value=TRUE)


## hivtest
## hivtest_pos

## --- time series
## === estimates
fn <- ddir + 'TB_burden_countries_2021-10-14.csv'
E <- fread(fn)
setkey(E,iso3,year)

## save
fn <- wdir + 'E.Rdata'
save(E,file=fn)


##total new/relapse by time/country, also retreatment not relapse
(NT <- N[,.(iso3,year,c_newinc,ret_nrel,## hivtest,hivtest_pos,
            newrel_art,newrel_hivpos
            )])

setkey(NT,iso3)

## treatment outcomes
fn <- ddir + 'TB_outcomes_2021-10-14.csv'
NO <- fread(fn)
grep('coh',names(NO),value=TRUE)

NTO <- NO[,.(iso3,year,newrel_coh,newrel_died,
             tbhiv_coh,tbhiv_died)]


## join
NT <- merge(NT,NTO,by=c('iso3','year'),all=TRUE)

NT <- merge(NT,E[,.(iso3,year,e_pop_num)],by=c('iso3','year'),all.x = TRUE) #oh dear, only upto 1980 TODO

setkey(NT,iso3,year)

NT['ZWE'][year>2015] #check

## save
fn <- wdir + 'NT.Rdata'
save(NT,file=fn)


## ----- HIV statistics
fn <- wdir + 'NT.Rdata'
load(fn)

## compute a meta-analytic fraction as a target:
TMP <- NT[iso3 %in% cz]
TMP <- TMP[!is.na(newrel_hivpos) & !is.na(c_newinc),
           .(iso3,year,c_newinc,newrel_hivpos)]

HF <- TMP[,{
  mam <- rma.glmm(measure = "PLO", #  binomial w/ logit link
                  xi = newrel_hivpos,     # numerator
                  ni = c_newinc,          # denominator
                  data = .SD)
  list(logitH=as.numeric(mam$b),logitH2=mam$tau2)
},
by=iso3]




fn <- wdir + 'HF.Rdata'
save(HF,file=fn)


## --- age/sex
## old age disaggregation
d2a <- c('04','514','1524','2534','3544','4554','5564','65')
agz <- c('0-4','5-14','15-24','25-34','35-44','45-54','55-64','65+')
hash <- data.table(sag=d2a,acat=agz)
oldage <- c(outer(c('sp','sn','ep'),c('m','f'),
                  function(x,y)paste(x,y,sep='_')))
oldage <- c(outer(oldage,d2a,paste0))
oldage <- paste0('new_',oldage)

## new age disaggregation
newage <- c(outer(c('m','f'),d2a,paste0))
newage <- paste0('newrel_',newage)
newnewage <- c('newrel_f59','newrel_f1014',
               'newrel_f1519','newrel_f2024',
               'newrel_m59','newrel_m1014',
               'newrel_m1519','newrel_m2024')

ioa <- c('iso3','year',oldage) #add iso3/year
ina <- c('iso3','year',newage) #add iso3/year
inna <- c('iso3','year',newnewage) #add iso3/year

## old age groups
nao <- N[,..ioa]
nao <- melt(nao,id=c('iso3','year'))
nao[,tota:=sum(value,na.rm=TRUE),by=.(iso3,year)]
nao <- nao[tota>0]
nao[is.na(value),value:=0]
nao[,variable:=gsub('new_','',variable)]
nao[,c('type','sag'):=tstrsplit(variable,split='_')]
nao[,sex:=ifelse(grepl('m',sag),'M','F')]
nao[,sag:=gsub('m|f','',sag)]
nao <- merge(nao,hash,by='sag',all.x=TRUE)
natype <- nao[,.(iso3,year,sex,acat,value,tota)] #keep for type

nao <- natype[,.(value=sum(value)),
              by=.(iso3,year,sex,acat,tota)] #sum over type

## save these
fn <- wdir + 'natype.Rdata'
save(natype,file=fn)
fn <- wdir + 'nao.Rdata'
save(nao,file=fn)


## newnew age groups
nna <- N[,..inna]
nna <- melt(nna,id=c('iso3','year'))
nna[,variable:=gsub('newrel_','',variable)]
nna[,table(is.na(value),year)]
nna <- nna[year>2018]
nna[,variable:=gsub('59','514',variable)] #old vars
nna[,variable:=gsub('1014','514',variable)]
nna[,variable:=gsub('1519','1524',variable)] #old vars
nna[,variable:=gsub('2024','1524',variable)]
nna[,unique(variable)]
nna <- nna[,.(value=sum(value,na.rm=TRUE)),by=.(iso3,year,variable)]

## new age groups
nan <- N[,..ina]
nan <- melt(nan,id=c('iso3','year'))
nan[,variable:=gsub('newrel_','',variable)]
nan <- merge(nan,nna[,.(iso3,year,variable,newvalue=value)],
             by=c('iso3','year','variable'),all.x=TRUE)
nan[,table(is.na(value),is.na(newvalue))]
nan[is.na(value) & !is.na(newvalue),value:=newvalue] #swap fine if used
nan[,newvalue:=NULL]
nan[,tota:=sum(value,na.rm=TRUE),by=.(iso3,year)]
nan <- nan[tota>0]
nan[is.na(value),value:=0]
nan[,sex:=ifelse(grepl('m',variable),'M','F')]
nan[,variable:=gsub('m|f','',variable)]
nan <- merge(nan,hash,by.x='variable',by.y='sag',all.x=TRUE)
nan <- nan[,.(iso3,year,sex,acat,value,tota)] #keep

fn <- wdir + 'nan.Rdata'
save(nan,file=fn)

## check old/new are disjoint:
intersect(nao[,unique(paste0(iso3,year))],
          nan[,unique(paste0(iso3,year))])

## join
NAS <- rbind(nao,nan)
NAS$sex <- factor(NAS$sex,levels=c('M','F'))
NAS$acat <- factor(NAS$acat,levels=agz,ordered = TRUE)
setkey(NAS,iso3,year)

NAS['ZWE'][year==2015] #check

## save
fn <- wdir + 'NAS.Rdata'
save(NAS,file=fn)


## check diffs
cz <- c('ETH','KEN','LSO','MOZ','MWI','NGA','SWZ','TZA','UGA','ZAF','ZMB','ZWE')
## library(glue)
## gh <- function(x)glue(here(x))

## for(cniso in cz){
##   tmp <- NT[cniso][!is.na(c_newinc),.(year,c_newinc)]
##   tmp2 <- NAS[iso3==cniso]
##   PP <- ggplot(unique(tmp2[,.(year,tota)]),aes(year,tota))+
##     geom_point()+
##     geom_point(data=tmp,aes(year,c_newinc),col=2)
##   ggsave(PP,file=gh('fitdata/agenote_{cniso}.pdf'),w=7,h=5)
## }

## === prevalence
fn <- ddir + 'notpublic/prev.rda'
load(fn)
P <- prev
P[,unique(age.group)]
P[,acat:=gsub('_','-',age.group)]
P[,acat:=gsub('plus','+',acat)]
P[,age.group:=NULL]
P[,unique(sex)]
P[sex=='a',sex:='all']
P[sex=='f',sex:='F']
P[sex=='m',sex:='M']

fn <- wdir + 'P.Rdata'
save(P,file=fn)

pcs <- P[,unique(iso3)] #only interested in prevalence survey countries


## --- age patterns
fn <- wdir + 'P.Rdata'
load(file=fn)


logit <- function(x)log(x/(1-x))
ilogit <- function(x)1/(1+exp(-x))

P <- P[iso3 %in% cz]
P <- P[sex=='all' & acat!='all']
P <- P[case.type=='b' & acat!='15+']
P$acat <- factor(P$acat,levels=c('15-24','25-34','35-44',
                                 '45-54','55-64','65+'),
                 ordered = TRUE)
P <- P[order(iso3,acat)]
P[,m:=(1e-5*prev)]
P[,v:=(((1e-5*prev.hi)-(1e-5*prev.lo))/(2*1.96))^2]

tmp <- P[acat=='15-24']
tmp <- tmp[,.(iso3,refm=m,refv=v)]

P <- merge(P,tmp,by='iso3',all.x = TRUE)


P[,RR:=m/refm]
P[,RR.sd:=RR * sqrt(m^2/v + refm^2/refv)]

P <- P[acat!='15-24']
PA <- P[,.(iso3,acat,RR,RR.sd)]

fn <- wdir + 'PA.Rdata'
save(PA,file=fn)



## === mortality?
fn <- '~/Dropbox/Holocron/TBMkids/VRwork/out/D10M.Rdata'
load(fn)
D10M <- D10M[,.(deaths=sum(deaths)),
             by=.(Admin1,SubDiv,Year,Sex,Age,Frmat,IM_Frmat,iso3)]

D10M <- D10M[iso3 %in% pcs]

summary(D10M)
D10M[,c('Admin1','SubDiv'):=NULL]

D10M[,unique(Age)]

VRtot <- D10M[,.(totd=sum(deaths)),by=.(iso3,Year)]

D10M[,sex:=c('M','F')[Sex]]
D10M <- D10M[!is.na(sex)]

## relabel ages
D10M[Age %in% c('0 day','0-6 days','1-6 days','7-27 days',
                '28-365 days','0-365 days',
                '0','1','2','3','4'),acat:='0-4']
D10M[Age %in% c('5-9','10-14'),acat:='5-14']
D10M[Age %in% c('15-19','20-24'),acat:='15-24']
D10M[Age%in%c('25-29','30-34','25-34'),acat:='25-34']
D10M[Age%in%c('35-39','40-44','35-44'),acat:='35-44']
D10M[Age%in%c('45-49','50-54','45-54'),acat:='45-54']
D10M[Age%in%c('55-59','60-64','55-64'),acat:='55-64']
D10M[Age%in%c('65-74','65-69','70-74','75-79','80-84',
              '85-89','90-94','95+','85+','65+','75+'),
     acat:='65+']
D10M[Age=='Unknown',acat:=NA]
D10M <- D10M[!is.na(acat)]

VR <- D10M[,.(deaths=sum(deaths)),by=.(iso3,year=Year,sex,acat)]
VR <- merge(VR,VRtot[,.(iso3,year=Year,totd)],by=c('iso3','year'),all=TRUE)

## save
fn <- wdir + 'VR.Rdata'
save(VR,file=fn)

