## preparing relevant mixing data
library(here)
library(data.table)

## read data
f <- '/home/pjd/Downloads/contacts/synthetic_contacts_2020.csv'
CD <- fread(f)
## --- mixing data
## contact data
CD <- CD[location_contact=="all" & setting == 'overall',
         .(iso3=iso3c,age_contactor,age_contactee=age_cotactee,
           ctx=mean_number_of_contacts)]
CD[age_contactor=='75+',age_contactor:='75-79']
CD[age_contactee=='75+',age_contactee:='75-79']
tmp1 <- CD[age_contactee=='75-79']
tmp2 <- CD[age_contactor=='75-79']
tmp3 <- CD[age_contactor=='75-79' & age_contactee=='75-79']
tmp1[,age_contactee:='80+']
tmp2[,age_contactor:='80+']
tmp3[,c('age_contactor','age_contactee'):='80+']
CD <- rbindlist(list(CD,tmp1,tmp2,tmp3))
CD[,age_contactee:=gsub(' to ','-',age_contactee)]
CD[,age_contactor:=gsub(' to ','-',age_contactor)]
## ordered
lvls <- paste(seq(from=0,to=75,by=5),seq(from=4,to=79,by=5),sep='-')
lvls <- c(lvls,'80+')
CD[,age_contactee:=factor(age_contactee,levels=lvls,ordered=TRUE)]
CD[,age_contactor:=factor(age_contactor,levels=lvls,ordered=TRUE)]

## countries of interest
cz <- c('ETH','KEN','LSO','MOZ','MWI','NGA','SWZ',
        'TZA','UGA','ZAF','ZMB','ZWE')

CML <- list()
for(cn in cz){
  tmp <- dcast(CD[iso3==cn],
               age_contactor ~ age_contactee,
               value.var = 'ctx')
  tmp <- as.matrix(tmp[,..lvls])
  rownames(tmp) <- lvls
  ## tots <- rowSums(tmp)
  CML[[cn]] <- tmp/mean(tmp)
}


CML[['ZWE']]
CML[['ETH']]

save(CML,file=here('indata/CML.Rdata'))
