## CD4/IRR heatmap plots
library(here)
library(glue)
library(data.table)
library(ggplot2)
library(viridis)

## for name data
load(here('indata/XD.Rdata'))           #data for AIM

## ranges for parameters to consider
## hrz <- seq(.1,.8,by=.1)
hrz <- c(0.1,0.3,0.5,0.7)
alph <- seq(.25,.45,by=.02)

## compute as array as in dynamics
IRR <- array(1,dim=c(length(alph),length(hrz),
                     length(XD[['nmhiv']]),length(XD[['nmart']])),
             dimnames=list(alpha=paste0(alph),HR=paste0(hrz),
                           'CD4'=XD[['nmhiv']],'ART'=XD[['nmart']]))

IRR[,,">=500",] <- exp((1e3-750)*alph/100)
IRR[,,"350-499",] <- exp((1e3-(350+500)/2)*alph/100)
IRR[,,"250-349",] <- exp((1e3-(250+350)/2)*alph/100)
IRR[,,"200-249",] <- exp((1e3-(200+250)/2)*alph/100)
IRR[,,"100-199",] <- exp((1e3-(100+200)/2)*alph/100)
IRR[,,"50-99",] <- exp((1e3-(50+100)/2)*alph/100)
IRR[,,"<50",] <- exp((1e3-(50)/2)*alph/100)
for(j in 1:length(hrz)){
  HR <- hrz[j]
  IRR[,j,,"[0,6)m"] <- IRR[,j,,"[0,6)m"]*(1-(1-HR)*.25)
  IRR[,j,,"[7,12)m"] <- IRR[,j,,"[7,12)m"]*(1-(1-HR)*.75)
  IRR[,j,,"[12,Inf)m"] <- IRR[,j,,"[12,Inf)m"]*HR
}

## for different alpha (rows)
## different duration on ART (columns)
## panels of CD4 (x-axis) vs rho (y-axis)
IRR <- as.data.table(IRR)
names(IRR)[5] <- 'IRR'
IRR$alpha <- as.numeric(IRR$alpha)
IRR$CD4 <- factor(IRR$CD4,levels=XD[['nmhiv']],ordered = TRUE)
IRR$ART <- factor(IRR$ART,levels=XD[['nmart']],ordered = TRUE)
IRR$HR <- factor(IRR$HR,levels=hrz,ordered = TRUE)

GP <- ggplot(IRR[CD4!='hiv-ve'],
       aes(CD4,alpha,fill=IRR))+
  geom_tile()+
  scale_fill_viridis(option = 'plasma')+
  facet_grid(HR~ART)+
  ggtitle('Time on ART in months (m)')+
  xlab('CD4 cell count category')+
  ylab('Parameter for IRR CD4-dependence (\u03c1)')+
  scale_y_continuous(sec.axis = sec_axis(~.,name='Hazard ratio for TB from established ART'))+
  theme(legend.position='left',
        plot.title = element_text(hjust = 0.5,size=11),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)
        )

ggsave(GP,file=here('plots/resfigs/CD4heatmap.pdf'),device = cairo_pdf,w=11,h=14)
ggsave(GP,file=here('plots/resfigs/CD4heatmap.png'),w=14,h=11)
