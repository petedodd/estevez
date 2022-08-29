## to do MR data afterwards
## read in arguments
## R --slave --vanilla --args < cornerplots.R ETH
args <- commandArgs(trailingOnly = TRUE)
cn <- as.character(args[1])
cat('--- ',cn,' ---\n')


useunicode <- TRUE

library(here)
library(glue)
library(data.table)
library(ggplot2)
## library(ggthemes)
## library(scales)
library(ggpubr)
library(GGally)

## see dataplots.R
load(file=here('data/parmkey.Rdata'))
der <- c(1,11,2, #tx
         4,20,5,3,16, #nh
         6:7,       #hiv
         8,17:19,   #cdr
         9:10,      #dur
         12:15      #cfr
         )


## ============== inference plot utils ==============

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

MakeInfPlots <- function(cniso){
  ## load data
  sminf <- glue(here('data/sampe_{cniso}.Rdata'))
  load(sminf)
  if(useunicode){
    setnames(sampe,old=parmkey$model,new=parmkey$writeup)
    setcolorder(sampe,der)
  }
  ## pairs plot
  PP <- ppairs(sampe)
  ppfn <- glue(here('plots/corner/C_{cniso}.pdf'))

  ## save
  ggsave(PP,file=ppfn,w=10,h=10,device = cairo_pdf)
}


## ============== inference outputs ==============

## make folder if absent
fn <- glue(here('plots/corner'))
if(!file.exists(fn)) dir.create(fn)

## inference
MakeInfPlots(cn)


