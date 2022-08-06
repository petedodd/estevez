library(here)
library(ggplot2)
library(ggthemes)
library(scales)
library(data.table)
library(ggdist)

## === load data
load(here('fitdata/NT.Rdata'))
load(here('fitdata/NAS.Rdata'))
load(here('fitdata/E.Rdata'))
load(here('fitdata/P.Rdata'))
load(here('fitdata/VR.Rdata'))


## === utilities
## graph utilities
absci_10 <- function(x) {
  x <- abs(x)
  parse(text=gsub("e", "%*%10^", scales::scientific_format()(x)))
}

absspace <- function(x,...) {
  format(abs(x), ..., big.mark=" ",scientific = FALSE, trim = TRUE)
}
rot45 <- theme(axis.text.x = element_text(angle = 45, hjust = 1))
