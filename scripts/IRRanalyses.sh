#! /bin/bash
## ETH KEN LSO MOZ MWI NGA
# SWZ TZA UGA ZAF ZMB ZWE

# batch 1
R --slave --vanilla --args < IRRanalysis.R SWZ & R --slave --vanilla --args < IRRanalysis.R TZA & R --slave --vanilla --args < IRRanalysis.R UGA & R --slave --vanilla --args < IRRanalysis.R ZAF & R --slave --vanilla --args < IRRanalysis.R ZMB &R --slave --vanilla --args < IRRanalysis.R ZWE

# batch 2
R --slave --vanilla --args < IRRanalysis.R ETH & R --slave --vanilla --args < IRRanalysis.R KEN & R --slave --vanilla --args < IRRanalysis.R LSO & R --slave --vanilla --args < IRRanalysis.R MOZ & R --slave --vanilla --args < IRRanalysis.R MWI &R --slave --vanilla --args < IRRanalysis.R NGA


# joining
# R --slave --vanilla --args < IRRanalysis.R ALL



