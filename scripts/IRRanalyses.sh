bash/#! /bin
## ETH KEN LSO MOZ MWI NGA
# SWZ TZA UGA ZAF ZMB ZWE

# NOTE edit to batch as per number of cores available
# batch 1
R --slave --vanilla --args < ../R/IRRanalysis.R SWZ & R --slave --vanilla --args < ../R/IRRanalysis.R TZA & R --slave --vanilla --args < ../R/IRRanalysis.R UGA & R --slave --vanilla --args < ../R/IRRanalysis.R ZAF & R --slave --vanilla --args < ../R/IRRanalysis.R ZMB &R --slave --vanilla --args < ../R/IRRanalysis.R ZWE & R --slave --vanilla --args < ../R/IRRanalysis.R ETH & R --slave --vanilla --args < ../R/IRRanalysis.R KEN & R --slave --vanilla --args < ../R/IRRanalysis.R LSO & R --slave --vanilla --args < ../R/IRRanalysis.R MOZ & R --slave --vanilla --args < ../R/IRRanalysis.R MWI &R --slave --vanilla --args < ../R/IRRanalysis.R NGA


# joining
# R --slave --vanilla --args < IRRanalysis.R ALL



