## handles arguments in same way for inference.R and multirunner.R
## R --slave --vanilla --args < inference.R ZWE
args <- commandArgs(trailingOnly = TRUE)
cniso <- as.character(args[1])
cat('cniso  = ',cniso,'...\n')

sensitivity.analysis <- ''
if(length(args)>1){
  sensitivity.analysis <- as.character(args[2])
  cat('---!! running sensitivity analysis = ',sensitivity.analysis,' !!---\n')
}
