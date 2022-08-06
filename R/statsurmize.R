library(here)
library(data.table)


fz <- dir(path=here('scripts'),pattern = "err",full.names = TRUE) #to change

P <- W <- WT <- L <- list()
for(i in 1:length(fz)){
  ## read
  D <- scan(fz[i],what='string',quiet=TRUE)
  ## progress
  prog <- D[length(D)-2]
  prog <- eval(parse(text=prog))*1e2
  P[[i]] <- prog
  ## time to go
  togo <- D[length(D)-1]
  togo <- gsub(",","",togo)
  togo <- strsplit(togo,"<")[[1]][2]
  togot <- as.POSIXct(togo,format="%H:%M:%S")
  W[[i]] <- togot
  WT[[i]] <- togo
  ## duration per iteration
  sperit <- D[length(D)]
  sperit <- gsub("s/it]","",sperit)
  sperit <- as.numeric(sperit)
  L[[i]] <- sperit
}

P <- unlist(P)
L <- unlist(L)
WT <- unlist(WT)
W <- unlist(W)

cat('Max progress = ',max(P),'% at ', which.max(P),'\n' )
cat('Min progress = ',min(P),'% at ', which.min(P),'\n' )
cat('Slowest = ',max(L),'s/it at ', which.max(L),'\n' )
cat('Fastest = ',min(L),'s/it at ', which.min(L),'\n' )
cat('Max est wait = ',WT[which.max(W)],' at ', which.max(W),'\n' )
cat('Min est wait = ',WT[which.min(W)],' at ', which.min(W),'\n' )
