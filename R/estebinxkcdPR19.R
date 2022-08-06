## EQUATIONS FOR CALCULATING RECENT INCIDENCE IN 2019

## === incidence
## X
rIXF[,] <- pv[i] * rXLR[i,j] #fast
rIXS[,] <- eps * rXLL[i,j]   #slow
rIXR[,] <- rel * rXR[i,j]    #relapse
## H
rIHF[,1] <- irrhm[i] * pv[i] * rHLR[i,1] #fast
rIHS[,1] <- irrhm[i] * eps * rHLL[i,1]   #slow
rIHR[,1] <- irrhm[i] * rel * rHR[i,1]    #relapse
rIHF[,2] <- irrhf[i] * pv[i] * rHLR[i,2] #fast
rIHS[,2] <- irrhf[i] * eps * rHLL[i,2]   #slow
rIHR[,2] <- irrhf[i] * rel * rHR[i,2]    #relapse
## A
rIAF[,1] <- irram[i] * pv[i] * rALR[i,1] #fast
rIAS[,1] <- irram[i] * eps * rALL[i,1]   #slow
rIAR[,1] <- irram[i] * rel * rAR[i,1]    #relapse
rIAF[,2] <- irraf[i] * pv[i] * rALR[i,2] #fast
rIAS[,2] <- irraf[i] * eps * rALL[i,2]   #slow
rIAR[,2] <- irraf[i] * rel * rAR[i,2]    #relapse

## all incidence
rIX[,] <- rIXF[i,j] + rIXS[i,j] + rIXR[i,j]
rIH[,] <- rIHF[i,j] + rIHS[i,j] + rIHR[i,j]
rIA[,] <- rIAF[i,j] + rIAS[i,j] + rIAR[i,j]

## ======= DYNAMICS ===========

## top hat function to multiply all inflows
tophat <- (if (t>=2017 && t <= 2019) 1 else 0)


## --- ODEs

## &c: demography/TB/HIV&ART
## -- Mtb INFECTED RECENT
deriv(rXLR[,]) <- rageXLR[i,j] - (muxb[i,j]+r)*rXLR[i,j] +
  (ariv[i]*XU[i,j] + v*ariv[i]*(XLL[i,j]+XR[i,j])) * tophat -
  arig*rXLR[i,j] - rIXF[i,j] -
  hpc[i,j]*rXLR[i,j]
deriv(rHLR[,]) <- rageHLR[i,j] - (muxb[i,j]+r+muhb[i,j])*rHLR[i,j] +
  (ariv[i]*HU[i,j] + v*ariv[i]*(HLL[i,j]+HR[i,j])) * tophat -
  arig*rHLR[i,j] - rIHF[i,j] +
  hpc[i,j]*rXLR[i,j] - apc[i,j]*rHLR[i,j]
deriv(rALR[,]) <- rageALR[i,j] - (muxb[i,j]+r+muab[i,j])*rALR[i,j] +
  (ariv[i]*AU[i,j] + v*ariv[i]*(ALL[i,j]+AR[i,j])) * tophat -
  arig*rALR[i,j] - rIAF[i,j] +
  apc[i,j]*rHLR[i,j]

## -- Mtb INFECTED DISTANT
deriv(rXLL[,]) <- rageXLL[i,j] - (muxb[i,j]+r)*rXLL[i,j] +
  arig*rXLR[i,j] - v*ariv[i]*rXLL[i,j]+
  ((1-cfrn)*rDX[i,j]/dx[i]) -
  rIXS[i,j] -
  hpc[i,j]*rXLL[i,j]
deriv(rHLL[,]) <- rageHLL[i,j] - (muxb[i,j]+r+muhb[i,j])*rHLL[i,j] +
  arig*rHLR[i,j] - v*ariv[i]*rHLL[i,j]+
  ((1-cfrpn)*rDH[i,j]/drnH) -
  rIHS[i,j] +
  hpc[i,j]*rXLL[i,j] - apc[i,j]*rHLL[i,j]
deriv(rALL[,]) <- rageALL[i,j] - (muxb[i,j]+r+muab[i,j])*rALL[i,j] +
  arig*rALR[i,j] -v*ariv[i]*rALL[i,j]+
  ((1-cfrpp)*rDA[i,j]/drnA) -
  rIAS[i,j] +
  apc[i,j]*rHLL[i,j]

## -- RECOVERED
deriv(rXR[,]) <- rageXR[i,j] - (r+muxb[i,j])*rXR[i,j] +
  2*(1-txf)*rTX[i,j] -
  rIXR[i,j] - v*ariv[i]*rXR[i,j] -
  hpc[i,j]*rXR[i,j]
deriv(rHR[,]) <- rageHR[i,j] - (r+muxb[i,j]+muhb[i,j])*rHR[i,j] +
  2*(1-txf)*rTH[i,j] -
  rIHR[i,j] - v*ariv[i]*rHR[i,j] +
  hpc[i,j]*rXR[i,j] - apc[i,j]*rHR[i,j]
deriv(rAR[,]) <- rageAR[i,j] - (r+muxb[i,j]+muab[i,j])*rAR[i,j] +
  2*(1-txf)*rTA[i,j] -
  rIAR[i,j] - v*ariv[i]*rAR[i,j] +
  apc[i,j]*rHR[i,j]

## -- TB DISEASE
deriv(rDX[,]) <- rageDX[i,j] - (muxb[i,j]+r)*rDX[i,j] - #demo
  hpc[i,j]*rDX[i,j] +                                       #HIV/ART
  rIX[i,j] - (1+cor[i]*CDR/(1-CDR)) * rDX[i,j] / dx[i]       #inc/end
deriv(rDH[,]) <- rageDH[i,j] - (muxb[i,j]+r+muhb[i,j])*rDH[i,j]+ #demo
  hpc[i,j]*rTX[i,j] - apc[i,j]*rDH[i,j]+               #HIV/ART
  rIH[i,j] - rDH[i,j] / (drnH * (1-CDR))               #inc/end
deriv(rDA[,]) <- rageDA[i,j] - (muxb[i,j]+r+muab[i,j])*rDA[i,j]+ #demo
  apc[i,j]*rDH[i,j]+                                             #HIV/ART
  rIA[i,j] - rDA[i,j] / (drnA * (1-CDR))                         #inc/end

## -- TB TREATMENT (6 mo durn)
deriv(rTX[,]) <- rageTX[i,j] - (muxb[i,j]+r)*rTX[i,j] - #demo
  hpc[i,j]*rTX[i,j] +                                       #HIV/ART
  CDR * rDX[i,j] * cor[i] / (dx[i] * (1-CDR))- ## rNoX[i,j]
  2 * rTX[i,j]                                   #tx S/F
deriv(rTH[,]) <- rageTH[i,j] - (muxb[i,j]+r+muhb[i,j])*rTH[i,j]+ #demo
  hpc[i,j]*rTX[i,j] - apc[i,j]*rTH[i,j]+                #HIV/ART
  CDR * rDH[i,j] * cor[i] / (drnH * (1-CDR))-## rNoH[i,j]
  2 * rTH[i,j]                                   #tx S/F
deriv(rTA[,]) <- rageTA[i,j] - (muxb[i,j]+r+muab[i,j])*rTA[i,j]+ #demo
  apc[i,j]*rTH[i,j]+                                               #HIV/ART
  CDR * rDA[i,j] * cor[i] / (drnA * (1-CDR))-## rNoA[i,j]
  2 * rTA[i,j]                                   #tx S/F


## === SIMPLER AGEING
rageXLL[1,1:2] <- 0
rageXLL[2:nage,1:2] <- r*rXLL[i-1,j]
rageHLL[1,1:2] <- 0
rageHLL[2:nage,1:2] <- r*rHLL[i-1,j]
rageALL[1,1:2] <- 0
rageALL[2:nage,1:2] <- r*rALL[i-1,j]

rageXLR[1,1:2] <- 0
rageXLR[2:nage,1:2] <- r*rXLR[i-1,j]
rageHLR[1,1:2] <- 0
rageHLR[2:nage,1:2] <- r*rHLR[i-1,j]
rageALR[1,1:2] <- 0
rageALR[2:nage,1:2] <- r*rALR[i-1,j]

rageXR[1,1:2] <- 0
rageXR[2:nage,1:2] <- r*rXR[i-1,j]
rageHR[1,1:2] <- 0
rageHR[2:nage,1:2] <- r*rHR[i-1,j]
rageAR[1,1:2] <- 0
rageAR[2:nage,1:2] <- r*rAR[i-1,j]

rageDX[1,1:2] <- 0
rageDX[2:nage,1:2] <- r*rDX[i-1,j]
rageDH[1,1:2] <- 0
rageDH[2:nage,1:2] <- r*rDH[i-1,j]
rageDA[1,1:2] <- 0
rageDA[2:nage,1:2] <- r*rDA[i-1,j]

rageTX[1,1:2] <- 0
rageTX[2:nage,1:2] <- r*rTX[i-1,j]
rageTH[1,1:2] <- 0
rageTH[2:nage,1:2] <- r*rTH[i-1,j]
rageTA[1,1:2] <- 0
rageTA[2:nage,1:2] <- r*rTA[i-1,j]


## === INITIAL STATES
initial(rXLR[,]) <- 0
initial(rHLR[,]) <- 0
initial(rALR[,]) <- 0

initial(rXLL[,]) <- 0
initial(rHLL[,]) <- 0
initial(rALL[,]) <- 0

initial(rXR[,]) <- 0
initial(rHR[,]) <- 0
initial(rAR[,]) <- 0

initial(rDX[,]) <- 0
initial(rDH[,]) <- 0
initial(rDA[,]) <- 0

initial(rTX[,]) <- 0
initial(rTH[,]) <- 0
initial(rTA[,]) <- 0


## === DIMENSIONS
dim(rXLR) <- c(nage,2)
dim(rHLR) <- c(nage,2)
dim(rALR) <- c(nage,2)

dim(rXLL) <- c(nage,2)
dim(rHLL) <- c(nage,2)
dim(rALL) <- c(nage,2)

dim(rXR) <- c(nage,2)
dim(rHR) <- c(nage,2)
dim(rAR) <- c(nage,2)

dim(rDX) <- c(nage,2)
dim(rDH) <- c(nage,2)
dim(rDA) <- c(nage,2)

dim(rTX) <- c(nage,2)
dim(rTH) <- c(nage,2)
dim(rTA) <- c(nage,2)

dim(rIXF) <- c(nage,2)
dim(rIXS) <- c(nage,2)
dim(rIXR) <- c(nage,2)

dim(rIHF) <- c(nage,2)
dim(rIHS) <- c(nage,2)
dim(rIHR) <- c(nage,2)

dim(rIAF) <- c(nage,2)
dim(rIAS) <- c(nage,2)
dim(rIAR) <- c(nage,2)

dim(rIX) <- c(nage,2)
dim(rIH) <- c(nage,2)
dim(rIA) <- c(nage,2)

output(rIX) <- TRUE
output(rIH) <- TRUE
output(rIA) <- TRUE


dim(rageXLL) <- c(nage,2)
dim(rageHLL) <- c(nage,2)
dim(rageALL) <- c(nage,2)
dim(rageXLR) <- c(nage,2)
dim(rageHLR) <- c(nage,2)
dim(rageALR) <- c(nage,2)
dim(rageXR) <- c(nage,2)
dim(rageHR) <- c(nage,2)
dim(rageAR) <- c(nage,2)

dim(rageDX) <- c(nage,2)
dim(rageDH) <- c(nage,2)
dim(rageDA) <- c(nage,2)
dim(rageTX) <- c(nage,2)
dim(rageTH) <- c(nage,2)
dim(rageTA) <- c(nage,2)

