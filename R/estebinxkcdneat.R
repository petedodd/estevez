## === INPUT DATA & PARAMETERS
nage <- 17                            #no age cats 80+ & 5 year
r <- 0.2                              #ageing rate
tscale <- user()
popinitM[] <- user()                  #initial state
popinitF[] <- user()
BF[] <- user()                        #births
BM[] <- user()
CDRdata[] <- user()                      #CDR targets
hitD[] <- user()                        #HIV/ART targets
hatD[] <- user()
rra[,] <- user()
ttp[] <- user()                       #data times
ttq[] <- user()                       #data times
## artbyage[] <- user()                  #art starts by age (per cap)
hsrD[] <- user()
hivbyage[,] <- user()                 #HIV by age/sex data
am[] <- user()                        #age midpoints
MM[,] <- user()                        #mixing matrix

## --- TB parameters
drnX <- user() #durn HIV-ve TB
drnH <- user() #durn HIV+ve/ART-ve TB
drnA <- user() #durn HIV+ve/ART+ve TB
ari0 <- user() #initial state parm
pp <- user()      #recent progression _rate_
arig <- user() #aging _rate_ from fast latent R -> L
eps <- user()     #slow progression _rate_
v <- user() #protection from LTBI - HR
bet <- user() #transmission parm
rel <- user() #from recovered
output(CDR) <- TRUE
## CFRs
txf <- user(0.05) #CFR on ATT
cfrn <- user(0.5) #CFR HIV-ve no ATT
cfrpn <- user(0.95) #CFR HIV+ve/ART-ve no ATT
cfrpp <- user(0.6)  #CFR HIV+ve/ART+ve no ATT
OR04 <- user()                         #OR CDR u5
OR514 <- user()                        #OR CDR 514
pp04 <- user()                         #u5 progn

## vector parameters for paediatric sector
pv[1] <- pp04
pv[2:nage] <- pp
dx[1:3] <- 0.5 #6 month duration
dx[4:nage] <- drnX
cor[1] <- OR04 #CDR ORs
cor[2:3] <- OR514
cor[4:nage] <- 1
relinf[1:3] <- 0 #relative infectiousness
relinf[4:nage] <- 1

## dealing with sexes more easily
BB[1,1] <- bzm
BB[1,2] <- bzf
BB[2:nage,] <- 0


## --- inpnut data interpolation
## mortality
muxfD[,] <- user()
muxmD[,] <- user()
muhfD[,] <- user()
muhmD[,] <- user()
muafD[,] <- user()
muamD[,] <- user()
irrhmD[,] <- user()
irrhfD[,] <- user()
irramD[,] <- user()
irrafD[,] <- user()
muxf2[] <- interpolate(ttq,muxfD,'linear') #ch
muxm2[] <- interpolate(ttq,muxmD,'linear') #ch
muhf[] <- interpolate(ttp,muhfD,'linear')
muhm[] <- interpolate(ttp,muhmD,'linear')
muaf[] <- interpolate(ttp,muafD,'linear')
muam[] <- interpolate(ttp,muamD,'linear')
irrhm[] <- interpolate(ttq,irrhmD,'linear')
irrhf[] <- interpolate(ttq,irrhfD,'linear')
irram[] <- interpolate(ttq,irramD,'linear')
irraf[] <- interpolate(ttq,irrafD,'linear')
## CDR
CDR <- interpolate(ttq,CDRdata,'linear')
## incidence & IRR
hit <- interpolate(ttq,hitD,'linear')
## hit2 <- interpolate(ttq,hitD2,'linear')
hsr <- interpolate(ttq,hsrD,'linear')
hat <- interpolate(ttq,hatD,'linear')
## births
bzf <- interpolate(ttq,BF,'linear')
bzm <- interpolate(ttq,BM,'linear')


## --- state key ---
## S = U
## E = LR
## L = LL
## I = D
## T = T
## R = R
## -----------------

## === SIMPLER AGEING
ageXU[1,1:2] <- 0
ageXU[2:nage,1:2] <- r*XU[i-1,j]
ageHU[1,1:2] <- 0
ageHU[2:nage,1:2] <- r*HU[i-1,j]
ageAU[1,1:2] <- 0
ageAU[2:nage,1:2] <- r*AU[i-1,j]

ageXLL[1,1:2] <- 0
ageXLL[2:nage,1:2] <- r*XLL[i-1,j]
ageHLL[1,1:2] <- 0
ageHLL[2:nage,1:2] <- r*HLL[i-1,j]
ageALL[1,1:2] <- 0
ageALL[2:nage,1:2] <- r*ALL[i-1,j]

ageXLR[1,1:2] <- 0
ageXLR[2:nage,1:2] <- r*XLR[i-1,j]
ageHLR[1,1:2] <- 0
ageHLR[2:nage,1:2] <- r*HLR[i-1,j]
ageALR[1,1:2] <- 0
ageALR[2:nage,1:2] <- r*ALR[i-1,j]

ageXR[1,1:2] <- 0
ageXR[2:nage,1:2] <- r*XR[i-1,j]
ageHR[1,1:2] <- 0
ageHR[2:nage,1:2] <- r*HR[i-1,j]
ageAR[1,1:2] <- 0
ageAR[2:nage,1:2] <- r*AR[i-1,j]

ageDX[1,1:2] <- 0
ageDX[2:nage,1:2] <- r*DX[i-1,j]
ageDH[1,1:2] <- 0
ageDH[2:nage,1:2] <- r*DH[i-1,j]
ageDA[1,1:2] <- 0
ageDA[2:nage,1:2] <- r*DA[i-1,j]

ageTX[1,1:2] <- 0
ageTX[2:nage,1:2] <- r*TX[i-1,j]
ageTH[1,1:2] <- 0
ageTH[2:nage,1:2] <- r*TH[i-1,j]
ageTA[1,1:2] <- 0
ageTA[2:nage,1:2] <- r*TA[i-1,j]



## === INITIAL STATE
## initial prevalence (utility in below)
friL[] <- (1-exp(-ari0*am[i]))
initF <- (rel + eps + v * ari0) / arig #0.2
dinit[] <- 2 * eps * drnX * (1-CDRdata[1])




## E:L balance as (mu + eps + v*ari) / arig - E equation
initial(XU[1:nage,1]) <- popinitM[i]*(1-friL[i])
initial(XU[1:nage,2]) <- popinitF[i]*(1-friL[i])
initial(XLL[1:nage,1]) <- popinitM[i]* friL[i] * (1.0 - initF - 4*dinit[i]/3)
initial(XLL[1:nage,2]) <- popinitF[i]* friL[i] * (1.0 - initF - 4*dinit[i]/3)
initial(XLR[1:nage,1]) <- popinitM[i]* friL[i] * initF/2 ## (0.02 + eps + v * ari0) / (2*arig)
initial(XLR[1:nage,2]) <- popinitF[i]* friL[i] * initF/2## (0.02 + eps + v * ari0) / (2*arig)
initial(XR[1:nage,1]) <- popinitM[i]* friL[i] * initF/2## (0.02 + eps + v * ari0) / (2*arig)
initial(XR[1:nage,2]) <- popinitF[i]* friL[i] * initF/2## (0.02 + eps + v * ari0) / (2*arig)


initial(HU[1:nage,1:2]) <- 1e-10
initial(AU[1:nage,1:2]) <- 1e-10
initial(HLR[1:nage,1:2]) <- 1e-10
initial(ALR[1:nage,1:2]) <- 1e-10
initial(HLL[1:nage,1:2]) <- 1e-10
initial(ALL[1:nage,1:2]) <- 1e-10
initial(HR[1:nage,1]) <- 1e-10
initial(HR[1:nage,2]) <- 1e-10
initial(AR[1:nage,1]) <- 1e-10
initial(AR[1:nage,2]) <- 1e-10

## tb disease
initial(DX[,1]) <- dinit[i] * popinitM[i] * friL[i]
initial(DX[,2]) <- dinit[i] * popinitF[i] * friL[i]
initial(DH[1:nage,1:2]) <- 1e-10
initial(DA[1:nage,1:2]) <- 1e-10

## tb treatment
initial(TX[,1]) <- dinit[i] * popinitM[i] * friL[i]/3
initial(TX[,1]) <- dinit[i] * popinitF[i] * friL[i]/3
initial(TH[1:nage,1:2]) <- 1e-10
initial(TA[1:nage,1:2]) <- 1e-10


## === DYNAMICS

## --- populations
Na[,] <- AU[i,j] + ALR[i,j] + ALL[i,j] + AR[i,j] + DA[i,j] + TA[i,j]
Nh[,] <- HU[i,j] + HLR[i,j] + HLL[i,j] + HR[i,j] + DH[i,j] + TH[i,j]
Nx[,] <- XU[i,j] + XLR[i,j] + XLL[i,j] + XR[i,j] + DX[i,j] + TX[i,j]
N[,] <- Nx[i,j] + Nh[i,j] + Na[i,j]
poptot <- sum(N)
pop15plus <- sum(N[4:nage,1:2])
pop1549 <- sum(N[4:10,1:2])
poprec <- sum(XR) + sum(HR) + sum(AR) #recovered
poph1549 <- sum(Nh[4:10,1:2]) + sum(Na[4:10,1:2]) #HIV 15-49
poph <- sum(Nh[4:nage,1:2]) + sum(Na[4:nage,1:2]) #HIV total
popx <- sum(Nx[4:nage,1:2])  #HIV-eligible
pophae <- sum(Nh[4:nage,1:2]) #ART-eligible
popart <- sum(Na[4:nage,1:2])


hivi <- (1)*max((hit*pop1549 - poph1549),0)/tscale #pursuit
## (pa*H)' = pa'*H + pa*H' = pa'*H + pa*ph'*H
arti <- max(hat*poph - popart,0)/tscale #pursuit
hivbyageN[,] <- hivbyage[i,j] * N[i,j]
rraN[,] <- rra[i,j] * (HU[i,j] + HLR[i,j] + HLL[i,j] + HR[i,j])
## rraN[,] <- H[i,j]                     #test with even assignment
raz[,] <- arti * rraN[i,j] / (sum(rraN) + 1e-15) #ART initiations
hfz[,1] <- hivi * (1-hsr) * hivbyageN[i,1]/ (sum(hivbyageN[1:nage,1])+1e-15) #HIV infections
hfz[,2] <- hivi * (hsr) * hivbyageN[i,2] / (sum(hivbyageN[1:nage,2])+1e-15)
hpc[,] <- 1*hfz[i,j]/popx #per capita HIV infection
apc[,] <- 1*raz[i,j]/(pophae+1e-15) #per capita ART starts


## --- processes

## calculated remaining mortality rate
##  muT * N = muU * U + muHt * H + muAt * A
##  muT * N = muU * N + muH * H + muA * A
##  muU = (muT * N - muH * H - muA * A) / N
muxm[1:nage] <- ((muxm2[i]-r)*N[i,1]-muhm[i]*Nh[i,1]-muam[i]*Na[i,1]
  - muX[i,1] - muH[i,1] - muA[i,1])/N[i,1]
muxf[1:nage] <- ((muxf2[i]-r)*N[i,2]-muhf[i]*Nh[i,2]-muaf[i]*Na[i,2]
  - muX[i,2] - muH[i,2] - muA[i,2])/N[i,2]


## new version
##  muT * N = muU * N + muH * H + muA * A + muX + muH + muA

## 2-sex objects
muxb[,1] <- muxm[i]
muxb[,2] <- muxf[i]
muhb[1:nage,1] <- muhm[i]
muhb[1:nage,2] <- muhf[i]
muab[1:nage,1] <- muam[i]
muab[1:nage,2] <- muaf[i]

## TB outcomes
## a = 1 / d = rate out with no detection
## ATT : no-ATT = CDR : (1-CDR)
## total rate w/detection = a + b; CDR=b/(a+b) -> b = CDR/(1-CDR)* a -> a+b = 1/(1-CDR) / d
## 1 / ((1-CDR)*d) = rate out with detection
## notes ~ b * D = (CDR/(1-CDR)) * D / d
## -> non-notes ~  D / d
## deaths = cfr * D / d; self-cure = (1-cfr) * D / d
## w OR, total rate = (1+OR*CDR/(1-CDR))/d

## TB deaths
muX[,] <- 2 * txf * TX[i,j] + cfrn  * DX[i,j] / dx[i]
muH[,] <- 2 * txf * TH[i,j] + cfrpn * DH[i,j] / drnH
muA[,] <- 2 * txf * TA[i,j] + cfrpp * DA[i,j] / drnA
mortX <- sum(muX)
mortH <- sum(muH)
mortA <- sum(muA)
mort <- mortX + mortH + mortA
mortFH <- sum(muH[4:nage,1:2]) / (sum(muX[4:nage,1:2]) + sum(muH[4:nage,1:2]) + sum(muA[4:nage,1:2]))

## notifications
NoX[,] <- CDR * DX[i,j] * cor[i] / (dx[i] * (1-CDR))
NoH[,] <- CDR * DH[i,j] * cor[i] / (drnH * (1-CDR))
NoA[,] <- CDR * DA[i,j] * cor[i] / (drnA * (1-CDR))
NoN[,] <- NoX[i,j] + NoH[i,j] + NoA[i,j]
Ntot <- sum(NoN)
NtotH <- sum(NoH) + sum(NoA)
NtotA <- sum(NoA)

## incidence
## X
IXF[,] <- pv[i] * XLR[i,j] #fast
IXS[,] <- eps * XLL[i,j]   #slow
IXR[,] <- rel * XR[i,j]    #relapse
## H
IHF[,1] <- irrhm[i] * pv[i] * HLR[i,j] #fast
IHS[,1] <- irrhm[i] * eps * HLL[i,j]   #slow
IHR[,1] <- irrhm[i] * rel * HR[i,j]    #relapse
IHF[,2] <- irrhf[i] * pv[i] * HLR[i,j] #fast
IHS[,2] <- irrhf[i] * eps * HLL[i,j]   #slow
IHR[,2] <- irrhf[i] * rel * HR[i,j]    #relapse
## A
IAF[,1] <- irram[i] * pv[i] * ALR[i,j] #fast
IAS[,1] <- irram[i] * eps * ALL[i,j]   #slow
IAR[,1] <- irram[i] * rel * AR[i,j]    #relapse
IAF[,2] <- irraf[i] * pv[i] * ALR[i,j] #fast
IAS[,2] <- irraf[i] * eps * ALL[i,j]   #slow
IAR[,2] <- irraf[i] * rel * AR[i,j]    #relapse

## all incidence
IX[,] <- IXF[i,j] + IXS[i,j] + IXR[i,j]
IH[,] <- IHF[i,j] + IHS[i,j] + IHR[i,j]
IA[,] <- IAF[i,j] + IAS[i,j] + IAR[i,j]
Inn <- sum(IX)
Ipn <- sum(IH)
Ipp <- sum(IA)
Itot <- Inn + Ipn + Ipp
IreactX <- eps * (sum(XLL[,])) / (sum(IX[,])+1e-15) #fraction of X-inc from ractivation
IrelX <- rel * (sum(XR[,])) / (sum(IX[,])+1e-15) #fraction of X-inc from ractivation

## prevalence
prevtot15plus <- sum(DX[4:nage,1:2]) + sum(DA[4:nage,1:2]) + sum(DH[4:nage,1:2])
prevTtot15plus <- sum(TX[4:nage,1:2]) + sum(TA[4:nage,1:2]) + sum(TH[4:nage,1:2])
prev1524 <- (sum(DX[4:5,1:2]) + sum(DA[4:5,1:2]) + sum(DH[4:5,1:2]))/(sum(N[4:5,1:2])+1e-10)
prev2534 <- (sum(DX[6:7,1:2]) + sum(DA[6:7,1:2]) + sum(DH[6:7,1:2]))/(sum(N[6:7,1:2])+1e-10)
prev3544 <- (sum(DX[8:9,1:2]) + sum(DA[8:9,1:2]) + sum(DH[8:9,1:2]))/(sum(N[8:9,1:2])+1e-10)
prev4554 <- (sum(DX[10:11,1:2]) + sum(DA[10:11,1:2]) + sum(DH[10:11,1:2]))/(sum(N[10:11,1:2])+1e-10)
prev5564 <- (sum(DX[12:13,1:2]) + sum(DA[12:13,1:2]) + sum(DH[12:13,1:2]))/(sum(N[12:13,1:2])+1e-10)
prev65pl <- (sum(DX[14:nage,1:2]) + sum(DA[14:nage,1:2]) + sum(DH[14:nage,1:2]))/(sum(N[14:nage,1:2])+1e-10)


## transmission rule
## ari <- bet * (sum(DX)  + sum(DH) + sum(DA)) / sum(N)
## arim[,] <- bet * (1/(17)) * (sum(DX[j,])  + sum(DH[j,]) + sum(DA[j,])) / (sum(N[j,])+1e-15) #random mixing
arim[,] <- bet * MM[i,j] * relinf[j] * (sum(DX[j,])  + sum(DH[j,]) + sum(DA[j,])) / (sum(N[j,])+1e-15)
ariv[] <- sum(arim[i,])
ariw[] <- ariv[i]*sum(N[i,])      #weighted by pop group (summed over sex) to get rate of infection in gp i
aris <- sum(ariw)/poptot #scalar summary for population = mean ari for popn
## infections from groups
arimf[,] <- arim[i,j]*sum(N[i,]) #number of infections in i from j
arif[] <- sum(arimf[,i])         #number of infections from i in total

## --- which processes to output:
## incidence
output(Itot) <- TRUE
output(Inn) <- TRUE
output(Ipn) <- TRUE
output(Ipp) <- TRUE
output(IX) <- TRUE
output(IH) <- TRUE
output(IA) <- TRUE
## notifs
output(Ntot) <- TRUE
output(NtotH) <- TRUE
output(NtotA) <- TRUE
output(NoN) <- TRUE
## prev & pop
output(prevtot15plus) <- TRUE
output(prevTtot15plus) <- TRUE
output(pop15plus) <- TRUE
output(prev1524) <- TRUE
output(prev2534) <- TRUE
output(prev3544) <- TRUE
output(prev4554) <- TRUE
output(prev5564) <- TRUE
output(prev65pl) <- TRUE


output(poprec) <- TRUE
output(poptot) <- TRUE
output(popart) <- TRUE
output(poph) <- TRUE
output(poph1549) <- TRUE
output(pop1549) <- TRUE
## output(N) <- TRUE
## mortality
output(mort) <- TRUE
output(mortH) <- TRUE
output(mortA) <- TRUE
output(mortFH) <- TRUE
output(IreactX) <- TRUE
output(IrelX) <- TRUE

## these in PR
## output(aris) <- TRUE
## output(ariv) <- TRUE
## output(arif) <- TRUE


## --- ODEs
## demography; TB; HIV/ART
## UNINFECTED
deriv(XU[,]) <- BB[i,j] + ageXU[i,j] - (muxb[i,j]+r)*XU[i,j] - #demo
  ariv[i]*XU[i,j] -                                        #TB infection
  hpc[i,j]*XU[i,j]                                         #HIV/ART
deriv(HU[,]) <- ageHU[i,j] - (muxb[i,j]+r+muhb[i,j])*HU[i,j]- #demo
  ariv[i]*HU[i,j]+                                                #TB infection
  hpc[i,j]*XU[i,j] - apc[i,j]*HU[i,j]                             #HIV/ART
deriv(AU[,]) <- ageAU[i,j] - (muxb[i,j]+r+muab[i,j])*AU[i,j]- #demo
  ariv[i]*AU[i,j]+                                                #TB infection
  apc[i,j]*HU[i,j]                                                #HIV/ART

## &c: demography/TB/HIV&ART
## -- Mtb INFECTED RECENT
deriv(XLR[,]) <- ageXLR[i,j] - (muxb[i,j]+r)*XLR[i,j] +
  ariv[i]*XU[i,j] + v*ariv[i]*(XLL[i,j]+XR[i,j]) - arig*XLR[i,j] - IXF[i,j] -
  hpc[i,j]*XLR[i,j]
deriv(HLR[,]) <- ageHLR[i,j] - (muxb[i,j]+r+muhb[i,j])*HLR[i,j] +
  ariv[i]*HU[i,j] + v*ariv[i]*(HLL[i,j]+HR[i,j]) - arig*HLR[i,j] - IHF[i,j] +
  hpc[i,j]*XLR[i,j] - apc[i,j]*HLR[i,j]
deriv(ALR[,]) <- ageALR[i,j] - (muxb[i,j]+r+muab[i,j])*ALR[i,j] +
  ariv[i]*AU[i,j] + v*ariv[i]*(ALL[i,j]+AR[i,j]) - arig*ALR[i,j] - IAF[i,j] +
  apc[i,j]*HLR[i,j]

## -- Mtb INFECTED DISTANT
deriv(XLL[,]) <- ageXLL[i,j] - (muxb[i,j]+r)*XLL[i,j] +
  arig*XLR[i,j] - v*ariv[i]*XLL[i,j] - IXS[i,j] + (1-cfrn)*DX[i,j]/dx[i] -
  hpc[i,j]*XLL[i,j]
deriv(HLL[,]) <- ageHLL[i,j] - (muxb[i,j]+r+muhb[i,j])*HLL[i,j] +
  arig*HLR[i,j] - v*ariv[i]*HLL[i,j] - IHS[i,j] + (1-cfrpn)*DH[i,j]/drnH +
  hpc[i,j]*XLL[i,j] - apc[i,j]*HLL[i,j]
deriv(ALL[,]) <- ageALL[i,j] - (muxb[i,j]+r+muab[i,j])*ALL[i,j] +
  arig*ALR[i,j] - v*ariv[i]*ALL[i,j] - IAS[i,j] + (1-cfrpp)*DA[i,j]/drnA +
  apc[i,j]*HLL[i,j]


## -- TB DISEASE
deriv(DX[,]) <- ageDX[i,j] - (muxb[i,j]+r)*DX[i,j] - #demo
  hpc[i,j]*DX[i,j] +                                       #HIV/ART
  IX[i,j] - (1+cor[i]*CDR/(1-CDR)) * DX[i,j] / dx[i]       #inc/end
deriv(DH[,]) <- ageDH[i,j] - (muxb[i,j]+r+muhb[i,j])*DH[i,j]+ #demo
  hpc[i,j]*TX[i,j] - apc[i,j]*DH[i,j]+                             #HIV/ART
  IH[i,j] - DH[i,j] / (drnH * (1-CDR))                             #inc/end
deriv(DA[,]) <- ageDA[i,j] - (muxb[i,j]+r+muab[i,j])*DA[i,j]+ #demo
  apc[i,j]*DH[i,j]+                                               #HIV/ART
  IA[i,j] - DA[i,j] / (drnA * (1-CDR))                            #inc/end

## -- TB TREATMENT (6 mo durn)
deriv(TX[,]) <- ageTX[i,j] - (muxb[i,j]+r)*TX[i,j] - #demo
  hpc[i,j]*TX[i,j] +                                       #HIV/ART
  NoX[i,j] - 2 * TX[i,j]                                   #tx S/F
deriv(TH[,]) <- ageTH[i,j] - (muxb[i,j]+r+muhb[i,j])*TH[i,j]+ #demo
  hpc[i,j]*TX[i,j] - apc[i,j]*TH[i,j]+                             #HIV/ART
  NoH[i,j] - 2 * TH[i,j]                                   #tx S/F
deriv(TA[,]) <- ageTA[i,j] - (muxb[i,j]+r+muab[i,j])*TA[i,j]+ #demo
  apc[i,j]*TH[i,j]+                                               #HIV/ART
  NoA[i,j] - 2 * TA[i,j]                                   #tx S/F


## -- RECOVERED
deriv(XR[,]) <- ageXR[i,j] - (r+muxb[i,j])*XR[i,j] +
  2*(1-txf) *TX[i,j] - IXR[i,j] - v*ariv[i]*XR[i,j] -
  hpc[i,j]*XR[i,j]
deriv(HR[,]) <- ageHR[i,j] - (r+muxb[i,j]+muhb[i,j])*HR[i,j] +
  2*(1-txf)*TH[i,j] - IHR[i,j] - v*ariv[i]*HR[i,j] +
  hpc[i,j]*XR[i,j] - apc[i,j]*HR[i,j]
deriv(AR[,]) <- ageAR[i,j] - (r+muxb[i,j]+muab[i,j])*AR[i,j] +
  2*(1-txf)*TA[i,j] - IAR[i,j] - v*ariv[i]*AR[i,j] +
  apc[i,j]*HR[i,j]



## === DIMENSIONS
dim(ttp) <- user()                    #note need this before length() use
dim(ttq) <- user()                    #note need this before length() use
lttp <- length(ttp)
lttq <- length(ttq)
dim(popinitF) <- c(nage)
dim(popinitM) <- c(nage)
dim(BF) <- lttq
dim(BM) <- lttq
dim(hitD) <- lttq
dim(CDRdata) <- lttq
dim(hatD) <- lttq
dim(hsrD) <- lttq
dim(muxmD) <- c(lttq,nage)            #ch
dim(muxfD) <- c(lttq,nage)            #ch
dim(irrhmD) <- c(lttq,nage)
dim(irrhfD) <- c(lttq,nage)
dim(irrafD) <- c(lttq,nage)
dim(irramD) <- c(lttq,nage)
dim(muhmD) <- c(lttp,nage)
dim(muhfD) <- c(lttp,nage)
dim(muamD) <- c(lttp,nage)
dim(muafD) <- c(lttp,nage)
dim(muxm) <- c(nage)
dim(muxf) <- c(nage)
dim(muxm2) <- c(nage)
dim(muxf2) <- c(nage)
dim(muhm) <- c(nage)
dim(muhf) <- c(nage)
dim(muam) <- c(nage)
dim(muaf) <- c(nage)
dim(irrhm) <- c(nage)
dim(irrhf) <- c(nage)
dim(irram) <- c(nage)
dim(irraf) <- c(nage)
dim(XU) <- c(nage,2)
dim(HU) <- c(nage,2)
dim(AU) <- c(nage,2)
dim(XLR) <- c(nage,2)
dim(HLR) <- c(nage,2)
dim(ALR) <- c(nage,2)
dim(XLL) <- c(nage,2)
dim(HLL) <- c(nage,2)
dim(ALL) <- c(nage,2)
dim(XR) <- c(nage,2)
dim(HR) <- c(nage,2)
dim(AR) <- c(nage,2)
dim(DX) <- c(nage,2)
dim(DH) <- c(nage,2)
dim(DA) <- c(nage,2)
dim(TX) <- c(nage,2)
dim(TH) <- c(nage,2)
dim(TA) <- c(nage,2)
dim(IX) <- c(nage,2)
dim(IH) <- c(nage,2)
dim(IA) <- c(nage,2)
dim(NoX) <- c(nage,2)
dim(NoH) <- c(nage,2)
dim(NoA) <- c(nage,2)
dim(muX) <- c(nage,2)
dim(muH) <- c(nage,2)
dim(muA) <- c(nage,2)
dim(N) <- c(nage,2)
dim(Nx) <- c(nage,2)
dim(Nh) <- c(nage,2)
dim(Na) <- c(nage,2)
dim(hfz) <- c(nage,2)
dim(hivbyage) <- c(nage,2)
dim(hivbyageN) <- c(nage,2)
dim(rra) <- c(nage,2)
dim(rraN) <- c(nage,2)
dim(raz) <- c(nage,2)
dim(hpc) <- c(nage,2)
dim(apc) <- c(nage,2)
dim(am) <- c(nage)
dim(ariv) <- c(nage)
dim(ariw) <- c(nage)
dim(arif) <- c(nage)
dim(arim) <- c(nage,nage)
dim(arimf) <- c(nage,nage)
dim(MM) <- c(nage,nage)
dim(pv) <- c(nage)
dim(dx) <- c(nage)
dim(cor) <- c(nage)
dim(dinit) <- c(nage)
dim(friL) <- c(nage)
## tests
dim(NoN) <- c(nage,2)
## ## combo dims
dim(IXF) <- c(nage,2)
dim(IXS) <- c(nage,2)
dim(IXR) <- c(nage,2)
dim(IHF) <- c(nage,2)
dim(IHS) <- c(nage,2)
dim(IHR) <- c(nage,2)
dim(IAF) <- c(nage,2)
dim(IAS) <- c(nage,2)
dim(IAR) <- c(nage,2)
dim(muxb) <- c(nage,2)
dim(muhb) <- c(nage,2)
dim(muab) <- c(nage,2)
dim(BB) <- c(nage,2)
## ageing dims
dim(ageXU) <- c(nage,2)
dim(ageHU) <- c(nage,2)
dim(ageAU) <- c(nage,2)
dim(ageXLL) <- c(nage,2)
dim(ageHLL) <- c(nage,2)
dim(ageALL) <- c(nage,2)
dim(ageXLR) <- c(nage,2)
dim(ageHLR) <- c(nage,2)
dim(ageALR) <- c(nage,2)
dim(ageXR) <- c(nage,2)
dim(ageHR) <- c(nage,2)
dim(ageAR) <- c(nage,2)
dim(ageDX) <- c(nage,2)
dim(ageDH) <- c(nage,2)
dim(ageDA) <- c(nage,2)
dim(ageTX) <- c(nage,2)
dim(ageTH) <- c(nage,2)
dim(ageTA) <- c(nage,2)
dim(relinf) <- c(nage)


