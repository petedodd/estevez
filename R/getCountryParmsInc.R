##' A function to prepare country specific parameteres for use in demographic/AIM model
##'
##' Content to be written
##' @title Parameter preparation for AIM
##' @param iso 
##' @return a list of parameters. See Details.
##' @author Pete Dodd
getCountryAimParmsInc <- function(iso,eyear=2015){

  ## usual treatment ---------
  ans <- list()
  ans$iso3 <- iso
  ans$pp <- ans$mh <-
    ans$ma <- ans$cdinit <- array(0,
                                  dim=XD[['dz']],
                                  dimnames=list(XD[['nmsex']],
                                                XD[['nmage']],
                                                XD[['nmhiv']],
                                                XD[['nmart']]))
  for( sx in XD[['nmsex']] ){
    for(cd4 in XD[['nmhiv']][-1]){               #CD4
      for(ag in XD[['nmagg']]){
        ## print(c(ans$Region,sx,ag,cd4))
        if(cd4!='<50')
          ans$pp[sx,XD[['aglst']][[ag]],cd4,1] <- XD[[iso]]$pp[cd4,paste0(sx,ag)]
        ## mortality no ART
        ans$mh[sx,XD[['aglst']][[ag]],cd4,1] <- XD[[iso]]$muna[cd4,paste0(sx,ag)]
        if(TRUE) ans$mh[sx,XD[['aglst']][[ag]],cd4,1] <-  ma(ans$mh[sx,XD[['aglst']][[ag]],cd4,1],n=6)

        ## mortality on ART
        ans$ma[sx,XD[['aglst']][[ag]],cd4,"[0,6)m"] <- XD[[iso]]$mua1[cd4,paste0(sx,ag)]
        ans$ma[sx,XD[['aglst']][[ag]],cd4,"[7,12)m"] <- XD[[iso]]$mua2[cd4,paste0(sx,ag)]
        ans$ma[sx,XD[['aglst']][[ag]],cd4,"[12,Inf)m"] <- XD[[iso]]$mua3[cd4,paste0(sx,ag)]

        ## CD4 init
        ans$cdinit[sx,XD[['aglst']][[ag]],cd4,1] <- XD[[iso]]$ni[cd4,paste0(sx,ag)]/1e2
      }
    }
  }

  ## age ratio, sex ratio
  ans$SR <- XD[['SR']]
  ans$AR <- XD[['AR']]

  ## ART mortality ratio
  ans$rr1 <- XD[[iso]]$rr1
  ans$rr2 <- XD[[iso]]$rr2

  ## eligibility NOTE not used
  ans$elig <- data.table(year=XD[[iso]]$y1,cd4=XD[[iso]]$th) #top category eligible in nmhiv
  ans$alc <- XD[[iso]]$alc        #allocation parm

  ## HIV incidence, prevalence etc. 
  ans$Adata <- data.table(
    iso3=iso,
    year=XD[[iso]]$y2,
    mlhiv=XD[[iso]]$mh,
    flhiv=XD[[iso]]$fh,
    artm=XD[[iso]]$ma,
    artf=XD[[iso]]$fa,
    ltfu=XD[[iso]]$ltfu
  )
  ans$Adata[,plhiv:=(mlhiv+flhiv)]
  ans$Adata[,art:=(artm*mlhiv+artf*flhiv)/(plhiv)]
  ans$Hdata <- data.table(iso3=iso,
                          year=XD[[iso]]$y3,
                          hivi=XD[[iso]]$hi, #1549
                          hp80=XD[[iso]]$hp80,
                          hi80=XD[[iso]]$hi80,
                          hm80=XD[[iso]]$hm80,
                          hp80pc=XD[[iso]]$hp80pc*1e2,
                          hi80pc=XD[[iso]]$hi80pc*1e2,
                          hm80pc=XD[[iso]]$hm80pc*1e2
                          )

    ## UNPD data
  ans$UPD <- XD[[iso]]$UPD
  ## births
  ans$MB <- XD[[iso]]$UPD$MB
  ans$FB <- XD[[iso]]$UPD$FB

  ## ---- new data pre-computed for use in dynamics -----
  ans$tstep <- .1
  ans$syear <- 1970
  ans$eyear <- eyear
  ans$ntime <- ceiling((ans$eyear-ans$syear)/ans$tstep) + 1
  ans$yrz <- seq(from=ans$syear,to=ans$eyear,by=ans$tstep)
  ans$npy <- ceiling(1/ans$tstep)                       #times per year
  ans$yrbrks <- as.numeric(unlist(
    strsplit(as.character(ans$MB[,years]),
             split='-'))[seq(from=1,to=2*nrow(ans$MB),by=2)]) #year breaks
  ans$xyr <- ans$Adata[,max(year)]
  ans$syr <- ans$Adata[art>0,min(year)]
  ans$asy <- ans$syr
  ans$myr <- ans$syr
  ans$Adata <- ans$Adata[year>=ans$myr]
  ans$ady <- ans$Adata[year==ans$myr+1,1e-2*art]
  ans$y1 <- XD[[iso]]$y1

  tmp <- getAIMtargetsInc(ans)           #precomputed targets and other data
  ans <- c(ans,tmp)

  ## ======= return ==========
  ans
}


ma <- function(x,n=5){y <- filter(x,rep(1/n,n), sides=2); y[is.na(y)] <- x[is.na(y)]; y}
artsmooth <- function(D,n=3) as.numeric(ma(D$acov,n))


##' Precomputing targets for dynamics
##' @param L 
##' @return list of HIV/ART targets
##' @author Pete Dodd
getAIMtargetsInc <- function(L){
    list2env(L,envir=environment())

    ## quantities: HIV/ART; M/F births; mortality(!); HIV SR; migration
    yrzs <- unique(floor(yrz))
    nsmll <- length(yrzs)
    nbm <- nbf <- srv <- rep(0,nsmll)
    hivtarget <- hivtarget2 <- arttarget <- arttarget2 <- rep(0,nsmll)
    migrates <- mortrates <- array(0,dim=c(nsmll,XD[['dz']][1:2])) #yrzs x M/F x age
    j <- 0
    for(i in 1:ntime){
      if(!(i-1)%%(1/tstep)){           #NB assumes 1/tstep is integer
        j <- j+1
        yr <- floor(yrz[i])
        ## HIV & ART targets
        yr2 <- yr + 1
        if(yr2>eyear)
          yr2 <- yr
        hivtarget[j] <- Hdata[year==yr,hp80pc*1e-2]   #HIV prevalence in those age <= 80
        hivtarget2[j] <- Hdata[year==yr2,hp80pc*1e-2] #new
        AD <- Adata[,.(year,acov=art)]             #switch to new NOTE
        AD[,acv:=artsmooth(AD)]
        if(yr >= myr ){ #period with data
          if(yr>xyr)  #over end
            arttarget[j] <- AD[year==xyr,1e-2*acv]
          else
            arttarget[j] <- AD[year==yr,1e-2*acv]
          if( yr+1 > xyr ){
            tmp <- AD[year==xyr,1e-2*acv]
            if(length(tmp)!=1)stop(paste0('art prob ',yr,' ',
                                          yr2,' ',length(tmp)))
            arttarget2[j] <- tmp
            adf <- (AD[year==xyr,1e-2*acv]-AD[year==xyr-1,1e-2*acv])
            arttarget2[j] <- arttarget2[j] + adf*(yr+1-xyr)
            arttarget[j] <- arttarget[j] + adf*(yr-xyr)
          } else{
            arttarget2[j] <- AD[year==yr+1,1e-2*acv]
          }
        } else if(yr+1 > asy){      #early extrapolation period
          arttarget2[j] <- (yr+1-asy)*ady
          arttarget[j] <- max( (yr-asy)*ady, 0)
        }

        ## births
        who <- sum(yr>yrbrks)       #which base element
        yover <- yr - yrbrks[who]   #how far to next 5 yrs
        p <- yover/5
        nbm[j] <- (1-p)*MB[who,value] + p*MB[who+1,value] #boys
        nbf[j] <- (1-p)*FB[who,value] + p*FB[who+1,value] #girls
        ## migration rate
        migrates[j,1,] <- UPD$migr[year==yr & sex==1,value]
        migrates[j,2,] <- UPD$migr[year==yr & sex==2,value]
        ## mortality
        mortrates[j,1,] <- UPD$lfts[year==yr & sex==1 & age<81,Sx]^tstep #tstep mu
        mortrates[j,2,] <- UPD$lfts[year==yr & sex==2 & age<81,Sx]^tstep #tstep mu

        ## hiv disaggregation
        hsyear <- syear; if(is.na(hsyear)) hsyear <- 1975
        ypost <- floor(yrz[i] - hsyear) #years after HIV start
        sryr <- min(31,max(1,ypost))    #in between 
        sr <- 1                         #default in case not there
        if('SR' %in% names(L))
          sr <- L$SR[sryr]                            #sex ratio HIV
        sr <- sr/(1+sr)                             #fraction male HIV
        srv[j] <- sr
      }
    }

    return(list(hivtarget=hivtarget,hivtarget2=hivtarget2,
                arttarget=arttarget,arttarget2=arttarget2,
                nbmv=nbm,nbfv=nbf,srv=srv,
                migrates=migrates,mortrates=mortrates))

}

