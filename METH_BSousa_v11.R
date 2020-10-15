
# Script to run MeTHODICAL.
# Author: Bruno Sousa
# Version v11
#
# Some comments are left to guide its usage.
#

# rngBen <- 2:2
# rngCost <- 2:4
# myRNGBEN <- 1:2
# myRNGCOST <- 3:5
# myWrngBEN <- 1:1
# myWrngCost <- 2:4


#
# Auxiliary function to apply weights in the matrix with the diverse criteria
#
aux_step03_Weight_of_Norm_matrix <-function(inMatrix, inWeight){

  outMat <- inMatrix
  ncols<-dim(inMatrix)[2]
  nrows<-dim(inMatrix)[1]

  #iMcalc2 <- matrix(inMatrix[,2:ncols], nrow=nrows, ncol=ncols-1, byrow=TRUE)
  iMcalc <- matrix(inMatrix[,2:ncols], nrow=nrows, ncol=ncols-1, byrow=FALSE)

  nC <- dim(iMcalc)[2] # Convert weigths into a matrix in order to perform Mult
  aMmul <- matrix(rep(inWeight, times=nrows), ncol=nC, byrow=TRUE)

  auxCalc <- iMcalc * aMmul

  outMat[,2:ncols]<-auxCalc

  return(outMat)
}




#
#
#
#
METH_METHODICALv11 <- function(iM, iVecWei, rngBEN=1:5, rngCost=6, rngWBen=1:4, rngWCost=5, MeTHBeta=0.5, MeTHOmega=0.5, itry=1 ){

  MINSUM_TOPSIS <- 1e-99 # To avoid divisions by zero in normalization

  mBen_Criteria <- iM[,rngBEN]
  mCost_Criteria <- iM[,c(1,rngCost)] # Put also the id in the matrix of costs.

  vBen_weight <- iVecWei[rngWBen]
  vCost_weight <- iVecWei[rngWCost]

  applyWeight <- TRUE
  applyNorm <- TRUE


  norm_method <- "vector"

  fAuxPos <- function(iN){

    if (iN < 0) {
      iN <- iN * (-1)
    }
    return(iN)
  }

  fAux_Positive_Dif <- function(i1, i2){
    adif <- (i1 - i2)
    #if (adif < 0) {
    #  adif <- adif * (-1)
    #}
    adif<-fAuxPos(adif)
    return(adif)
  }

  #
  # Determine RScore
  #
  calcRscorev10 <- function(inM, alpha=0.5, omega=0.5){
    nrow <- dim(inM)[1]
    outMat <- inM
    outMat <- outMat[,-3]

    for (nR in 1:nrow){
      #outMat[nR,2] <-  sqrt(alpha*(inM[nR,2]) + omega*(inM[nR,3])   ) # no square
      outMat[nR,2] <-  sqrt((inM[nR,2]) + (inM[nR,3])   ) # no square
      outMat[nR,2] <- outMat[nR,2] #+ auxSD
    }

    return(outMat)
  }


  #
  # Example of normalization mechanism, others can be employed
  #
  aux_step02_VEC_Normalization <- function(inMatrix){

    #internal function to help in normalization
    fSum <- function(i){
      sumrow <- MINSUM_TOPSIS
      aSu <- sum(i^2) + sumrow
      return(aSu)
    }

    outMat <- inMatrix
    ncols<-dim(inMatrix)[2]
    nrows<-dim(inMatrix)[1]

    iMcalc <- matrix(inMatrix[,2:ncols], nrow=nrows, ncol=ncols-1, byrow=FALSE)
    #print(iMcalc)

    #Apply Sum by col
    auxSum <- apply(iMcalc, 2, FUN=fSum)
    #print(auxSum)

    nC <- dim(iMcalc)[2]
    aMSum <- matrix(rep(auxSum, times=nrows), ncol=nC, byrow=TRUE)
    auxCalc <- iMcalc / sqrt(aMSum)
    #print(auxCalc)

    outMat[,2:ncols]<-auxCalc
    return(outMat)

  }


  calcDist_METHODICALv10b<- function(inMatrix, inIdeal, inMin, inMax, BenCriteria=TRUE, ignFirstCol=TRUE){

    ncols<-dim(inMatrix)[2]
    nrows<-dim(inMatrix)[1]

    if (ignFirstCol){
      iMcalc <- matrix(inMatrix[,2:ncols], nrow=nrows, ncol=ncols-1, byrow=FALSE)

    }else{
      iMcalc <- inMatrix

    }

    fCalcAux <- function(iMcalc){
      auxMean <- mean(iMcalc)
      auxSD   <- var(iMcalc) # MAIN difference in MeTH_v10
      #auxSD   <- sd(iMcalc) # MAIN difference in MeTH_v10

      if (BENC){
        SDMean <- auxMean + auxSD
      }else{
        SDMean <- auxMean - auxSD
      }
      return(SDMean)
    }

    outMat <- matrix(ncol=2,nrow=nrows)
    outMat[,1] <- inMatrix[,1]

    nColMat <- dim(iMcalc)[2]
    ncolIdeal <- dim(inIdeal)[2]

    #}
    # Just to be sure
    stopifnot(nColMat == ncolIdeal)

    BENC <<- BenCriteria

    SDMean <- apply(iMcalc, 2, FUN=fCalcAux)

    mSDMean <- matrix(SDMean, byrow=TRUE, nrow=dim(iMcalc)[1], ncol=dim(iMcalc)[2])
    mIdeal  <- matrix(inIdeal, byrow=TRUE, nrow=dim(iMcalc)[1], ncol=dim(iMcalc)[2])
    mMax  <- matrix(inMax, byrow=TRUE, nrow=dim(iMcalc)[1], ncol=dim(iMcalc)[2])
    mMin  <- matrix(inMin, byrow=TRUE, nrow=dim(iMcalc)[1], ncol=dim(iMcalc)[2])

    #Perform Operation
    #auxTETA <- -1 *  (mSDMean - iMcalc)
    #Yaux <- (iMcalc - mIdeal)^2
    #auxDiffDeno <- (mMax - mMin)^2 #originally without sqrt
    #auxDeno <-  (auxTETA + auxDiffDeno)


    #auxTETA <- -1 *  (mSDMean - iMcalc)
    #Yaux <- (iMcalc - mIdeal)^2 + (mIdeal - mSDMean)^2
    #auxDiffDeno <- (mMax - mMin)^2 #originally without sqrt
    #auxDeno <-  (auxDiffDeno - auxTETA)

    Yaux <- (iMcalc - mIdeal)^2
    auxDeno <-  abs(mIdeal - mSDMean) + 0.001


    auxCalcM <- (Yaux /  auxDeno )
    auxCalcM[is.na(auxCalcM)]<-0

    auxCalc <- apply(auxCalcM, 1, sum)

    outMat[,2] <- auxCalc

    return(outMat)
  }

  ncolMB <- dim(mBen_Criteria)[2]
  MeTHTOPsisBenefits <- as.matrix(mBen_Criteria )

  ncolMC <- dim(mCost_Criteria)[2]
  MeTHTOPsisCosts <- as.matrix(mCost_Criteria )

  MeTHWeiBenTOP <- vBen_weight
  MeTHWeiCostTOP <- vCost_weight

  #
  # Step 01
  #

  # apply Normalization
  if (applyNorm==TRUE){

    if (norm_method =="vector"){
      MeTHTOPsisBenefits <- aux_step02_VEC_Normalization(MeTHTOPsisBenefits)
      MeTHTOPsisCosts    <- aux_step02_VEC_Normalization(MeTHTOPsisCosts)
    }

    # TODO: (if needed) Other normalization methods can be implemented
    # Add them here according to your needs.

  }


  #Apply Weighting
  if (applyWeight==TRUE){
    MeTHTOPsisBenefits <-aux_step03_Weight_of_Norm_matrix(MeTHTOPsisBenefits, MeTHWeiBenTOP)
    MeTHTOPsisCosts    <-aux_step03_Weight_of_Norm_matrix(MeTHTOPsisCosts, MeTHWeiCostTOP)
  }


  #
  # Step 02
  #
  # Retrieve Maximum and Minimum
  #
  rngBBen <- seq(from=2, by=1, to=dim(MeTHTOPsisBenefits)[2])
  rngCCost <- seq(from=2, by=1, to=dim(MeTHTOPsisCosts)[2])

  MeTHMaxIdeal<- as.numeric(apply(as.matrix(MeTHTOPsisBenefits[,rngBBen]), 2, max))
  MeTHBenMin  <- as.numeric(apply(as.matrix(MeTHTOPsisBenefits[,rngBBen]), 2, min))
  MeTHBenMax  <- MeTHMaxIdeal

  MeTHMinIdeal<- as.numeric(apply(as.matrix(MeTHTOPsisCosts[,rngCCost]), 2, min))
  MeTHCostMax <- as.numeric(apply(as.matrix(MeTHTOPsisCosts[,rngCCost]), 2, max))
  MeTHCostMin <- MeTHMinIdeal

  MeTHBenDist<-calcDist_METHODICALv10b(MeTHTOPsisBenefits,MeTHMaxIdeal, MeTHBenMin, MeTHBenMax, BenCriteria=TRUE)
  MeTHCostDist<-calcDist_METHODICALv10b(MeTHTOPsisCosts,MeTHMinIdeal, MeTHCostMin, MeTHCostMax, BenCriteria=FALSE)


  #
  # Step 04
  #
  inMatForR <- MeTHBenDist
  inMatForR <- cbind(inMatForR, MeTHCostDist[,2])

  MeTHRScore<-calcRscorev10(inMatForR, MeTHBeta, MeTHOmega)

  orderRet  <-  MeTHRScore[order(MeTHRScore[,2]),]

  return(orderRet)

}
