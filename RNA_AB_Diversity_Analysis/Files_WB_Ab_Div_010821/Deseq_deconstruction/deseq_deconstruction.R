## Analysis using DEseq2, controlling for AB diversity ##

library("DESeq2")

#==============================================================================
# Preparing data
#==============================================================================

## Loading Data ##

data <- read.csv("rna_data_raw_counts.csv", sep=',')
meta_data <- read.csv("meta_data_with_Qvalues_250721.csv")

## 'Standard' Deseq2 analysis, comparing old v young ##

row.names(data) <- data$ensgene
data <- data[-1]

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = meta_data,
                              design = ~ age +Q4.00)

#==============================================================================
# Auxiliary functions
#==============================================================================

sanitizeRowRanges <- function(object) {
  if (is.null(mcols(mcols(object)))) {
    mcols(mcols(object)) <- DataFrame(type=rep("input",ncol(mcols(object))),
                                      description=character(ncol(mcols(object))))
  }
  class(mcols(mcols(object))$type) <- "character"
  class(mcols(mcols(object))$description) <- "character"
  mcols(mcols(object))$type[ is.na(mcols(mcols(object))$type) ] <- ""
  mcols(mcols(object))$description[ is.na(mcols(mcols(object))$description) ] <- ""
  object
}

designAndArgChecker <- function(object, betaPrior) {
  termsOrder <- attr(terms.formula(design(object)),"order")
  hasIntercept <- attr(terms(design(object)),"intercept") == 1
  interactionPresent <- any(termsOrder > 1)
  if (betaPrior & !hasIntercept) {
    stop("betaPrior=TRUE can only be used if the design has an intercept.
  if specifying + 0 in the design formula, use betaPrior=FALSE")
  }
  if (betaPrior & interactionPresent) {
    stop("betaPrior=FALSE should be used for designs with interactions")
  }
  
  if (!betaPrior) {
    mm <- stats::model.matrix(design(object), data=as.data.frame(colData(object)))
    q <- qr(mm)
    if (q$rank < ncol(mm))
      stop("full model matrix is less than full rank")
  }
  
  design <- design(object)
  designVars <- all.vars(design)
  if (length(designVars) > 0) {
    if (any(sapply(designVars, function(v) any(is.na(colData(object)[[v]]))))) {
      stop("variables in the design formula cannot have NA values")
    }
    designFactors <- designVars[sapply(designVars, function(v) is(colData(object)[[v]], "factor"))]
    if (length(designFactors) > 0 && any(sapply(designFactors,function(v) any(table(colData(object)[[v]]) == 0)))) {
      stop("factors in design formula must have samples for each level.
  this error can arise when subsetting a DESeqDataSet, in which
  all the samples for one or more levels of a factor in the design were removed.
  if this was intentional, use droplevels() to remove these levels, e.g.:
  dds$condition <- droplevels(dds$condition)
")
    }
    if (any(sapply(designVars, function(v) is(colData(object)[[v]], "ordered")))) {
      stop("the design formula contains an ordered factor. The internal steps
do not work on ordered factors as a formula. Instead you should provide a matrix to
the 'design' slot or to the 'full' argument of DESeq(), constructed using model.matrix.")
    }
  }
}

checkFullRank <- function(modelMatrix) {
  if (qr(modelMatrix)$rank < ncol(modelMatrix)) {
    if (any(apply(modelMatrix, 2, function(col) all(col == 0)))) {
      stop("the model matrix is not full rank, so the model cannot be fit as specified.
  Levels or combinations of levels without any samples have resulted in
  column(s) of zeros in the model matrix.
  Please read the vignette section 'Model matrix not full rank':
  vignette('DESeq2')")
    } else {
      stop("the model matrix is not full rank, so the model cannot be fit as specified.
  One or more variables or interaction terms in the design formula are linear
  combinations of the others and must be removed.
  Please read the vignette section 'Model matrix not full rank':
  vignette('DESeq2')")
    }
  }
}

getBaseMeansAndVariances <- function(object) {
  cts.norm <- counts(object,normalized=TRUE)
  if ("weights" %in% assayNames(object)) {
    wts <- assays(object)[["weights"]]
    cts.norm <- wts * cts.norm
  }
  meanVarZero <- DataFrame(baseMean = unname(rowMeans(cts.norm)),
                           baseVar = unname(rowVars(cts.norm)),
                           allZero = unname(rowSums(counts(object)) == 0))
  mcols(meanVarZero) <- DataFrame(type = rep("intermediate",ncol(meanVarZero)),
                                  description = c("mean of normalized counts for all samples",
                                                  "variance of normalized counts for all samples",
                                                  "all counts for a gene are zero"))
  if (all(c("baseMean","baseVar","allZero") %in% names(mcols(object)))) {
    mcols(object)[c("baseMean","baseVar","allZero")] <- meanVarZero
  } else {
    mcols(object) <- cbind(mcols(object),meanVarZero)
  }
  return(object)
}

getModelMatrix <- function(object) {
  if (is(design(object), "matrix")) {
    design(object)
  } else if (is(design(object), "formula")) {
    stats::model.matrix.default(design(object), data=as.data.frame(colData(object)))
  }
}

linearModelMu <- function(y, x) {
  qrx <- qr(x)    
  Q <- qr.Q(qrx)  
  Rinv <- solve(qr.R(qrx))
  # old code:
  # hatmatrix <- x %*% Rinv %*% t(Q)
  # t(hatmatrix %*% t(y))
  # Wolfgang Huber's rewrite is up to 2 orders of magnitude faster (Sept 2018):
  (y %*% Q) %*% t(x %*% Rinv)
}

roughDispEstimate <- function(y, x) {
  
  # must be positive
  mu <- linearModelMu(y, x)
  mu <- matrix(pmax(1, mu), ncol=ncol(mu))
  
  m <- nrow(x)
  p <- ncol(x)
  
  # an alternate rough estimator with higher mean squared or absolute error
  # (rowSums( (y - mu)^2/(mu * (m - p)) ) - 1)/rowMeans(mu)
  
  # rough disp ests will be adjusted up to minDisp later
  est <- rowSums( ((y - mu)^2 - mu) / mu^2 ) / (m - p)
  pmax(est, 0)
}

momentsDispEstimate <- function(object) {
  xim <- if (!is.null(normalizationFactors(object))) {
    mean(1/colMeans(normalizationFactors(object)))
  } else {
    mean(1/sizeFactors(object))
  }
  bv <- mcols(object)$baseVar
  bm <- mcols(object)$baseMean
  (bv - xim*bm)/bm^2
}

modelMatrixGroups <- function(x) {
  factor(unname(apply(x, 1, paste0, collapse="__")))
}

# convenience function to make more descriptive names
# for factor variables
renameModelMatrixColumns <- function(data, design) {
  data <- as.data.frame(data)
  designVars <- all.vars(design)
  designVarsClass <- sapply(designVars, function(v) class(data[[v]]))
  factorVars <- designVars[designVarsClass == "factor"]
  colNamesFrom <- make.names(do.call(c,lapply(factorVars, function(v) paste0(v,levels(data[[v]])[-1]))))
  colNamesTo <- make.names(do.call(c,lapply(factorVars, function(v) paste0(v,"_",levels(data[[v]])[-1],"_vs_",levels(data[[v]])[1]))))
  data.frame(from=colNamesFrom,to=colNamesTo,stringsAsFactors=FALSE)
}

# convenience function for building results tables
# out of a list and filling in NA rows
buildDataFrameWithNARows <- function(resultsList, NArows) {
  lengths <- sapply(resultsList,length)
  if (!all(lengths == lengths[1])) {
    stop("lengths of vectors in resultsList must be equal")
  }
  if (sum(!NArows) != lengths[1]) {
    stop("number of non-NA rows must be equal to lengths of vectors in resultsList")
  }
  if (sum(NArows) == 0) {
    return(DataFrame(resultsList))
  }
  dfFull <- DataFrame(lapply(resultsList, function(x) vector(mode(x), length(NArows))))
  dfFull[NArows,] <- NA
  dfFull[!NArows,] <- DataFrame(resultsList)
  dfFull
}

# convenience function for building larger matrices
# by filling in NA rows
buildMatrixWithNARows <- function(m, NARows) {
  mFull <- matrix(NA, ncol=ncol(m), nrow=length(NARows))
  mFull[!NARows,] <- m
  mFull
}

# simple function to return a matrix of size factors
# or normalization factors
getSizeOrNormFactors <- function(object) {
  if (!is.null(normalizationFactors(object))) {
    return(normalizationFactors(object))
  } else { 
    return(matrix(rep(sizeFactors(object),each=nrow(object)),
                  ncol=ncol(object)))
  }
}

getAndCheckWeights <- function(object, modelMatrix, weightThreshold=1e-2) {
  if ("weights" %in% assayNames(object)) {
    useWeights <- TRUE
    weights <- unname(assays(object)[["weights"]])
    stopifnot(all(weights >= 0))
    weights <- weights / apply(weights, 1, max)
    # some code for testing whether still full rank
    # only performed once per analysis, by setting object attribute
    if (is.null(attr(object, "weightsOK"))) {
      m <- ncol(modelMatrix)
      full.rank <- qr(modelMatrix)$rank == m
      weights.ok <- logical(nrow(weights))
      # most designs are full rank with current version of DESeq2
      if (full.rank) {
        for (i in seq_len(nrow(weights))) {
          # note: downweighting of samples very low will still be full rank
          # so this test is kind of minimally in play -- good for checking
          # the user input however, e.g. all zero weights for a gene
          test1 <- qr(weights[i,] * modelMatrix)$rank == m
          # we test that it will be possible to calculate the CR term
          # following subsetting based on weightThreshold
          mm.sub <- modelMatrix[weights[i,] > weightThreshold,,drop=FALSE]
          mm.sub <- mm.sub[,colSums(abs(mm.sub)) > 0,drop=FALSE]
          test2 <- qr(mm.sub)$rank == ncol(mm.sub)
          weights.ok[i] <- test1 & test2
        }
      } else {
        # model matrix is not full rank (backwards compatibility for betaPrior=TRUE)
        # just check zero columns
        weights.ok <- rep(TRUE, nrow(weights))
        for (j in seq_len(ncol(modelMatrix))) {
          num.zero <- colSums(t(weights) * modelMatrix[,j] == 0)
          weights.ok <- weights.ok & (num.zero != nrow(modelMatrix))
        }
      }
      # instead of giving an error, switch allZero to TRUE for the problem rows
      if (!all(weights.ok)) {
        mcols(object)$allZero[!weights.ok] <- TRUE
        weightsDF <- DataFrame(weightsFail = !weights.ok)
        mcols(weightsDF) <- DataFrame(type="intermediate",
                                      description="weights fail to allow parameter estimation")
        mcols(object) <- cbind(mcols(object), weightsDF)
        warning(paste("for", sum(!weights.ok),
                      "row(s), the weights as supplied won't allow parameter estimation, producing a
  degenerate design matrix. These rows have been flagged in mcols(dds)$weightsFail
  and treated as if the row contained all zeros (mcols(dds)$allZero set to TRUE).
  If you are blocking for donors/organisms, consider design = ~0+donor+condition."))
      }
    }
    attr(object, "weightsOK") <- TRUE
  } else {
    useWeights <- FALSE
    weights <- matrix(1, nrow=nrow(object), ncol=ncol(object))
  }
  list(object=object,weights=weights,useWeights=useWeights)
}

fitBetaWrapper <- function (ySEXP, xSEXP, nfSEXP, alpha_hatSEXP, contrastSEXP,
                            beta_matSEXP, lambdaSEXP, weightsSEXP, useWeightsSEXP,
                            tolSEXP, maxitSEXP, useQRSEXP, minmuSEXP) {
  if ( missing(contrastSEXP) ) {
    # contrast is not required, just give 1,0,0,...
    contrastSEXP <- c(1,rep(0,ncol(xSEXP)-1))
  }
  # test for any NAs in arguments
  arg.names <- names(formals(fitBetaWrapper))
  na.test <- sapply(mget(arg.names), function(x) any(is.na(x)))
  if (any(na.test)) stop(paste("in call to fitBeta, the following arguments contain NA:",
                               paste(arg.names[na.test],collapse=", ")))
  
  fitBeta(ySEXP=ySEXP, xSEXP=xSEXP, nfSEXP=nfSEXP, alpha_hatSEXP=alpha_hatSEXP,
          contrastSEXP=contrastSEXP, beta_matSEXP=beta_matSEXP,
          lambdaSEXP=lambdaSEXP, weightsSEXP=weightsSEXP, useWeightsSEXP=useWeightsSEXP,
          tolSEXP=tolSEXP, maxitSEXP=maxitSEXP, useQRSEXP=useQRSEXP, minmuSEXP=minmuSEXP)
}

#==============================================================================
# Define DESeq args
#==============================================================================

object_raw <- dds
test <- "Wald"
fitType <- "parametric"
sfType <- "ratio"
betaPrior <- FALSE
full = design(object_raw)
# reduced = NULL # LRT test only
quiet = TRUE
minReplicatesForReplace = 7
modelMatrixType = NULL
useT = FALSE
minmu = 0.5
parallel = FALSE
BPPARAM = bpparam()
betaTol = 1e-08
maxit = 100
useOptim = TRUE
useQR = TRUE
alphaInit=NULL
niter=1
minDisp=1e-8
linearMu=NULL
#test#
useWeights <- FALSE
type = c("DESeq2", "glmGamPoi")

#==============================================================================
# Initial bookkeeping and data checks
#==============================================================================

modelAsFormula <- !is.matrix(full) & is(design(object_raw), "formula")
object <- sanitizeRowRanges(object_raw)

if (modelAsFormula) {
  designAndArgChecker(object, betaPrior)
  modelMatrix <- NULL
} else {
  checkFullRank(full)
  modelMatrix <- full
}
attr(object, "betaPrior") <- betaPrior
stopifnot(length(parallel) == 1 & is.logical(parallel))

#==============================================================================
# Estimate size factors
#==============================================================================

if (!is.null(sizeFactors(object)) || !is.null(normalizationFactors(object))) {
  if (!quiet) {
    if (!is.null(normalizationFactors(object))) {
      message("using pre-existing normalization factors")
    } else {
      message("using pre-existing size factors")
    }
  }
} else {
  object <- estimateSizeFactors(object, type = sfType, quiet = TRUE)
}

#==============================================================================
# De-constructing estimate size factors - nevermind. Doesn't use model. ##
#==============================================================================

#==============================================================================
# Estimate dispersions (for normalisation)
#==============================================================================

object <- estimateDispersions(object, fitType = fitType,
                              quiet = quiet, modelMatrix = modelMatrix, minmu = minmu)

#==============================================================================
# De-constructing estimate dispersions - this is the first time model is used #
#==============================================================================

#==========================================
## Define model matrix
#==========================================

if (is.null(modelMatrix)) {
  modelMatrix <- getModelMatrix(object) 
}

#==========================================
## takes object and gets base mean and variances of object ##
# An internally used function to calculate the row means and variances
# from the normalized counts, which requires that \code{\link{estimateSizeFactors}}
# has already been called.  Adds these and a logical column if the row sums
# are zero to the mcols of the object.
#==========================================
object <- getBaseMeansAndVariances(object)

#==========================================
## only continue on rows with non-zero row mean ##
## this doesn't change anything in our dataset, so we can drop ##
## but we run here to keep object name as objectNZ ##
#==========================================
#objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]
objectNZ <- object
#==========================================

if (is.null(alphaInit)) {
  # this rough dispersion estimate (alpha_hat)
  # is for estimating mu
  # and for the initial starting point for line search
  roughDisp <- roughDispEstimate(y = counts(objectNZ,normalized=TRUE),
                                 x = modelMatrix)
  momentsDisp <- momentsDispEstimate(objectNZ)
  alpha_hat <- pmin(roughDisp, momentsDisp)
} else {
  if (length(alphaInit) == 1) {
    alpha_hat <- rep(alphaInit, nrow(objectNZ))
  } else {
    stopifnot(length(alphaInit) == nrow(objectNZ))
    alpha_hat <- alphaInit
  }
}

# bound the rough estimated alpha between minDisp and maxDisp for numeric stability
maxDisp <- max(10, ncol(object))
alpha_hat <- alpha_hat_new <- alpha_init <- pmin(pmax(minDisp, alpha_hat), maxDisp)

stopifnot(length(niter) == 1 & niter > 0)

#==========================================

# use a linear model to estimate the expected counts
# if the number of groups according to the model matrix
# is equal to the number of columns

# what does this do exactly? #

#==========================================

if (is.null(linearMu)) {
  modelMatrixGroups <- modelMatrixGroups(modelMatrix)
  linearMu <- nlevels(modelMatrixGroups) == ncol(modelMatrix)
  # also check for weights (then can't do linear mu)
  if (useWeights) {
    linearMu <- FALSE
  }
}


#==============================================================================
# Fit binomial GLM function
#==============================================================================

fitNbinomGLMs <- function(object, modelMatrix=NULL, modelFormula, alpha_hat, lambda,
                          renameCols=TRUE, betaTol=1e-8, maxit=100, useOptim=TRUE,
                          useQR=TRUE, forceOptim=FALSE, warnNonposVar=TRUE, minmu=0.5,
                          type = c("DESeq2", "glmGamPoi")) {
  type <- match.arg(type, c("DESeq2", "glmGamPoi"))
  
  if (missing(modelFormula)) {
    modelFormula <- design(object)
  }
  if (is.null(modelMatrix)) {
    modelAsFormula <- TRUE
    modelMatrix <- stats::model.matrix.default(modelFormula, data=as.data.frame(colData(object)))
  } else {
    modelAsFormula <- FALSE
  }
  
  stopifnot(all(colSums(abs(modelMatrix)) > 0))
  
  # rename columns, for use as columns in DataFrame
  # and to emphasize the reference level comparison
  modelMatrixNames <- colnames(modelMatrix)
  modelMatrixNames[modelMatrixNames == "(Intercept)"] <- "Intercept"
  modelMatrixNames <- make.names(modelMatrixNames)
  
  if (renameCols) {
    convertNames <- renameModelMatrixColumns(colData(object),
                                             modelFormula)
    convertNames <- convertNames[convertNames$from %in% modelMatrixNames,,drop=FALSE]
    modelMatrixNames[match(convertNames$from, modelMatrixNames)] <- convertNames$to
  }
  colnames(modelMatrix) <- modelMatrixNames
  
  normalizationFactors <- getSizeOrNormFactors(object)
  
  if (missing(alpha_hat)) {
    alpha_hat <- dispersions(object)
  }
  
  if (length(alpha_hat) != nrow(object)) {
    stop("alpha_hat needs to be the same length as nrows(object)")
  }
  
  # set a wide prior for all coefficients
  if (missing(lambda)) {
    lambda <- rep(1e-6, ncol(modelMatrix))
  }
  
  # use weights if they are present in assays(object)
  wlist <- getAndCheckWeights(object, modelMatrix)
  weights <- wlist$weights
  useWeights <- wlist$useWeights
  
  if(type == "glmGamPoi"){
    stopifnot("type = 'glmGamPoi' cannot handle weights" = ! useWeights,
              "type = 'glmGamPoi' does not support NA's in alpha_hat" = all(! is.na(alpha_hat))) 
    gp_res <- glmGamPoi::glm_gp(counts(object), design = modelMatrix,
                                size_factors = FALSE, offset = log(normalizationFactors),
                                overdispersion = alpha_hat, verbose = FALSE)
    logLikeMat <- dnbinom(counts(object), mu=gp_res$Mu, size=1/alpha_hat, log=TRUE)
    logLike <- rowSums(logLikeMat)
    res <- list(logLike = logLike, betaConv =  rep(TRUE, nrow(object)), betaMatrix = gp_res$Beta / log(2),
                betaSE = NULL, mu = gp_res$Mu, betaIter = rep(NA,nrow(object)),
                modelMatrix=modelMatrix, 
                nterms=ncol(modelMatrix), hat_diagonals = NULL)
    return(res)
  }
  
  # bypass the beta fitting if the model formula is only intercept and
  # the prior variance is large (1e6)
  # i.e., LRT with reduced ~ 1 and no beta prior
  justIntercept <- if (modelAsFormula) {
    modelFormula == formula(~ 1)
  } else {
    ncol(modelMatrix) == 1 & all(modelMatrix == 1)
  }
  if (justIntercept & all(lambda <= 1e-6)) {
    alpha <- alpha_hat
    betaConv <- rep(TRUE, nrow(object))
    betaIter <- rep(1,nrow(object))
    betaMatrix <- if (useWeights) {
      matrix(log2(rowSums(weights*counts(object, normalized=TRUE))
                  /rowSums(weights)),ncol=1)
    } else {
      matrix(log2(rowMeans(counts(object, normalized=TRUE))),ncol=1)
    }
    mu <- normalizationFactors * as.numeric(2^betaMatrix)
    logLikeMat <- dnbinom(counts(object), mu=mu, size=1/alpha, log=TRUE)
    logLike <- if (useWeights) {
      rowSums(weights*logLikeMat)
    } else {
      rowSums(logLikeMat)
    }
    modelMatrix <- stats::model.matrix.default(~ 1, data=as.data.frame(colData(object)))
    colnames(modelMatrix) <- modelMatrixNames <- "Intercept"
    w <- if (useWeights) {
      weights * (mu^-1 + alpha)^-1
    } else {
      (mu^-1 + alpha)^-1
    }
    xtwx <- rowSums(w)
    sigma <- xtwx^-1
    betaSE <- matrix(log2(exp(1)) * sqrt(sigma),ncol=1)      
    hat_diagonals <- w * xtwx^-1;
    res <- list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
                betaSE = betaSE, mu = mu, betaIter = betaIter,
                modelMatrix=modelMatrix, 
                nterms=1, hat_diagonals=hat_diagonals)
    return(res)
  }
  
  qrx <- qr(modelMatrix)
  # if full rank, estimate initial betas for IRLS below
  if (qrx$rank == ncol(modelMatrix)) {
    Q <- qr.Q(qrx)
    R <- qr.R(qrx)
    y <- t(log(counts(object,normalized=TRUE) + .1))
    beta_mat <- t(solve(R, t(Q) %*% y))
  } else {
    if ("Intercept" %in% modelMatrixNames) {
      beta_mat <- matrix(0, ncol=ncol(modelMatrix), nrow=nrow(object))
      # use the natural log as fitBeta occurs in the natural log scale
      logBaseMean <- log(rowMeans(counts(object,normalized=TRUE)))
      beta_mat[,which(modelMatrixNames == "Intercept")] <- logBaseMean
    } else {
      beta_mat <- matrix(1, ncol=ncol(modelMatrix), nrow=nrow(object))
    }
  }
  
  # here we convert from the log2 scale of the betas
  # and the beta prior variance to the log scale
  # used in fitBeta.
  # so we divide by the square of the
  # conversion factor, log(2)
  lambdaNatLogScale <- lambda / log(2)^2
  
  betaRes <- fitBetaWrapper(ySEXP = counts(object), xSEXP = modelMatrix,
                            nfSEXP = normalizationFactors,
                            alpha_hatSEXP = alpha_hat,
                            beta_matSEXP = beta_mat,
                            lambdaSEXP = lambdaNatLogScale,
                            weightsSEXP = weights,
                            useWeightsSEXP = useWeights,
                            tolSEXP = betaTol, maxitSEXP = maxit,
                            useQRSEXP=useQR, minmuSEXP=minmu)
  
  # Note on deviance: the 'deviance' calculated in fitBeta() (C++)
  # is not returned in mcols(object)$deviance. instead, we calculate
  # the log likelihood below and use -2 * logLike.
  # (reason is that we have other ways of estimating beta:
  # above intercept code, and below optim code)
  
  mu <- normalizationFactors * t(exp(modelMatrix %*% t(betaRes$beta_mat)))
  dispersionVector <- rep(dispersions(object), times=ncol(object))
  logLike <- nbinomLogLike(counts(object), mu, dispersions(object), weights, useWeights)
  
  # test for stability
  rowStable <- apply(betaRes$beta_mat,1,function(row) sum(is.na(row))) == 0
  
  # test for positive variances
  rowVarPositive <- apply(betaRes$beta_var_mat,1,function(row) sum(row <= 0)) == 0
  
  # test for convergence, stability and positive variances
  betaConv <- betaRes$iter < maxit
  
  # here we transform the betaMatrix and betaSE to a log2 scale
  betaMatrix <- log2(exp(1))*betaRes$beta_mat
  colnames(betaMatrix) <- modelMatrixNames
  colnames(modelMatrix) <- modelMatrixNames
  # warn below regarding these rows with negative variance
  betaSE <- log2(exp(1))*sqrt(pmax(betaRes$beta_var_mat,0))
  colnames(betaSE) <- paste0("SE_",modelMatrixNames)
  
  # switch based on whether we should also use optim
  # on rows which did not converge
  rowsForOptim <- if (useOptim) {
    which(!betaConv | !rowStable | !rowVarPositive)
  } else {
    which(!rowStable | !rowVarPositive)
  }
  
  if (forceOptim) {
    rowsForOptim <- seq_along(betaConv)
  }
  
  if (length(rowsForOptim) > 0) {
    # we use optim if didn't reach convergence with the IRLS code
    resOptim <- fitNbinomGLMsOptim(object,modelMatrix,lambda,
                                   rowsForOptim,rowStable,
                                   normalizationFactors,alpha_hat,
                                   weights,useWeights,
                                   betaMatrix,betaSE,betaConv,
                                   beta_mat,
                                   mu,logLike,minmu=minmu)
    betaMatrix <- resOptim$betaMatrix
    betaSE <- resOptim$betaSE
    betaConv <- resOptim$betaConv
    mu <- resOptim$mu
    logLike <- resOptim$logLike
  }
  
  stopifnot(!any(is.na(betaSE)))
  nNonposVar <- sum(rowSums(betaSE == 0) > 0)
  if (warnNonposVar & nNonposVar > 0) warning(nNonposVar,"rows had non-positive estimates of variance for coefficients")
  
  list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
       betaSE = betaSE, mu = mu, betaIter = betaRes$iter, modelMatrix=modelMatrix, 
       nterms=ncol(modelMatrix), hat_diagonals=betaRes$hat_diagonals)
}

# this function calls fitNbinomGLMs() twice:
# 1 - without the beta prior, in order to calculate the
#     beta prior variance and hat matrix
# 2 - again but with the prior in order to get beta matrix and standard errors
fitGLMsWithPrior <- function(object, betaTol, maxit, useOptim, useQR, betaPriorVar, modelMatrix=NULL, minmu=0.5) {
  
  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]
  modelMatrixType <- attr(object, "modelMatrixType")
  
  if (missing(betaPriorVar) | !(all(c("mu","H") %in% assayNames(objectNZ)))) {
    
    # stop unless modelMatrix was NOT supplied, the code below all works
    # by building model matrices using the formula, doesn't work with incoming model matrices
    stopifnot(is.null(modelMatrix))
    
    # fit the negative binomial GLM without a prior,
    # used to construct the prior variances
    # and for the hat matrix diagonals for calculating Cook's distance
    fit <- fitNbinomGLMs(objectNZ,
                         betaTol=betaTol, maxit=maxit,
                         useOptim=useOptim, useQR=useQR,
                         renameCols = (modelMatrixType == "standard"),
                         minmu=minmu)
    modelMatrix <- fit$modelMatrix
    modelMatrixNames <- colnames(modelMatrix)
    H <- fit$hat_diagonal
    betaMatrix <- fit$betaMatrix
    mu <- fit$mu
    
    modelMatrixNames[modelMatrixNames == "(Intercept)"] <- "Intercept"
    modelMatrixNames <- make.names(modelMatrixNames)
    colnames(betaMatrix) <- modelMatrixNames
    
    # save the MLE log fold changes for addMLE argument of results
    convertNames <- renameModelMatrixColumns(colData(object),
                                             design(objectNZ))
    convertNames <- convertNames[convertNames$from %in% modelMatrixNames,,drop=FALSE]
    modelMatrixNames[match(convertNames$from, modelMatrixNames)] <- convertNames$to
    mleBetaMatrix <- fit$betaMatrix
    colnames(mleBetaMatrix) <- paste0("MLE_",modelMatrixNames)
    
    # store for use in estimateBetaPriorVar below
    mcols(objectNZ) <- cbind(mcols(objectNZ), DataFrame(mleBetaMatrix))
  } else {
    # we can skip the first MLE fit because the
    # beta prior variance and hat matrix diagonals were provided
    if (is.null(modelMatrix)) {
      modelMatrix <- getModelMatrix(object)
    }
    H <- assays(objectNZ)[["H"]]
    mu <- assays(objectNZ)[["mu"]]
    mleBetaMatrix <- as.matrix(mcols(objectNZ)[,grep("MLE_",names(mcols(objectNZ))),drop=FALSE])
  }
  
  if (missing(betaPriorVar)) {
    betaPriorVar <- estimateBetaPriorVar(objectNZ, modelMatrix=modelMatrix)
  } else {
    # else we are provided the prior variance:
    # check if the lambda is the correct length
    # given the design formula
    if (modelMatrixType == "expanded") {
      modelMatrix <- makeExpandedModelMatrix(objectNZ)
    }
    p <- ncol(modelMatrix)
    if (length(betaPriorVar) != p) {
      stop(paste("betaPriorVar should have length",p,"to match:",paste(colnames(modelMatrix),collapse=", ")))
    }
  }
  
  # refit the negative binomial GLM with a prior on betas
  if (any(betaPriorVar == 0)) {
    stop("beta prior variances are equal to zero for some variables")
  }
  lambda <- 1/betaPriorVar
  
  if (modelMatrixType == "standard") {
    fit <- fitNbinomGLMs(objectNZ, lambda=lambda,
                         betaTol=betaTol, maxit=maxit,
                         useOptim=useOptim, useQR=useQR,
                         minmu=minmu)
    modelMatrix <- fit$modelMatrix
  } else if (modelMatrixType == "expanded") {
    modelMatrix <- makeExpandedModelMatrix(objectNZ)
    fit <- fitNbinomGLMs(objectNZ, lambda=lambda,
                         betaTol=betaTol, maxit=maxit,
                         useOptim=useOptim, useQR=useQR,
                         modelMatrix=modelMatrix, renameCols=FALSE,
                         minmu=minmu)
  } else if (modelMatrixType == "user-supplied") {
    fit <- fitNbinomGLMs(objectNZ, lambda=lambda,
                         betaTol=betaTol, maxit=maxit,
                         useOptim=useOptim, useQR=useQR,
                         modelMatrix=modelMatrix, renameCols=FALSE,
                         minmu=minmu)
  }
  
  res <- list(fit=fit, H=H, betaPriorVar=betaPriorVar, mu=mu,
              modelMatrix=modelMatrix, mleBetaMatrix=mleBetaMatrix)
  res
}

#==========================================
### complicated function ### 
#==========================================

alpha_hatSEXP <- alpha_hat
alpha_hat

# below, iterate between mean and dispersion estimation (niter) times
fitidx <- rep(TRUE,nrow(objectNZ))
mu <- matrix(0, nrow=nrow(objectNZ), ncol=ncol(objectNZ))
dispIter <- numeric(nrow(objectNZ))
# bound the estimated count by 'minmu'
# this helps make the fitting more robust,
# because 1/mu occurs in the weights for the NB GLM
for (iter in seq_len(niter)) {
  if (!linearMu) {
    fit <- fitNbinomGLMs(objectNZ[fitidx,,drop=FALSE],
                         alpha_hat=alpha_hat[fitidx],
                         modelMatrix=modelMatrix, type=type)
    fitMu <- fit$mu
  } else {
    fitMu <- linearModelMuNormalized(objectNZ[fitidx,,drop=FALSE],
                                     modelMatrix)
  }
  fitMu[fitMu < minmu] <- minmu
  mu[fitidx,] <- fitMu
  
  
  # use of kappa_0 in backtracking search
  # initial proposal = log(alpha) + kappa_0 * deriv. of log lik. w.r.t. log(alpha)
  # use log(minDisp/10) to stop if dispersions going to -infinity
  if (type == "DESeq2") {
    dispRes <- fitDispWrapper(ySEXP = counts(objectNZ)[fitidx,,drop=FALSE],
                              xSEXP = modelMatrix,
                              mu_hatSEXP = fitMu,
                              log_alphaSEXP = log(alpha_hat)[fitidx],
                              log_alpha_prior_meanSEXP = log(alpha_hat)[fitidx],
                              log_alpha_prior_sigmasqSEXP = 1, min_log_alphaSEXP = log(minDisp/10),
                              kappa_0SEXP = kappa_0, tolSEXP = dispTol,
                              maxitSEXP = maxit, usePriorSEXP = FALSE,
                              weightsSEXP = weights,
                              useWeightsSEXP = useWeights,
                              weightThresholdSEXP = weightThreshold,
                              useCRSEXP = useCR)
    
    dispIter[fitidx] <- dispRes$iter
    alpha_hat_new[fitidx] <- pmin(exp(dispRes$log_alpha), maxDisp)
    last_lp <- dispRes$last_lp
    initial_lp <- dispRes$initial_lp
    # only rerun those rows which moved
  } else if (type == "glmGamPoi") {
    if (!requireNamespace("glmGamPoi", quietly=TRUE)) {
      stop("type='glmGamPoi' requires installing the Bioconductor package 'glmGamPoi'")
    }
    if (!quiet) message("using 'glmGamPoi' as fitType. If used in published research, please cite:
    Ahlmann-Eltze, C., Huber, W. (2020) glmGamPoi: Fitting Gamma-Poisson
    Generalized Linear Models on Single Cell Count Data. Bioinformatics.
    https://doi.org/10.1093/bioinformatics/btaa1009")
    Counts <- counts(objectNZ)
    initial_lp <- vapply(which(fitidx), function(idx){
      sum(dnbinom(Counts[idx, ], mu = fitMu[idx, ], size = 1 / alpha_hat[idx], log = TRUE))
    }, FUN.VALUE = 0.0)
    dispersion_fits <- glmGamPoi::overdispersion_mle(Counts[fitidx, ], mean = fitMu[fitidx, ],
                                                     model_matrix = modelMatrix, verbose = ! quiet)
    dispIter[fitidx] <- dispersion_fits$iterations
    alpha_hat_new[fitidx] <- pmin(dispersion_fits$estimate, maxDisp)
    last_lp <- vapply(which(fitidx), function(idx){
      sum(dnbinom(Counts[idx, ], mu = fitMu[idx, ], size = 1 / alpha_hat_new[idx], log = TRUE))
    }, FUN.VALUE = 0.0)
  }
  fitidx <- abs(log(alpha_hat_new) - log(alpha_hat)) > .05
  alpha_hat <- alpha_hat_new
  if (sum(fitidx) == 0) break
}
# dont accept moves if the log posterior did not
# increase by more than one millionth,
# and set the small estimates to the minimum dispersion
dispGeneEst <- alpha_hat
if (niter == 1) {
  noIncrease <- last_lp < initial_lp + abs(initial_lp)/1e6
  dispGeneEst[which(noIncrease)] <- alpha_init[which(noIncrease)]
}
# didn't reach the maxmium and iterated more than once
dispGeneEstConv <- dispIter < maxit & !(dispIter == 1)

# if lacking convergence from fitDisp() (C++)...
refitDisp <- !dispGeneEstConv & dispGeneEst > minDisp*10
if (sum(refitDisp) > 0) {
  dispGrid <- fitDispGridWrapper(y = counts(objectNZ)[refitDisp,,drop=FALSE],
                                 x = modelMatrix,
                                 mu = mu[refitDisp,,drop=FALSE],
                                 logAlphaPriorMean = rep(0,sum(refitDisp)),
                                 logAlphaPriorSigmaSq = 1, usePrior = FALSE,
                                 weightsSEXP = weights[refitDisp,,drop=FALSE],
                                 useWeightsSEXP = useWeights,
                                 weightThresholdSEXP = weightThreshold,
                                 useCRSEXP = useCR)
  dispGeneEst[refitDisp] <- dispGrid
}
dispGeneEst <- pmin(pmax(dispGeneEst, minDisp), maxDisp)

dispDataFrame <- buildDataFrameWithNARows(list(dispGeneEst=dispGeneEst,
                                               dispGeneIter=dispIter),
                                          mcols(object)$allZero)
mcols(dispDataFrame) <- DataFrame(type=rep("intermediate",ncol(dispDataFrame)),
                                  description=c("gene-wise estimates of dispersion",
                                                "number of iterations for gene-wise"))
mcols(object) <- cbind(mcols(object), dispDataFrame)
assays(object, withDimnames=FALSE)[["mu"]] <- buildMatrixWithNARows(mu, mcols(object)$allZero)

object

#==============================================================================
# Prepare to run GLM fit
#==============================================================================

object <- nbinomWaldTest(object, betaPrior = betaPrior,
                         quiet = quiet, modelMatrix = modelMatrix, modelMatrixType = modelMatrixType, 
                         useT = useT, minmu = minmu)

# object <- sanitizeRowRanges(object)
# 
# if ("results" %in% mcols(mcols(object))$type) {
#   object <- removeResults(object)
# }
# 
# if (is.null(mcols(object)$allZero)) {
#   object <- getBaseMeansAndVariances(object)
# }
# 
# objectNZ <- object[!mcols(object)$allZero, , drop = FALSE] # Dropping something or other
# 
# 
# 
# if (is.null(modelMatrix)) {
#   modelAsFormula <- TRUE
#   termsOrder <- attr(terms.formula(design(object)), "order")
#   interactionPresent <- any(termsOrder > 1)
#   if (missing(betaPrior)) {
#     betaPrior <- FALSE
#   }
#   designAndArgChecker(object, betaPrior)
#   stopifnot(is.logical(betaPrior))
#   blindDesign <- design(object) == formula(~1)
#   if (blindDesign) {
#     betaPrior <- FALSE
#   }
#   if (missing(modelMatrixType) || is.null(modelMatrixType)) {
#     modelMatrixType <- if (betaPrior) {
#       "expanded"
#     }
#     else {
#       "standard"
#     }
#   }
#   if (modelMatrixType == "expanded" & !betaPrior) {
#     stop("expanded model matrices require a beta prior")
#   }
#   attr(object, "modelMatrixType") <- modelMatrixType
#   hasIntercept <- attr(terms(design(object)), "intercept") == 
#     1
#   renameCols <- hasIntercept
# }
