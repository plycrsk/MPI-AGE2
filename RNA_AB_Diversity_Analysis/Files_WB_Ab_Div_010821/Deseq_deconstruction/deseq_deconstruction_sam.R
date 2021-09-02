## Analysis using DEseq2, controlling for AB diversity ##

library("DESeq2")
library("tidyverse")

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
                              design = ~ Q4.00 + age)

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
# Estimate dispersions (for normalisation)
#==============================================================================

object <- estimateDispersions(object, fitType = fitType,
                              quiet = quiet, modelMatrix = modelMatrix, minmu = minmu)

#==============================================================================
# Prepare to run GLM fit
#==============================================================================

#object <- nbinomWaldTest(object, betaPrior = betaPrior,
#                         quiet = quiet, modelMatrix = modelMatrix, modelMatrixType = modelMatrixType, 
#                         useT = useT, minmu = minmu)

if (is.null(modelMatrix)) {
  modelAsFormula <- TRUE
  termsOrder <- attr(terms.formula(design(object)),"order")
  interactionPresent <- any(termsOrder > 1)
  if (missing(betaPrior)) {
    betaPrior <- FALSE
  }
  
  # run some tests common to DESeq, nbinomWaldTest, nbinomLRT
  designAndArgChecker(object, betaPrior)
  
  # what kind of model matrix to use
  stopifnot(is.logical(betaPrior))
  blindDesign <- design(object) == formula(~ 1)
  if (blindDesign) {
    betaPrior <- FALSE
  }
  if (missing(modelMatrixType) || is.null(modelMatrixType)) {
    modelMatrixType <- if (betaPrior) {
      "expanded"
    } else {
      "standard"
    }
  }
  if (modelMatrixType == "expanded" & !betaPrior) {
    stop("expanded model matrices require a beta prior")
  }
  # store modelMatrixType so it can be accessed by estimateBetaPriorVar
  attr(object, "modelMatrixType") <- modelMatrixType
  hasIntercept <- attr(terms(design(object)),"intercept") == 1
  renameCols <- hasIntercept
} else {
  # modelMatrix is not NULL, user-supplied
  if (missing(betaPrior)) {
    betaPrior <- FALSE
  }
  if (betaPrior) {
    if (missing(betaPriorVar)) stop("user-supplied model matrix with betaPrior=TRUE requires supplying betaPriorVar")
  }
  modelAsFormula <- FALSE
  attr(object, "modelMatrixType") <- "user-supplied"
  renameCols <- FALSE
}

# only continue on the rows with non-zero row mean
objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]
#weights <- weights[!mcols(object)$allZero,,drop=FALSE]

#==============================================================================
# Defining GLM fitting function
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

#==============================================================================
# Convenience function, required for GLM fit 
#==============================================================================


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

#==============================================================================
# Running GLM fitting function
#==============================================================================

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

# breaking out the optim backup code from fitNbinomGLMs
fitNbinomGLMsOptim <- function(object,modelMatrix,lambda,
                               rowsForOptim,rowStable,
                               normalizationFactors,alpha_hat,
                               weights,useWeights,
                               betaMatrix,betaSE,betaConv,
                               beta_mat,
                               mu,logLike,minmu=0.5) {
  x <- modelMatrix
  lambdaNatLogScale <- lambda / log(2)^2
  large <- 30
  for (row in rowsForOptim) {
    betaRow <- if (rowStable[row] & all(abs(betaMatrix[row,]) < large)) {
      betaMatrix[row,]
    } else {
      beta_mat[row,]
    }
    nf <- normalizationFactors[row,]
    k <- counts(object)[row,]
    alpha <- alpha_hat[row]
    objectiveFn <- function(p) {
      mu_row <- as.numeric(nf * 2^(x %*% p))
      logLikeVector <- dnbinom(k,mu=mu_row,size=1/alpha,log=TRUE)
      logLike <- if (useWeights) {
        sum(weights[row,] * logLikeVector)
      } else {
        sum(logLikeVector)
      }
      logPrior <- sum(dnorm(p,0,sqrt(1/lambda),log=TRUE))
      negLogPost <- -1 * (logLike + logPrior)
      if (is.finite(negLogPost)) negLogPost else 10^300
    }
    o <- optim(betaRow, objectiveFn, method="L-BFGS-B",lower=-large, upper=large)
    ridge <- if (length(lambdaNatLogScale) > 1) {
      diag(lambdaNatLogScale)
    } else {
      as.matrix(lambdaNatLogScale,ncol=1)
    }
    # if we converged, change betaConv to TRUE
    if (o$convergence == 0) {
      betaConv[row] <- TRUE
    }
    # with or without convergence, store the estimate from optim
    betaMatrix[row,] <- o$par
    # calculate the standard errors
    mu_row <- as.numeric(nf * 2^(x %*% o$par))
    # store the new mu vector
    mu[row,] <- mu_row
    mu_row[mu_row < minmu] <- minmu
    w <- if (useWeights) {
      diag((mu_row^-1 + alpha)^-1)
    } else {
      diag(weights[row,] * (mu_row^-1 + alpha)^-1)
    }
    xtwx <- t(x) %*% w %*% x
    xtwxRidgeInv <- solve(xtwx + ridge)
    sigma <- xtwxRidgeInv %*% xtwx %*% xtwxRidgeInv
    # warn below regarding these rows with negative variance
    betaSE[row,] <- log2(exp(1)) * sqrt(pmax(diag(sigma),0))
    logLikeVector <- dnbinom(k,mu=mu_row,size=1/alpha,log=TRUE)
    logLike[row] <- if (useWeights) {
      sum(weights[row,] * logLikeVector)
    } else {
      sum(logLikeVector)
    }
  }
  return(list(betaMatrix=betaMatrix,betaSE=betaSE,
              betaConv=betaConv,mu=mu,logLike=logLike))
}


# convenience function for testing the log likelihood
# for a count matrix, mu matrix and vector disp
nbinomLogLike <- function(counts, mu, disp, weights, useWeights) {
  if (is.null(disp)) return(NULL)
  if (useWeights) {
    rowSums(weights * matrix(dnbinom(counts,mu=mu,size=1/disp,
                                     log=TRUE),ncol=ncol(counts)))
  } else {
    rowSums(matrix(dnbinom(counts,mu=mu,size=1/disp,
                           log=TRUE),ncol=ncol(counts)))    
  }
}

fitDisp <- function(ySEXP, xSEXP, mu_hatSEXP, log_alphaSEXP, log_alpha_prior_meanSEXP, log_alpha_prior_sigmasqSEXP, min_log_alphaSEXP, kappa_0SEXP, tolSEXP, maxitSEXP, usePriorSEXP, weightsSEXP, useWeightsSEXP, weightThresholdSEXP, useCRSEXP) {
  .Call('_DESeq2_fitDisp', PACKAGE = 'DESeq2', ySEXP, xSEXP, mu_hatSEXP, log_alphaSEXP, log_alpha_prior_meanSEXP, log_alpha_prior_sigmasqSEXP, min_log_alphaSEXP, kappa_0SEXP, tolSEXP, maxitSEXP, usePriorSEXP, weightsSEXP, useWeightsSEXP, weightThresholdSEXP, useCRSEXP)
}

fitBeta <- function(ySEXP, xSEXP, nfSEXP, alpha_hatSEXP, contrastSEXP, beta_matSEXP, lambdaSEXP, weightsSEXP, useWeightsSEXP, tolSEXP, maxitSEXP, useQRSEXP, minmuSEXP) {
  .Call('_DESeq2_fitBeta', PACKAGE = 'DESeq2', ySEXP, xSEXP, nfSEXP, alpha_hatSEXP, contrastSEXP, beta_matSEXP, lambdaSEXP, weightsSEXP, useWeightsSEXP, tolSEXP, maxitSEXP, useQRSEXP, minmuSEXP)
}

fitDispGrid <- function(ySEXP, xSEXP, mu_hatSEXP, disp_gridSEXP, log_alpha_prior_meanSEXP, log_alpha_prior_sigmasqSEXP, usePriorSEXP, weightsSEXP, useWeightsSEXP, weightThresholdSEXP, useCRSEXP) {
  .Call('_DESeq2_fitDispGrid', PACKAGE = 'DESeq2', ySEXP, xSEXP, mu_hatSEXP, disp_gridSEXP, log_alpha_prior_meanSEXP, log_alpha_prior_sigmasqSEXP, usePriorSEXP, weightsSEXP, useWeightsSEXP, weightThresholdSEXP, useCRSEXP)
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


# fit the negative binomial GLM without a prior
# (in actuality a very wide prior with standard deviation 1e3 on log2 fold changes)
fit <- fitNbinomGLMs(objectNZ,
                     betaTol=betaTol, maxit=maxit,
                     useOptim=useOptim, useQR=useQR,
                     renameCols=renameCols,
                     modelMatrix=modelMatrix,
                     minmu=minmu)
H <- fit$hat_diagonals
mu <- fit$mu
modelMatrix <- fit$modelMatrix
modelMatrixNames <- fit$modelMatrixNames
# record the wide prior variance which was used in fitting
betaPriorVar <- rep(1e6, ncol(fit$modelMatrix))

fit$betaMatrix
fit$betaSE
fit$modelMatrix
fit$betaConv

model_fit_output <- fit$betaMatrix

fit$

model_fit_output <- as.data.frame(model_fit_output)

write_csv(model_fit_output, "Model_fit_output_beta_matrix_050821.csv")

colnames(modelMatrix)

# convenience function for building larger matrices
# by filling in NA rows
buildMatrixWithNARows <- function(m, NARows) {
  mFull <- matrix(NA, ncol=ncol(m), nrow=length(NARows))
  mFull[!NARows,] <- m
  mFull
}

# store 'mu' and 'H', the hat matrix diagonals
dimnames(mu) <- NULL
assays(objectNZ, withDimnames=FALSE)[["mu"]] <- mu
assays(object, withDimnames=FALSE)[["mu"]] <- buildMatrixWithNARows(mu, mcols(object)$allZero)
dimnames(H) <- NULL
assays(objectNZ, withDimnames=FALSE)[["H"]] <- H
assays(object, withDimnames=FALSE)[["H"]] <- buildMatrixWithNARows(H, mcols(object)$allZero)

# store the prior variance directly as an attribute
# of the DESeqDataSet object, so it can be pulled later by
# the results function (necessary for setting max Cook's distance)
attr(object,"betaPrior") <- betaPrior
attr(object,"betaPriorVar") <- betaPriorVar
attr(object,"modelMatrix") <- modelMatrix
attr(object,"test") <- "Wald"

#convenience function 

getModelMatrix <- function(object) {
  if (is(design(object), "matrix")) {
    design(object)
  } else if (is(design(object), "formula")) {
    stats::model.matrix.default(design(object), data=as.data.frame(colData(object)))
  }
}

calculateCooksDistance <- function(object, H, modelMatrix) {
  p <- ncol(modelMatrix)
  dispersions <- robustMethodOfMomentsDisp(object, modelMatrix)
  V <- assays(object)[["mu"]] + dispersions * assays(object)[["mu"]]^2
  PearsonResSq <- (counts(object) - assays(object)[["mu"]])^2 / V
  cooks <- PearsonResSq / p  * H / (1 - H)^2
  cooks
}

# calculate a robust method of moments dispersion,
# in order to estimate the dispersion excluding
# individual outlier counts which would raise the variance estimate
robustMethodOfMomentsDisp <- function(object, modelMatrix) {
  cnts <- counts(object,normalized=TRUE)
  # if there are 3 or more replicates in any cell
  threeOrMore <- nOrMoreInCell(modelMatrix,n=3)
  v <- if (any(threeOrMore)) {
    cells <- apply(modelMatrix,1,paste0,collapse="")
    cells <- unname(factor(cells,levels=unique(cells)))
    levels(cells) <- seq_along(levels(cells))
    levelsThreeOrMore <- levels(cells)[table(cells) >= 3]
    idx <- cells %in% levelsThreeOrMore
    cntsSub <- cnts[,idx,drop=FALSE]
    cellsSub <- factor(cells[idx])
    trimmedCellVariance(cntsSub, cellsSub)
  } else {
    trimmedVariance(cnts)
  }
  m <- rowMeans(cnts)
  alpha <- ( v - m ) / m^2
  # cannot use the typical minDisp = 1e-8 here or else all counts in the same
  # group as the outlier count will get an extreme Cook's distance
  minDisp <- 0.04
  alpha <- pmax(alpha, minDisp)
  alpha
}


# for each sample in the model matrix,
# are there n or more replicates in the same cell
# (including that sample)
# so for a 2 x 3 comparison, the returned vector for n = 3 is:
# FALSE, FALSE, TRUE, TRUE, TRUE
nOrMoreInCell <- function(modelMatrix, n){
  row_hash <- apply(modelMatrix, 1, paste0, collapse = "_")
  hash_table <- table(row_hash)
  numEqual <- as.vector(unname(hash_table[row_hash]))
  numEqual >= n
}

trimmedVariance <- function(x) {
  rm <-  apply(x,1,mean,trim=1/8)
  sqerror <- (x - rm)^2
  # scale due to trimming of large squares
  1.51 * apply(sqerror,1,mean,trim=1/8)
}

# calculate Cook's distance
dispModelMatrix <- if (modelAsFormula) {
  getModelMatrix(object)
} else {
  modelMatrix
}
attr(object,"dispModelMatrix") <- dispModelMatrix
cooks <- calculateCooksDistance(objectNZ, H, dispModelMatrix)

# this function breaks out the logic for calculating the max Cook's distance:
# the samples over which max Cook's distance is calculated:
#
# Cook's distance is considered for those samples with 3 or more replicates per cell
#
# if m == p or there are no samples over which to calculate max Cook's, then give NA
recordMaxCooks <- function(design, colData, modelMatrix, cooks, numRow) {
  samplesForCooks <- nOrMoreInCell(modelMatrix, n=3)
  p <- ncol(modelMatrix)
  m <- nrow(modelMatrix)
  maxCooks <- if ((m > p) & any(samplesForCooks)) {
    apply(cooks[,samplesForCooks,drop=FALSE], 1, max)
  } else {
    rep(NA, numRow)
  }
  maxCooks
}

# record maximum Cook's
maxCooks <- recordMaxCooks(design(object), colData(object), dispModelMatrix, cooks, nrow(objectNZ))

# store Cook's distance for each sample
assays(object, withDimnames=FALSE)[["cooks"]] <- buildMatrixWithNARows(cooks, mcols(object)$allZero)

# add betas, standard errors and Wald p-values to the object
modelMatrixNames <- colnames(modelMatrix)
betaMatrix <- fit$betaMatrix
colnames(betaMatrix) <- modelMatrixNames
betaSE <- fit$betaSE
colnames(betaSE) <- paste0("SE_",modelMatrixNames)
WaldStatistic <- betaMatrix/betaSE
colnames(WaldStatistic) <- paste0("WaldStatistic_",modelMatrixNames)

modelMatrix
## t distribution for p-values ##
#################################

if (useT) {
  # if the `df` was provided to nbinomWaldTest...
  if (!missing(df)) {
    stopifnot(length(df) == 1 | length(df) == nrow(object))
    if (length(df) == 1) {
      df <- rep(df, nrow(objectNZ))
    } else {
      # the `WaldStatistic` vector is along nonzero rows of `object`
      df <- df[!mcols(object)$allZero]
    }
  } else {
    # df was missing, so compute it from the number of samples (w.r.t. weights)
    # and the number of coefficients
    if ("weights" %in% assayNames(object)) {
      # this checks that weights are OK and normalizes to have rowMax == 1
      # (although this has already happened earlier in estDispGeneEst and estDispMAP...
      wlist <- getAndCheckWeights(objectNZ, dispModelMatrix)
      num.samps <- rowSums(wlist$weights)
    } else {
      num.samps <- rep(ncol(object), nrow(objectNZ))
    }
    df <- num.samps - ncol(dispModelMatrix)
  }
  df <- ifelse(df > 0, df, NA)
  stopifnot(length(df) == nrow(WaldStatistic))
  # use a t distribution to calculate the p-value
  WaldPvalue <- 2*pt(abs(WaldStatistic),df=df,lower.tail=FALSE)
} else {
  # the default DESeq2 p-value: use the standard Normal
  WaldPvalue <- 2*pnorm(abs(WaldStatistic),lower.tail=FALSE)
}
colnames(WaldPvalue) <- paste0("WaldPvalue_",modelMatrixNames)

betaConv <- fit$betaConv

if (any(!betaConv)) {
  if (!quiet) message(paste(sum(!betaConv),"rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest"))
}

mleBetas <- if (betaPrior) {
  matrixToList(mleBetaMatrix)
} else {
  NULL
}

# convenience function for breaking up matrices
# by column and preserving column names
matrixToList <- function(m) {
  l <- split(m, col(m))
  names(l) <- colnames(m)
  l
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

# if useT need to add the t degrees of freedom to the end of resultsList
tDFList <- if (useT) list(tDegreesFreedom=df) else NULL

resultsList <- c(matrixToList(betaMatrix),
                 matrixToList(betaSE),
                 mleBetas,
                 matrixToList(WaldStatistic),
                 matrixToList(WaldPvalue),
                 list(betaConv = betaConv,
                      betaIter = fit$betaIter,
                      deviance = -2 * fit$logLike,
                      maxCooks = maxCooks),
                 tDFList)

WaldResults <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)

modelMatrixNamesSpaces <- gsub("_"," ",modelMatrixNames)

lfcType <- if (attr(object,"betaPrior")) "MAP" else "MLE"
coefInfo <- paste(paste0("log2 fold change (",lfcType,"):"),modelMatrixNamesSpaces)
seInfo <- paste("standard error:",modelMatrixNamesSpaces)
mleInfo <- if (betaPrior) {
  gsub("_"," ",colnames(mleBetaMatrix))
} else {
  NULL
}
statInfo <- paste("Wald statistic:",modelMatrixNamesSpaces)
pvalInfo <- paste("Wald test p-value:",modelMatrixNamesSpaces)

tDFDescription <- if (useT) "t degrees of freedom for Wald test" else NULL  
mcolsWaldResults <- DataFrame(type = rep("results",ncol(WaldResults)),
                              description = c(coefInfo, seInfo, mleInfo, statInfo, pvalInfo,
                                              "convergence of betas",
                                              "iterations for betas",
                                              "deviance for the fitted model",
                                              "maximum Cook's distance for row",
                                              tDFDescription))

mcols(WaldResults) <- mcolsWaldResults

mcols(object) <- cbind(mcols(object),WaldResults)

all_data <- mcols(object)

write.csv(all_data, "all_data_Deseq2_050821.csv")

pvalInfo

mcols(object)

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
