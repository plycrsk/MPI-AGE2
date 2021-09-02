function (object, test = c("Wald", "LRT"), fitType = c("parametric", 
                                                       "local", "mean", "glmGamPoi"), sfType = c("ratio", "poscounts", 
                                                                                                 "iterate"), betaPrior, full = design(object), reduced, quiet = FALSE, 
          minReplicatesForReplace = 7, modelMatrixType, useT = FALSE, 
          minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5, parallel = FALSE, 
          BPPARAM = bpparam()) 
{
  stopifnot(is(object, "DESeqDataSet"))
  test <- match.arg(test, choices = c("Wald", "LRT"))
  fitType <- match.arg(fitType, choices = c("parametric", "local", 
                                            "mean", "glmGamPoi"))
  dispersionEstimator <- if (fitType == "glmGamPoi") {
    "glmGamPoi"
  }
  else {
    "DESeq2"
  }
  if (fitType == "glmGamPoi") {
    minReplicatesForReplace <- Inf
    if (parallel) {
      warning("parallelization of DESeq() is not implemented for fitType='glmGamPoi'")
    }
  }
  sfType <- match.arg(sfType, choices = c("ratio", "poscounts", 
                                          "iterate"))
  stopifnot(is.logical(quiet))
  stopifnot(is.numeric(minReplicatesForReplace))
  stopifnot(is.logical(parallel))
  modelAsFormula <- !is.matrix(full) & is(design(object), "formula")
  if (missing(betaPrior)) {
    betaPrior <- FALSE
  }
  else {
    stopifnot(is.logical(betaPrior))
  }
  object <- sanitizeRowRanges(object)
  if (test == "LRT") {
    if (missing(reduced)) {
      stop("likelihood ratio test requires a 'reduced' design, see ?DESeq")
    }
    if (betaPrior) {
      stop("test='LRT' does not support use of LFC shrinkage, use betaPrior=FALSE")
    }
    if (!missing(modelMatrixType) && modelMatrixType == "expanded") {
      stop("test='LRT' does not support use of expanded model matrix")
    }
    if (is.matrix(full) | is.matrix(reduced)) {
      if (!(is.matrix(full) & is.matrix(reduced))) {
        stop("if one of 'full' and 'reduced' is a matrix, the other must be also a matrix")
      }
    }
    if (modelAsFormula) {
      checkLRT(full, reduced)
    }
    else {
      checkFullRank(full)
      checkFullRank(reduced)
      if (ncol(full) <= ncol(reduced)) {
        stop("the number of columns of 'full' should be more than the number of columns of 'reduced'")
      }
    }
  }
  if (test == "Wald" & !missing(reduced)) {
    stop("'reduced' ignored when test='Wald'")
  }
  if (dispersionEstimator == "glmGamPoi" && test == "Wald") {
    warning("glmGamPoi dispersion estimator should be used in combination with a LRT and not a Wald test.", 
            call. = FALSE)
  }
  if (modelAsFormula) {
    designAndArgChecker(object, betaPrior)
    if (design(object) == formula(~1)) {
      warning("the design is ~ 1 (just an intercept). is this intended?")
    }
    if (full != design(object)) {
      stop("'full' specified as formula should equal design(object)")
    }
    modelMatrix <- NULL
  }
  else {
    if (!quiet) 
      message("using supplied model matrix")
    if (betaPrior == TRUE) {
      stop("betaPrior=TRUE is not supported for user-provided model matrices")
    }
    checkFullRank(full)
    modelMatrix <- full
  }
  attr(object, "betaPrior") <- betaPrior
  stopifnot(length(parallel) == 1 & is.logical(parallel))
  if (!is.null(sizeFactors(object)) || !is.null(normalizationFactors(object))) {
    if (!quiet) {
      if (!is.null(normalizationFactors(object))) {
        message("using pre-existing normalization factors")
      }
      else {
        message("using pre-existing size factors")
      }
    }
  }
  else {
    if (!quiet) 
      message("estimating size factors")
    object <- estimateSizeFactors(object, type = sfType, 
                                  quiet = quiet)
  }
  if (!parallel) {
    if (!quiet) 
      message("estimating dispersions")
    object <- estimateDispersions(object, fitType = fitType, 
                                  quiet = quiet, modelMatrix = modelMatrix, minmu = minmu)
    if (!quiet) 
      message("fitting model and testing")
    if (test == "Wald") {
      object <- nbinomWaldTest(object, betaPrior = betaPrior, 
                               quiet = quiet, modelMatrix = modelMatrix, modelMatrixType = modelMatrixType, 
                               useT = useT, minmu = minmu)
    }
    else if (test == "LRT") {
      object <- nbinomLRT(object, full = full, reduced = reduced, 
                          quiet = quiet, minmu = minmu, type = dispersionEstimator)
    }
  }
  else if (parallel) {
    if (!missing(modelMatrixType)) {
      if (betaPrior) 
        stopifnot(modelMatrixType == "expanded")
    }
    object <- DESeqParallel(object, test = test, fitType = fitType, 
                            betaPrior = betaPrior, full = full, reduced = reduced, 
                            quiet = quiet, modelMatrix = modelMatrix, useT = useT, 
                            minmu = minmu, BPPARAM = BPPARAM)
  }
  sufficientReps <- any(nOrMoreInCell(attr(object, "modelMatrix"), 
                                      minReplicatesForReplace))
  if (sufficientReps) {
    object <- refitWithoutOutliers(object, test = test, betaPrior = betaPrior, 
                                   full = full, reduced = reduced, quiet = quiet, minReplicatesForReplace = minReplicatesForReplace, 
                                   modelMatrix = modelMatrix, modelMatrixType = modelMatrixType)
  }
  metadata(object)[["version"]] <- packageVersion("DESeq2")
  object
}
<bytecode: 0x563a5df79ad0>
  <environment: namespace:DESeq2>