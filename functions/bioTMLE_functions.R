# =============================================================================
# functions/bioTMLE_functions.R
#
# This script's input/output paths are assembled at runtime, so they
# cannot be listed here without executing it.
#
# Header generated from the code itself; it makes no claim about method.
# See README.md for run order and results/ARTIFACT_MANIFEST.csv for the
# exhibit each script feeds.
# =============================================================================




biomarkertmle_imic <- function(adjust_biomarker=FALSE, include_biomarkers_in_W=FALSE, se, varInt, 
                               normalized = TRUE, ngscounts = FALSE, bppar_type = BiocParallel::MulticoreParam(),
                               bppar_debug = FALSE, cv_folds = 1, 
                               g_lib = c("SL.mean", "SL.glm","SL.bayesglm"),
                               Q_lib = c("SL.mean", "SL.bayesglm", "SL.earth","SL.ranger"), ...){
  
  call <- match.call(expand.dots = TRUE)
  biotmle <- .biotmle(SummarizedExperiment(assays = list(expMeasures = assay(se)), 
                                           rowData = rowData(se), colData = colData(se)), call = call, 
                      tmleOut = tibble::as_tibble(matrix(NA, 10, 10), .name_repair = "minimal"), 
                      topTable = tibble::as_tibble(matrix(NA, 10, 10), .name_repair = "minimal"))
  
  BiocParallel::bpprogressbar(bppar_type) <- TRUE
  BiocParallel::register(bppar_type, default = TRUE)
  

  #Note: scaling now done outside of this function
  #if(!normalized){
  #   exp_normed <- limma::normalizeBetweenArrays(as.matrix(assay(se)),
  #                                               method = "scale")
  #   Y <- tibble::as_tibble(t(exp_normed), .name_repair = "minimal")
  # }else{
     Y <- tibble::as_tibble(t(as.matrix(assay(se))), .name_repair = "minimal")
   #}
  if(!all(apply(Y, 2, class) == "numeric")){
    stop("Warning - values in Y do not appear to be numeric.")
  }
  #A <- as.numeric(SummarizedExperiment::colData(se)[, varInt])
  A <- (SummarizedExperiment::colData(se)[, varInt])
  
  W <- tibble::as_tibble(SummarizedExperiment::colData(se)[, -varInt], .name_repair = "minimal")
  if (is.null(dim(W)[2])) {
    W <- as.numeric(rep(1, length(A)))
  }
  if(!all(is.numeric(apply(W, 2, class)))) {
    require(caret)
    W <- droplevels(W)
    if(length(nearZeroVar(W))>0 & ncol(W)>1){
      W<-W[,-nearZeroVar(W)]
    }
    W <- design_matrix(as.data.frame(W))
    W <- tibble::as_tibble(apply(W, 2, as.numeric), .name_repair = "minimal")
  }
  
  if(!bppar_debug){
    
    if(!adjust_biomarker){
      biomarkertmle_out <- BiocParallel::bplapply(Y[, seq_along(Y)], 
                                                  exp_biomarkertmle_imic, 
                                                  W = W, A = A, 
                                                  g_lib = g_lib, Q_lib = Q_lib, 
                                                  cv_folds = cv_folds)
    }else{
      biomarkertmle_out <- BiocParallel::bplapply(seq_along(Y), function(i){
        # Current column of Y being processed
        current_Y <- Y[, i, drop = FALSE]
        
        # Exclude the current column from Y for adjustment variables
        Y_excluded <- Y[, -i, drop = FALSE]
        
        # Combine W with the modified Y (excluding the current column)
        W_modified <- cbind(W, Y_excluded)
        
        # Call the exp_biomarkertmle function with the modified arguments
        exp_biomarkertmle_imic(Y = current_Y, W = W_modified, A = A, g_lib = g_lib, 
                               Q_lib = Q_lib, cv_folds = cv_folds)
      })
    }
    
    
  }else{
    

    if(!adjust_biomarker){
      biomarkertmle_out <- lapply(Y[, seq_along(Y)], exp_biomarkertmle_imic, 
                                  W = W, A = A, g_lib = g_lib, Q_lib = Q_lib, cv_folds = cv_folds)
      
    }else{
      
      biomarkertmle_out <- lapply(seq_along(Y), function(i) {
        # Current column of Y being processed
        current_Y <- Y[, i, drop = FALSE]
        
        # Exclude the current column from Y for adjustment variables
        Y_excluded <- Y[, -i, drop = FALSE]
        
        # Combine W with the modified Y (excluding the current column)
        W_modified <- cbind(W, Y_excluded)
        
        # Call the exp_biomarkertmle function with the modified arguments
        exp_biomarkertmle_imic(Y = current_Y, W = W_modified, A = A, g_lib = g_lib, 
                               Q_lib = Q_lib, cv_folds = cv_folds)
      })
    }
  }
  
  #extract lists
  biomarkertmle_params <- do.call(c, lapply(biomarkertmle_out, `[[`, "param"))
  biomarkertmle_eifs <- do.call(cbind.data.frame, lapply(biomarkertmle_out, `[[`, "eif"))
  biomarkertmle_res <- do.call(rbind.data.frame, lapply(biomarkertmle_out, `[[`, "res"))
  
  biotmle@ateOut <- as.numeric(biomarkertmle_params)
  
  biomarker_eifs <- t(as.matrix(biomarkertmle_eifs))
  colnames(biomarker_eifs) <- colnames(se)
  biotmle@tmleOut <- tibble::as_tibble(biomarker_eifs, .name_repair = "minimal")
  
  #format TMLE results
  biomarkertmle_res$biomarker <- rep(colnames(Y), each=2*length(unique(A))-1)
  rownames(biomarkertmle_res) = NULL
  #pvalue correct heterogeneity test Pval
  
  #try(biomarkertmle_res <- p.adjust_QEp(biomarkertmle_res))
  try(biomarkertmle_res <- p.adjust_chi_p(biomarkertmle_res))
  
  
  output = list(biotmle=biotmle, res=biomarkertmle_res)
  
  return(output)
}


# Yfull=Y
# Y=Yfull[,17]

exp_biomarkertmle_imic <- function(Y, A, W, g_lib, Q_lib, cv_folds, ...){
  
  if(any(class(Y) == "data.frame")){
    Y <- as.numeric(unlist(Y[, 1]))
  } 
  if(any(class(A) == "data.frame")){
    A <- as.numeric(unlist(A[, 1])) #convert to factor?
  }
  
  missingY <- is.na(Y)
  if(sum(missingY)>0){
    A <- A[!missingY]
    if(nrow(W)>0){
      if(ncol(W)>1){
        W <- W[!missingY,]
      }
      if(ncol(W)==1){
        W <- W[!missingY,1]
      }
    }
    Y <- Y[!missingY]
  }
  
  
  if(length(nzv(Y))==0){
    
    assertthat::assert_that(length(unique(A)) > 1)
    if(is.factor(A)){
      A<-droplevels(A)
      a_0 <- levels(A[!is.na(A)])
    }else{
      a_0 <- sort(unique(A[!is.na(A)]))
    }
    
    
    suppressWarnings(tmle_fit <- drtmle::drtmle(Y = Y, A = A,  family = stats::gaussian(),
                                                W = W, a_0 = a_0, SL_g = g_lib, SL_Q = Q_lib, cvFolds = cv_folds, 
                                                stratify = TRUE, guard = NULL, parallel = FALSE, use_future = FALSE))
    res=imic_drtmle_summary(tmle_fit, a_0)
    ate_tmle <- tmle_fit$tmle$est[seq_along(a_0)[-1]] - tmle_fit$tmle$est[1]
    eif_tmle_delta <- tmle_fit$ic$ic[, seq_along(a_0)[-1]] - tmle_fit$ic$ic[, 1]
    
    if(!is.vector(eif_tmle_delta)){
      param_out <- ate_tmle[length(ate_tmle)]
      eif_out <- eif_tmle_delta[, ncol(eif_tmle_delta)] + ate_tmle[length(ate_tmle)]
    }else{
      param_out <- ate_tmle
      eif_out <- eif_tmle_delta + ate_tmle
    }
    
    #add missingness back into the eif vector so they can be merged across outcomes
    if(sum(missingY)>0){
      # Step 3: Create a new vector of the same length as the original, filled with NAs
      eif_NA_vector <- rep(NA, length(missingY))
      
      # Step 4: Insert calculated values back into their original positions
      eif_NA_vector[!missingY] <- eif_out
      eif_out <- eif_NA_vector
    }else{
      out <- list(param = NA, eif = NA, res=NA)
    }
    
    assertthat::assert_that(is.vector(eif_out))
    out <- list(param = param_out, eif = eif_out, res=res)
  }else{
    out <- list(param = NA, eif = NA, res=NA)
  }
  return(out)
}




#function to get global p-val from test of heterogeneity, and summarize means and treatment contrasts
imic_drtmle_summary <- function(fit, a_0){
  

  # Calculate the variance-covariance matrix of the IC and the ATE
  IC_Ya=fit$ic_drtmle$ic
  K=ncol(IC_Ya)-1
  IC = matrix(NA, nrow=nrow(IC_Ya), ncol = K)
  ate = rep(NA, K)
  
    for(i in 1:K){
      IC[,i] = IC_Ya[,i+1]-IC_Ya[,1]
      ate[i] = fit$drtmle$est[i+1]-fit$drtmle$est[1]
    }
  
  # Step 4: Compute the empirical variance-covariance matrix
  vmat <- var(IC)/nrow(IC)
  
  # Calculate the inverse of the variance-covariance matrix
  inv_vmat <- solve(vmat)
  
  # Calculate the square root of the inverse matrix
  # Note: You might need to install the `Matrix` package for sqrtm function
  # install.packages("Matrix")
  #library(Matrix)
  require(pracma)
  sdIC <- sqrtm(inv_vmat)$B
  
  # Assuming that Zvec is already in the form where each row corresponds to a transformed IC component,
  # you can then compute the chi-square test statistic as the sum of the squared Z components
  # This step might vary depending on your specific transformation and data structure
  
  # Calculate the chi-square test statistic
  zstats <- ((sdIC%*%ate)^2)
  chi_square_statistic <- sum(zstats)
  
  chi_pval<-pchisq(chi_square_statistic, K,lower.tail=F)
  
  res_tab <- data.frame(contrast=a_0, ci(fit)$drtmle, measure="MN")
  
  if(length(unique(a_0))==4){
    res_ATE= bind_rows(
      data.frame( ci(fit, contrast =c(-1,1,0,0))$drtmle, measure="ATE"),
      data.frame( ci(fit, contrast =c(-1,0,1,0))$drtmle, measure="ATE"),
      data.frame( ci(fit, contrast =c(-1,0,0,1))$drtmle, measure="ATE")
    )
    res_ATE$contrast = rownames(res_ATE)
  }
  
  if(length(unique(a_0))==3){
    res_ATE= bind_rows(
      data.frame( ci(fit, contrast =c(-1,1,0))$drtmle, measure="ATE"),
      data.frame( ci(fit, contrast =c(-1,0,1))$drtmle, measure="ATE")
    )
    res_ATE$contrast = rownames(res_ATE)
  }
  
  if(length(unique(a_0))==2){
    res_ATE= bind_rows(
      data.frame( ci(fit, contrast =c(-1,1))$drtmle, measure="ATE"),
    )
    res_ATE$contrast = rownames(res_ATE)
  }
  
    res_tab$contrast=as.character(res_tab$contrast)
    res_ATE$contrast=as.character(res_ATE$contrast)
  res_tab <- bind_rows(res_tab, res_ATE)
  rownames(res_tab)=NULL
  #res_tab$QEp <- QEp
  res_tab$chi_pval <- chi_pval
  
  #clean up intervention labels
  res_tab$contrast <- gsub("E\\[Y\\((.*?)\\)\\]-E\\[Y\\(Control\\)\\]", "\\1", res_tab$contrast)
  
  return(res_tab)
}



p.adjust_chi_p <- function(res){
  resP <- res %>% distinct(chi_pval , biomarker)
  resP$chi_pval_adj <- p.adjust(resP$chi_pval , method="BH") 
  #resP[resP$QEp_adj<0.05,]
  res <- left_join(res, resP, by = c("chi_pval", "biomarker"))
  #res[res$QEp_adj<0.05,] 
  res <- res %>% arrange(chi_pval_adj)
  res$sig <- ifelse(res$chi_pval<0.05, 1, 0)
  res$sigFDR <- ifelse(res$chi_pval_adj<0.05, 1, 0)
  
  return(res)
}



run_bioTMLE <- function(d, Wvars, g_lib = c("SL.glm"),Q_lib = c("SL.glm"),
                        Yvars=c(all_milk_components$macro, all_milk_components$micro, all_milk_components$bvit),
                        scale=FALSE, cv.folds=1, bppar.type = BiocParallel::SerialParam(), bppar.debug=FALSE, adjust_biomarker=FALSE,
                        seed = 12345){

  # Fixed seed for reproducibility: the SuperLearner library (glmnet CV folds,
  # xgboost) is stochastic, which made borderline FDR calls (e.g., Mumta-LW
  # vitamin A at 2 mo, q ~ 0.05) flip between runs. Seeding per study x visit
  # group makes every run deterministic across machines (serial execution).
  set.seed(seed)

  d <- droplevels(d)
  confounderData <- d %>% select(all_of(Wvars))
  biomarkerData <- d %>% select(all_of(Yvars))
  
  #drop NZV columns
  if(colnames(confounderData)[2]!="dummy"){ #skip if unadjusted analysis
    #confounderData <- predict(preProcess(confounderData, method = c("nzv")),confounderData)
    
    if(length(nearZeroVar(confounderData))>0){
      confounderData<-confounderData[,-nearZeroVar(confounderData)]
    }
    ARM <- confounderData[,1]
    confounderData <- cbind(ARM, design_matrix(as.data.frame(confounderData[,-1])))
  }
  #biomarkerData <- predict(preProcess(biomarkerData, method = c("nzv")),biomarkerData)   %>% as.matrix()
  if(length(nearZeroVar(biomarkerData))>0){
    biomarkerData<-biomarkerData[,-nearZeroVar(biomarkerData)]
  }
  
  if(scale){
    biomarkerData = as.matrix(as.data.frame(scale_continuous(biomarkerData)))
  }
 
  d_se <- SummarizedExperiment(assays = list(hm_components = t(biomarkerData)),
                               colData = DataFrame(confounderData))
  
  biotmle_out <- biomarkertmle_imic(adjust_biomarker=adjust_biomarker,
                                    se = d_se,
                                    varInt = 1,
                                    normalized=scale,
                                    g_lib = g_lib,
                                    Q_lib = Q_lib,
                                    cv_folds = cv.folds, #temp
                                    bppar_debug=bppar.debug,
                                    #bppar_type = BiocParallel::MulticoreParam())
                                    #bppar_type = BiocParallel::SnowParam())
                                    bppar_type = bppar.type)  
  
  return(biotmle_out)
}


optWeight_imic <- function (Y, X, SL.library, family = "gaussian", CV.SuperLearner.V = 10, 
                            seed = 12345, whichAlgorithm = "SuperLearner", return.SuperLearner = TRUE, 
                            return.CV.SuperLearner = FALSE, return.IC = TRUE, parallel = FALSE, 
                            n.cores = parallel::detectCores(), ...) 
{
  n <- length(Y[, 1])
  J <- ncol(Y)
  Ymat <- data.matrix(Y)
  if (is.null(colnames(Ymat))) {
    colnames(Ymat) <- paste0("Y", 1:J)
  }
  CV.SuperLearner.list <- apply(Ymat, 2, function(y) {
    set.seed(seed)
    if (parallel) {
      options(mc.cores = n.cores)
    }
    fit <- SuperLearner::CV.SuperLearner(Y = y, X = X, SL.library = SL.library, 
                                         family = family, V = CV.SuperLearner.V, parallel = ifelse(parallel, 
                                                                                                   "multicore", "seq"), method = "method.NNLS")
  })
  psiHat.Pnv0 <- r2weight:::getPredictionsOnValidation(out = CV.SuperLearner.list, 
                                                       whichAlgorithm = whichAlgorithm)
  if (!is.matrix(psiHat.Pnv0)) {
    psiHat.Pnv0 <- matrix(psiHat.Pnv0, ncol = J)
  }
  univariateResults <- r2weight:::getUnivariateR2(Y = Ymat, psiHat.Pnv0 = psiHat.Pnv0, 
                                                  return.IC = return.IC)
  alpha_n <- r2weight:::alphaHat(Y = Ymat, psiHat.Pnv0 = psiHat.Pnv0)
  SuperLearner.list <- apply(Ymat, 2, function(y) {
    set.seed(seed)
    fit <- SuperLearner::SuperLearner(Y = y, X = X, SL.library = SL.library, 
                                      family = family, method = "method.NNLS")
  })
  out <- vector(mode = "list")
  out$SL.fits <- NULL
  if (return.SuperLearner) {
    out$SL.fits <- SuperLearner.list
  }
  out$SL.weights <- alpha_n
  out$SL.library <- SL.library
  out$CV.SL.fits <- NULL
  out$whichAlgorithm <- whichAlgorithm
  out$CV.SuperLearner.V <- CV.SuperLearner.V
  out$family <- family
  if (return.CV.SuperLearner) {
    out$CV.SL.fits <- CV.SuperLearner.list
  }
  out$univariateR2 <- univariateResults[colnames(Ymat)]
  out$IC <- NULL
  if (return.IC) {
    out$IC <- univariateResults[["IC"]]
  }
  out$MSE <- univariateResults[["MSE"]]
  out$Var <- univariateResults[["Var"]]
  out$Ynames <- colnames(Ymat)
  class(out) <- "optWeight"
  return(out)
}

