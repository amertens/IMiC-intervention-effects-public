

#Set up SL components
lrnr_mean <- make_learner(Lrnr_mean)

lrnr_glmnet <- Lrnr_glmnet$new()
random_forest <- Lrnr_randomForest$new()
glm_fast <- Lrnr_glm_fast$new()
nnls_lrnr <- Lrnr_nnls$new()

xgboost_lrnr <- Lrnr_xgboost$new()
ranger_lrnr <- Lrnr_ranger$new()
earth_lrnr <- Lrnr_earth$new()
hal_lrnr <- Lrnr_hal9001$new()
polyspline<-Lrnr_polspline$new()
screen_glmnet <- Lrnr_pkg_SuperLearner_screener$new("screen.glmnet")
screened_hal <- make_learner(Pipeline, screen_glmnet, hal_lrnr)



stack <- make_learner(
  Stack,
  lrnr_mean,
  glm_fast,
  #cor_glm,
  lrnr_glmnet,
  ranger_lrnr,
  xgboost_lrnr#,
  #cor_spline#,
  #screened_hal
)


metalearner <- make_learner(Lrnr_nnls)

sl <- make_learner(Lrnr_sl,
                   learners = stack,
                   metalearner = metalearner)




#wrapper function to fit and clean SL
fit_SuperLearner <- function(dat, outcome, covars, slmod=sl, CV=CV_setting, family="binomial"){
  
  covars=covars[covars %in% colnames(dat)]
  
  res <- dat %>% group_by(studyid, visit,contrast) %>% 
    do(fit=try(fit_SL(., outcome = outcome, family=family, covars=covars, slmod=slmod, CV=CV))) %>%
    ungroup() %>% do(clean_res(.))
  
  return(res)
}




fit_SL_fun <- function(dat,  
                       outcome = "haz",
                       id="subjid",
                       family="gaussian",
                       covars,
                       #fit_regression=T, 
                       slmod=sl,
                       CV=TRUE,
                       folds=5){
  
  dat <- dat %>% mutate(arm=ifelse(arm == "Control",0,1))
  
  Y_miss <- sum(is.na(dat[,outcome]))
  if(Y_miss>0){
    dat <- dat[!is.na(dat[,outcome]),]
    cat("Dropping ",Y_miss," missing outcome observations\n")
  }
  
  start_time <- Sys.time()
  full_covars <- covars
  #Drop missing variables or all-NA variables 
  missing_covars <- covars[!(covars %in% colnames(dat))]
  
  
  #Drop near-zero variance predictors
  not_cov <- dat[,!(colnames(dat) %in% covars)]
  cov <- dat[,covars]
  #impute missing
  cov <- cov %>% 
    do(impute_missing_values(., type = "standard", add_indicators = T, prefix = "missing_")$data) %>%
    as.data.frame()
  nzv_cols <- nearZeroVar(cov)
  dropped_covars <-  colnames(cov)[nzv_cols]
  #cat("Dropping for low variance: ", dropped_covars)
  if(length(nzv_cols) > 0){ 
    covars <- covars[!(covars %in% colnames(cov)[nzv_cols])]
    cov <- cov[, -nzv_cols]
  }
  not_cov <- data.frame(not_cov)
  cov <- data.frame(cov)
  dat <-cbind(not_cov,cov)
  dat <- data.table(dat)
  
  
  
  n= nrow(dat)
  y <- as.matrix(dat[, outcome, with=FALSE])
  colnames(y) <- NULL
  
  
  #parallelize
  if(.Platform$OS.type == "windows"){
    # multicore (fork-based) isn't available on Windows; use multisession instead.
    cpus_logical <- parallel::detectCores()
    plan(multisession, workers=floor(cpus_logical/2))
  } else {
    uname <- system("uname -a", intern = TRUE)
    os <- sub(" .*", "", uname)
    if(os=="Darwin"){
      cpus_logical <- as.numeric(system("sysctl -n hw.logicalcpu", intern = TRUE))
    } else if(os=="Linux"){
      cpus_logical <- system("lscpu | grep '^CPU(s)'", intern = TRUE)
      cpus_logical <- as.numeric(gsub("^.*:[[:blank:]]*","", cpus_logical))
    } else {
      stop("unsupported OS")
    }
    plan(multicore, workers=floor(cpus_logical/2))
  }
  
  # create the sl3 task
  set.seed(12345)
  SL_task <- make_sl3_Task(
    data = dat,
    covariates = covars,
    outcome = outcome,
    id="subjid",
    folds = make_folds(dat, fold_fun = folds_vfold, V = folds)
  )
  
  
  ##3. Fit the full model
  
  
  #get cross validated fit
  if(CV){
    suppressMessages(cv_sl <- make_learner(Lrnr_cv, slmod, full_fit = TRUE))
    #cv_glm <- make_learner(Lrnr_cv, glm_sl, full_fit = TRUE)
    suppressMessages(sl_fit <- cv_sl$train(SL_task))
    #glm_fit <- cv_glm$train(SL_task)
  }else{
    suppressMessages(sl_fit <- slmod$train(SL_task))
    #glm_fit <- glm_sl$train(SL_task)
  }
  
  #get outcome predictions
  yhat_full <- sl_fit$predict_fold(SL_task,"validation")
  #yhat_glm <- glm_fit$predict_fold(SL_task,"validation")
  
  #save residuals
  SL_residuals <- y-yhat_full
  #glm_residuals <- y-yhat_glm
  
  
  ##4. Fit the null model
  lrnr_cv_null <- make_learner(Lrnr_cv, make_learner(Lrnr_mean))
  fit_null <- lrnr_cv_null$train(SL_task)
  yhat_null <- fit_null$predict_fold(SL_task,"validation")
  
  if(family=="gaussian"){
    
    mse_full <- 1/n * sum((yhat_full-y)^2)
    mse_glm <- 1/n * sum((yhat_glm-y)^2)
    mse_null <- 1/n * sum((yhat_null-y)^2)
    
    
    
    ##5. Construct CI for mse_full/mse_null
    IC <- 1/mse_null * ((yhat_full-y)^2 - mse_full) - mse_full/(mse_null)^2 * ((yhat_null-y)^2- mse_null)
    IC_glm <- 1/mse_null * ((yhat_glm-y)^2 - mse_glm) - mse_glm/(mse_null)^2 * ((yhat_null-y)^2- mse_null)
    
    psi <- mse_full/mse_null
    se <- sqrt(var(IC)/n)
    CI_up <- psi + 1.96*se
    CI_low <- psi - 1.96*se
    R2 <- 1 - psi
    R2.ci1 <- 1 - CI_up
    R2.ci2 <- 1 - CI_low
    
    glm.psi <- mse_glm/mse_null
    glm.se <- sqrt(var(IC_glm)/n)
    glm.CI_up <- glm.psi + 1.96*glm.se
    glm.CI_low <- glm.psi - 1.96*glm.se
    glm.R2 <- 1 - glm.psi
    glm.R2.ci1 <- 1 - glm.CI_up
    glm.R2.ci2 <- 1 - glm.CI_low
    
    
    R2_full <- data.frame(R2=R2, R2.se=se, R2.ci1=R2.ci1, R2.ci2=R2.ci2)
    R2_glm <- data.frame(R2=glm.R2, R2.se=glm.se, R2.ci1=glm.R2.ci1, R2.ci2=glm.R2.ci2)
    
    ##6. add more results to pool R2
    
    ic_mse_full <- (yhat_full-y)^2 - mse_full
    ic_mse_null <- (yhat_null-y)^2 - mse_null
    
    se_mse_full <- sqrt(var(ic_mse_full)/n)
    se_mse_null <- sqrt(var(ic_mse_null)/n)
    
    mse_ic <- list(ic_mse_full=ic_mse_full,
                   ic_mse_null=ic_mse_null,
                   se_mse_full=se_mse_full,
                   se_mse_null=se_mse_null)
    
    perf_metrics <- list(mse_full = mse_full,
                         mse_glm = mse_glm,
                         mse_null = mse_null,
                         R2_full = R2_full,
                         R2_glm = R2_glm,
                         mse_ic=mse_ic
    )
    
    fold_index <- origami::folds2foldvec(SL_task$folds) 
    
  }else{
    fold_index <- origami::folds2foldvec(SL_task$folds) 
    
    outcome <- SL_task$Y
    
    auc.plotdf <- cvAUC(predictions=yhat_full, labels=outcome, folds=fold_index)
    auc.ci<-ci.cvAUC(predictions=yhat_full, labels=outcome, folds=fold_index, confidence=0.95)
    #glm.auc.ci<-ci.cvAUC(predictions=yhat_glm, labels=outcome, folds=fold_index, confidence=0.95)
    
    perf_metrics <- list(auc.plotdf = auc.plotdf,
                         auc.ci = auc.ci,
                         #glm.auc.ci = glm.auc.ci,
                         fold_index = fold_index)
  }
  
  
  end_time <- Sys.time()
  
  run_time <- end_time - start_time
  

  return(
    list(
      sl_fit = sl_fit,
      covars = list(specified_covars=full_covars,
                    missing_covars=missing_covars, 
                    dropped_covars=dropped_covars,
                    used_covars=covars),
      N_obs=n,
      yhat_full = yhat_full,
      #yhat_glm = yhat_glm,
      Y= SL_task$Y,
      fold_index=fold_index,
      perf_metrics=perf_metrics,
      SL_residuals = SL_residuals,
      #glm_residuals = glm_residuals,
      run_time=run_time,
      Y_miss=Y_miss,
      id=dat$subjid
    )
  )
  #}
}


#wrapper function to detect errors
fit_SL <- purrr::safely(fit_SL_fun)

clean_res <- function(res){
  
  clean_res <- res %>% #ungroup() %>%
    mutate(no_error = fit %>% purrr::map_lgl(.f = ~ is.null(.x$error)))
  return(clean_res)
  
}


#-------------------------------------------
# List extraction function
#-------------------------------------------

extract_results = function(x) x[["result"]]
extract_perf_metrics = function(x) x[["perf_metrics"]]
extract_R2_full = function(x) x[["R2_full"]]
extract_R2_glm = function(x) x[["R2_glm"]]
extract_N_obs = function(x) x[["N_obs"]]

extract_mse_full = function(x) x[["mse_full"]]
extract_mse_null = function(x) x[["mse_null"]]
extract_mse_ic = function(x) x[["mse_ic"]]
extract_se_mse_full = function(x) x[["se_mse_full"]]
extract_se_mse_null = function(x) x[["se_mse_null"]]

extract_yhat_full = function(x) x[["yhat_full"]]
extract_Y = function(x) x[["Y"]]
extract_fold_index = function(x) x[["fold_index"]]
extract_id = function(x) x[["id"]]

extract_AUC_full = function(x) x[["auc.ci"]]
extract_cvAUC = function(x) x[["cvAUC"]]
extract_ci.AUC = function(x) x[["ci"]]


#-------------------------------------------
# Extract AUC
#-------------------------------------------
extract_AUC <- function(res){
  results <- Map(extract_results, res$fit)
  N_obs <- Map(extract_N_obs, results)
  perf_metrics <- Map(extract_perf_metrics, results)
  AUC_full <- Map(extract_AUC_full, perf_metrics)
  cvAUC <- Map(extract_cvAUC, AUC_full)
  ci.AUC <- Map(extract_ci.AUC, AUC_full)
  names(ci.AUC) <- letters[1:length(ci.AUC)]
  ci.AUC <- t(bind_rows(ci.AUC))

  AUC_full <- data.frame(studyid=res$studyid,visit =res$visit ,arm =res$contrast, N_obs=unlist(N_obs, use.names=FALSE), cvAUC=unlist(cvAUC), ci.lb=ci.AUC[,1],  ci.ub=ci.AUC[,2])
  rownames(AUC_full) <- NULL
  
  return(AUC_full)
}

