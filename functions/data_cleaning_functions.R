

#get labels from SAS file
makeVlist <- function(dta) { 
  labels <- sapply(dta, function(x) attr(x, "label"))
  labs <- tibble(name = names(labels),
                 label = labels)
  labs$label<- as.character(labs$label)
  return(labs)
}



assetPCA<-function(ret, reorder=F){
 
  #Select assets
  ret<-as.data.frame(ret) 
  id<-subset(ret, select=c("subjido")) #drop subjectid
  ret<-ret[,which(!(colnames(ret) %in% c("subjido")))]
  
  #drop rows with no asset data
  id<-id[rowSums(is.na(ret[,2:ncol(ret)])) != ncol(ret)-1,]  
  ret<-ret[rowSums(is.na(ret[,2:ncol(ret)])) != ncol(ret)-1,]  
  
  #Drop assets with great missingness
  for(i in 1:ncol(ret)){
    cat(colnames(ret)[i],"\n")
    print(table(is.na(ret[,i])))
    print(class((ret[,i])))
  }
  
  # #Set missingness to zero
  # table(is.na(ret))
  # for(i in 1:ncol(ret)){
  #   ret[,i]<-as.character(ret[,i])
  #   ret[is.na(ret[,i]),i]<-"miss"
  #   ret[,i]<-as.factor(ret[,i])
  #   
  #}
  table(is.na(ret))
  
  #if income is present, median impute
  if("INCTOT" %in% colnames(ret)){
    ret[["INCTOT"]] <- as.numeric(as.character(ret[["INCTOT"]]))
    ret[["INCTOT"]][is.na( ret[["INCTOT"]])] <- median(ret[["INCTOT"]], na.rm=T)
  }
  
  # #Remove columns with almost no variance
  # if(length(nearZeroVar(ret))>0){
  #   ret<-ret[,-nearZeroVar(ret)]
  # }
  # 
  #Convert factors into indicators
  ret<-droplevels(ret)
  ret<-design_matrix(ret)
  # if(length(nearZeroVar(ret))>0){
  #   ret<-ret[,-nearZeroVar(ret)]
  # }
  
  #Set missingness to zero
  table(is.na(ret))
  ret[is.na(ret)]<-0
  table(is.na(ret))
  
  # #Remove columns with almost no variance
  # if(length(nearZeroVar(ret))>0){
  #   ret<-ret[,-nearZeroVar(ret)]
  # }
  
  ## Convert the data into matrix ##
  ret<-as.matrix(ret)
  
  ##Computing the principal component using eigenvalue decomposition ##
  princ.return <- princomp(ret) 
  
  ## To get the first principal component in a variable ##
  load <- loadings(princ.return)[,1]   
  
  pr.cp <- ret %*% load  ## Matrix multiplication of the input data with the loading for the 1st PC gives us the 1st PC in matrix form. 
  
  hhwealth <- scale(as.numeric(pr.cp)) ## Gives us the scaled 1st PC in numeric form in pr.
  
  #Create 4-level household weath index
  quartiles<-quantile(hhwealth, probs=seq(0, 1, 0.25))
  print(quartiles)
  ret<-as.data.frame(ret)
  ret$hhwealth <- hhwealth
  ret$hhwealth_quart<-rep(1, nrow(ret))
  ret$hhwealth_quart[hhwealth>=quartiles[2]]<-2
  ret$hhwealth_quart[hhwealth>=quartiles[3]]<-3
  ret$hhwealth_quart[hhwealth>=quartiles[4]]<-4
  if(length(unique(ret$hhwealth_quart))==3){ #make sure there are 4 levels if the 0 level is high
    quartiles<-quantile(hhwealth, probs=seq(0, 1, 0.2))
    print(quartiles)
    ret<-as.data.frame(ret)
    ret$hhwealth_quart<-rep(1, nrow(ret))
    ret$hhwealth_quart[hhwealth>=quartiles[3]]<-2
    ret$hhwealth_quart[hhwealth>=quartiles[4]]<-3
    ret$hhwealth_quart[hhwealth>=quartiles[5]]<-4
  }
  table(ret$hhwealth_quart)
  ret$hhwealth_quart<-factor(ret$hhwealth_quart)
  
  if(reorder==T){
    levels(ret$hhwealth_quart)<-c("Wealth Q4","Wealth Q3","Wealth Q2","Wealth Q1")
    ret$hhwealth_quart<-factor(ret$hhwealth_quart, levels=c("Wealth Q1", "Wealth Q2","Wealth Q3","Wealth Q4"))
  }else{
    levels(ret$hhwealth_quart)<-c("Wealth Q1", "Wealth Q2","Wealth Q3","Wealth Q4")
  }
  
  #Table assets by pca quartile to identify wealth/poverty levels
  d<-data.frame(id, ret)
  wealth.tab <- d %>% subset(., select=-c(id)) %>%
    group_by(hhwealth_quart) %>%
    summarise_all(funs(mean)) %>% as.data.frame()
  print(wealth.tab)
  
  #Save just the wealth data
  pca.wealth <- d %>% subset(select=c(id, hhwealth, hhwealth_quart)) %>% rename(subjido=id)
  pca.wealth$hhwealth<-as.numeric(pca.wealth$hhwealth)
  #pca.wealth$SUBJID<-as.numeric(as.character(pca.wealth$SUBJID))
  
  # d <-dfull %>% subset(., select=c("subjido"))
  # #d$SUBJID<-as.numeric(as.character(d$SUBJID))
  # d<-left_join(d, pca.wealth, by=c("subjido"))
  return(pca.wealth)
}



#function to convert all asset variables to indicators
design_matrix <- function(W){
  if (class(W) != "matrix" & class(W) != "data.frame") {
    W <- data.frame(W)
    if (is.null(ncol(W)) | ncol(W) == 0) {
      stop("Something is wrong with W.\nTo be safe, please try specifying it as class=data.frame.")
    }
  }
  ncolW <- ncol(W)
  flist <- numeric()
  for (i in 1:ncolW) {
    if (class(W[, i]) != "factor") {
      next
    }
    else {
      flist <- c(flist, i)
      W[, i] <- factor(W[, i])
      mm <- model.matrix(~-1 + W[, i])
      mW <- mm[, -c(1)]
      levs <- gsub(" ", "", levels(W[, i]))[-c(1)]
      if (length(levs) < 2) 
        mW <- matrix(mW, ncol = 1)
      colnames(mW) <- paste(names(W)[i], levs, sep = "")
      W <- data.frame(W, mW)
    }
  }
  if (length(flist) > 0) {
    W <- subset(W, select = -c(flist))
  }
  return(W)
}



#Results extraction function

extract_res = function(x) x[["res"]]

extract_bioTMLE_results <- function(res_list, single_group=FALSE){
  # FDR (pval_adj / sigFDR) is applied by study x visit (studytime) within each outcome
  # group, matching the pre-specified scheme reported in the manuscript Methods.
  # pval_adj_global retains the more conservative correction pooled across study and visit
  # within each outcome group (kept only for sensitivity / change-tracking; not the primary).

  if(single_group){
    res <- Map(extract_res, res_list$res)
    res <- rbindlist(res, idcol='studytime') %>% 
      as.data.frame()
    
    res$pval <- ci_to_pvalue(cil = res$cil, 
                 ciu = res$ciu)
    res = res %>% group_by(studytime, measure) %>% mutate(pval_adj=p.adjust(pval , method="BH") ) %>% group_by(measure) %>% mutate(pval_adj_global=p.adjust(pval , method="BH") ) %>% ungroup()
    
    res$sig = 1*(res$pval < 0.05)
    res$sigFDR = 1*(res$pval_adj < 0.05)
    
  }else{
    res_primary <- Map(extract_res, res_list$res_primary$res)
    res_primary <- rbindlist(res_primary, idcol='studytime') %>% 
      mutate(outcome_group='primary') %>% 
      as.data.frame()
    
    res_primary$pval <- ci_to_pvalue(cil = res_primary$cil, 
                             ciu = res_primary$ciu)
    res_primary = res_primary %>% group_by(studytime, measure) %>% mutate(pval_adj=p.adjust(pval , method="BH") ) %>% group_by(measure) %>% mutate(pval_adj_global=p.adjust(pval , method="BH") ) %>% ungroup()
    

    res_secondary <- Map(extract_res, res_list$res_secondary$res)
    res_secondary <- rbindlist(res_secondary, idcol='studytime') %>% 
      mutate(outcome_group='secondary') %>% 
      as.data.frame()
    
    res_secondary$pval <- ci_to_pvalue(cil = res_secondary$cil, 
                                     ciu = res_secondary$ciu)
    res_secondary = res_secondary %>% group_by(studytime, measure) %>% mutate(pval_adj=p.adjust(pval , method="BH") ) %>% group_by(measure) %>% mutate(pval_adj_global=p.adjust(pval , method="BH") ) %>% ungroup()

    res_tertiary <- Map(extract_res, res_list$res_tertiary$res)
    res_tertiary <- rbindlist(res_tertiary, idcol='studytime') %>% 
      mutate(outcome_group='tertiary') %>% 
      as.data.frame()
    
    res_tertiary$pval <- ci_to_pvalue(cil = res_tertiary$cil, 
                                       ciu = res_tertiary$ciu)
    res_tertiary = res_tertiary %>% group_by(studytime, measure) %>% mutate(pval_adj=p.adjust(pval , method="BH") ) %>% group_by(measure) %>% mutate(pval_adj_global=p.adjust(pval , method="BH") ) %>% ungroup()
    
    res= bind_rows(res_primary, res_secondary, res_tertiary)
    res$chi_sig = res$sig
    res$chi_sigFDR = res$sigFDR
    res$sig = 1*(res$pval < 0.05)
    res$sigFDR = 1*(res$pval_adj < 0.05)

  }

  
  res <- res %>% mutate(
    biomarker= case_when(
      biomarker=="carbohydrate" ~ "Total Carbohydrate",
      biomarker==biomarker ~ biomarker))
  
  unique(res$studytime)
  res <- res %>% mutate(
    studytime=as.character(studytime),
    studytime = case_when(
      studytime=="Misame-1" ~ "Misame (14-21 days)",
      studytime=="Misame-2" ~ "Misame (1-2 mo.)",
      studytime=="Misame-3" ~ "Misame (3-4 mo.)",
      studytime=="Vital-40" ~ "Vital (1.5 mo.)",
      studytime=="Vital-56" ~ "Vital (2 mo.)",
      studytime=="Elicit-1" ~ "Elicit (1 mo.)",
      studytime=="Elicit-5" ~ "Elicit (5 mo.)",
      studytime=="Misame-pooled" ~ "Misame-pooled",
      studytime=="Vital-pooled" ~ "Vital-pooled",
      studytime=="Elicit-pooled" ~ "Elicit-pooled",
    ),
    studytime=factor(studytime, levels=c("Misame (14-21 days)", "Misame (1-2 mo.)", "Misame (3-4 mo.)",
                                         "Vital (1.5 mo.)", "Vital (2 mo.)", "Elicit (1 mo.)", "Elicit (5 mo.)",
                                         "Misame-pooled","Vital-pooled","Elicit-pooled"))
  ) %>% droplevels()
  
  
  if(grepl("pooled",res$studytime[1])){
    res$study <- str_split(res$studytime, '-', simplify = T)[,1]
    res$visit <- str_split(res$studytime, '-', simplify = T)[,2]
  }else{
    res$study <- str_split(res$studytime, ' \\(', simplify = T)[,1]
    res$visit <- str_split(res$studytime, ' \\(', simplify = T)[,ncol(str_split(res$studytime, ' \\(', simplify = T))]
    res$visit <- gsub("\\)","",res$visit)
  }
  
  res <- clean_biomarker_labels(res)
  
  return(res)
}



# Function to calculate p-value from confidence interval to get individual intervention effects
ci_to_pvalue <- function(cil, ciu, null_value = 0, conf_level = 0.95) {
  # Calculate the estimate (midpoint of the CI)
  estimate <- (cil + ciu) / 2
  
  # Calculate standard error
  z_critical <- qnorm(1 - (1 - conf_level) / 2)  # For 95% CI, z_critical = 1.96
  se <- (ciu - cil) / (2 * z_critical)
  
  # Calculate z-statistic
  z <- (estimate - null_value) / se
  
  # Calculate two-sided p-value
  p_value <- 2 * pnorm(-abs(z))
  
  return(p_value)
}

#-------------------------------------------------------------------------------
# biomarker labels
#-------------------------------------------------------------------------------

clean_biomarker_labels <- function(res){
  
  labels <- read.csv(paste0(here::here(),"/metadata/Milk_Component_Spec_IMiC_V02.csv")) %>%
    mutate(biomarker=tolower(VarName)) %>% select(biomarker, VARLABEL, VARDESCRIPTION, VarSubCategory) %>% 
    rename("label"="VARLABEL", "category"="VarSubCategory", "description"="VARDESCRIPTION") %>% 
    distinct(biomarker, .keep_all = TRUE)
  labels$biomarker[labels$biomarker=="ca"] <- "ca_bio"
  labels$category_raw <- labels$category
  #add in missing biomarkers
  extra_labels <-data.frame(biomarker=c("Total Carbohydrate", 
                                        "as", "ca", "cr", "cu",
                                        "fe", "fgf.21", "iga", "k", "mg",
                                        "mn", "mo", "na", "p", "se",
                                        "trp_bio", "zn"), 
                            
                            label=c("Total Carbohydrate","Arsenic", "Calcium", "Chromium", "Copper",
                                    "Iron", "fgf.21", "iga", "Potassium", "Magnesium",
                                    "Manganese", "Molybdenum", "Sodium", "Phosphorus", "Selenium",
                                    "trp_bio", "Zinc"), 
                            category=c("Macronutrient", "Micronutrient", "Micronutrient", "Micronutrient", "Micronutrient",
                                       "Micronutrient", "fgf.21", "iga", "Micronutrient", "Micronutrient",
                                       "Micronutrient", "Micronutrient", "Micronutrient", "Micronutrient", "Micronutrient",
                                       "trp_bio", "Micronutrient")) 
  
  labels <- bind_rows(labels, extra_labels)
  
  labels <- labels %>% mutate(
    category= case_when(
      biomarker %in% all_milk_components$macro ~ "Macronutrient",
      biomarker %in% all_milk_components$micro ~ "Micronutrient",
      category== category ~  category)
  ) 
  
  
  labels$category[labels$category %in% c("Individual HMO")] <- "Other individual HMO"
  labels$category[labels$category %in% c("B3","B3 (catabolite)","B3 (precursor)")] <- "B3 or related"
  labels$category[labels$category %in% c("B12","B5","B7")] <- "Other B vitamins"
  labels$category[labels$category %in% c("trp_bio","Vitamins and Cofactors","Nucleobases and Related","Hormones and Related",
                                         "Alkaloids","Amine Oxides","Cresols","Carbohydrates and Related")] <- "Other metabolomics"
  # fgf.21/iga: the manuscript's Methods text names both as example "selected bioactive
  # proteins" alongside secretory IgA, calprotectin, leptin, insulin -- category them to
  # match (they previously fell through to a placeholder self-named category "fgf.21"/"iga"
  # from the extra_labels patch above, then got miscoded to "Other individual HMO", which
  # contradicted the manuscript text).
  labels$category[labels$biomarker %in% c("secretor","sum_nmol.ml",  "sia_nmol.ml", "insulin",
                                          "leptin","calprotectin","fuc_nmol.ml", "fgf.21", "iga")] <- "Bioactive"
  # fsh/lh (reproductive hormones) and evenness/diversity (microbiome alpha-diversity) are
  # intentionally NOT bucketed into "Bioactive": FSH/LH are never described as an outcome
  # anywhere in the manuscript, and diversity/evenness are explicitly "exploratory" outcomes
  # in the manuscript's own outcome hierarchy, not secondary/bioactive. Leaving their
  # category as whatever the source metadata gives (blank/NA for all four, since
  # Milk_Component_Spec_IMiC_V02.csv doesn't classify them) correctly excludes them from
  # any primary/secondary-outcomes table (e.g. Table S7) built by filtering on category.
  
  res <- left_join(res, labels, by="biomarker") %>% 
    mutate(label_f=ifelse(is.na(label), biomarker, label))
  
  #clean up biomarker labels for plotting
  res$biomarker <- str_to_sentence(res$biomarker)
  #res$biomarker <- str_replace_all(res$biomarker, "\\_", " ")
  #res$biomarker <- str_replace_all(res$biomarker, "\\.", " ")
  #res$biomarker <- str_replace_all(res$biomarker, "\\-", " ")
  res$biomarker <- str_replace_all(res$biomarker, "nmol ml", "conc.")
  res$biomarker <- str_replace_all(res$biomarker, " pct", " %")
  
  return(res)
  
}




scale_continuous <- function(df) {
  # Identify columns to scale (numeric and not binary)
  cols_to_scale <- sapply(df, function(col) {
    is.numeric(col) && !all(col %in% c(0, 1, NA))
  })
  
  # Scale the identified columns
  df[cols_to_scale] <- lapply(df[cols_to_scale], scale)
  
  return(df)
}



convert_high_missing_to_indicators <- function(df, threshold = 0.5) {
  # Calculate missing percentage for each column
  missing_pct <- sapply(df, function(x) sum(is.na(x)) / length(x))
  
  # Identify columns with missing percentage > threshold
  high_missing_cols <- names(missing_pct)[missing_pct > threshold]
  
  # Convert high missing columns to indicators
  for (col in high_missing_cols) {
    df[[col]] <- as.numeric(!is.na(df[[col]]))
  }
  
  # Print summary of what was converted
  if (length(high_missing_cols) > 0) {
    cat("Converted", length(high_missing_cols), "columns to indicators:\n")
    for (col in high_missing_cols) {
      cat("  -", col, "(", round(missing_pct[col] * 100, 1), "% missing)\n")
    }
  } else {
    cat("No columns had >", threshold * 100, "% missing values.\n")
  }
  
  return(df)
}
