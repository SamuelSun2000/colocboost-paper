library(tidyverse)
library(moloc)
run_moloc <- function(file, actual_pheno_number){
    rds = readRDS(file)                 
    X = list()
    Y = list()
    variant = list()
    for(i in 1:actual_pheno_number){
      X[[i]] = rds$X[, , drop = FALSE]
      Y[[i]] = as.matrix(rds$Y[, i, drop = FALSE])
      variant[[i]] = rds$variant[[i]]
    }          
    tb <- list()
    cnt = 1
    for (i in 1:actual_pheno_number){
        betas = c()
        sebetas = c()
        MAF = c()
        for(mm in 1:ncol(X[[i]])){
            rr <- susieR::univariate_regression(X[[i]][,mm], Y[[i]])
            betas[mm] = rr$betahat
            sebetas[mm] <- rr$sebetahat
            MAF[mm] = min(sum(X[[i]][,mm])/(2*nrow(X[[i]])), 1-sum(X[[i]][,mm])/(2*nrow(X[[i]])))
        }

        tb[[cnt]] = tibble(
            SNP = 1:ncol(X[[i]]),
            BETA = betas,
            SE = sebetas,
            N = nrow(X[[i]]),
            MAF = MAF
        )
        cnt = cnt+1

     }
    options(scipen = 1, digits = 3)
                     
    t1 <- Sys.time()
    moloc <- moloc_test(listData = tb)
    t2 <- Sys.time()
                     
    all_var = unlist(variant)
    coloc_var_index = as.numeric(names(which(table(all_var)>= 2)))
    
    # make a list recording what traits the coloc_vairants are in.
    true_coloc_var = coloc_var_index

    true_coloc_trait = list()

    # Iterate through the list
    for (variant_index in 1:length(coloc_var_index)){
        temp_vec = c()
        for(i in 1:length(variant)){
            if(coloc_var_index[variant_index] %in% variant[[i]]){
                temp_vec = c(temp_vec, i)
            }
        }
        true_coloc_trait[[variant_index]] = temp_vec

        }
                     
    index = which.max(moloc$priors_lkl_ppa$PPA)
    hypothesis_name = rownames(moloc$priors_lkl_ppa)[index]
    coloc_traits = unlist(strsplit(hypothesis_name, ","))
    trait_num = max(nchar(coloc_traits))
    coloc_trait = coloc_traits[which.max(nchar(coloc_traits))]
    coloc_trait_index = c()
    for (item in unlist(strsplit(coloc_trait, NULL))){
        temp_index = as.integer(charToRaw(item)) - as.integer(charToRaw("a")) + 1
        coloc_trait_index = c(coloc_trait_index, temp_index)
    }

    index2= which(rownames(moloc$best_snp) == coloc_trait)
    coloc_snp = moloc$best_snp$best.snp.coloc[index2]
    moloc_report = moloc
    moloc_report$file = file
    moloc_report$true_variant = true_coloc_var
    moloc_report$true_trait = true_coloc_trait
    # if only find one trait colocalized, then it's just a fine map, no coloc evidence found, therefore should return NULL as result
    if(length(coloc_trait_index) >= 2){
    moloc_report$coloc_set = list(as.numeric(coloc_snp))
    moloc_report$coloc_trait = list(coloc_trait_index)
    }else{
    moloc_report$coloc_set = NULL
    moloc_report$coloc_trait = NULL
    }
    moloc_report$time = (t2 - t1)
    return(moloc_report)
}