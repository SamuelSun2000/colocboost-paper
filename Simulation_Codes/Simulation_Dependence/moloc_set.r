library(tidyverse)
library(moloc)
library(susieR)
for(file in list.files("/home/xc2270/COLOCBoost/AD_target_ColocBoost/Code/v10/", full.names = T)) {source(file)}
run_moloc <- function(file, actual_pheno_number, setting = "normal"){
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
        moloc <- moloc_test(listData = tb, save.SNP.info = TRUE)
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

        index2= which(rownames(moloc$best_snp[[1]]) == coloc_trait)
        colname =  colnames(moloc$best_snp[[2]])[index2]  
        row_index <- moloc$best_snp[[2]] %>% as_tibble() %>%
                    arrange(desc(!!sym(colname)))%>%
                    mutate(row_number = row_number()) %>%  
                    mutate(cumulative_sum = cumsum(!!sym(colname))) %>%
                    filter(cumulative_sum >= 0.95) %>% pull(row_number) %>% min()

# set based result
        snps = moloc$best_snp[[2]] %>% as_tibble() %>%
                        mutate(row_number = row_number())  %>%
                        arrange(desc(!!sym(colname)))%>% pull(row_number)
        snps = snps[1:row_index]

    moloc_report = moloc
    moloc_report$file = file
    moloc_report$true_variant = true_coloc_var
    moloc_report$true_trait = true_coloc_trait
    # if only find one trait colocalized, then it's just a fine map, no coloc evidence found, therefore should return NULL as result
    if(length(coloc_trait_index) >= 2){
    moloc_report$coloc_set = list(snps)
    moloc_report$coloc_trait = list(coloc_trait_index)
    }else{
    moloc_report$coloc_set = NULL
    moloc_report$coloc_trait = NULL
    }
    purity_check = lapply(moloc_report$coloc_set, function(x) (get_purity(x, X= X[[1]]))[1])
    purity_index = which(unlist(purity_check) > 0.5)     
    if(length(purity_index) > 0){
        moloc_report$coloc_set = moloc_report$coloc_set[purity_index]
        moloc_report$coloc_trait = moloc_report$coloc_trait[purity_index]
    }else{
        print("All filtered out")
        moloc_report$coloc_set = NULL
        moloc_report$coloc_trait = NULL
    }
                          
    moloc_report$time = (t2 - t1)
    moloc_report$purity = purity_check
    return(moloc_report)
}