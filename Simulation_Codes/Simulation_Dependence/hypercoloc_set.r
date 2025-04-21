### Step 2: colocalize

library(hyprcoloc)
library(susieR)
for(file in list.files("/home/xc2270/COLOCBoost/AD_target_ColocBoost/Code/v10/", full.names = T)) {source(file)}

run_hypercoloc = function(file, actual_pheno_number, setting = "normal"){
    library(hyprcoloc)
    rds = readRDS(file)
    X = list()
    Y = list()
    variant = list()
    for(i in 1:actual_pheno_number){
      X[[i]] = rds$X[, , drop = FALSE]
      Y[[i]] = as.matrix(rds$Y[, i, drop = FALSE])
      variant[[i]] = rds$variant[[i]]
    }
    betas = sebetas = matrix(NA, nrow = ncol(X[[1]]), ncol = actual_pheno_number) ## calculate beta and se(beta)
    for (i in 1:actual_pheno_number){
        for(mm in 1:ncol(X[[1]])){
            rr <- susieR::univariate_regression(scale(X[[i]][,mm]), scale(Y[[i]]))
            betas[mm,i] = rr$betahat
            sebetas[mm, i] <- rr$sebetahat
        }
    }
    traits <- c(1:actual_pheno_number)
    rsid <- paste0("snp", c(1:ncol(X[[1]])))
    colnames(betas) <- colnames(sebetas) <- traits
    rownames(betas) <- rownames(sebetas) <- rsid
    source("~/cloud_colocalization/simulation_code/simulate_linreg.R")
    source("~/cloud_colocalization/simulation_code/misc.R")
    
    t1 <- Sys.time()
    res_hyprcoloc <- hyprcoloc(betas, sebetas, trait.names=traits, snp.id=rsid, snpscores=TRUE) 
    t2 <- Sys.time()
    #hyprcoloc_csets <- cred.sets(res_hyprcoloc_snpscores, value = 0.95)
    
    
    tb = res_hyprcoloc$results
    index = which(tb$traits != "None")
    tb = tb[index,,drop = FALSE]
    res_hyprcoloc$results = tb
    num = nrow(res_hyprcoloc$results)
    if(setting == "normal"){
        if(num ==0){
          causal_snp = NULL
          causal_trait = NULL
        }else{
        hyprcoloc_csets <- cred.sets(res_hyprcoloc, value = 0.95)
        causal_snp = list()
        causal_trait = list()
        for (i in 1:num){
            causal_snp[[i]] = as.numeric(gsub("[^0-9]", "", names(hyprcoloc_csets[[i]])))
            causal_trait[[i]] =as.numeric(unlist(strsplit(res_hyprcoloc$results$traits[i], ", ")))

          }
        }
    
    
                     
        ### find true coloc_var and coloc_trait, apply for all methods. ###
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
        ### find true coloc_var and coloc_trait, apply for all methods. ###

        hypercoloc_report = res_hyprcoloc 
        hypercoloc_report$true_variant = true_coloc_var
        hypercoloc_report$true_trait = true_coloc_trait
        hypercoloc_report$coloc_set = causal_snp    
        hypercoloc_report$coloc_trait = causal_trait
        hypercoloc_report$file = file
        hypercoloc_report$time = (t2 - t1)
        }else{
            if(num ==0){
          causal_snp = NULL
            }else{
            hyprcoloc_csets <- cred.sets(res_hyprcoloc, value = 0.95)
            causal_snp = list()
            for (i in 1:num){
                causal_snp[[i]] = as.numeric(gsub("[^0-9]", "", names(hyprcoloc_csets[[i]])))
               

              }
        }
        hypercoloc_report = res_hyprcoloc 
        hypercoloc_report$coloc_set = causal_snp    
        hypercoloc_report$file = file
        hypercoloc_report$time = (t2 - t1)
    
    }
    purity_check = lapply(hypercoloc_report$coloc_set, function(x) (get_purity(x, X= X[[1]]))[1])
    purity_index = which(unlist(purity_check) > 0.5)     
    if(length(purity_index) > 0){
        hypercoloc_report$coloc_set = hypercoloc_report$coloc_set[purity_index]
        hypercoloc_report$coloc_trait = hypercoloc_report$coloc_trait[purity_index]
    }else{
        print("All filtered out")
        hypercoloc_report$coloc_set = NULL
        hypercoloc_report$coloc_trait = NULL
    }
    hypercoloc_report$purity = purity_check
    return(hypercoloc_report)
}
