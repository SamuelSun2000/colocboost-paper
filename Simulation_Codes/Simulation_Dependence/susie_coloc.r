for(file in list.files("/home/xc2270/COLOCBoost/AD_target_ColocBoost/Code/v10/", full.names = T)) {source(file)}
run_susie_coloc <- function(file, actual_pheno_number = 2){
    library("susieR")
    library("coloc")
    library("tidyverse")
    rds_list <- list()
    
    # step1: extract specific number of traits
    cnt = 1
    # read 2 traits as group, if we are going to analyse 3 traits-colocalization, we
    # select the first 3 traits X-Y pari in this group
    rds = readRDS(file)                 
    X = list()
    Y = list()
    variant = list()
    for(i in 1:actual_pheno_number){
      X[[i]] = rds$X[, , drop = FALSE]
      Y[[i]] = as.matrix(rds$Y[, i, drop = FALSE])
      variant[[i]] = rds$variant[[i]]
    }   
    Xe = X[[1]]
    Xp = X[[2]]
    Ye = Y[[1]]
    Yp = Y[[2]]           
    beta_p = varbeta_p = rep(0, ncol(Xp))
    beta_e = varbeta_e = rep(0, ncol(Xe)) ## calculate beta and varbeta
    for(mm in 1:ncol(Xe)){
      rr <- susieR::univariate_regression(Yp, Xp[,mm])
      beta_p[mm] = rr$betahat
      varbeta_p[mm] <- rr$sebetahat^2
      rr <- susieR::univariate_regression(Ye, Xe[,mm])
      beta_e[mm] = rr$betahat
      varbeta_e[mm] <- rr$sebetahat^2
    }
    LDp = cor(Xp)
    LDe = cor(Xe)
    snp = paste0("snp", c(1:ncol(Xe)))
    colnames(LDe) = rownames(LDe) = colnames(LDp) = rownames(LDp) = snp
    dp = list("beta" = beta_p, "varbeta" = varbeta_p,
              "LD" = LDp, "snp" = snp, "N" = nrow(Xp),
              "type" = "quant", "sdY" = sd(Yp))
    de = list("beta" = beta_e, "varbeta" = varbeta_e,
              "LD" = LDe, "snp" = snp, "N" = nrow(Xe),
              "type" = "quant", "sdY" = sd(Ye))
    t1 <- Sys.time()
    out_p <- runsusie(dp)
    out_e <- runsusie(de)
    out_coloc = coloc.susie(out_e, out_p)
    t2 <- Sys.time()
               
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
               
    ### summary 
               
df <- out_coloc$summary

if(is.null(df)){
    out_coloc = list()
    sets = NULL
    coloc_trait_index = NULL
}else{# step 1: in summary, get which row has PP.H4 >= 0.8
    index <- which(df$PP.H4.abf > df$PP.H1.abf & df$PP.H4.abf > df$PP.H2.abf & 
    df$PP.H4.abf > df$PP.H3.abf)
    if(length(index) >0){
        # step 2: results: each column stands for one row in summary rows
        tb <- out_coloc$results
        number = 1
        # because there can be multiple rows that has H4 as max so add a loop
            for (i in index) {
        # row 2 correspond to column 3 because the first column is SNP so + 1
                new_index = i + 1
                colname = names(tb)[new_index]
        # !!sym(colname) can pass variable as column names, desc order that columns
        # use cumsum to get threshold of 0.95
                row_index <- tb %>% as_tibble() %>%
                arrange(desc(!!sym(colname)))%>%
                mutate(row_number = row_number()) %>%  
                mutate(cumulative_sum = cumsum(!!sym(colname))) %>%
                filter(cumulative_sum >= 0.95) %>% pull(row_number) %>% min()
                snps = tb %>% as_tibble() %>%
                arrange(desc(!!sym(colname)))%>% pull(snp)
         # get "sets" from susie_coloc
                if (number == 1){
                    sets = list(as.numeric(gsub("snp", "", snps[1:row_index])))
                }else{
         # if more than 2 sets, combine them together in a vector
                    new_list = as.numeric(gsub("snp", "", snps[1:row_index]))
                    sets = c(sets, list(new_list))
                }
                number = number + 1

            }
        }else{
        sets = NULL
        }
    coloc_trait_index = list()
    if(length(sets) == 0){
        coloc_trait_index = NULL
    }else{
            for(j in 1:length(sets)){
                coloc_trait_index[[j]] = c(1,2)
            }
    }
}
    susie_coloc_report = out_coloc
    susie_coloc_report$file = file
    susie_coloc_report$true_variant = true_coloc_var
    susie_coloc_report$true_trait = true_coloc_trait
    susie_coloc_report$coloc_set = sets
    susie_coloc_report$coloc_trait = coloc_trait_index
    susie_coloc_report$time = (t2 - t1)
    
    purity_check = lapply(susie_coloc_report$coloc_set, function(x) (get_purity(x, X= X[[1]]))[1])
    purity_index = which(unlist(purity_check) > 0.5)     
    if(length(purity_index) > 0){
        susie_coloc_report$coloc_set = susie_coloc_report$coloc_set[purity_index]
        susie_coloc_report$coloc_trait = susie_coloc_report$coloc_trait[purity_index]
    }else{
        print("All filtered out")
        susie_coloc_report$coloc_set = NULL
        susie_coloc_report$coloc_trait = NULL
    }
    susie_coloc_report$purity = purity_check
    
    return(susie_coloc_report)
}
