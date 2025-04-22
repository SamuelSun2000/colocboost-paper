library(tidyverse)
variant_vector = function(true_variant, coloc_set){
    var_vec = numeric(length(true_variant))
    for(i in c(1:length(true_variant))){
        for (j in c(1:length(coloc_set))){
            if(true_variant[i] %in% unlist(coloc_set[[j]])){
                var_vec[i] = 1
            }else{
                var_vec[i] = var_vec[i]
            }
        }
    }
    return(var_vec)

}

set_vector =  function(true_variant, coloc_set){
    set_vec = numeric(length(coloc_set))
    for(i in c(1:length(true_variant))){
         for (j in c(1:length(coloc_set))){
                 inter = intersect(unlist(coloc_set[[j]]),
                                  true_variant)
             if(length(inter > 0)){
                set_vec[j] = 1 
                 }else{
                 set_vec[j] = 0
                 }
             }
    return(set_vec)
    }
}

set_size = function(coloc_set){
    size_vec = numeric(length(coloc_set))
    for (j in c(1:length(coloc_set))){
    size_vec[j] = length(unlist(coloc_set[j]))
    }
    return(size_vec)
}

fineboost_summary = function(true_variant, coloc_set){
    return(list(
    match_variant = sum(variant_vector(true_variant, coloc_set)),
    true_trait_number = sum(set_vector(true_variant, coloc_set)),
    total_variant_number = length(true_variant),
    total_set_number = length(coloc_set)
    ))
}
    