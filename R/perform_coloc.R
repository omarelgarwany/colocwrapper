perform_coloc <- function(summstat_list,line_idx) {
  snp_names <- summstats_list[[line_idx]][['snps']]
  snp_positions <- as.integer(sapply(str_split(snp_names,pattern=":"),'[',2))
  
  #D1
  trait_1_idx <- 'T1'
  trait_1 <- summstats_list[[line_idx]][['traits']][[trait_1_idx]]
  type_D1 <-   summstats_list[[line_idx]][['trait_types']][[trait_1_idx]]
  
  beta_D1 <- summstats_list[[line_idx]][['beta']][[trait_1_idx]]
  varbeta_D1 <- summstats_list[[line_idx]][['se']][[trait_1_idx]]^2
  MAF_D1 <-summstats_list[[line_idx]][['maf']][[trait_1_idx]]
  N_D1 <- summstats_list[[line_idx]][['N']][[trait_1_idx]]
  
  
  D1 <- list(beta=setNames(beta_D1, snp_names),varbeta=setNames(varbeta_D1, snp_names),type=type_D1,snp=snp_names,position=snp_positions)
  if(type_D1=='quant') {D1[['N']] <- N_D1 ; D1[['MAF']] <- MAF_D1}
  
  
  #D2
  trait_2_idx <- 'T2'
  trait_2 <- summstats_list[[line_idx]][['traits']][[trait_2_idx]]
  type_D2 <-   summstats_list[[line_idx]][['trait_types']][[trait_2_idx]]
  
  beta_D2 <- summstats_list[[line_idx]][['beta']][[trait_2_idx]]
  varbeta_D2 <- summstats_list[[line_idx]][['se']][[trait_2_idx]]^2
  MAF_D2 <-summstats_list[[line_idx]][['maf']][[trait_2_idx]]
  N_D2 <- summstats_list[[line_idx]][['N']][[trait_2_idx]]
  
  
  D2 <- list(beta=setNames(beta_D2, snp_names),varbeta=setNames(varbeta_D2, snp_names),type=type_D2,snp=snp_names,position=snp_positions)
  
  if(type_D2=='quant') {D2[['N']] <- N_D2 ; D2[['MAF']] <- MAF_D2} 
  
  
  
  res <- coloc.abf(dataset1=D1,dataset2=D2)
  annotated_res <- res$summary %>% t() %>% as.data.frame()
  annotated_res$region <- region
  return(annotated_res)
  
}

perform_hyprcoloc <- function(summstat_list,line_idx,tr) {
  res <- hyprcoloc(summstats_list[[line_idx]][['beta']] %>% as.matrix(), summstats_list[[line_idx]][['se']] %>% as.matrix(), trait.names=names(summstat_list[[line_idx]][['traits']]), snp.id=summstat_list[[line_idx]][['snps']],bb.alg	=F)
  annotated_res=res$results
  annotated_res$nsnp <- nrow(summstats_list[[line_idx]][['beta']])
  annotated_res$line_idx <- line_idx
  return(annotated_res)
  
}
