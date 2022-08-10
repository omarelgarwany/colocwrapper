determine_num_traits <- function(config_dat) {
  num_traits <- (dim(config_dat)[[2]]-1)/2
  return(num_traits)
}
process_yaml <- function(yaml_dat,num_traits) {

  
  suffixes <- paste0('T',1:num_traits)
  tr_val_cols <-setNames(yaml_dat[['tr_val_cols']],suffixes)
  tr_id_cols <- setNames(yaml_dat[['tr_id_cols']],suffixes)
  tr_types <- setNames(yaml_dat[['type']],suffixes)
  
  yaml_meta <- list(suffixes=suffixes,
                    tr_val_cols=setNames(yaml_dat[['tr_val_cols']],suffixes),
                    tr_id_cols=setNames(yaml_dat[['tr_id_cols']],suffixes),
                    tr_types=setNames(yaml_dat[['type']],suffixes)
                    )
  return(yaml_meta)
}

generate_cmds <- function(config_dat,num_traits,tabix_binary) {
  cmds <- data.frame(row.names = 1:nrow(config_dat),stringsAsFactors = F)
  for(i in 1:num_traits) {
    f_col <- paste0('V',3+((i-1)*2))
    new_col <- data.frame(paste(tabix_binary,config_dat[[f_col]],config_dat[['V1']]),stringsAsFactors = F)
    cmds <- cbind(cmds,new_col)
    
  }
  colnames(cmds) <- paste0('V',1:num_traits)
  return(cmds)
}

get_all_regions <- function(config_dat){
  return(config_dat[['V1']])
}

collect_summstats <- function(config_f,yaml_f,tabix_binary='/software/team152/oe2/bin/tabix') {
  config_dat <- read.csv(config_f,sep='\t',header=F,stringsAsFactors = F)
  num_traits <- determine_num_traits(config_dat)
  
  yaml_dat <- config::get(file=yaml_f)
  yaml_meta <- process_yaml(yaml_dat,num_traits)
  suffixes <- yaml_meta[['suffixes']]
  tr_val_cols <- yaml_meta[['tr_val_cols']]
  tr_id_cols <- yaml_meta[['tr_id_cols']]
  tr_types <- yaml_meta[['tr_types']]
  
  regions <- get_all_regions(config_dat)
  
  #Building cmds dataframe#
  cmds <- generate_cmds(config_dat,num_traits,tabix_binary)
  
  
  collected_summstat_list <- list()
  print(nrow(cmds))
  collected_summstat_list[['num_lines']] <- nrow(cmds)
  
  
  for (i in 1:nrow(cmds)) {
    line_idx <- as.character(i)
    collected_summstat_list[[line_idx]] <- list()
    region <- regions[[i]]
    skip_region <- F
    tr_ph_list <- c()
    # i <- 1
    #################################
    #START OF INNER LOOP OVER TRAITS#
    #################################
    for (j in 1:num_traits) {
      
      # j <- 2
      suffix <- suffixes[[j]]
      tr_type <- tr_types[[suffix]]
      
      tr_cmd <- cmds[i,paste0('V',j)]
      
      tr_ph <- config_dat[i,paste0('V',2*j)]
      tr_ph_list[suffix] <- tr_ph
      
      tr_dat_txt <- system(tr_cmd,intern=T)
      print(paste0('Running ', tr_cmd, '...') )
      
      
      
      
      
      if (length(tr_dat_txt) == 0) {
        print(paste0('SKIPPED: no trait data read for ',tr_ph,'...') )
        skip_region <- T
        break
      }
      
      
      
      tr_dat <- read.table(text=tr_dat_txt,sep='\t',header=F,stringsAsFactors = F)
      current_tr_id_cols <- paste0('V',tr_id_cols[[suffix]])
      tr_dat <- tr_dat %>% unite('ph',all_of(current_tr_id_cols),remove=F,sep='|')
      print(tr_ph)
      tr_dat <- tr_dat %>% filter( ph %in% all_of(tr_ph) )

      if(dim(tr_dat)[1] == 0) {
        print(paste0('SKIPPED: no data read  or satisfying criteria for ',tr_ph,'...'))
        skip_region <- T
        break
      }
      
      
      #Extracting column names. We add N/MAF if it's a quantitative trait
      snp_chr_col <- paste0('V',tr_val_cols[[suffix]][[1]])
      snp_pos_col <- paste0('V',tr_val_cols[[suffix]][[2]])
      beta_col <- paste0('V',tr_val_cols[[suffix]][[3]])
      se_col <- paste0('V',tr_val_cols[[suffix]][[4]])
      pval_col <- paste0('V',tr_val_cols[[suffix]][[5]])
      
      if(tr_type=='quant') {
        N_col <- paste0('V',tr_val_cols[[suffix]][[6]])
        

        
        maf_col <- paste0('V',tr_val_cols[[suffix]][[7]])
      }
      
      #Building summstats dataframe depending on type
      if (tr_type=='quant') {
        N <- ifelse( startsWith(N_col,'N:') , as.integer( str_split(N_col,':')[[1]][2] ), tr_dat[[N_col]])
        maf=tr_dat[[maf_col]]
      } else {
        N=-1
        maf=-1
      }
      
      tr_summstat <- data.frame(snp=paste0(tr_dat[[snp_chr_col]],':',tr_dat[[snp_pos_col]]),beta=tr_dat[[beta_col]],se=tr_dat[[se_col]],pval=tr_dat[[pval_col]],N=N,maf=maf)
      colnames(tr_summstat) <- c(c('snp'),paste(c('beta','se','pval','N','maf'),j,sep='_'))
      # print(tr_summstat %>% head(20))
      #Merged summstats
      if (j > 1) {
        all_summstat <- tr_summstat %>% inner_join(all_summstat,by='snp')
        
      } else {
        all_summstat <- tr_summstat
      }
      
    }
    ###############################
    #END OF INNER LOOP OVER TRAITS#
    ###############################
    if(skip_region) {
      print(paste0('Region couldnt be processed for command: ',tr_cmd))
      collected_summstat_list[[i]] <- NA
      next
    }
    print(paste0("Joined summstats: ",dim(all_summstat)[1]) )
    
    if(dim(all_summstat)[1] == 0) {
      print('SKIPPED: no joint summstats data...')
      collected_summstat_list[[i]] <- NA
      next
    }
    
    se_cols <- all_summstat[paste('se',1:num_traits,sep='_')] %>% names()
    beta_cols <- all_summstat[paste('beta',1:num_traits,sep='_')] %>% names()
    pval_cols <- all_summstat[paste('pval',1:num_traits,sep='_')] %>% names()
    N_cols <- all_summstat[paste('N',1:num_traits,sep='_')] %>% names()
    maf_cols <- all_summstat[paste('maf',1:num_traits,sep='_')] %>% names()
    
    #Removing all duplicated entries
    all_summstat <- all_summstat %>% filter(! (snp %in% all_summstat$snp[all_summstat$snp %>% duplicated()]))
    #Removing NA's
    all_summstat <- all_summstat %>% drop_na()
    #Removing SE=0
    all_summstat <- all_summstat %>% filter(if_all(all_of(se_cols), ~ . > 0)) %>% filter(if_all(all_of(se_cols), ~ . != Inf))  %>% filter(if_all(all_of(beta_cols), ~ . != Inf))
    if(dim(all_summstat)[1] < 2) {
      print('SKIPPED: filtering join summstats resulted in data with less than 2 SNPs...')
      next
    }
    
    #
    all_summstat$snp <- as.character(all_summstat$snp)
    rownames(all_summstat) <- all_summstat$snp
    
    #The columns are always named T1,..,Tn to avoid conflicts in trait names
    new_cols <- suffixes
    #
    beta_dat <- all_summstat[beta_cols] %>% dplyr::rename(setNames(beta_cols,new_cols))
    se_dat <- all_summstat[se_cols] %>% dplyr::rename(setNames(se_cols,new_cols))
    pval_dat <- all_summstat[pval_cols] %>% dplyr::rename(setNames(pval_cols,new_cols))
    N_dat <- all_summstat[N_cols] %>% dplyr::rename(setNames(N_cols,new_cols))
    maf_dat <- all_summstat[maf_cols] %>% dplyr::rename(setNames(maf_cols,new_cols))
    #
    
    collected_summstat_list[[line_idx]][['beta']] <- beta_dat
    collected_summstat_list[[line_idx]][['se']] <- se_dat
    collected_summstat_list[[line_idx]][['pval']] <- pval_dat
    collected_summstat_list[[line_idx]][['N']] <- N_dat
    collected_summstat_list[[line_idx]][['maf']] <- maf_dat
    collected_summstat_list[[line_idx]][['snps']] <- all_summstat$snp
    collected_summstat_list[[line_idx]][['region']] <- region
    collected_summstat_list[[line_idx]][['trait_types']] <- tr_types
    collected_summstat_list[[line_idx]][['traits']] <- tr_ph_list
    collected_summstat_list[[line_idx]][['num_traits']] <- num_traits

  }
  


  return(collected_summstat_list) 
}


#####
