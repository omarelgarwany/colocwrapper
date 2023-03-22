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


collect_summstats <- function(config_f,yaml_f,sample_size_f,tabix_binary='/software/team152/oe2/bin/tabix') {
  sample_size_dat <- read.csv(sample_size_f,sep='\t',header=F)
  
  
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
      
      # tr_cmd <- cmds[i,paste0('V',j)]
      
      
      tr_ph <- config_dat[i,paste0('V',2*j)]
      tr_f <- config_dat[i,paste0('V',(2*j)+1)]
      
      tr_cmd <- paste(tabix_binary, tr_f, region)
      
      print(paste0('Running ', tr_cmd, '...') )
      tr_dat_txt <- system(tr_cmd,intern=T)
      
      #If looking up by the given region fails try adding or omitting "chr". it's annoying how different summstats files chromosomes are coded differently
      if(length(tr_dat_txt) == 0) {
        
        attempt_chr_region <- ifelse(grepl("^chr", region ), str_sub(region,4), paste0("chr", region))
        tr_cmd <- paste(tabix_binary, tr_f, attempt_chr_region)
        tr_dat_txt <- system(tr_cmd,intern=T)
      }
      
      #Applying some checks
      if (length(tr_dat_txt) == 0) {
        print(paste0('SKIPPED: no trait data read for ',tr_ph,'...') )
        skip_region <- T
        collected_summstat_list[[line_idx]][['skipped']] <- skip_region
        break
      }
      
      #Creating data frame from read text and checking if any column names specific
      tr_dat <- read.table(text=tr_dat_txt,sep='\t',header=F,stringsAsFactors = F)
      

      
      #This check indicates that even  after attempting numeric conversion it failed (definitely characters)
      any_non_numeric_cols <-  (as.numeric(tr_val_cols[[suffix]]) %>% suppressWarnings() %>% is.na() %>% any())
      if (any_non_numeric_cols) {
        tr_f_header <- read.csv(tr_f,sep='\t',nrows=1,header=F,stringsAsFactors = F) 
        
        new_col_names <- setNames(names(tr_dat),tr_f_header)
        tr_dat <- tr_dat %>% dplyr::rename(all_of(new_col_names))

      }
      
      #Extracting column names. We add N/MAF if it's a quantitative trait. We also check if any provided column name is non-numeric. if any column name is non-numeric we load the header
      snp_chr_col <- ifelse(tr_val_cols[[suffix]][[1]] %>% as.numeric() %>% suppressWarnings() %>% is.na(),  tr_val_cols[[suffix]][[1]], as.numeric(tr_val_cols[[suffix]][[1]]))
      snp_pos_col <- ifelse(tr_val_cols[[suffix]][[2]] %>% as.numeric() %>% suppressWarnings() %>% is.na(),  tr_val_cols[[suffix]][[2]],as.numeric(tr_val_cols[[suffix]][[2]]))
      beta_col <- ifelse(tr_val_cols[[suffix]][[3]] %>% as.numeric() %>% suppressWarnings() %>% is.na(),  tr_val_cols[[suffix]][[3]],as.numeric(tr_val_cols[[suffix]][[3]]))
      se_col <- ifelse(tr_val_cols[[suffix]][[4]] %>% as.numeric() %>% suppressWarnings() %>% is.na(),  tr_val_cols[[suffix]][[4]],as.numeric(tr_val_cols[[suffix]][[4]]))
      pval_col <- ifelse(tr_val_cols[[suffix]][[5]] %>% as.numeric() %>% suppressWarnings() %>% is.na(),  tr_val_cols[[suffix]][[5]],as.numeric(tr_val_cols[[suffix]][[5]]))


      
      
      #Correcting chromosome names in snp name
      tr_dat[[snp_chr_col]] <- ifelse(grepl("^chr", as.character(tr_dat[[snp_chr_col]]) ), as.character(tr_dat[[snp_chr_col]]), paste0("chr", as.character(tr_dat[[snp_chr_col]]) ))
 
      
      
      #If it's a quant trait, then sometimes tabix a certain region might give you more than one phenotype
      if(tr_type=='quant') {
        tr_dat <- tr_dat %>% unite('ph',all_of(tr_id_cols[[suffix]]),remove=F,sep='|')
        tr_dat <- tr_dat %>% filter( ph %in% all_of(tr_ph) )
      }
      
      #MAF columns
      if(tr_type=='quant') { 
        maf_col <- ifelse(tr_val_cols[[suffix]][[6]] %>% as.numeric() %>% suppressWarnings() %>% is.na(),  tr_val_cols[[suffix]][[6]],as.numeric(tr_val_cols[[suffix]][[6]]))
        maf=tr_dat[[maf_col]] 
      } else { 
          maf <- -1 
      }
      
      #N column. if user provided a column name use it if not look for sample size file
      if(length(tr_val_cols[[suffix]])==7) {
        N_col <- ifelse(tr_val_cols[[suffix]][[7]] %>% as.numeric() %>% suppressWarnings() %>% is.na(),  tr_val_cols[[suffix]][[7]],as.numeric(tr_val_cols[[suffix]][[7]]))
        N <- first(as.numeric(tr_dat[[N_col]]))
      } else {
        if(tr_type=='quant') {
          N_info <- sample_size_dat[sample_size_dat['V1']==tr_f,]

          if(dim(N_info)[[1]] == 0) {
            print(paste0('SKIPPED: no sample size information matching file: ',tr_f,'...'))
            skip_region <- T
            collected_summstat_list[[line_idx]][['skipped']] <- skip_region
            break 
          }
          N <- first(N_info[['V2']])
          
        } else {
          #We don't care if it's not quant trait...
          N=-1
        }        
      }

      
      


      
      
      #Applying some checks
      if(dim(tr_dat)[1] == 0) {
        print(paste0('SKIPPED: no data read  or satisfying criteria for ',tr_ph,'...'))
        skip_region <- T
        collected_summstat_list[[line_idx]][['skipped']] <- skip_region
        break
      }
      

      
      tr_summstat <- data.frame(snp=paste0(tr_dat[[snp_chr_col]],':',tr_dat[[snp_pos_col]]),beta=tr_dat[[beta_col]],se=tr_dat[[se_col]],pval=tr_dat[[pval_col]],N=N,maf=maf)
      colnames(tr_summstat) <- c(c('snp'),paste(c('beta','se','pval','N','maf'),j,sep='_'))

      #Merged summstats
      if (j > 1) {
        all_summstat <- tr_summstat %>% inner_join(all_summstat,by='snp')
        
      } else {
        all_summstat <- tr_summstat
      }
      
      #Keeping track of phenotypes names
      tr_ph_list[suffix] <- tr_ph
    }
    ###############################
    #END OF INNER LOOP OVER TRAITS#
    ###############################
    if(skip_region) {
      print(paste0('Region couldnt be processed for command: ',tr_cmd))
      # collected_summstat_list[[i]] <- NA
      skip_region <- T
      collected_summstat_list[[line_idx]][['skipped']] <- skip_region
      next
    }
    print(paste0("Joined summstats: ",dim(all_summstat)[1]) )
    
    if(dim(all_summstat)[1] == 0) {
      print('SKIPPED: no joint summstats data...')
      # collected_summstat_list[[i]] <- NA
      skip_region <- T
      collected_summstat_list[[line_idx]][['skipped']] <- skip_region
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
      skip_region <- T
      collected_summstat_list[[line_idx]][['skipped']] <- skip_region
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
    collected_summstat_list[[line_idx]][['skipped']] <- skip_region
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
