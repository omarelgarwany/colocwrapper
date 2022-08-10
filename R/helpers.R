read_tabix <- function(pheno_file,region,tabix_binary='/software/team152/oe2/bin/tabix') {
  #Load data
  tr_cmd <- paste(tabix_binary,pheno_file,region)
  tr_dat_txt <- system(tr_cmd,intern=T)
  if (length(tr_dat_txt) == 0) {
    warning(paste0('No trait data read for ',region,' from ',pheno_file,'. Returning NA...') )
    return(NA)
  }
  tr_dat <- read.table(text=tr_dat_txt,sep='\t',header=F,stringsAsFactors = F) 
  return(tr_dat)
}

create_window_nomax <- function(pos_c,window) {
  s_c <- c()
  e_c <- c()
  for (pos in pos_c) {
    s <- max(pos-(window/2),0)
    e <- s+window    
    
    s_c <- c(s_c,s)
    e_c <- c(e_c,e)
  }

  return(list(s_c,e_c))
}