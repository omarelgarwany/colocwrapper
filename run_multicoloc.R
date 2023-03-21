library(coloc)
# library(hyprcoloc)
library(tidyverse)
library(config)
library(optparse)



test_run <- F

if(test_run == F) {
  
  
  option_list <- list(
    make_option(c("-c", "--config-file"), type="character",
                help="Configuration file in the format specified on Github", metavar="FILE"),
    make_option(c("-y", "--yaml-file"), type="character",
                help="YAML file in the format specified on Github", metavar="FILE"),
    make_option(c("-o", "--output-file"), type="character",
                help="Output file name [default= %default]", metavar="FILE"),
    make_option(c("-s", "--sample-sizes"), type="character",
                help="Sample size info file name [default= %default]", metavar="FILE"),
    make_option(c("-t", "--coloc-type"), type="character",default='coloc',
                help="Which type of colocalization programme? Currently available COLOC and HyprColoc [default= %default]")
  ); 
  
  parser <- OptionParser()
  opt <- parse_args(OptionParser(option_list=option_list))
  
  config_f <- opt[['config-file']]
  yaml_f <- opt[['yaml-file']]
  output_f <- opt[['output-file']]
  coloc_type <- opt[['coloc-type']]
  sample_size_f <- opt[['sample-sizes']]
  
  available_colocs <- c('coloc','hyprcoloc')
  
  if(!any(available_colocs==coloc_type)) {
    stop('coloc type has to be either coloc or hyprcoloc')
  }
} else {
  config_f <- '/nfs/team152/oe2/sqtl/scripts/dissect_locus/configs/LACC1_eqtl.IBD.txt'
  yaml_f <- '/nfs/team152/oe2/sqtl/scripts/dissect_locus/yamls/macromapeqtl_gwas.yaml'
  output_f <- '/nfs/team152/oe2/sqtl/scripts/colocwrapper/example/out.o'
  coloc_type <- 'coloc'
}
source('/nfs/users/nfs_o/oe2/oe2_packages/colocwrapper/R/perform_coloc.R')
source('/nfs/users/nfs_o/oe2/oe2_packages/colocwrapper/R/collect_summstats.R')

print('===================')
print(paste0('Config file: ',config_f))
print(paste0('YAML file: ',yaml_f))
print(paste0('Sample size file: ',sample_size_f))
print(paste0('Output file: ',output_f))
print('===================')
summstats_list <- collect_summstats(config_f,yaml_f,sample_size_f)

coloc_res <- data.frame()
for(line in 1:summstats_list[['num_lines']]) {
  line_idx <-as.character(line)
  print(paste0('Doing colocalization for region in line: ',line_idx))
  # print(summstats_list)
  if (summstats_list[[line_idx]][['skipped']]==F) {
    region <- summstats_list[[line_idx]][['region']]
    #Formatting output
    if (coloc_type=='coloc') {
      #COLOC
      annotated_res <- perform_coloc(summstats_list,as.character(line_idx))

    } else if (coloc_type=='hyprcoloc') {
      #HYPRCOLOC
      annotated_res <- perform_hyprcoloc(summstats_list,as.character(line_idx))
      
    }
    
    num_traits <- summstats_list[[line_idx]][['num_traits']]
    for (j in 1:num_traits) {
      trait_id <- paste0('T',j)
      annotated_res[[trait_id]] <- summstats_list[[line_idx]][['traits']][trait_id]
    }    
    annotated_res <- annotated_res %>% mutate(n=line_idx)
    coloc_res <- rbind(coloc_res,annotated_res)
  } 
  
  
}

#print(coloc_res)
write.table(coloc_res,file = output_f,quote=F,col.names = F,row.names = F,sep='\t')


