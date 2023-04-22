import pandas as pd
import sys
import os
import pandas as pd
import glob



coloc_fs_prefix=sys.argv[1]
coloc_fs_suffix=sys.argv[2]
config_fs_prefix=sys.argv[3]
config_fs_suffix=sys.argv[4]
num_chunks=int(sys.argv[5])
output_f=sys.argv[6]

coloc_fs=[coloc_fs_prefix+str(i)+coloc_fs_suffix for i in range(1,num_chunks+1)]
config_fs=[config_fs_prefix+str(i)+config_fs_suffix for i in range(1,num_chunks+1)]




all_merged_dat=pd.DataFrame()

for coloc_f,config_f in zip(coloc_fs,config_fs):
    if (os.path.exists(coloc_f)) and (os.path.exists(config_f)):
        if (os.stat(coloc_f).st_size != 0) and (os.stat(config_f).st_size != 0):
            config_dat=pd.read_csv(config_f,sep='\t',header=None)
            config_dat['_line_idx']=[i for i in range(1,config_dat.shape[0]+1)]
            coloc_dat=pd.read_csv(coloc_f,sep='\t',header=None)

            merged_dat=coloc_dat.merge(config_dat,left_on=[9],right_on=['_line_idx'],how='left')
            all_merged_dat=pd.concat([all_merged_dat,merged_dat],ignore_index=True).drop(columns=['_line_idx'])

all_merged_dat.to_csv(output_f,sep='\t',header=False,index=False)
