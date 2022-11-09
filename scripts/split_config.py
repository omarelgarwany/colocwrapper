import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Split configuration file into chunks')

parser.add_argument('--config-file', dest='conf_f',required=True, help='Configuration file')
parser.add_argument('--num-chunks', dest='num_chunks',type=int,required=True, help='number of chunks')
parser.add_argument('--chunk-files-list', dest='chunkoutput_fs_list',required=True ,help='List of chunked config files')
args = parser.parse_args()

conf_f=args.conf_f
num_chunks=args.num_chunks
chunkoutput_fs_list=args.chunkoutput_fs_list


chunkoutput_fs=(pd.read_csv(chunkoutput_fs_list,sep='\t',header=None))[0].tolist()



conf_dat=pd.read_csv(conf_f,sep='\t',header=None)
conf_dat_chunked=np.array_split(conf_dat, num_chunks)
chunks=[i for i in range(num_chunks)]

for i,df in enumerate(conf_dat_chunked):
  df.to_csv(chunkoutput_fs[i],sep='\t',header=False,index=False)
  print("Written chunk {0} to file: {1}".format(i+1,chunkoutput_fs[i]))
