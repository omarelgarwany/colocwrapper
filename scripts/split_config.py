import numpy as np
import pandas as pd
import sys
conf_f=sys.argv[1]
num_chunks=sys.argv[2]
chunkoutput_fs=sys.argv[3]

conf_dat=pd.read_csv(conf_f,sep='\t',header=None)
conf_dat_chunked=np.array_split(conf_dat, num_chunks)
chunks=[i for i in range(num_chunks)]

for i,df in enumerate(conf_dat_chunked):
  df.to_csv(chunkoutput_fs[i],sep='\t',header=False,index=False)
  print("Written chunk {0} to file: {1}".format(i,chunkoutput_fs[i]))
