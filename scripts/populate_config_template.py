# import pandas as pd
import argparse
import string

my_parser = argparse.ArgumentParser(description='Populate config file placeholders')
my_parser.add_argument('-c','--config-template',type=str,required=True,dest='config_f',help='Config template file')
my_parser.add_argument('-s','--substitution-file',type=str,required=True,dest='subs_f',help='Substitution file')
my_parser.add_argument('-o','--output-file',type=str,required=True,dest='output_f',help='Output file')

args = my_parser.parse_args()
config_f=args.config_f
subs_f=args.subs_f
output_f=args.output_f

populated_fields=[]
# config_f='template.config'
# subs_f='example.subs'

with open(config_f) as f:
    template = f.readline().rstrip()
    f.close()

with open(subs_f) as f:
    lines = [line.rstrip().split(',') for line in f]
    for subs_l in lines:
        a = template.format(*subs_l)
        populated_fields.append("{0}{1}".format(a,'\n'))
# print(populated_fields)
f = open(output_f, "w")
f.writelines(populated_fields)
f.close()
