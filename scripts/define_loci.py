import pandas as pd
import subprocess as sp
import io
import sys
import re

def add_chr(dat,col):
    dat.loc[~(dat[col].astype(str).str.startswith('chr')),col]='chr'+dat.loc[~(dat[col].astype(str).str.startswith('chr')),col].astype(str)
    return dat

def remove_chr(dat,col):
    dat.loc[dat[col].astype(str).str.startswith('chr'),col]=dat.loc[dat[col].astype(str).str.startswith('chr'),col].astype(str).str.slice(3)
    return dat

def get_chr_pos_from_snp_name(variants):
    chr=[x.split('_')[0].split(':')[0] for x in variants]
    pos=[x.split('_')[0].split(':')[1] for x in variants]
    return chr,pos

#snp_to_pos_map has to be a dataframe with columns named: variant, chr, pos
def get_chr_pos_from_rsid(variants,extra_args):
    snp_to_pos_map=extra_args
    chr_pos=pd.DataFrame({'variant': variants}).merge(snp_to_pos_map,how='left')
    return chr_pos['chr'].tolist(),chr_pos['pos'].tolist()

def parse_tabix(f,region,fix_chr=False,remove_chr=False,header=0,sep='\t'):
    if fix_chr:
        region=region if region.startswith('chr') else 'chr'+region
    if remove_chr:
        region=region[3:] if region.startswith('chr') else region

    # o=sp.run(['tabix',f,region],stdout=sp.PIPE)
    # tsv = io.BytesIO(o.stdout).getvalue()
    tsv=sp.check_output(['tabix',f,region]).decode(sys.stdout.encoding)
    if tsv.strip()=='':
        return pd.DataFrame()
    dat=pd.read_csv(io.StringIO(tsv),sep=sep,header=None)

    if header==0:
        h=pd.read_csv(f,sep=sep,header=header,nrows=1)
        dat=dat.rename(columns={old_col:new_col for old_col,new_col in zip(dat.columns.values, h.columns.values)})

    return dat
def snp_pos_to_region(snp_pos,padding=500000):
    chr,pos=[str(x) for x in snp_pos.split(':')]
    start_pos=int(pos)-int(padding)
    end_pos=int(pos)+int(padding)
    region=chr+":"+str(start_pos)+"-"+str(end_pos)
    return region
def region_to_chr_s_e(region):
    out=re.split('[:-]',region)
    if(len(out)==3):
        return chr,int(s),int(e)
    chr,s=out
    return chr,int(s)

def chr_s_e_to_region(chr,s,e=None):
    if e != None:
        chr=chr if str(chr).startswith('chr') else 'chr{0}'.format(chr)
        return '{0}:{1}-{2}'.format(chr,str(s),str(e))
    return '{0}:{1}'.format(chr,str(s))

def window_left_boundary(s,p,chrom_size):
    if s-p < 0:
        l=0

    elif s+p > chrom_size:
        l=chrom_size-1-(2*p)

    else:
        l=s-p

    return int(l)
def window_right_boundary(s,p,chrom_size):
    if s-p < 0:
        r=2*p
    elif s+p > chrom_size:
        r=chrom_size-1
    else:
        r=s+p
    return int(r)

def pad_region(x,chrpos_col,p,chrom_size_dict):
    # print(x[chrpos_col])
    chr=x[chrpos_col].split(':')[0]
    pos=int(x[chrpos_col].split(':')[1])

    lt_boundary=window_left_boundary(pos,p,chrom_size_dict[chr])
    rt_boundary=window_right_boundary(pos,p,chrom_size_dict[chr])

    region=chr_s_e_to_region(chr,lt_boundary,rt_boundary)

    return region
def save_summstat_chunk(df,f):
    df.to_csv(f,sep='\t',index=False)
def outside_region(snp_coord,region):
    snp_chr,snp_pos=snp_coord.split(':')
    region_chr,region_s,region_e=re.split('[:-]',region)
    snp_pos,region_s,region_e=[int(snp_pos),int(region_s),int(region_e)]
    inside=((snp_chr==region_chr) and (snp_pos >= region_s) and (snp_pos <= region_e))
    return not inside


import sys
import pandas as pd




thresh=0.00000005
thresh=0.8


#USAGE   screen -DR hyprcoloc.gwas_loci ; bsub -I -R "select[mem>16000] rusage[mem=16000]" -M16000 -J "get_sigloci" -m "modern_hardware"
# source /nfs/team152/oe2/sqtl/scripts/configs.sh; source /nfs/team152/oe2/sqtl/scripts/paths.sh ; for gw in "${gwas_names[@]}"; do j_c="python /nfs/team152/oe2/sqtl/scripts/hyprcoloc/lookup_locus.py /lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/gwas/${gw}.formatted.sorted.tsv.gz ${GENCODE_GENEID_FILE} /lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/hyprcoloc/${gw}.coloc 1000000 /lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/gwas/${gw}.sigloci" ; j_o="/lustre/scratch123/hgi/projects/macromapsqtl/oe2/job_logs/lookup_locus/lookup_locus_${gw}.o" ; j_n="lookup_locus_${gw}" ; j_e="/lustre/scratch123/hgi/projects/macromapsqtl/oe2/job_logs/lookup_locus/lookup_locus_${gw}.e" ; bsub -R "select[mem>16000] rusage[mem=16000]" -M16000 -o "$j_o" -e "$j_e" -J "$j_n" -m "modern_hardware" "$j_c"; done;

#FOR one example: gw=IBD; python ./lookup_locus.py "/lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/gwas/${gw}.formatted.sorted.tsv.gz" "$GENCODE_GENEID_FILE" "/lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/hyprcoloc/${gw}.coloc" 1000000 "/lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/gwas/${gw}.sigloci";


#To collect all sig loci by disease
#source /nfs/team152/oe2/sqtl/scripts/configs.sh; for gw in "${gwas_names[@]}";do cat /lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/gwas/"$gw".sigloci | awk -v gw="${gw}" 'BEGIN{FS="\t"; OFS="\t"}{print $1,$2,gw}'; done > /lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/gwas/all.sigloci
gwas_f=sys.argv[1]
gene_lookup_f=sys.argv[2]
padding=int(sys.argv[3])
chrom_size_f=sys.argv[4]
header=sys.argv[5]
output_f=sys.argv[6]

p = padding/2

chr_col='CHROM'
pos_col="GENPOS"
pval_col="PVAL"

chrom_size_dat=pd.read_csv(chrom_size_f,sep='\t',header=None)


if(header=='T'):
    gwas_dat=pd.read_csv(gwas_f,sep='\t',nrows=3)
else:
    gwas_dat=pd.read_csv(gwas_f,sep='\t',header=None,nrows=3)

gwas_dat=add_chr(gwas_dat,chr_col)
# chrom_size_dat=remove_chr(chrom_size_dat,0)

chrom_size_dict={chr:sz for chr,sz in zip(chrom_size_dat[0].tolist(),chrom_size_dat[1].tolist())}
gwas_hits=gwas_dat.loc[gwas_dat[pval_col] <= thresh,:].sort_values(by=[pval_col],ascending=True)

#Gene data
tss_dat=pd.read_csv(gene_lookup_f,sep='\t',header=None)
tss_regions=[region_to_chr_s_e(i) for i in tss_dat[1].tolist()]
tss_dat['_chr']=[i[0] for i in tss_regions]
tss_dat['_s']=[i[1] for i in tss_regions]

all_hit_genes_dat=pd.DataFrame()
if gwas_hits.shape[0] > 0:
    final_regions=[]
    lead_snps=[]
    while gwas_hits.shape[0] > 0:
        lowest_pval_snp_chr,lowest_pval_snp_pos=gwas_hits.sort_values(by=[pval_col],ascending=True).iloc[0,:][[chr_col,pos_col]]
        lowest_pval_snp='{0}:{1}'.format(lowest_pval_snp_chr,lowest_pval_snp_pos)
        lt_boundary = window_left_boundary(lowest_pval_snp_pos,p,chrom_size_dict['{0}'.format(lowest_pval_snp_chr)])
        rt_boundary = window_right_boundary(lowest_pval_snp_pos,p,chrom_size_dict['{0}'.format(lowest_pval_snp_chr)])

        lowest_pval_region=chr_s_e_to_region(lowest_pval_snp_chr,lt_boundary,rt_boundary)



        gwas_hits['_snp_chrpos']=gwas_hits[chr_col].astype(str)+":"+gwas_hits[pos_col].astype(str)
        gwas_hits['_outside']=gwas_hits.agg(lambda x: outside_region(x['_snp_chrpos'],lowest_pval_region),axis=1)

        gwas_hits=gwas_hits.loc[gwas_hits['_outside']==True,:].drop(['_outside','_snp_chrpos'],axis=1)

        #Getting gene information
        hit_genes_dat=tss_dat.loc[(tss_dat['_chr']==lowest_pval_snp_chr) & (tss_dat['_s']>=lt_boundary) & (tss_dat['_s']<=rt_boundary),:]
        hit_genes_dat['region']=lowest_pval_region
        hit_genes_dat['lead_snp']=lowest_pval_snp
        hit_genes_dat=hit_genes_dat.drop(columns=['_s','_chr'])
        hit_genes_dat=hit_genes_dat.loc[:,['region','lead_snp',0,1]]
        all_hit_genes_dat=pd.concat([all_hit_genes_dat,hit_genes_dat],ignore_index=True)
        # final_regions.append(lowest_pval_region)
        # lead_snps.append(lowest_pval_snp)
        print('Finished locus: {0}..'.format(lowest_pval_region))

loci_dat=all_hit_genes_dat.loc[:,['region','lead_snp']].drop_duplicates()

# loci_dat=pd.DataFrame({0:final_regions, 1:lead_snps})
# gws_regions=[region_to_chr_s_e(i) for i in loci_dat[0].tolist()]
# loci_dat['_chr']=[i[0] for i in gws_regions]
# loci_dat['_s']=[i[1] for i in gws_regions]
# loci_dat['_e']=[i[2] for i in gws_regions]




print('Outputting to: {0}..'.format(output_f))
loci_dat.to_csv(output_f+'.sigloci',sep='\t',index=False,header=False)
all_hit_genes_dat.to_csv(output_f+'.siglocigenes',sep='\t',index=False,header=False)
