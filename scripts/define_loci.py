import pandas as pd
import sys
import re
import pandas as pd

#Adds "chr" to chromsome name. This may need to change in future versions
def add_chr(dat,col):
    dat.loc[~(dat[col].astype(str).str.startswith('chr')),col]='chr'+dat.loc[~(dat[col].astype(str).str.startswith('chr')),col].astype(str)
    return dat

#Takes a string chr:s[-e] and converts them to chr, and s [and e]
def region_to_chr_s_e(region):
    out=re.split('[:-]',region)
    if(len(out)==3):
        return chr,int(s),int(e)
    chr,s=out
    return chr,int(s)

#Takes chr, s, and e and converts them into a string chr:s-e
def chr_s_e_to_region(chr,s,e=None):
    if e != None:
        chr=chr if str(chr).startswith('chr') else 'chr{0}'.format(chr)
        return '{0}:{1}-{2}'.format(chr,str(s),str(e))
    return '{0}:{1}'.format(chr,str(s))

#Determines left-end boundary of a region s given padding p. Also requires chromosome size chrom_size
def window_left_boundary(s,p,chrom_size):
    if s-p < 0:
        l=0

    elif s+p > chrom_size:
        l=chrom_size-1-(2*p)

    else:
        l=s-p

    return int(l)

#Determines right-end boundary of a region s given padding p. Also requires chromosome size chrom_size
def window_right_boundary(s,p,chrom_size):
    if s-p < 0:
        r=2*p
    elif s+p > chrom_size:
        r=chrom_size-1
    else:
        r=s+p
    return int(r)

#Defines variants outside the position (usually defined by top SNP)
def outside_region(snp_coord,region):
    snp_chr,snp_pos=snp_coord.split(':')
    region_chr,region_s,region_e=re.split('[:-]',region)
    snp_pos,region_s,region_e=[int(snp_pos),int(region_s),int(region_e)]
    inside=((snp_chr==region_chr) and (snp_pos >= region_s) and (snp_pos <= region_e))
    return not inside










#USAGE   screen -DR hyprcoloc.gwas_loci ; bsub -I -R "select[mem>16000] rusage[mem=16000]" -M16000 -J "get_sigloci" -m "modern_hardware"
# source /nfs/team152/oe2/sqtl/scripts/configs.sh; source /nfs/team152/oe2/sqtl/scripts/paths.sh ; for gw in "${gwas_names[@]}"; do j_c="python /nfs/team152/oe2/sqtl/scripts/hyprcoloc/lookup_locus.py /lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/gwas/${gw}.formatted.sorted.tsv.gz ${GENCODE_GENEID_FILE} /lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/hyprcoloc/${gw}.coloc 1000000 /lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/gwas/${gw}.sigloci" ; j_o="/lustre/scratch123/hgi/projects/macromapsqtl/oe2/job_logs/lookup_locus/lookup_locus_${gw}.o" ; j_n="lookup_locus_${gw}" ; j_e="/lustre/scratch123/hgi/projects/macromapsqtl/oe2/job_logs/lookup_locus/lookup_locus_${gw}.e" ; bsub -R "select[mem>16000] rusage[mem=16000]" -M16000 -o "$j_o" -e "$j_e" -J "$j_n" -m "modern_hardware" "$j_c"; done;

#FOR one example: gw=IBD; python ./lookup_locus.py "/lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/gwas/${gw}.formatted.sorted.tsv.gz" "$GENCODE_GENEID_FILE" "/lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/hyprcoloc/${gw}.coloc" 1000000 "/lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/gwas/${gw}.sigloci";


#To collect all sig loci by disease
#source /nfs/team152/oe2/sqtl/scripts/configs.sh; for gw in "${gwas_names[@]}";do cat /lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/gwas/"$gw".sigloci | awk -v gw="${gw}" 'BEGIN{FS="\t"; OFS="\t"}{print $1,$2,gw}'; done > /lustre/scratch123/hgi/projects/macromapsqtl/oe2/output/gwas/all.sigloci

gwas_f=sys.argv[1]
gene_lookup_f=sys.argv[2]
chrom_size_f=sys.argv[3]
window=int(sys.argv[4])
header=sys.argv[5]
output_f=sys.argv[6]

#Defining relevant column names
thresh=0.00000005
chr_col='CHROM'
pos_col="GENPOS"
pval_col="PVAL"


p = window/2


#Reading/formatting GWAS data
if(header=='T'):
    gwas_dat=pd.read_csv(gwas_f,sep='\t')
else:
    gwas_dat=pd.read_csv(gwas_f,sep='\t',header=None)
gwas_dat=add_chr(gwas_dat,chr_col)

#Loading and formatting chromosome size data
chrom_size_dat=pd.read_csv(chrom_size_f,sep='\t',header=None)
chrom_size_dict={chr:sz for chr,sz in zip(chrom_size_dat[0].tolist(),chrom_size_dat[1].tolist())}

#Filtering GWAS hits
gwas_hits=gwas_dat.loc[gwas_dat[pval_col] <= thresh,:].sort_values(by=[pval_col],ascending=True)

#Formatting gene_tss data
tss_dat=pd.read_csv(gene_lookup_f,sep='\t',header=None)
tss_regions=[region_to_chr_s_e(i) for i in tss_dat[1].tolist()]
tss_dat['_chr']=[i[0] for i in tss_regions]
tss_dat['_s']=[i[1] for i in tss_regions]

#Storing genes data
all_hit_genes_dat=pd.DataFrame()

#Checking if there are any hits
if gwas_hits.shape[0] > 0:
    final_regions=[]
    lead_snps=[]
    #Loop that keeps going until there are no more GWAS hits left
    while gwas_hits.shape[0] > 0:

        #Defining lowest p-value hits and region around them defined by padding
        lowest_pval_snp_chr,lowest_pval_snp_pos=gwas_hits.sort_values(by=[pval_col],ascending=True).iloc[0,:][[chr_col,pos_col]]
        lowest_pval_snp='{0}:{1}'.format(lowest_pval_snp_chr,lowest_pval_snp_pos)
        lt_boundary = window_left_boundary(lowest_pval_snp_pos,p,chrom_size_dict[lowest_pval_snp_chr])
        rt_boundary = window_right_boundary(lowest_pval_snp_pos,p,chrom_size_dict[lowest_pval_snp_chr])
        lowest_pval_region=chr_s_e_to_region(lowest_pval_snp_chr,lt_boundary,rt_boundary)

        #Keeping only hits that are inside the region
        gwas_hits['_snp_chrpos']=gwas_hits[chr_col].astype(str)+":"+gwas_hits[pos_col].astype(str)
        gwas_hits['_outside']=gwas_hits.agg(lambda x: outside_region(x['_snp_chrpos'],lowest_pval_region),axis=1)
        gwas_hits=gwas_hits.loc[gwas_hits['_outside']==True,:].drop(['_outside','_snp_chrpos'],axis=1)

        #Calculating relevant information for the genes and hits
        hit_genes_dat=tss_dat.loc[(tss_dat['_chr']==lowest_pval_snp_chr) & (tss_dat['_s']>=lt_boundary) & (tss_dat['_s']<=rt_boundary),:]
        tss_lt_boundaries = [window_left_boundary(tss_pos,p,chrom_size_dict[tss_chrom]) for tss_chrom,tss_pos in zip(hit_genes_dat['_chr'],hit_genes_dat['_s'])]
        tss_rt_boundaries = [window_right_boundary(tss_pos,p,chrom_size_dict[tss_chrom]) for tss_chrom,tss_pos in zip(hit_genes_dat['_chr'],hit_genes_dat['_s'])]
        tss_regions=['{0}:{1}-{2}'.format(lowest_pval_snp_chr,tss_lt_boundary,tss_rt_boundry) for tss_lt_boundary,tss_rt_boundry in zip(tss_lt_boundaries,tss_rt_boundaries) ]

        #Adding relevant columns
        hit_genes_dat['tss_region']=tss_regions
        hit_genes_dat['leadsnp_region']=lowest_pval_region
        hit_genes_dat['leadsnp']=lowest_pval_snp

        #Removing irrelevant columns and renaming columns
        hit_genes_dat=hit_genes_dat.drop(columns=['_s','_chr'])
        hit_genes_dat=hit_genes_dat.loc[:,['leadsnp_region','leadsnp',0,1,'tss_region']].rename(columns={0:'gene_id',1:'gene_tss'})

        #Adding to bigger dataframe
        all_hit_genes_dat=pd.concat([all_hit_genes_dat,hit_genes_dat],ignore_index=True)


        print('Finished locus: {0}..'.format(lowest_pval_region))

#Loci only data (no gene info)
loci_dat=all_hit_genes_dat.loc[:,['leadsnp_region','leadsnp']].drop_duplicates()


print('Found {0} loci and {1} genes'.format(loci_dat.shape[0],all_hit_genes_dat.shape[0]))
print('Outputting to: {0}..'.format(output_f))

#Outputting two files {output_prefix}.sigloci and {output_prefix}.siglocigenes
loci_dat.to_csv(output_f+'.sigloci',sep='\t',index=False)
all_hit_genes_dat.to_csv(output_f+'.siglocigenes',sep='\t',index=False)
