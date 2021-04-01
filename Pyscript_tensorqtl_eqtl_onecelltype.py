# salloc --time=2:0:0 --ntasks=1 --cpus-per-task=40  --mem=96Gb --gres=gpu:1 --account=def-awadalla-ac
# cd projects/ctb-awadalla/GROUP/sc_rnaseq/tensorqtl/

# module load StdEnv/2020 gcc/9.3.0
# module load python
# module load arrow/2.0.0

# source env_tensorqtl/bin/activate

# python

import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print('PyTorch {}'.format(torch.__version__))
print('Pandas {}'.format(pd.__version__))

# define paths to data
plink_prefix_path = '/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/data/geno/ATACseq_imputed_322_expression_Fst001_ncancer_Jan2021_eqtl_Mono'
expression_bed = '/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/data/expression/ATACseq_imputed_322_jan2021_expression_AveCounts_nocancer_Bcells.bed.gz'
covariates_file = '/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/data/covariates/ATACseq_imputed_322_Bcells_nocancer_covariates_SexAgeBatch.txt'
prefix = 'ATACseq_Bcells'
#expression_gene_bed = '/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/data/geno/ATACseq_imputed_344_expression_AveCounts_bulk.bed.gz'


#expression_bed = '/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/data/geno/ATACseq_imputed_344_expression_AveCounts_bulk.bed.gz'

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
#phenotype_gene_df, phenotype_gene_pos_df = tensorqtl.read_phenotype_bed(expression_gene_bed)
covariates_df = covariates_df.astype('float64')

#phenotype_pos_df = phenotype_pos_df.astype({'chr': 'object'}).dtypes
#phenotype_pos_df = phenotype_pos_df.astype({'tss': 'float64'}).dtypes

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

phenotype_pos_df["chr"]= phenotype_pos_df["chr"].astype(str)
#phenotype_pos_df["tss"]= phenotype_pos_df["tss"].astype(float)
#phenotype_df = phenotype_df.astype(float)


#cis.map_nominal(genotype_df, variant_df,
#                phenotype_df.loc[phenotype_pos_df['chr']=='1'],
#                phenotype_pos_df.loc[phenotype_pos_df['chr']=='1'],
#                prefix, covariates_df=covariates_df)
                
cis.map_nominal(genotype_df, variant_df,
                phenotype_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'])],
                phenotype_pos_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'])],
                prefix, covariates_df=covariates_df)
                
                

                
                
                
#cis.map_nominal(genotype_df, variant_df,phenotype_df,phenotype_pos_df,prefix, covariates_df=covariates_df)
                

                
#cis.map_nominal(genotype_df, variant_df,
#                phenotype_df,phenotype_pos_df,prefix, covariates_df=covariates_df)
                

# save results
pairs_df_1 = pd.read_parquet('{}.cis_qtl_pairs.1.parquet'.format(prefix))
pairs_df_1.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr1_nominal.txt',sep='\t')

pairs_df_2 = pd.read_parquet('{}.cis_qtl_pairs.2.parquet'.format(prefix))
pairs_df_2.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr2_nominal.txt',sep='\t')

pairs_df_3 = pd.read_parquet('{}.cis_qtl_pairs.3.parquet'.format(prefix))
pairs_df_3.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr3_nominal.txt',sep='\t')

pairs_df_4 = pd.read_parquet('{}.cis_qtl_pairs.4.parquet'.format(prefix))
pairs_df_4.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr4_nominal.txt',sep='\t')

pairs_df_5 = pd.read_parquet('{}.cis_qtl_pairs.5.parquet'.format(prefix))
pairs_df_5.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr5_nominal.txt',sep='\t')

pairs_df_6 = pd.read_parquet('{}.cis_qtl_pairs.6.parquet'.format(prefix))
pairs_df_6.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr6_nominal.txt',sep='\t')

pairs_df_7 = pd.read_parquet('{}.cis_qtl_pairs.7.parquet'.format(prefix))
pairs_df_7.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr7_nominal.txt',sep='\t')

pairs_df_8 = pd.read_parquet('{}.cis_qtl_pairs.8.parquet'.format(prefix))
pairs_df_8.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr8_nominal.txt',sep='\t')

pairs_df_9 = pd.read_parquet('{}.cis_qtl_pairs.9.parquet'.format(prefix))
pairs_df_9.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr9_nominal.txt',sep='\t')

pairs_df_10 = pd.read_parquet('{}.cis_qtl_pairs.10.parquet'.format(prefix))
pairs_df_10.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr10_nominal.txt',sep='\t')

pairs_df_11 = pd.read_parquet('{}.cis_qtl_pairs.11.parquet'.format(prefix))
pairs_df_11.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr11_nominal.txt',sep='\t')

pairs_df_12 = pd.read_parquet('{}.cis_qtl_pairs.12.parquet'.format(prefix))
pairs_df_12.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr12_nominal.txt',sep='\t')

pairs_df_13 = pd.read_parquet('{}.cis_qtl_pairs.13.parquet'.format(prefix))
pairs_df_13.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr13_nominal.txt',sep='\t')

pairs_df_14 = pd.read_parquet('{}.cis_qtl_pairs.14.parquet'.format(prefix))
pairs_df_14.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr14_nominal.txt',sep='\t')

pairs_df_15 = pd.read_parquet('{}.cis_qtl_pairs.15.parquet'.format(prefix))
pairs_df_15.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr15_nominal.txt',sep='\t')

pairs_df_16 = pd.read_parquet('{}.cis_qtl_pairs.16.parquet'.format(prefix))
pairs_df_16.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr16_nominal.txt',sep='\t')

pairs_df_17 = pd.read_parquet('{}.cis_qtl_pairs.17.parquet'.format(prefix))
pairs_df_17.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr17_nominal.txt',sep='\t')

pairs_df_18 = pd.read_parquet('{}.cis_qtl_pairs.18.parquet'.format(prefix))
pairs_df_18.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr18_nominal.txt',sep='\t')

pairs_df_19 = pd.read_parquet('{}.cis_qtl_pairs.19.parquet'.format(prefix))
pairs_df_19.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr19_nominal.txt',sep='\t')

pairs_df_20 = pd.read_parquet('{}.cis_qtl_pairs.20.parquet'.format(prefix))
pairs_df_20.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr20_nominal.txt',sep='\t')

pairs_df_21 = pd.read_parquet('{}.cis_qtl_pairs.21.parquet'.format(prefix))
pairs_df_21.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr21_nominal.txt',sep='\t')

pairs_df_22 = pd.read_parquet('{}.cis_qtl_pairs.22.parquet'.format(prefix))
pairs_df_22.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_chr22_nominal.txt',sep='\t')


### FIX THIS!!!

#for i in 1:22:
 #pairs_df_+i = pd.read_parquet('{}.cis_qtl_pairs.'+i+'.parquet'.format(prefix))
    
#for i in 1:22:
    #pairs_df_+i.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/chromatin/cis_caqtl_chr'+i+'_nominal.txt',sep='\t')


# all genes
# cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df)

# genes on chr18
#cis_df = cis.map_cis(genotype_df, variant_df,
#                     phenotype_df.loc[phenotype_pos_df['chr']=='1'],
#                     phenotype_pos_df.loc[phenotype_pos_df['chr']=='1'],
#                     covariates_df=covariates_df, seed=123456)
                     
cis_df = cis.map_cis(genotype_df, variant_df,
                phenotype_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'])],
                phenotype_pos_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'])],
                 covariates_df=covariates_df, seed=123456)
                 
cis_df.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Tcis_eqtl_Bcells_permutations.txt',sep='\t')

                
trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df, batch_size=10000,
                           return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05)

trans_df.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Ttrans_eqtl_Bcells_nominal.txt',sep='\t')


#cis_df.head()
                

