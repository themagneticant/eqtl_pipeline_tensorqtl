#salloc --time=2:0:0 --ntasks=1 --cpus-per-task=40  --mem=96Gb --gres=gpu:1 --account=def-awadalla-ac
#cd projects/ctb-awadalla/GROUP/sc_rnaseq/tensorqtl/

 #module load StdEnv/2020 gcc/9.3.0
 #module load python
 #module load arrow/2.0.0

 #source env_tensorqtl/bin/activate

 #python
import sys
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print('PyTorch {}'.format(torch.__version__))
print('Pandas {}'.format(pd.__version__))

# define paths to data

#prefix = 'Bcells'
prefix = (sys.argv[1])
prefix1 = prefix + '_canonical'
suffix1 = prefix + '_modelPCAgeSexBatch_canonical'

prefix2 = prefix + '_interactionAge'
suffix2 = prefix + '_modelPCAgeSexBatch_interactionAge'

prefix22 = prefix + '_interactionIRS'
suffix22 = prefix + '_modelPCAgeSexBatch_interactionIRS'


prefix3 = prefix + '_interactionIRS_young'
suffix3 = prefix + '_modelPCAgeSexBatch_interactionIRS_young'

prefix4 = prefix + '_interactionIRS_old'
suffix4 = prefix + '_modelPCAgeSexBatch_interactionIRS_old'

plink_prefix_path = '/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/data/geno/all_samples/ATACseq_imputed_geno_HW10m16_'+ prefix
expression_bed = '/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/data/expression/ATACseq_imputed_expr_'+ prefix +'.bed.gz'
covariates_file1 = '/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/data/covariates/ATACseq_imputed_covariates_'+ prefix +'_SexAgeBatch.txt'
covariates_file2 = '/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/data/covariates/ATACseq_imputed_covariates_'+ prefix +'_SexAgeBatchIRS.txt'

#interaction_s = pd.Series(data=covariates_df['SDC_AGE_CALC'], index=covariates_df.index.values)


# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df1 = pd.read_csv(covariates_file1, sep='\t', index_col=0).T
#phenotype_gene_df, phenotype_gene_pos_df = tensorqtl.read_phenotype_bed(expression_gene_bed)
covariates_df1 = covariates_df1.astype('float64')
#interaction_s = pd.Series(data=covariates_df['IRSnoAge'], index=covariates_df.index.values)
covariates_df2 = pd.read_csv(covariates_file2, sep='\t', index_col=0).T
covariates_df2 = covariates_df2.astype('float64')


young = covariates_df2['SDC_AGE_CALC'] < 46
covariates_df2_young = covariates_df2[young]

old = covariates_df2['SDC_AGE_CALC'] > 64
covariates_df2_old = covariates_df2[old]

interaction_s_young = pd.Series(data=covariates_df2_young['IRSnoAge'], index=covariates_df2_young.index.values)
interaction_s_old = pd.Series(data=covariates_df2_old['IRSnoAge'], index=covariates_df2_old.index.values)

interaction_age = pd.Series(data=covariates_df1['SDC_AGE_CALC'], index=covariates_df1.index.values)
interaction_IRS = pd.Series(data=covariates_df2['IRSnoAge'], index=covariates_df2.index.values)


# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

phenotype_pos_df["chr"]= phenotype_pos_df["chr"].astype(str)


##########################
# Nominal eqtls - canonical
##########################

#cis.map_nominal(genotype_df, variant_df,
#                phenotype_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'])],
#                phenotype_pos_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'])],
#                prefix1, covariates_df=covariates_df1)
                

                
##########################
# save results
##########################

#d = {}
#for i in range(1, 22):
#    d["pairs_df_{0}".format(i)] = pd.read_parquet('{p}.cis_qtl_pairs.{it}.parquet'.format(p=prefix1, it=i))
#    d["pairs_df_{0}".format(i)].to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/filters/canonical/cis_nominal_chr{it}_{p}.txt'.format(p=prefix1, it=i),sep='\t')

             
             
##########################
# permutations
##########################

#
#cis_df = cis.map_cis(genotype_df, variant_df,
#                phenotype_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'])],
#                phenotype_pos_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'])],
#                 covariates_df=covariates_df1, seed=123456)
##tensorqtl.calculate_qvalues(cis_df, qvalue_lambda=0.85)
#
#cis_df.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/filters/canonical/cis_eqtl_permutations' + suffix1 + '.txt',sep='\t')
##
##

##########################
# Trans eqtls
##########################
#
#trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df1, batch_size=10000,
#                           return_sparse=True, pval_threshold=1e-5, maf_threshold=0.05)
#
#trans_df.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/filters/trans/trans_eqtl_nominal_' + suffix1 + '.txt',sep='\t')

##########################
# Independant eqtls as in Gtex paper
# qvalue is not working - need to install r2py on beluga and it keeps failing
##########################

#indep_df = cis.map_independent(genotype_df, variant_df, cis_df,
#                              phenotype_df, phenotype_pos_df, covariates_df)


#indep_df.to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/Cis_independant_eqtl_nominal_' + suffix + '.txt',sep='\t')
#cis_df.head()
                
##########################
# interaction
##########################

######## AGE #############

cis.map_nominal(genotype_df, variant_df,
                phenotype_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']),covariates_df1.index],
                phenotype_pos_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'])],
                prefix2, covariates_df=covariates_df1,interaction_s=interaction_age, maf_threshold_interaction=0.05,run_eigenmt=True)

# save

d = {}
for i in range(1, 22):
    d["pairs_df_{0}".format(i)] = pd.read_parquet('{p}.cis_qtl_pairs.{it}.parquet'.format(p=prefix2, it=i))
    d["pairs_df_{0}".format(i)].to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/filters/interaction/Age/cis_nominal_chr{it}_{p}.txt'.format(p=prefix2, it=i),sep='\t')


####################################
####################################
# IRS #
#interaction_IRS
cis.map_nominal(genotype_df, variant_df,
                phenotype_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']),covariates_df2.index],
                phenotype_pos_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'])],
                prefix22, covariates_df=covariates_df2,interaction_s=interaction_IRS, maf_threshold_interaction=0.05,run_eigenmt=True)

# save

d = {}
for i in range(1, 22):
    d["pairs_df_{0}".format(i)] = pd.read_parquet('{p}.cis_qtl_pairs.{it}.parquet'.format(p=prefix22, it=i))
    d["pairs_df_{0}".format(i)].to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/filters/interaction/IRS/cis_nominal_chr{it}_{p}.txt'.format(p=prefix22, it=i),sep='\t')

####################################
####################################
# Subgroups
####################################

#old

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path,select_samples=covariates_df2_old.index)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

phenotype_pos_df["chr"]= phenotype_pos_df["chr"].astype(str)


cis.map_nominal(genotype_df, variant_df,
                phenotype_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']),covariates_df2_old.index],
                phenotype_pos_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'])],
                prefix4, covariates_df=covariates_df2_old,interaction_s=interaction_s_old, maf_threshold_interaction=0.05,run_eigenmt=True)

# save

d = {}
for i in range(1, 22):
    d["pairs_df_{0}".format(i)] = pd.read_parquet('{p}.cis_qtl_pairs.{it}.parquet'.format(p=prefix4, it=i))
    d["pairs_df_{0}".format(i)].to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/filters/interaction/IRS_old/cis_nominal_chr{it}_{p}.txt'.format(p=prefix4, it=i),sep='\t')
    
###############################
#young

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path,select_samples=covariates_df2_young.index)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

phenotype_pos_df["chr"]= phenotype_pos_df["chr"].astype(str)


cis.map_nominal(genotype_df, variant_df,
                phenotype_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']),covariates_df2_young.index],
                phenotype_pos_df.loc[phenotype_pos_df['chr'].isin(['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'])],
                prefix3, covariates_df=covariates_df2_young,interaction_s=interaction_s_young, maf_threshold_interaction=0.05,run_eigenmt=True)

# save

d = {}
for i in range(1, 22):
    d["pairs_df_{0}".format(i)] = pd.read_parquet('{p}.cis_qtl_pairs.{it}.parquet'.format(p=prefix3, it=i))
    d["pairs_df_{0}".format(i)].to_csv('/lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/results/expression/filters/interaction/IRS_young/cis_nominal_chr{it}_{p}.txt'.format(p=prefix3, it=i),sep='\t')

