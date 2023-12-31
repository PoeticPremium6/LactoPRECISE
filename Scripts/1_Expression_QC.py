#Expression Quality Control (Part 1)¶
#This is a template notebook for performing preliminary quality
#control on your organism's expression data.
#Setup
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from os import path
sns.set_style('ticks')

### Inputs
#Enter path of log-TPM, MultiQC, and metadata files here
logTPM_file = path.join('..','data','raw_data','log_tpm_final.csv') # Enter log-TPM filename here
multiqc_file = path.join('..','data','raw_data','multiqc_stats_final.tsv') # Enter MultiQC stats filename here
metadata_file = path.join('..','data','raw_data','L_reuteri_metadata_final.tsv') # Enter metadata filename here

#Load Expression Data
DF_log_tpm = pd.read_csv(logTPM_file,index_col=0).fillna(0)
print('Number of genes:',DF_log_tpm.shape[0])
print('Number of samples:',DF_log_tpm.shape[1])
DF_log_tpm.head()

### Load QC data
#There may be some datasets that failed along the processing pipeline,
#so the number of samples with QC data may be higher than
#the number of samples with expression data.
DF_qc_stats = pd.read_csv(multiqc_file,index_col=0, sep='\t')
print('Number of samples with QC data:',DF_qc_stats.shape[0])

DF_qc_stats.fillna(0,inplace=True)
DF_qc_stats.head()

#Load Metadata
DF_metadata = pd.read_csv(metadata_file,index_col=0,sep='\t')
print('Number of samples with metadata:',DF_metadata.shape[0])
DF_metadata.head()

#Remove extra sample rows
#Ensure that metadata and qc_stats
# data contain all log_tpm sample information.

assert(set(DF_log_tpm.columns) - set(DF_metadata.index) == set())
assert(set(DF_log_tpm.columns) - set(DF_qc_stats.index) == set())

DF_metadata = DF_metadata.loc[DF_log_tpm.columns]
DF_qc_stats = DF_qc_stats.loc[DF_log_tpm.columns]

## Check QC statistics
#FASTQC Quality Control
fastqc_cols = ['per_base_sequence_quality',
       'per_tile_sequence_quality', 'per_sequence_quality_scores',
       'per_base_sequence_content', 'per_sequence_gc_content',
       'per_base_n_content', 'sequence_length_distribution',
       'sequence_duplication_levels', 'overrepresented_sequences',
       'adapter_content']

DF_fastqc = DF_qc_stats[fastqc_cols]
ax = sns.heatmap(DF_fastqc.replace('pass',1).replace('warn',0).replace('fail',-1),
            cmap='RdYlBu',vmax=1.3,vmin=-1.3)
cbar = ax.collections[0].colorbar
cbar.set_ticks([-1,0,1])
cbar.set_ticklabels(['fail','warn','pass'])

#The following four categories are the most important:
#per_base_sequence_quality
#per_sequence_quality_scores
#per_base_n_content
#adapter_content
#If a sample does not pass any of these four categories, discard the sample.
fastqc_fail_cols = ['per_base_sequence_quality','per_sequence_quality_scores','per_base_n_content','adapter_content']
DF_failed_fastqc = DF_fastqc[fastqc_fail_cols][(DF_fastqc[fastqc_fail_cols] != 'pass').any(axis=1)]
DF_failed_fastqc[fastqc_fail_cols]
#Mark passes samples
DF_metadata['passed_fastqc'] = ~DF_metadata.index.isin(DF_failed_fastqc.index)

### Number of aligned reads
#The following histogram shows how many reads map to coding sequences (i.e. mRNA).
#Too few aligned reads reduces the sensitivity of the resulting data.
min_mrna_reads = 500000 # Minimum number of reads mapped to mRNA (500,000)
fig,ax = plt.subplots()
ax.hist(DF_qc_stats['Assigned']/1e6,bins=50,alpha=0.8)
ymin,ymax = ax.get_ylim()
ax.vlines(min_mrna_reads/1e6,ymin,ymax,color='r')
ax.set_ylim((ymin,ymax))
ax.set_xlabel('# Reads (M)',fontsize=18)
ax.set_ylabel('# Samples',fontsize=18)
ax.set_title('Number of reads mapped to CDS',fontsize=18)

#Identify samples with poor read depth
DF_failed_mrna = DF_qc_stats[DF_qc_stats['Assigned'] < min_mrna_reads].sort_values('Assigned')
DF_failed_mrna

#Mark samples that passed.
DF_metadata['passed_reads_mapped_to_CDS'] = ~DF_metadata.index.isin(DF_failed_mrna.index)

##Examine Global Correlations
#Only examine data that passed the first two steps.
metadata_passed_step2 = DF_metadata[DF_metadata[['passed_fastqc','passed_reads_mapped_to_CDS']].all(axis=1)]
DF_log_tpm_passed_step2 = DF_log_tpm[metadata_passed_step2.index]

#A clustermap is a great way to visualize the global correlations between one sample and all others. The global_clustering function uses hierarchical clustering to identify specific clusters in the clustermap. The optional arguments are:
#threshold: Threshold used to extract clusters from the hierarchy. To increase the number of clusters, decrease the value of threshold.
#To decrease the number of clusters, increase the value of threshold (default: 0.3)
#figsize: A tuple describing the length and width of the final clustermap. A larger figsize can make x and y-axis labels clearer.
#xticklabels: Show NCBI SRA accession numbers on the x-axis
#yticklabels: Show NCBI SRA accession numbers on the y-axis

import scipy.cluster.hierarchy as sch
import matplotlib.patches as patches

def global_clustering(data, threshold=0.3, xticklabels=False, yticklabels=False, figsize=(9,9)):

    # Retrieve clusters using fcluster
    corr = data.corr()
    corr.fillna(0,inplace=True)
    dist = sch.distance.pdist(corr)
    link = sch.linkage(dist, method='complete')
    clst = pd.DataFrame(index=data.columns)
    clst['cluster'] = sch.fcluster(link, threshold * dist.max(), 'distance')

    # Get colors for each cluster
    cm = plt.cm.get_cmap('tab20')
    cluster_colors = dict(zip(clst.cluster.unique(), cm.colors))
    clst['color'] = clst.cluster.map(cluster_colors)

    print('Number of cluster: ', len(cluster_colors))

    legend_items = [patches.Patch(color=c, label=l) for l,c in cluster_colors.items()]

    sns.set(rc={'figure.facecolor':'white'})

    clst_map = sns.clustermap(data.corr(),
                              figsize=figsize,
                              row_linkage=link,
                              col_linkage=link,
                              col_colors=clst.color,
                              yticklabels=yticklabels,
                              xticklabels=xticklabels,
                              vmin=0,
                              vmax=1)

    legend = clst_map.ax_heatmap.legend(loc='upper left',
                                        bbox_to_anchor=(1.01,0.85),
                                        handles=legend_items,
                                        frameon=True)

    legend.set_title(title='Clusters',prop={'size':10})

    return clst['cluster']
clusters = global_clustering(DF_log_tpm_passed_step2)

#Print out all clusters
with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(clusters)
#Select clusters to remove
remove_clusters = [1]
remove_clusters = [3]
passed_global_corr = clusters[~clusters.isin(remove_clusters)].index
#The following code can be adapted to see the NCBI SRA accession for samples in each cluster.
clusters[clusters == 1-2]
#Re-cluster samples to ensure all outliers were removed.
DF_log_tpm_passed_step3 = DF_log_tpm[passed_global_corr]
clusters = global_clustering(DF_log_tpm_passed_step3)

#Once you are satisfied with your dataset, mark the samples that passed the global correlation
DF_metadata['passed_global_correlation'] = DF_metadata.index.isin(passed_global_corr)
DF_metadata.head()

##Remove Failed Samples
qc_columns = ['passed_fastqc',
              'passed_reads_mapped_to_CDS',
              'passed_global_correlation']
pass_qc = DF_metadata[qc_columns].all(axis=1)
DF_metadata_passed = DF_metadata[pass_qc]
_,_,pcts = plt.pie(pass_qc.value_counts().reindex([False,True]),
        labels = ['Failed','Passed'],
        colors=['tab:red','tab:blue'],
        autopct='%.0f%%',textprops={'size':16});

#Save current metadata
metadata_all_qc_file = path.join('..', 'data', 'interim', 'metadata_qc_part1_all.tsv') # Enter filename for full metadata QC file
metadata_qc_file = path.join('..', 'data', 'interim', 'metadata_qc_part1.tsv') # Enter filename for metadata QC file with only passing datasets
DF_metadata.to_csv(metadata_all_qc_file, sep='\t')
DF_metadata_passed.to_csv(metadata_qc_file, sep='\t')


# Expression Quality Control (Part 2)
#This is a template notebook for performing the final quality control on your organism's expression data. ' \
#'This requires a curated metadata sheet.

import itertools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from os import path
from scipy import stats
from tqdm.notebook import tqdm
sns.set_style('ticks')

#Input files from curation QC
logTPM_file = path.join('..','data','raw_data','log_tpm_final.csv') # Enter log-TPM filename here
all_metadata_file = path.join('..','data','interim','metadata_qc_part1_all.tsv') # Enter full metadata filename here
metadata_file = path.join('..','data','interim','metadata_qc_curated.tsv') # Enter curated metadata filename here

##Load Expression Data
DF_log_tpm = pd.read_csv(logTPM_file,index_col=0).fillna(0)
print('Number of genes:',DF_log_tpm.shape[0])
print('Number of samples:',DF_log_tpm.shape[1])
DF_log_tpm.head()

#Load Metadata
DF_metadata = pd.read_csv(metadata_file,index_col=0,sep='\t')
print('Number of samples with curated metadata:',DF_metadata.shape[0])
DF_metadata.head()
DF_metadata_all = pd.read_csv(all_metadata_file,index_col=0,sep='\t')

#Remove samples due to poor metadata
#After curation, some samples either
#did not have enough replicates or metadata to warrant inclusion in this database.
DF_metadata_passed_step4 = DF_metadata[~DF_metadata.skip.fillna(False)].copy()
print('New number of samples with curated metadata:',DF_metadata_passed_step4.shape[0])
DF_metadata_passed_step4.head()
#Check curation
#Since manual curation is error-prone, we want to make sure that all samples have labels for their project and condition. In addition, there should only be one reference condition in each project, and it should be in the project itself.

#Any samples that fail these checks will be printed below.
assert(DF_metadata_passed_step4.project.notnull().all())
assert(DF_metadata_passed_step4.condition.notnull().all())

for name,group in DF_metadata_passed_step4.groupby('project'):
    ref_cond = group.reference_condition.unique()

    # Ensure that there is only one reference condition per project
    if not len(ref_cond) == 1:
        print('Multiple reference conditions for:, name')

    # Ensure the reference condition is in fact in the project
    ref_cond = ref_cond[0]
    if not ref_cond in group.condition.tolist():
        print('Reference condition not in project:', name)

#Next, make a new column called full_name that gives every experimental condition a unique,
# human-readable identifier.
DF_metadata_passed_step4['full_name'] = DF_metadata_passed_step4['project'].str.cat(DF_metadata_passed_step4['condition'],sep=':')

#Remove samples with only one replicate
#Find sample names which have at least two replicates
counts = DF_metadata_passed_step4.full_name.value_counts()
keep_samples = counts[counts >= 2].index
print(keep_samples[:5])

#Keep only these samples
DF_metadata_passed_step4 = DF_metadata_passed_step4[DF_metadata_passed_step4.full_name.isin(keep_samples)]
print('New number of samples with curated metadata:',DF_metadata_passed_step4.shape[0])
DF_metadata_passed_step4.head()

### Save this information to the full metadata dataframe
DF_metadata_all['passed_curation'] = DF_metadata_all.index.isin(DF_metadata_passed_step4.index)
#Check correlations between replicates
#Remove failed data from log_tpm files
DF_log_tpm = DF_log_tpm[DF_metadata_passed_step4.index]
### Compute Pearson R Score
#Biological replicates should have a Pearson R correlation above 0.95.
#For samples with more than 2 replicates, the replicates must have R >= 0.95 with at least one other replicate or it will be dropped.
#The correlation threshold can be changed below:
rcutoff = 0.95


#The following code computes correlations between all samples and collects correlations
#between replicates and non-replicates.
rep_corrs = {}
rand_corrs = {}
num_comparisons = len(DF_metadata_passed_step4)*(len(DF_metadata_passed_step4)-1)/2
for exp1,exp2 in tqdm(itertools.combinations(DF_metadata_passed_step4.index,2),total=num_comparisons):
    if DF_metadata_passed_step4.loc[exp1,'full_name'] == DF_metadata_passed_step4.loc[exp2,'full_name']:
        rep_corrs[(exp1,exp2)] = stats.pearsonr(DF_log_tpm[exp1],DF_log_tpm[exp2])[0]
    else:
        rand_corrs[(exp1,exp2)] = stats.pearsonr(DF_log_tpm[exp1],DF_log_tpm[exp2])[0]

#Correlations can be plotted on a histogram

fig,ax = plt.subplots(figsize=(5,5))
ax2 = ax.twinx()
ax2.hist(rep_corrs.values(),bins=50,range=(0.2,1),alpha=0.8,color='green',linewidth=0)
ax.hist(rand_corrs.values(),bins=50,range=(0.2,1),alpha=0.8,color='blue',linewidth=0)
ax.set_title('Pearson R correlation between experiments',fontsize=14)
ax.set_xlabel('Pearson R correlation',fontsize=14)
ax.set_ylabel('Different Conditions',fontsize=14)
ax2.set_ylabel('Known Replicates',fontsize=14)

med_corr = np.median([v for k,v in rep_corrs.items()])
print('Median Pearson R between replicates: {:.2f}'.format(med_corr))

#Remove samples without any high-correlation replicates
dissimilar = []
for idx, grp in DF_metadata_passed_step4.groupby('full_name'):
    ident = np.identity(len(grp))
    corrs = (DF_log_tpm[grp.index].corr() - ident).max()
    dissimilar.extend(corrs[corrs<rcutoff].index)

# Save this information in both the original metadata dataframe and the new metadata dataframe
DF_metadata_all['passed_replicate_correlations'] = ~DF_metadata_all.index.isin(dissimilar)
DF_metadata_passed_step4['passed_replicate_correlations'] = ~DF_metadata_passed_step4.index.isin(dissimilar)

DF_metadata_final = DF_metadata_passed_step4[DF_metadata_passed_step4['passed_replicate_correlations']]
print('# Samples that passed replicate correlations:',len(DF_metadata_final))
## Check that reference conditions still exist
#If a reference condition was removed due to poor replicate correlations,
# a new reference condition needs to be defined.
#Again, any samples that fail these checks will be printed below.

project_exprs = []
for name,group in DF_metadata_final.groupby('project'):
    # Get reference condition
    ref_cond = group.reference_condition.iloc[0]

    # Ensure the reference condition is still in the project
    if ref_cond not in group.condition.tolist():
        print('Reference condition missing from:', name)

    # Check that each project has at least two conditions (a reference and at least one test condition)
    if len(group.condition.unique()) <= 1:
        print('Only one condition in:', name)
#If necessary, choose a new condition for failed projects and re-run notebook.
## Normalize dataset to reference conditions
DF_log_tpm_final = DF_log_tpm[DF_metadata_final.index]

project_exprs = []
for name,group in DF_metadata_final.groupby('project'):

    # Get reference condition
    ref_cond = group.reference_condition.iloc[0]

    # Get reference condition sample ids
    ref_samples = group[group.condition == ref_cond].index

    # Get reference condition expression
    ref_expr = DF_log_tpm_final[ref_samples].mean(axis=1)

    # Subtract reference expression from project
    project_exprs.append(DF_log_tpm_final[group.index].sub(ref_expr,axis=0))

DF_log_tpm_norm = pd.concat(project_exprs,axis=1)

#Save final datasets
logTPM_qc_file = path.join('..','data','processed_data','log_tpm_20230223.csv')
logTPM_norm_file = path.join('..','data','processed_data','log_tpm_norm_20230223.csv')
final_metadata_file = path.join('..','data','processed_data','metadata_20230223.tsv')
final_metadata_all_file = path.join('..','data','interim','metadata_qc_part2_all_20230223.tsv')

DF_log_tpm_final.to_csv(logTPM_qc_file)
DF_log_tpm_norm.to_csv(logTPM_norm_file)
DF_metadata_final.to_csv(final_metadata_file, sep='\t')
DF_metadata_all.to_csv(final_metadata_all_file, sep='\t')
