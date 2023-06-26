#######FIGURE 5 ###########
#########ArgR_Bitomics#####
#ArgR has been suggested to have some interesting biological functions
import sys
sys.path.append('bitome2')
from bitome.core import Bitome
from feature_functions import *
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sparse
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_validate, RandomizedSearchCV, LeaveOneOut, train_test_split, KFold
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn.utils.class_weight import compute_class_weight
from sklearn.utils import resample
from sklearn.preprocessing import StandardScaler
import copy
import random
import statistics
import matplotlib.patheffects as path_effects

#Let's plot the gene weights for both L. reuteri & Ecoli
#Let's plot the gene weights for both L. reuteri & Ecoli

fig, ax = plt.subplots(1, 2, figsize=(15, 5))  # increased figure size for additional space

# Plot the first figure
plot_gene_weights(ecoli_data, 'ArgR', show_labels=False, ax=ax[0])
title = ax[0].set_title("Ecoli ArgR")
title.set_path_effects([path_effects.withStroke(linewidth=1.5, foreground='black')])
ax[0].tick_params(axis='both', which='major', labelsize=10)
leg1 = ax[0].legend(fontsize='small', bbox_to_anchor=(1.05, 1), loc='upper left')  # shift the legend to right of plot

# Plot the second figure
plot_gene_weights(ica_data, 'ArgR', show_labels=False, ax=ax[1])
title = ax[1].set_title("Lreuteri ArgR")
title.set_path_effects([path_effects.withStroke(linewidth=1.5, foreground='black')])
ax[1].tick_params(axis='both', which='major', labelsize=10)
leg2 = ax[1].legend(fontsize='x-small', bbox_to_anchor=(1.05, 1), loc='upper left')  # shift the legend to right of plot

# Adjust space between plots and increase left margin
plt.subplots_adjust(left=0.05, wspace=0.8, right=0.75)

# Save the figure
plt.savefig("Shared_ArgR.png", dpi=300, bbox_inches='tight')



#Add L. reuteri dataframe
gene_table = pd.read_csv('gene_table.csv')
df = gene_table[['locus_tag','accession','start','end','strand']]
df

# Grouping genes into operons
# Load gene information into a dataframe
gene_data = df
d=[]
gene_clusters=[]
accession =[]

# Set a distance threshold (in base pairs)
distance_threshold = 500

# Group genes by accession ID
accession_groups = gene_data.groupby('accession')

# Iterate over each accession group
for accession_id, accession_group in accession_groups:
    # Sort genes by start position
    sorted_genes = accession_group.sort_values('start')
    
    # Initialize variables to track operon membership
    current_operon = []
    last_gene_end = None
    
    # Iterate over genes in the sorted list
    for _, gene in sorted_genes.iterrows():
        # Check if this gene is close enough to the last gene to be in the same operon
        if last_gene_end is not None and \
           gene['start'] - last_gene_end <= distance_threshold and \
           gene['strand'] == last_gene_strand:
            # Add the gene to the current operon
            current_operon.append(gene['locus_tag'])
        else:
            # Start a new operon
            if current_operon:
                d.append([accession_id,[", ".join(current_operon)]])
                gene_clusters.append(",".join(current_operon).split(","))
                accession.append(accession_id)
                print(f'Operon for {accession_id}: {",".join(current_operon)}')
            current_operon = [gene['locus_tag']]
        
        # Update tracking variables
        last_gene_end = gene['end']
        last_gene_strand = gene['strand']
    
    # Print the last operon
    if current_operon:
        d.append([accession_id,[", ".join(current_operon)]])
        gene_clusters.append(",".join(current_operon).split(","))
        accession.append(accession_id)
        print(f'Operon for {accession_id}: {", ".join(current_operon)}')

# gene_clusters is a nested list of gene names for each operon/cluster

# create a dictionary to map locus tags to their corresponding start, end, and strand values
locus_dict = dict(zip(gene_table['locus_tag'], gene_table[['start', 'end', 'strand']].values))

# create an empty list to store the start, end, and strand values for each cluster
cluster_values = []

for cluster in gene_clusters:
    # get the start, end, and strand values for each gene in the cluster
    values = [locus_dict[gene_table.loc[gene_table['locus_tag'] == gene, 'locus_tag'].iloc[0]] for gene in cluster]
    
    # find the smallest start and largest end values based on the strand value
    # positive strand - operon start site was determined based on the start position of the gene with the lowest start position
    # negative strand - operon start site was determined based on the end position of the gene with the highest end position 
    start, end = None, None
    for value in values:
        if value[2] == '+':
            if start is None or value[0] < start:
                start = value[0]
        elif value[2] == '-':
            if end is None or value[1] > end:
                end = value[1]
    
    # append the start, end, and strand values for the cluster to cluster_values
    cluster_values.append([start, end, values[0][2]])


from itertools import chain
result = [[value]*len(cluster) for value, cluster in zip(cluster_values, gene_clusters)]
cluster_list = [[i+1]*len(sublist) for i, sublist in enumerate(result)]
cluster_list = list(chain(*cluster_list))
df_1 = pd.DataFrame(list(chain(*gene_clusters)))
df_2 = pd.DataFrame(list(chain(*result)))
merged = pd.concat([df_1, df_2], axis=1)
merged.columns = ['locus_tag','start','end','strand']
merged_df = pd.merge(merged, df, on="locus_tag", how="inner")
merged_df = merged_df[['locus_tag','start_x','end_x','strand_x','accession']]
merged_df


# ArgR motif scores
tf_pwm_db = rpwm('motif_pwm_db.txt')
argr_pwm_old = tf_pwm_db['ArgR']
argr_pssm = {base: [pos_dict[base] for pos_dict in argr_pwm_old] for base in 'ATCG'}
score_argr=[]
for i,j,k,l in zip(merged_df['accession'],merged_df['strand_x'], merged_df['start_x'],merged_df['end_x']):
    bitome_fasta = Bitome(i+'.fasta')
    if j == '+':
        if int(k) > 300:
            score_argr.append(int(bitome_fasta.motif_search(k-300, k+50, 1, argr_pssm, n_best_matches=1)['log_odds']))
        else:
            score_argr.append(0)
    elif j =='-':
        score_argr.append(int(bitome_fasta.motif_search(l-50, l+300, -1, argr_pssm, n_best_matches=1)['log_odds']))
merged_df['ArgR'] = score_argr

#creating target labels
target_labels=[]

#target_argr can be modified to have either the former locus_tags or the updated locus_tags (after re-running ICA)
#this version has the updated locus_tags
target_argr =['LMB90_RS07505','LMB90_RS07510','LMB90_RS08635','LMB90_RS03500','LMB90_RS03505','LMB90_RS05535','LMB90_RS03525','LMB90_RS03530','LMB90_RS10070','LMB90_RS10075','LMB90_RS09710','LMB90_RS04610','LMB90_RS04615','LMB90_RS04690','LMB90_RS04695','LMB90_RS04700']
for i in merged_df['locus_tag']:
    if i in target_argr:
        target_labels.append(1)
    else:
        target_labels.append(0)
merged_df['target_labels'] = target_labels
target_label={}
target_label['ArgR'] = target_labels

# Building the ML model

df_score = merged_df[['ArgR']]
X_eng = df_score.values
X_eng = StandardScaler().fit_transform(np.array(X_eng).reshape(-1, 1))
X_eng.shape

df_score.set_index(merged_df['locus_tag'], inplace=True)
df_score_stdd = pd.DataFrame(data = X_eng, columns=df_score.columns)

subset_dict = {}
subset_dict['ArgR'] = target_argr
search_cluster_df = pd.DataFrame(data=cluster_list, columns=['cluster'], index=df_score.index)
IM_promoters = {}
for IM in subset_dict:
    gene_list = [gene for gene in subset_dict[IM] if gene in list(df_score.index)]
    cluster_df = search_cluster_df.loc[gene_list, 'cluster'] 
    promoter_list = []
    for cluster_num in pd.unique(cluster_df):
        each_promoter_list = list(cluster_df[cluster_df==cluster_num].index)
        promoter_list.append(each_promoter_list)
    # ArgR iModulon genes grouped into potential promoters/operons
    IM_promoters[IM] = promoter_list
# Grouping all genes into promoters/operons
cluster_df = search_cluster_df['cluster'] 
All_promoters = []
for cluster_num in pd.unique(cluster_df):
    each_promoter_list = list(cluster_df[cluster_df==cluster_num].index)
    All_promoters.append(each_promoter_list)
All_promoters_df = pd.DataFrame(columns=['first_gene', 'genes'], index=range(len(All_promoters)))
ind = 0
for promoter in All_promoters:
    All_promoters_df.iloc[ind, 0] = promoter[0]
    All_promoters_df.iloc[ind, 1] = promoter
    ind += 1 
All_promoters_df.to_csv('gene_groupby_promoter.csv')

models_to_try = {
    
#     'SVM': LinearSVC(
#     penalty = 'l1',
#     class_weight='balanced',
#     dual=False,
#     random_state=42,
#     verbose=0,
#     ),
    
    'LR': LogisticRegression(
        penalty='elasticnet',
        solver='saga',
        class_weight='balanced',
        l1_ratio = 0.5,
        random_state=42
    )
    
#     'RF': RandomForestClassifier(
#         class_weight='balanced',
#         random_state=42,
#         verbose=0
#     )
}

# from imblearn.under_sampling import ClusterCentroids, CondensedNearestNeighbour
from imblearn.over_sampling import SMOTE, ADASYN, RandomOverSampler, KMeansSMOTE, SVMSMOTE
from imblearn.combine import SMOTETomek
from sklearn.decomposition import PCA
import seaborn as sns

def resample_nfold(oversampler, num_neighbors=5, SMOTE_only = False, exact_repeats_boost=0.1,
                features_table=df_score, target_labels=target_label
                #, undersampler=None
                 ,  final_ratio=1, N_CV=5, random_state=4
                ):
    '''
    oversampler: selected oversampling method from imblearn
    num_neighbor: Number of k-nearest neighbor used to perform SMOTE oversampling
    final_ratio: If undersample=True, the ratio of the n_negative_aftersampling to n_negative_beforesampling.
    
    '''
    
    result_df = pd.DataFrame(
        columns=['model', 'im', 'train_auc', 'test_auc']
    )

    for model_name, model in models_to_try.items():

        for IM_name, y_labels in target_labels.items():          
            # skip the promoters that are too small
            if len(IM_promoters[IM_name])<N_CV:
                continue
            
            #index_num = 0
            print(f'{model_name}: {IM_name}')

            temp_eng_features = features_table.copy()
            temp_eng_features['label'] = y_labels        

            ingroup_set = IM_promoters[IM_name]
            outgroup_set = [l for l in All_promoters if not any(item in l for item in sum(IM_promoters[IM_name], [])) ]

            kf = KFold(n_splits=N_CV)

            train_scores = []
            test_scores = []
            
            for (out_index, in_index) in zip(kf.split(outgroup_set), kf.split(ingroup_set)):
                
                test_promoters_out = sum([outgroup_set[i] for i in out_index[1]], [])          
                train_promoters_out = sum([outgroup_set[i] for i in out_index[0]], [])
                
                test_promoters_in = sum([ingroup_set[i] for i in in_index[1]], [])          
                train_promoters_in = sum([ingroup_set[i] for i in in_index[0]], [])
                
                test_promoters = test_promoters_out + test_promoters_in
                train_promoters = train_promoters_out + train_promoters_in
                
                train_X = temp_eng_features.loc[train_promoters].iloc[:, :-1]
                test_X = temp_eng_features.loc[test_promoters].iloc[:, :-1]

                train_Y = temp_eng_features.loc[train_promoters, 'label']
                test_Y = temp_eng_features.loc[test_promoters, 'label']
            
            
                if SMOTE_only == True:
                    # Only use RandomOverSampling to boost the samples number over num_neighbors+1
                    if len([l for l in train_Y if l==1])<=num_neighbors+1:
                        
                        exact_repeats_boost = (num_neighbors*2)/len([l for l in train_Y if l==0])                
                        rs = RandomOverSampler(sampling_strategy = exact_repeats_boost,random_state=random_state)
                        train_X, train_Y = rs.fit_resample(train_X, train_Y)
                    
                else:
                    # Use RandomOverSampling to boost the positive samples number to a certain pos/neg ratio.
                    rs = RandomOverSampler(sampling_strategy = exact_repeats_boost,random_state=random_state)
                    train_X, train_Y = rs.fit_resample(train_X, train_Y)             
                
                if exact_repeats_boost != 1:
                    # Oversampling by SMOTE
                    if oversampler == SMOTETomek:
                        over = oversampler(sampling_strategy = final_ratio, smote = SMOTE(k_neighbors = num_neighbors),
                                           random_state=random_state)
                        train_X, train_Y = over.fit_resample(train_X, train_Y)
                        
                    else:
                        over = oversampler(sampling_strategy = final_ratio, k_neighbors = num_neighbors,
                                           random_state=random_state)
                        train_X, train_Y = over.fit_resample(train_X, train_Y)
                    

                model.fit(train_X, train_Y)

                if model_name == 'SVM':
                    train_auc = roc_auc_score(train_Y, model.decision_function(train_X))
                    test_auc = roc_auc_score(test_Y, model.decision_function(test_X))

                else:
                    train_auc = roc_auc_score(train_Y, model.predict_proba(train_X)[:,1])
                    test_auc = roc_auc_score(test_Y, model.predict_proba(test_X)[:,1])

                train_scores.append(train_auc)
                test_scores.append(test_auc)
                
                
            sub_result_df = pd.DataFrame(
            data={
                'model': [model_name] * N_CV,
                'im': [IM_name] * N_CV,
                'train_auc': train_scores,
                'test_auc': test_scores
            }

            )

            result_df = pd.concat([result_df,sub_result_df], axis=0)
        

    return result_df

result = resample_nfold(oversampler=SMOTETomek, SMOTE_only = True, num_neighbors=5, target_labels=target_label, 
                              random_state=4)

result_LR = result[result['model']=='LR']
plt.figure( figsize=(5,5))
sns.boxplot(x='im', y='test_auc', data=result_LR)
plt.axhline(y= 0.5, color = 'r', ls = '--')
plt.axhline(y= 0.8, color = 'r', ls = '--')
plt.ylabel('AUC_ROC score')
plt.xlabel('ArgR')
plt.savefig("ArgR_AUC_ROC_score.pdf", format="pdf", bbox_inches="tight")

import matplotlib
y = df_score['ArgR']
x = target_label['ArgR']
for i in range(len(x)):
    if x[i] == 0:
        x[i] = 'Random'
    else:
        x[i] = 'ArgR iModulon'
sns.boxplot(x=x,y=y)

sns.boxplot(x=x,y=df_score['ArgR'])
plt.ylabel("Motif score", fontsize=10)
plt.savefig("ArgR_motif_score_updated.pdf", format="pdf", bbox_inches="tight")
