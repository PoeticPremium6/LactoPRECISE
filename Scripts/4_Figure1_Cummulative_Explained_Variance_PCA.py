#Figure 1#
#Cummulative Explained Variance & PCA
from sklearn.decomposition import PCA
import numpy as np
#Load Data
DF_metadata = path.join('Metadata_Final.tsv') # Enter curated metadata filename here
DF_metadata = pd.read_csv(DF_metadata,index_col=0,sep='\t')

DF_log_tpm = pd.read_csv(os.path.join('log_tpm_final.csv'),index_col=0)
#First, compute principal components
pca = PCA()
DF_weights = pd.DataFrame(pca.fit_transform(DF_log_tpm.T),index=DF_log_tpm.columns)
DF_components = pd.DataFrame(pca.components_.T,index=DF_log_tpm.index)
#Next plot the cumulative explained variance
# Set the explained variance threshold
var_cutoff = 0.99
fig,ax = plt.subplots(figsize=(5,3.5))
pca_var = np.cumsum(pca.explained_variance_ratio_)
ax.plot(pca_var)
dims = np.where(pca_var > var_cutoff)[0][0] + 1
ax.vlines(dims,0,1,linestyles='dotted')
ax.hlines(var_cutoff,0,len(DF_log_tpm.columns),linestyles='dotted')
ax.set_ylim(0,1)
ax.set_xlim(0,len(DF_log_tpm.columns))
ax.set_ylabel('Fraction of Explained Variance',fontsize=12)
ax.set_xlabel('Number of Dimensions',fontsize=12)
ax.set_title('Cumulative Explained Variance',fontsize=16)
print('Number of dimensions for 99% of variance:',dims)
plt.savefig("Cummulative_Explained_Variance.png", dpi=300, bbox_inches='tight')

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.cm as cm

# Compute correlation matrix
correlation_matrix = DF_log_tpm.T.corr()

#### Perform PCA
pca = PCA()
DF_weights = pd.DataFrame(pca.fit_transform(DF_log_tpm.T), index=DF_log_tpm.columns)

# Get unique conditions and assign unique colors
conditions = DF_metadata['condition'].unique()
num_conditions = len(conditions)
colors = cm.tab20(np.linspace(0, 1, num_conditions, endpoint=False))

# Plot PCA in 3D
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

for i, condition in enumerate(conditions):
    group = DF_metadata.loc[DF_metadata['condition'] == condition]
    idx = DF_log_tpm.loc[:, group.index.tolist()].columns.tolist()
    ax.scatter(DF_weights.loc[idx, 0], DF_weights.loc[idx, 1], DF_weights.loc[idx, 2], label=condition, alpha=0.8, color=colors[i])

ax.set_xlabel('Component 1: %.1f%%' % (pca.explained_variance_ratio_[0] * 100), fontsize=14)
ax.set_ylabel('Component 2: %.1f%%' % (pca.explained_variance_ratio_[1] * 100), fontsize=14)
ax.set_zlabel('Component 3: %.1f%%' % (pca.explained_variance_ratio_[2] * 100), fontsize=14)
ax.set_title('Principal Component Plot', fontsize=18)
plt.legend(bbox_to_anchor=(1.15, 1), fontsize=18, ncol=2)

plt.savefig("PCA_3D.png", dpi=300, bbox_inches='tight')


#ICA Summary
#Printout inputs for Treemap
tree_df = ica_data.imodulon_table[['imodulon_size', 'function']].reset_index()
tree_df['index'] = tree_df['index'].str.replace('Uncharacterized', 'Unc')
tree_df = tree_df.rename(columns={'index': 'name'})
tree_df.to_csv('tree_df.csv')

# Read the data
tree_df <- read.csv('tree_df.csv', row.names = 1)

#Let's go to R & create the treemap, since I like it better there
#library(treemap)
#library(RColorBrewer)
# Read the data
#tree_df <- read.csv('tree_df.csv', row.names = 1)

# Define color palette
#color_palette <- brewer.pal(8, "Pastel1")

# Create the treemap
#treemap(tree_df,
#        index=c('function.', 'name'),
#        vSize='imodulon_size',
#        palette=color_palette,
#        border.col=c("white","black"),
#        aspRatio=4.5,
#        fontsize.labels = c(12, 10),  # Adjust font size for groups and subgroups
#        align.labels = list(c("center", "center"), c("left", "top")),  # Align labels for groups and subgroups
#        fontcolor.labels = c("black", "black"),  # Font color for groups and subgroups
#        fontface.labels = c(2, 1)  # Font type for groups and subgroups
#)
# plot the variance explained of indiv iModulons
from pymodulon import io, core, plotting, util
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from matplotlib.ticker import MultipleLocator

rec_var = {}
for k in ica_data.imodulon_names:
    rec_var[k] = explained_variance(ica_data,imodulons=k)
df_rec_var = pd.Series(rec_var)
df_rec_var = df_rec_var.sort_values(ascending=False)
df_rec_var.head(100)

# Calculate explained variance for each iModulon
ic_var = [util.explained_variance(ica_data, imodulons=imodulon) for imodulon in ica_data.imodulon_names]

# Sort and assign names to the explained variance values
ic_var_name = dict(zip(ic_var, ica_data.imodulon_names))
ic_var_sum = np.insert(np.cumsum(sorted(ic_var, reverse=True)), 0, 0)

sns.reset_orig()
fig, ax = plt.subplots(figsize=(8, 6))
clr = sns.diverging_palette(250, 15, s=75, l=40, n=9, center="dark")
c1, c2 = clr[0], clr[-1]

ax2 = ax.twinx()
ax2.bar(np.arange(len(ic_var)), sorted(ic_var), color=c2)
ax2.set_ylim(0, 0.085)

ax.plot(range(len(ic_var_sum)), ic_var_sum, label="Independent Components", color=c2, zorder=0, alpha=0.9)
ax.set_ylim(0, 1.1)

# Get PCA & explained variance
pca = PCA().fit(ica_data.X.T)
pc_var = np.insert(np.cumsum(pca.explained_variance_ratio_), 0, 0)
pc_var = pc_var[:len(ic_var)]

# Plot PCs
ax.plot(range(len(ic_var)), pc_var, label="Principal Components", color=c1)

# Set tick parameters
ax.tick_params(which='both', axis='both', direction='in')
ax2.tick_params(which='both', axis='both', direction='in')
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))

ax.set_xticks(range(len(ic_var)))
ax.set_xticklabels(ica_data.imodulon_names, fontsize=8, rotation=45, ha='right')
ax.set_yticks([])
ax.set_yticklabels([])

ax.legend(loc=2, frameon=False)

ax.set_zorder(1)  # default zorder is 0 for ax1 and ax2
ax.patch.set_visible(False)  # prevents ax1 from hiding ax2

ax.set_xlabel('iModulon', fontsize=10)
ax.set_ylabel('Cumulative explained variance', fontsize=10)
ax2.set_ylabel('Individual explained variance', fontsize=10)

# Remove random words in the middle of the graph
ax2.text(11, 0.05, '', zorder=0, fontsize=10)
ax2.text(10.5, 0.035, '', fontsize=10)
ax2.text(11, 0.029, '', fontsize=10)

# Save the plot as a PNG file
plt.savefig("iModulon_Variance.png", dpi=300, bbox_inches='tight')

#Top 10 Variable iModulons
explained_var = explained_variance(ica_data)
plt.figure(figsize=(16, 8))  # Increase the width of the figure
plt.barh(range(10, 0, -1), df_rec_var.head(10), tick_label=df_rec_var.head(10).index)
plt.xlabel('Fraction of Explained Variance', fontsize=18, fontweight='bold')
new_labels[1] = 'Isoprenoid Biosyns. & AR'  # Specify the shortened name for the second from the top label
plt.yticks(range(10, 0, -1), labels=new_labels, fontsize=30, fontweight='bold')  # Set font size and weight for y-axis labels
plt.xticks(fontsize=14, fontweight='bold')  # Set font size and weight for x-axis labels
plt.title('Top 10 Explained Variance', fontsize=18, fontweight='bold')  # Set title font weight
plt.tight_layout()
plt.savefig("/Users/josspa/iMODULON/Manuscript/Figures/Figure_1/top10_Variance.png", dpi=300, bbox_inches='tight')
plt.show()
df_rec_var.head(10)
