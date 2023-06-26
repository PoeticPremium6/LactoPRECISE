#Metabolic Clustering
#Let's check out how similar these iModulons are
select_imods = ica_data.A.loc[['T-box(Met)', 'Lysine_Riboswitch', 'Stress/Homeostasis', 'recBCD', 'FMN-Box (Riboswitch)',
                               'RP-Operon-1', 'ArgR','FabT', 'Nar', 'Sugar Fermentation Regulator', 'Intestinal Carbohydrate Utilization',
                               'TPP Riboswitch', 'CcpA', 'YxeR', 'Fe-S Cluster Metabolism', 'RP-Operon-2', 'Zur', 'DeoR', 'Oxidative Stress Response',
                               'Isoprenoid Biosynthesis and Adaptive Response', 'PerR', 'PurR', 'PyrR', 'Stress Response and Proteostasis',
                               'uncharacterized-1', 'uncharacterized-2', 'uncharacterized-3', 'uncharacterized-4', 'uncharacterized-5',
                               'uncharacterized-6', 'uncharacterized-7', 'uncharacterized-8', 'uncharacterized-9', 'SG_1', 'SG_2']]
clusters = select_imods.T.corr()
g = sns.clustermap(clusters)
reorder = [clusters.index[i] for i in g.dendrogram_row.reordered_ind]
triu = np.triu(clusters.loc[reorder, reorder], k=1)
sns.set(font_scale=0.45)
# generate first clustermap
g = sns.clustermap(clusters)
plt.setp(g.ax_heatmap.get_xticklabels(), fontweight='bold')
plt.setp(g.ax_heatmap.get_yticklabels(), fontweight='bold')
plt.savefig("First_Clustering.png", dpi=300, bbox_inches='tight')
# clear the figure after saving it
plt.close()
# generate second clustermap
ax = sns.clustermap(clusters, cmap='vlag', annot=True)
plt.setp(ax.ax_heatmap.get_xticklabels(), fontweight='bold')
plt.setp(ax.ax_heatmap.get_yticklabels(), fontweight='bold')
plt.savefig("Second_Clustering.png", dpi=300, bbox_inches='tight')
# clear the figure after saving it
plt.close()

#Now let's check out which iModulons come together to form clusters that are related to biological function
# Call the function to generate the plots
cluster_activities(ica_data, show_best_clusters=True, n_best_clusters=5)
# Get all current figures
figures = [plt.figure(n) for n in plt.get_fignums()]
# Save each figure independently
for i, figure in enumerate(figures):
    figure.savefig(f"Best_Clustering_{i}.png", dpi=300, bbox_inches='tight')
