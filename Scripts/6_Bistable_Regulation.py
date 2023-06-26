########FIGURE 3###############
#####Bistable-Regulation#######
#Plot simple bar graph
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('/Users/josspa/iMODULON/Manuscript/Figure3_Input.csv')
#Negative#LAC_011/12 Pectin / LAC_064/65 woThiamin / LAC_105-106 Ecoli_CDM,
#Positive #LAC_019/20 Glucose_ph67/ LAC_021-23 Glucose_pH74/ LAC_086-87 Human Milk

conditions = ['Pectin', 'woThiamine', 'E.coli_CDM', 'ph67', 'ph74', 'Human_Milk']  # Modify with your desired conditions
# Extract the necessary data for plotting
sample_names = data.columns[1:]
imodulon_names = data.iloc[[0, 1], 0]
activities = data.iloc[[0, 1], 1:].values
# Filter the activities based on the specified conditions
filtered_activities = activities[:, [sample_names.tolist().index(cond) for cond in conditions]]
# Plot the bar plot
fig, ax = plt.subplots(figsize=(10, 6))
bar_width = 0.35
index = range(len(conditions))
# Plot the first iModulon
ax.bar(index, filtered_activities[0], bar_width, label=imodulon_names[0])
# Plot the second iModulon next to the first one
ax.bar([i + bar_width for i in index], filtered_activities[1], bar_width, label=imodulon_names[1])
# Customize the plot
ax.set_xlabel('Conditions')
ax.set_ylabel('Activity')
ax.set_title('Activities of iModulons')
ax.set_xticks([i + bar_width/2 for i in index])
ax.set_xticklabels(conditions)
ax.legend()
plt.tight_layout()
plt.show()
plt.savefig("/Users/josspa/iMODULON/Manuscript/Figures/Figure_3/Barplot.png", dpi=300, bbox_inches='tight')


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#Plot Heatmap
data = pd.read_csv('/Users/josspa/iMODULON/Manuscript/Figure3_Input.csv')
# Extract the necessary data for plotting
sample_names = data.columns[1:]
imodulon_names = data.iloc[[0, 1], 0]
activities = data.iloc[[0, 1], 1:].values
# Convert activity values to binary (positive or negative)
activities_binary = np.where(activities >= 0, -1, 1)  # Swap positive and negative values
# Create the colormap with blue and red colors
cmap_custom = plt.cm.colors.ListedColormap(['red', 'blue'])  # Swap colors: red for negative, blue for positive
# Create the heatmap
fig, ax = plt.subplots(figsize=(10, 6))
im = ax.imshow(activities_binary, cmap=cmap_custom)
# Customize the plot
ax.set_xlabel('Samples', fontsize=10)  # Adjust the font size for x-axis labels
ax.set_ylabel('iModulons', fontsize=12)  # Adjust the font size for y-axis labels
ax.set_title('Binary Heatmap: Positive/Negative Activity')
# Create a custom colorbar with blue and red colors and reduce the size
bounds = [-1, 0, 1]
norm = plt.Normalize(vmin=-1, vmax=1)
cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap_custom), ax=ax, boundaries=bounds, ticks=[-1, 1], shrink=0.25)
cbar.ax.tick_params(labelsize=8)  # Adjust the font size for the colorbar labels
# Set the tick labels for the colorbar as "Negative" and "Positive"
cbar.set_ticklabels(['Negative', 'Positive'])
# Rotate x-axis labels vertically and adjust font size
ax.set_xticks(np.arange(len(sample_names)))
ax.set_xticklabels(sample_names, rotation='vertical', fontsize=8)
# Set y-axis tick labels as iModulon names
ax.set_yticks(np.arange(len(imodulon_names)))
ax.set_yticklabels(imodulon_names, fontsize=10)
# Adjust the layout to accommodate the rotated x-axis labels
plt.tight_layout()
# Save the plot to a file
plt.savefig("/Users/josspa/iMODULON/Manuscript/Figures/Figure_3/Heatmap.png", dpi=300, bbox_inches='tight')

from pymodulon.plotting import *
#Let's look at outliers here
plot_gene_weights(ica_data,'FabT', show_labels=True)
#Compare iModulons
compare_activities(ica_data, 'FMN-Box (Riboswitch)','FabT', groups=groups, show_labels=True)
#We have several genes which may have a bidirectional correlation such as
#LAC_050, LAC_052, LAC_053, LAC_093 LAC_094, LAC_095,
#LAC_099,  LAC_101, LAC_102,LAC_107,
#LAC_111, LAC_112, LAC_115, LAC_116, LAC_143, LAC_144
#Now that we know who is there, highlight them
groups = {'LAC_050': 'pABA+glutamate+GTP', 'LAC_050': 'pABA+glutamate+GTP', 'LAC_051': 'pABA+glutamate+GTP', 'LAC_052': 'pABA+glutamate+GTP', 'LAC_053': 'pABA+glutamate+GTP',
          'LAC_093': 'MRS_43C', 'LAC_094': 'MRS_43C', 'LAC_095': 'MRS_43C', 'LAC_096': 'MRS_43C',
          'LAC_099': '0.5g_Cysteine', 'LAC_100': '0.5g_Cysteine', 'LAC_101': 'S. thermophilus', 'LAC_102': 'S. thermophilus',
          'LAC_107': 'E.coli_MRS', 'LAC_108': 'E.coli_MRS', 'LAC_111': '1/4MRS+3/4BHI_37C', 'LAC_112': '1/4MRS+3/4BHI_37C',
          'LAC_115': '1/4MRS+3/4BHI_43C', 'LAC_116': '1/4MRS+3/4BHI_43C',
          'LAC_143': 'GTP', 'LAC_144': 'GTP'}

compare_activities(ica_data, 'FMN-Box (Riboswitch)', 'FabT', groups=groups, show_labels=False)
# Adjust the position of the legend and set the font size
legend = plt.legend(loc='upper right')
for text in legend.get_texts():
    text.set_fontsize(7)  # Adjust the font size of the legend
# Save the plot to a file
plt.savefig("/Users/josspa/iMODULON/Manuscript/Figures/Figure_3/FMN_vs_FabT.png", dpi=300, bbox_inches='tight')

#Let's plot a DIMA of the bistable conditions
plot_dima(ica_data, [groups=groups],show_labels=True, threshold=9, table=True)

plot_dima(ica_data,['LAC_050', 'LAC_051', 'LAC_052', 'LAC_053', 'LAC_143', 'LAC_143'], []
          show_labels=False, adjust_labels= False,
          threshold=-10, table=False)


#Negative#LAC_011/12 Pectin / LAC_064/65 woThiamin / LAC_105-106 Ecoli_CDM,
#Positive #LAC_019/20 Glucose_ph67/ LAC_021-23 Glucose_pH74/ LAC_086-87 Human Milk
plot_dima(ica_data,['LAC_011', 'LAC_012', 'LAC_064', 'LAC_065', 'LAC_105', 'LAC_106'],
          ['LAC_019','LAC_020', 'LAC_021', 'LAC_022', 'LAC_023',
           'LAC_086', 'LAC_087'],
          show_labels=False, adjust_labels= False,
          threshold=20, table=False)
plt.savefig("/Users/josspa/iMODULON/Manuscript/Figures/Figure_3/DIMA_Raw.png", dpi=300, bbox_inches='tight')
