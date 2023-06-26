######FIGURE 4 #########
######Context-specific_ArgR####
plot_gene_weights(ica_data,'ArgR')
plot_activities(ica_data,'ArgR')

ArgR = ica_data.view_imodulon('ArgR')
ArgR.to_csv(path.join('..','data','interim','ArgR.csv'))

#Negative Gene weights
LMB90_RS03530  (-0.210557) (argH)
LMB90_RS03525 (-0.215007)
LMB90_RS03500 (-0.207930)
LMB90_RS03505 (-0.194097)
LMB90_RS07510  (-0.128690) (carA)
LMB90_RS07505 (-0.126191)
LMB90_RS10070 (-0.097815)
LMB90_RS05535 (-0.094347)
LMB90_RS10075 (-0.090404)
#Positive Gene Weights
LMB90_RS04610     (0.248290) (argF)
LMB90_RS04615     (0.248019) (arcC)
LMB90_RS04690     (0.185828) (arcA)
LMB90_RS08635     (0.153912)
LMB90_RS09710     (0.126116) (gtfA)
LMB90_RS06155     (0.124790)
LMB90_RS09705     (0.117866)
LMB90_RS04700     (0.115997)
LMB90_RS04695     (0.111116)
LMB90_RS06160     (0.111666)

ica_data.num2name('LMB90_RS06160')

groups = {'LAC_009': 'Glycerol', 'LAC_010': 'Glycerol', 'LAC_086': 'HumanMilk', 'LAC_087': 'HumanMilk', 'LAC_089': 'Fruit_Juice', 'LAC_090': 'Fruit_Juice',
          'LAC_103': 'L.lactis_cremoris_JS102', 'LAC_104': 'L.lactis_cremoris_JS102', 'LAC_101': 'S. thermophilus', 'LAC_102': 'S. thermophilus',
          'LAC_072': 'Piperacillin-tazobactam', 'LAC_073': 'Piperacillin-tazobactam', 'LAC_074': 'Ciprofloxacin', 'LAC_075': 'Ciprofloxacin',
          'LAC_145': 'Tetracycline', 'LAC_146': 'Tetracycline',
          'LAC_107': 'E. coli', 'LAC_108': 'E. coli', 'LAC_111': 'L.reuteri_37C', 'LAC_112': 'L.reuteri_37C',
          'LAC_115': 'L.reuteri_43C', 'LAC_116': 'L.reuteri_43C',
          'LAC_156': 'wo_Glutamate', 'LAC_139': 'woFe', 'LAC_140': 'woFe',
          'LAC_143': 'GTP', 'LAC_144': 'GTP', 'LAC_093': 'MRS_43C', 'LAC_094': 'MRS_43C', 'LAC_095': 'MRS_30C', 'LAC_096': 'MRS_30C',
          'LAC_099': '0.5g_Cysteine', 'LAC_100': '0.5g_Cysteine', 'LAC_155': 'wo_Glutamate', 'LAC_157': 'wo_Glutamate', 'LAC_158': 'wo_Glutamate',
          'LAC_141':'pABA', 'LAC_142': 'pABA'}

ica_data.num2name('LMB90_RS04700')
#Color options : https://matplotlib.org/stable/gallery/color/named_colors.html
import matplotlib.pyplot as plt

# The scatterplot function call
scatterplot(ica_data.X.T['LMB90_RS09710'],ica_data.A.T['ArgR'],
            fit_metric= 'pearson',fit_line='True',
            xlabel='gtfA expression', ylabel='ArgR iModulon Activity', show_labels=False,
            groups=groups, colors=['black', 'lightcoral', 'firebrick', 'red',
                                   'darkorange', 'khaki', 'yellow', 'olivedrab', 'chartreuse',
                                   'aquamarine', 'aqua', 'deepskyblue', 'dodgerblue',
                                   'cornflowerblue', 'royalblue', 'lavender', 'lime', 'blue',
                                   'rebeccapurple', 'darkorchid', 'plum', 'magenta', 'deeppink', 'crimson'])

# Get the current figure and axes
fig = plt.gcf()
ax = plt.gca()
# Increase x-axis label font size
ax.xaxis.label.set_size(20)
# Adjust the legend size
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
# Save the figure
plt.savefig("gtfA.png", dpi=300, bbox_inches='tight')

