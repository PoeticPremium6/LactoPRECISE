#######FIGURE 6 ###########
#########T3PKS_BGCs#####

#We just ran our genome against antiSMASH version 7
#RIPP-Like region was found but none of the genes that were predicted
#are in our iModulon model
#However T3PKS has several genes that are found in our model
#The following genes were identified by AntiSMASH for the T3PKS BGC:
#LMB90_RS06635 (Core Biosynthetic Gene)	hydroxymethylglutaryl-CoA synthase
#LMB90_RS06655	(Additional Biosynthetic Genes) 1-acyl-sn-glycerol-3-phosphate acyltransferase
#LMB90_RS06670	(Additional Biosynthetic Genes) D-2-hydroxyacid dehydrogenase
#Found in 'Sugar Fermentation Regulator'
#LMB90_RS06715 (rseP)	(Additional Biosynthetic Genes) RIP metalloprotease RseP
ica_data.num2name('LMB90_RS06670')
ica_data.imodulons_with('LMB90_RS06670')
#We can look at all amino acids
groups = {'LAC_149':'Uracil', 'LAC_150':'Uracil', 'LAC_151':'wo_Uracil', 'LAC_152':'wo_Uracil','LAC_155':'wo_Glutamate','LAC_156':'wo_Glutamate','LAC_157':'wo_Glutamate', 'LAC_158':'wo_Glutamate',
          'LAC_003': 'NH4Cl', 'LAC_004': 'NH4Cl', 'LAC_084': 'woAdenine', 'LAC_085': 'woAdenine',
          'LAC_133':'woTryptophan','LAC_134':'woTryptophan','LAC_135':'woHistidine','LAC_136':'woHistidine',
          'LAC_028':'woGlutamine','LAC_029':'woGlutamine', 'LAC_030':'woGlutamine', 'LAC_031':'woGlutamine', 'LAC_032':'1g_Cysteine', 'LAC_033':'1g_Cysteine', 'LAC_034':'1g_Cysteine',
          'LAC_024':'Glutamine', 'LAC_025':'Glutamine', 'LAC_026':'Glutamine', 'LAC_027':'Glutamine',
          'LAC_046':'Glutamate', 'LAC_047':'Glutamate', 'LAC_048':'Glutamate','LAC_049':'Glutamate',
          'LAC_099': '0.5g_Cysteine', 'LAC_100': '0.5g_Cysteine',
          'LAC_042':'Histidine', 'LAC_043':'Histidine', 'LAC_044':'Histidine', 'LAC_045':'Histidine',
          'LAC_131':'Tyrosine', 'LAC_132':'Tyrosine', 'LAC_129':'woCysteine', 'LAC_130':'woCysteine'}
#Or a subset
groups = {'LAC_003':'NH4Cl', 'LAC_004':'NH4Cl',
          'LAC_149':'Uracil', 'LAC_150':'Uracil','LAC_151':'wo_Uracil', 'LAC_152':'wo_Uracil',
          'LAC_155':'wo_Glutamate','LAC_156':'wo_Glutamate','LAC_157':'wo_Glutamate', 'LAC_158':'wo_Glutamate',
          'LAC_046':'Glutamate', 'LAC_047':'Glutamate', 'LAC_048':'Glutamate','LAC_049':'Glutamate',}
#Quick looks at conditions that activate secondary metabolites
plot_activities(ica_data,'Isoprenoid Biosynthesis and Adaptive Response',highlight=['wo_Glutamate', 'wo_Uracil', 'NH4Cl']);

#Plot the expression of our BGC genes against secondary metabolite iModulon

scatterplot(ica_data.X.T['LMB90_RS06655'],ica_data.A.T['Isoprenoid Biosynthesis and Adaptive Response'],
            fit_metric= 'pearson',fit_line='True',
            xlabel='', ylabel='',
            show_labels=False, groups=groups, colors=['black', 'rosybrown','lightcoral', 'firebrick', 'red',
                                                     'darkorange', 'khaki', 'yellow', 'olivedrab', 'chartreuse',
                                                     'aquamarine', 'aqua', 'deepskyblue', 'dodgerblue',
                                                     'pink', 'royalblue', 'lavender', 'lime', 'blue',
                                                     'rebeccapurple', 'darkorchid', 'plum', 'magenta', 'deeppink', 'crimson'])
plt.xlabel('LMB90_RS06655 expression', fontsize=14)
plt.ylabel('Isoprenoid Biosynthesis and Adaptive Response', fontsize=11.5)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig("BGC_LMB90_RS06655.png", dpi=300, bbox_inches='tight')

#It appears that we have an interesting relationship with nitrogen sources and secondary metabolite production
#Let's further compare activity between isoprenoid iModulon & amino acid control iModulons
compare_activities(ica_data,'T-box(Met)','Isoprenoid Biosynthesis and Adaptive Response',groups=groups, show_labels=False,colors=['black', 'rosybrown','lightcoral', 'firebrick', 'red',
                                                                                                                       'darkorange', 'khaki', 'yellow', 'olivedrab', 'chartreuse',
                                                                                                                       'aquamarine', 'aqua', 'deepskyblue', 'dodgerblue',
                                                                                                                       'pink', 'royalblue', 'lavender', 'lime', 'blue',
                                                                                                                       'rebeccapurple', 'darkorchid', 'plum', 'magenta', 'deeppink', 'crimson'])
plt.xlabel('T-box(Methionine)', fontsize=14)
plt.ylabel('Isoprenoid Biosynthesis and Adaptive Response', fontsize=11.5)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', prop={'size': 12})
plt.savefig("Methionine.png", dpi=300, bbox_inches='tight')
#Pretty interesting and strong correlations:
#Methione regulation has a Pearson's R of 0.49
#Lysine regulation has a Pearson's R of 0.66
#ArgR regulation has a Pearson's R of -0.51
#Uncharacterized-9 has a Pearson's R of 0.53
#This uncharacterized iModulon has a single gene for Ammonia transport
