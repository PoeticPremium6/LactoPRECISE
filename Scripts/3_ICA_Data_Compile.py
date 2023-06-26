from pymodulon.core import *
from pymodulon.plotting import *
from pymodulon.enrichment import *
from pymodulon.compare import *
from pymodulon.io import *
from os import path

#Load M,A,X matrix, look for any duplications
M = ('M.csv')
DF_M = pd.read_csv(M,index_col=0,sep=',')
print(DF_M)
DF_M.columns.duplicated().any()
A = ('A_copy.csv')
DF_A = pd.read_csv(A,index_col=0,sep=',')
print(DF_A)
DF_A = DF_A.reset_index(drop=True)
DF_A.columns.duplicated().any()


X = 'log_tpm_final.csv'
DF_X = pd.read_csv(X, index_col=0, sep=',')
DF_X.head()

#Let's add auxilirary info
gene_table = ('gene_table.csv')
DF_gene = pd.read_csv(gene_table,index_col=0,sep=';')
print(DF_gene)
df_metadata = ('Metadata_Final.tsv')
df_MetaData= pd.read_csv(df_metadata,index_col=0,sep='\t')
df_MetaData[['project','condition']].head()
print(df_MetaData.project.notnull().all())
print(df_MetaData.condition.notnull().all())

df_trn = (TRN.csv')
DF_TRN = pd.read_csv(df_trn,sep=',')
DF_TRN.head()
trn = DF_TRN
#Make IcaData Object
ica_data = IcaData(M = DF_M,
                   A = DF_A,
                   X = DF_X,
                   gene_table = DF_gene,
                   sample_table = df_MetaData,
                   trn = DF_TRN,
                   optimize_cutoff=True)
#Compute TRN enrichment
ica_data.compute_trn_enrichment()
# First search for regulator enrichments with 2 regulators
ica_data.compute_trn_enrichment(max_regs=2,save=True)

# Next, search for regulator enrichments with just one regulator. This will supersede the 2 regulator enrichments.
ica_data.compute_trn_enrichment(max_regs=1,save=True)
#Make a list of regulatory iModulons
regulatory_imodulons = ica_data.imodulon_table[ica_data.imodulon_table.regulator.notnull()]
print(len(ica_data.imodulon_table),'Total iModulons')
print(len(regulatory_imodulons),'Regulatory iModulons')
regulatory_imodulons
#If two iModulons have the same regulatory, they will be named 'Reg-1' &'Reg-2'
ica_data.rename_imodulons(regulatory_imodulons.regulator.to_dict())
ica_data.imodulon_table.head()
regulatory_imodulons = ica_data.imodulon_table[ica_data.imodulon_table.regulator.notnull()]

#Look for single gene iModulons
sg_imods = ica_data.find_single_gene_imodulons(save=True)
len(sg_imods)
for i,mod in enumerate(sg_imods):
    ica_data.rename_imodulons({mod:'SG_'+str(i+1)})
ica_data.imodulon_table[ica_data.imodulon_table.single_gene == True]

#Save iModulon Object
# Add iModulon sizes and explained variance
for im in ica_data.imodulon_names:
    ica_data.imodulon_table.loc[im,'imodulon_size'] = len(ica_data.view_imodulon(im))
    ica_data.imodulon_table.loc[im,'explained_variance'] = explained_variance(ica_data,imodulons=im)

#save what we have
save_to_json(ica_data, path.join('..','data','interim','lreu_ICA_06062023.json.gz'))
ica_data.imodulon_table.to_csv(path.join('..','data','interim','imodulon_table_raw_01062023.csv'))

#Manually Curate iModulon
ica_data = load_json_model('lreu_ICA_06062023.json.gz')
DF_enrichments = pd.read_csv(path.join('Functional_enrichments.csv'),index_col=0)

#Add iModulon Categories
for i,row in ica_data.imodulon_table.iterrows():
    if pd.notnull(row.regulator):
        ica_data.imodulon_table.loc[i, 'category'] = 'regulatory'
    elif pd.notnull(row.single_gene):
        ica_data.imodulon_table.loc[i, 'category'] = 'single_gene'
    else:
        ica_data.imodulon_table.loc[i, 'category'] = 'uncharacterized'

ica_data.imodulon_table.head()
