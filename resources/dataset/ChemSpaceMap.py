from typing import List

#Dataframe
import pandas as pd
import numpy as np

#Plot
import plotly.graph_objects as go

#molecules
from rdkit import Chem, DataStructs
from rdkit.Chem.rdchem import Mol

#PCA
from sklearn.decomposition import PCA

#t-SNE
from sklearn.manifold import TSNE

#umap
import umap

#Fingerprint definition
def _compute_single_fp_descriptor(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as E:
        return None

    if mol:
        fp = Chem.AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        return np.array(fp)
    
    return None
    
#Computes ecfp descriptors  
def compute_fp_descriptors(smiles_list):
    
    idx_to_keep = list()
    descriptors = list()
    for i, smiles in enumerate(smiles_list):
        fp = _compute_single_fp_descriptor(smiles)
        if fp is not None:
            idx_to_keep.append(i)
            descriptors.append(fp)

    return np.vstack(descriptors), idx_to_keep


#Plotting functions 2D
def plot_2d(component1, component2, title=None):
    
    fig = go.Figure(data = go.Scatter(
        x = component1,
        y = component2,
        mode = 'markers',
        marker = dict(
                        size = 20,
                        color = dataframe.pIC50, #set color equal to a variable
                        colorbar = dict(
                        title = "pIC50"),
                        colorscale = 'Rainbow', # one of plotly colorscales
                        showscale = True,
                        line_width = 1
                    )
                                    )
                    )
    
    fig.update_layout(margin = dict( l = 100, r = 100, b = 100, t = 100),
                      width = 2000,
                      height = 1200, 
                      title = {
                            'text': title,
                            'y':0.9,
                            'x':0.5,
                            'xanchor': 'center',
                            'yanchor': 'top',
                            'font_size' : 100
                              },
                      font = dict(
                                family = "Courier New, monospace",
                                size = 30,
                                color = "RebeccaPurple"
                                  ) 
                        )     
    
    fig.update_traces(customdata = dataframe.molecule_chembl_id,
                  hovertemplate = "%{customdata}<extra></extra>"
                  )
    plotly.offline.plot(fig, filename= title +'_2D.html')
    fig.show()
    
#Plotting functions 3D
def plot_3d(component1, component2, component3, title = None):
    fig = go.Figure(data = [go.Scatter3d(
                                        x = component1,
                                        y = component2,
                                        z = component3,
                                        mode = 'markers',
                                        marker = dict(
                                                    size = 10,
                                                    color = dataframe.pIC50,
                                                    colorbar = dict(
                                                    title = "pIC50"),# set color to an array/list of desired values
                                                    colorscale = 'Rainbow',# choose a colorscale
                                                    opacity = 1,
                                                    line_width  =1
                                                    )
                                        )
                           ]
                   )

    fig.update_layout(margin = dict( l = 100, r = 100, b = 100, t = 100),
                      width = 2000,
                      height = 1200, 
                      title = {
                            'text': title,
                            'y':0.9,
                            'x':0.5,
                            'xanchor': 'center',
                            'yanchor': 'top',
                            'font_size' : 100
                              },
                      font = dict(
                            family = "Courier New, monospace",
                            size = 18,
                            color = "RebeccaPurple"
                      ) 
                     )         
    fig.update_traces(customdata = dataframe.molecule_chembl_id,
                  hovertemplate = "%{customdata}<extra></extra>"
                  )
    plotly.offline.plot(fig, filename= title +'_3D.html')

    fig.show()
    
#open and modify data
dataframe = pd.read_csv('dataset.csv', index_col = 0)
dataframe = dataframe[dataframe.bioactivity_class != 'intermediate']
y = pd.get_dummies(dataframe['bioactivity_class'])['active']
 
# Compute desrciptors and keep track of which failed to featurize
fp_descriptors, keep_idx = compute_fp_descriptors(dataframe["canonical_smiles"])

# Only keep those that sucessfully featurized
dataframe = dataframe.iloc[keep_idx]

#embed with PCA
pca = PCA(n_components=3)

principalComponents = pca.fit_transform(fp_descriptors)

principal = pd.DataFrame(
                        data = principalComponents,
                        columns = ['PCA 1', 'PCA 2','PCA 3']
                        )


#embed with umap
umap_data = umap.UMAP(
                    random_state=794,
                    n_components=3
                    )

embedding = umap_data.fit_transform(fp_descriptors)

#embed with t-SNE
pca_tsne = PCA(n_components=270)

pca_result = pca_tsne.fit_transform(fp_descriptors)

tsne = TSNE(
            init = 'pca',
            random_state = 794,
            n_components = 3,
            verbose = 2,
            perplexity = 10,
            n_iter = 3000
            ).fit_transform(pca_result)
            
            
#Plot
#PCA
plot_2d(
        principalComponents[:, 0],
        principalComponents[:, 1],
        'PCA'
        )
        
plot_3d(
        principalComponents[:, 0],
        principalComponents[:, 1],
        principalComponents[:, 2],
        title = 'PCA'
        )

#umap
plot_2d(
        umap_data.embedding_[:, 0],
        umap_data.embedding_[:, 1],
        title = 'UMAP'
        )

plot_3d(
        umap_data.embedding_[:, 0],
        umap_data.embedding_[:, 1],
        umap_data.embedding_[:, 2],
        title = 'UMAP'
       )

#t-SNE
plot_2d(
        tsne[:, 0],
        tsne[:, 1],
        title = 't-SNE'
        )

plot_3d(
    tsne[:, 0],
    tsne[:, 1],
    tsne[:, 2],
    title = 't-SNE'
        )        