
# ChemSpaceMap

Use ChemSpaceMap to map molecular similarity in a chemical space in 2D and 3D with PCA, t-SNE and UMAP.

<p align="center">
<img align="center" style="width: 400px" src="https://github.com/carineribeirost/ChemSpaceMap/blob/main/resources/images/PCA_2D.png?"/>
</p>
<p align="center">
<img align="center" style="width: 500px" src="https://github.com/carineribeirost/ChemSpaceMap/blob/main/resources/images/PCA_3D.png?"/>
</p>

## Molecular Similarity in a Chemical Space

ChemSpaceMap creates a 2D and 3D visualization of a chemical dataset with **PCA, t-SNE** and **UMAP**, using a curated Tuberculosis
dataset from ChEMBL. 

By exploring the source code you will see

* How to install the relevant packages 
* Embed your SMILES as (extended-connectivity fingerprints) ECFPs
* Reduce high-dimensional vectors to 2 and 3 dimensions  
* How to represent chemical space in 2D and 3D

## Libraries used

* [Pandas](https://pandas.pydata.org/) - python package for easy and intuitive data manipulation and analysis

* [NumPy](https://numpy.org/) -  the fundamental package for array computing with Python

* [RDKit](https://www.rdkit.org/) - Open source toolkit for cheminformatics

* [Scikit-learn](https://scikit-learn.org/stable/) - Machine Learning in Python

* [UMAP](https://umap-learn.readthedocs.io/en/latest/) - Dimension reduction technique similar to t-SNE 

* [Plotly](https://plotly.com/) - provides graphing, analytics, and statistics tools, as well as scientific graphing libraries for Python, R and other languages.


Libraries were used in a [Conda3](https://docs.conda.io/en/latest/) environment using python 3.10.4

## Installation

Conda3: [Installation](https://docs.anaconda.com/anaconda/install/index.html)

pandas:
```
conda install -c anaconda pandas
```
numpy
```
conda install -c anaconda numpy
```
RDKit
```
conda install -c rdkit rdkit
```
scikit-learn
```
conda install -c anaconda scikit-learn
```

plotly
```
conda install -c anaconda plotly
```
umap
```
conda install -c conda-forge umap-learn
```
## How to run

* Download the code and unzip it on the desired directory

To run use the following command:

```
python ChemMapSpace.py
```
## Observations

ChemSpaceMap has been elaborated using 
as references the following articles and codes:

* [Dimensionality Reduction for Data Visualization: PCA vs TSNE vs UMAP vs LDA](https://towardsdatascience.com/dimensionality-reduction-for-data-visualization-pca-vs-tsne-vs-umap-be4aa7b1cb29)

* [Visualize chemical space as a grid](https://iwatobipen.wordpress.com/2019/08/27/visualize-chemical-space-as-a-grid-chemoinformatics-rdkit/)

* [Mapping Chemical Space with UMAP](https://blog.reverielabs.com/mapping-chemical-space-with-umap/)

## Authorship
* Author: **Carine Ribeiro** ([carineribeirost](https://github.com/carineribeirost))

Social preview original photo by **Carine Ribeiro** ([carineribeirost](https://github.com/carineribeirost))
