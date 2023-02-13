# Liver_Lung_PDAC
Analysis of immunofluorescence staining to detect replication stress foci

Analysis notebook [here](https://github.com/engjen/Liver_Lung_PDAC/blob/main/20221025_PDAC_pipeline_Link.ipynb) used to:

1. Detect nuclear foci of replication stress markers pRPA, gH2AX and RAD51.
2. Calculate mean foci per epithelial cell in primary PDAC tumors and link to organotropism (i.e. liver or lung metastasis)
3. Calculate fraction of multiplex IHC cell types per tissue and link to organotropism.

Image data and large files to run the first and third part of the code can be found [here](https://www.synapse.org/#!Synapse:syn51068458/files/).

If utilizing images, data or code, please cite our work: [Ongoing Replication Stress Response and New Clonal T Cell Development Discriminate Between Liver and Lung Recurrence Sites and Patient Outcomes in Pancreatic Ductal Adenocarcinoma](https://www.biorxiv.org/content/10.1101/2022.05.04.490552v1)


## Analysis environment

To run the main analysis, installing python3/miniconda, and enter the following in the terminal to set up an `analysis` environment. 

`conda create -n analysis`

`conda activate analysis`

`conda install seaborn pytables pandas ipykernel`

`conda install -c conda-forge jupyterlab matplotlib scikit-image tifffile statsmodels`

`pip install statannotations`

Finally, clone my repo for processing, visualization and analysis of multiplex imaging data

`git clone https://gitlab.com/engje/mplex_image.git`
