# Liver_Lung_PDAC
Analysis of patitnet survival, DNA alterations, gene expression, TCR sequencing and immunofluorescence staining versus organotropism.

## Code

- **Main analysis notebook** to produce figures in paper is found [here](https://github.com/engjen/Liver_Lung_PDAC/blob/main/20221025_PDAC_pipeline_Link.ipynb).
Contents:
1. Detect nuclear foci of replication stress markers pRPA, gH2AX and RAD51.
2. Calculate mean foci per epithelial cell in primary PDAC tumors and link to organotropism (i.e. liver or lung metastasis)
3. Calculate fraction of multiplex IHC cell types per tissue and link to organotropism.
4. Patient metadata. Primary vs met. DDR vs TMB.
5. CPH modeling CPH forest plots
6. gene expression analysis
7. TCR analysis TCR survival
8. GSVA violins GSEA bar plots

- Image data and large files to run the first and third part of the code can be found [here](https://www.synapse.org/#!Synapse:syn51068458/files/).

- **Additional analysis notebook** to load Adaptive TCR seq data and calculate TCR seq metrics found [here](https://github.com/engjen/Liver_Lung_PDAC/blob/main/20231222_TCR_seq_RNA_seq_data_processing.ipynb). 

## Citation

If utilizing images, data or code, please cite our work: [Ongoing Replication Stress Response and New Clonal T Cell Development Discriminate Between Liver and Lung Recurrence Sites and Patient Outcomes in Pancreatic Ductal Adenocarcinoma](https://www.biorxiv.org/content/10.1101/2022.05.04.490552v1)


## Analysis environment

### Python
To run the analysis notebooks, install [python3/miniconda](https://docs.conda.io/en/latest/miniconda.html) (installers for Windows, macOS and Linux), and enter the following in the terminal to set up an `analysis` environment. 

`conda create -n analysis`

`conda activate analysis`

`conda install seaborn pytables pandas ipykernel`

`conda install -c conda-forge jupyterlab matplotlib scikit-image tifffile statsmodels`

`pip install statannotations`

Finally, clone my repo for processing, visualization and analysis of multiplex imaging data

`git clone https://gitlab.com/engje/mplex_image.git`

### R

...

