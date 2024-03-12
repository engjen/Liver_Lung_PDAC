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

- **Figures made using R**
1.  Figure 2 A-C [here](https://github.com/engjen/Liver_Lung_PDAC/blob/main/Figure_2A-C.R)
2.  Figure 7D [here](https://github.com/engjen/Liver_Lung_PDAC/blob/main/Figure_7D.R)
3.  Supplemental figure 4A [here](https://github.com/engjen/Liver_Lung_PDAC/blob/main/Supplemental_Figure_4A.R)
4.  Supplemental figure 5 [here](https://github.com/engjen/Liver_Lung_PDAC/blob/main/Supplemental_Figure5A-D_scRNA.R)

- **Large files** including raw image data, single cell image features, and detailed Adaptive TCRseq and DNA sequence data can be found [here](https://www.synapse.org/#!Synapse:syn51068458/files/).

- **Additional analysis notebook** to load Adaptive TCR seq data and calculate TCR seq metrics found [here](https://github.com/engjen/Liver_Lung_PDAC/blob/main/20231222_TCR_seq_RNA_seq_data_processing.ipynb).
- **Immunarch code** to generate repertoire overlap found [here](https://github.com/engjen/Liver_Lung_PDAC/blob/main/TCR_analysis.R).

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

### R Packages and Versions

R version 4.1.2 was used with R packages DESeq2, GSVA, msigdbr, gplots, and ggplot.  GSEA was run in JAVA using the command line interface.

- immunarch (v0.9.0)
- ClusterProfiler (v4.6.2)
- immunedeconv (v2.1.0)
- enrichplot (v1.18.4) 
- Seurat (v4.3.0)
- pheatmap (v1.0.12) 
- MSigDB database (v7.5.1)
- FastQC (ver 0.11.8) 
- MultiQC (ver 1.7)
- trim-galore (ver 0.6.3)
- kallisto (ver 0.44.0)
- genome assembly GRCh38.p5 
- gencode (ver 24) 
- CNAtools
- GenVisR 

