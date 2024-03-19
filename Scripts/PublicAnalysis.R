### Public data and scripts downloads ### ======================================
#
### PurIST subtyping tool ###
# Download the zip file from the github repo here:
# https://github.com/naimurashid/PurIST
# (NOTE: Disclaimer: PurIST is patent pending and access/use is for not-for-profit research only.)
# Extract the archive to ScriptsPath. This should create a "PurIST-master" folder.
# Place this folder in the "Scripts" folder.

###  ICGC data ### 
# https://static-content.springer.com/esm/art%3A10.1038%2Fnature16965/MediaObjects/41586_2016_BFnature16965_MOESM271_ESM.xlsx
# Save the "PDAC normalised exp" tab as a tab-separated-file and place it in the
# "PublicDataUsed" folder. 
# 
# The survival data for the above was received as a personal communication and 
# stored in the file "ICGC_APGI_Survival_Metadata.tsv".  This is acknowledged in
# the publication. It should be in the "data" folder.
# 
###  TCGA data ### 
# Download full dataset from bioportal:
# https://cbioportal-datahub.s3.amazonaws.com/paad_tcga_pan_can_atlas_2018.tar.gz
# Extract and use: data_RNA_Seq_v2_expression_median.txt
# Go to: https://www.cbioportal.org/study/summary?id=paad_tcga_pan_can_atlas_2018
# and click "Download clinical data for the selected cases".  This should download
# "paad_tcga_pan_can_atlas_2018_clinical_data.tsv" for 184 samples.
# Place both of these files in the "PublicDataUsed" folder.
#
#= Setup paths and load scripts and data.  =====================================

library(XLConnect)
library(ROCit)
library(survival)
library(ggplot2)
library(ggfortify)
library("Cairo", character.only=TRUE)
library(ComplexHeatmap)
library("DESeq2")
library("edgeR")
library(GSEABase)
library(GSVA)
library(msigdbr)
library("fdrtool")
library(biomaRt)

BaseDir = paste(getwd(), "/", sep="")
ScriptsPath = paste(BaseDir, "Scripts/", sep="")
SourcePath = BaseDir
AnalysisOutputPath = paste(BaseDir, "Analysis/", sep="")

source(paste(ScriptsPath, "SupportFunctions.R", sep=""))
BCC_RNA_Data = LoadRnaSeqData(SourcePath, AnalysisOutputPath)

#= Do PurIST subtyping of PDAC samples.  ========================================

# Calculate PurIST subtype for all RNA-Seq samples.
PurIST_Result = ApplyPurIST(ScriptsPath = ScriptsPath, 
                            Data = BCC_RNA_Data$CollapsedTPM)

# Add to metadata.
BCC_RNA_Data$Metadata$PurIST_Subtype = PurIST_Result$PurIST_Subtype
BCC_RNA_Data$Metadata$PurIST_Detailed_Subtype = PurIST_Result$PurIST_Detailed_Subtype
BCC_RNA_Data$Metadata$PurIST_Score = PurIST_Result$PurIST_Score

#= Generate pORG/pSUB gene sets and GSVA scores. ===============================

BCC_RNA_Data = CalculateBCCpORGpSUBandGSVAscores(AnalysisOutputPath, BCC_RNA_Data)

Idx1 = match(x=BCC_RNA_Data$Metadata$Public_Specimen_ID, rownames(BCC_RNA_Data$GSVA_pORGpSUB_Primaries))
Idx2 = match(x=BCC_RNA_Data$Metadata$Public_Specimen_ID, rownames(BCC_RNA_Data$GSVA_pORGpSUB_Mets))
Idx3 = match(x=BCC_RNA_Data$Metadata$Public_Specimen_ID, rownames(BCC_RNA_Data$GSVA_pORGpSUB_All))
write.table(cbind(BCC_RNA_Data$Metadata,
                  pORG_Up_55_Primaries = BCC_RNA_Data$GSVA_pORGpSUB_Primaries[Idx1, "pORG_Up_55"],
                  pSUB_Up_51_Primaries = BCC_RNA_Data$GSVA_pORGpSUB_Primaries[Idx1, "pSUB_Up_51"],
                  pORG_Up_55_Mets = BCC_RNA_Data$GSVA_pORGpSUB_Mets[Idx2, "pORG_Up_55"],
                  pSUB_Up_51_Mets = BCC_RNA_Data$GSVA_pORGpSUB_Mets[Idx2, "pSUB_Up_51"],
                  pORG_Up_55_All = BCC_RNA_Data$GSVA_pORGpSUB_All[Idx3, "pORG_Up_55"],
                  pSUB_Up_51_All = BCC_RNA_Data$GSVA_pORGpSUB_All[Idx3, "pSUB_Up_51"]),
            paste(AnalysisOutputPath, "RNA_Specimen_Metadata_Plus.tsv", sep=""), 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
# Write tables used for Supplemental Dataset 5 GSVA hallmarks scores for 
# primaries, mets, and all.
write.table(cbind(Public_Specimen_ID = rownames(BCC_RNA_Data$GSVA_Hallmarks_Primaries),
                  BCC_RNA_Data$GSVA_Hallmarks_Primaries),
            paste(AnalysisOutputPath, "GSVA_Hallmarks_Primaries.tsv", sep=""), 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(cbind(Public_Specimen_ID = rownames(BCC_RNA_Data$GSVA_Hallmarks_Mets),
                  BCC_RNA_Data$GSVA_Hallmarks_Mets),
            paste(AnalysisOutputPath, "GSVA_Hallmarks_Mets.tsv", sep=""), 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(cbind(Public_Specimen_ID = rownames(BCC_RNA_Data$GSVA_Hallmarks_All),
                  BCC_RNA_Data$GSVA_Hallmarks_All),
            paste(AnalysisOutputPath, "GSVA_Hallmarks_All.tsv", sep=""), 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#= The following produces all the OncoPrint plots used. ========================
{
#= Load BCC xT panel alterations data and metadata =============================

BCC_DNA_Data = LoadDnaPanelData(SourcePath)
# Define a sortable ImpactLevel:
#unique(BCC_DNA_Data$Variants$Significance)
BCC_DNA_Data$Variants$ImpactLevel = NA
BCC_DNA_Data$Variants$ImpactLevel[BCC_DNA_Data$Variants$Significance == "Biologically relevant"] = 1
BCC_DNA_Data$Variants$ImpactLevel[BCC_DNA_Data$Variants$Significance == "VUS"] = 3

# Note that there are two versions of the panel (the second has just one
# additional gene).
unique(BCC_DNA_Data$Variants$Panel_Name)

#= HR_DDR_Genes ================================================================

# Add HR_DDR status to DNA data (using all called variants).
BCC_DNA_Data$Metadata$HR_DDR_Altered_All = FALSE
idx = which(BCC_DNA_Data$Variants$Gene %in% BCC_DNA_Data$HR_DDR_Genes)
length(idx)
BCC_DNA_Data$Metadata$HR_DDR_Altered_All[BCC_DNA_Data$Metadata$Public_Specimen_ID %in% 
                                           unique(BCC_DNA_Data$Variants$Public_Specimen_ID[idx])] = TRUE
length(which(BCC_DNA_Data$Metadata$HR_DDR_Altered_All == FALSE))
length(which(BCC_DNA_Data$Metadata$HR_DDR_Altered_All == TRUE))

# Add HR_DDR status to DNA data (using just variants called Biologically Relevant).
BCC_DNA_Data$Metadata$HR_DDR_Altered_BioRel = FALSE
idx = which(BCC_DNA_Data$Variants$Significance == "Biologically relevant" &
            BCC_DNA_Data$Variants$Gene %in% BCC_DNA_Data$HR_DDR_Genes)
length(idx)
BCC_DNA_Data$Metadata$HR_DDR_Altered_BioRel[BCC_DNA_Data$Metadata$Public_Specimen_ID %in% 
                                           unique(BCC_DNA_Data$Variants$Public_Specimen_ID[idx])] = TRUE
length(which(BCC_DNA_Data$Metadata$HR_DDR_Altered_BioRel == FALSE))
length(which(BCC_DNA_Data$Metadata$HR_DDR_Altered_BioRel == TRUE))

# Number of patients with xT DNA panels:
length(unique(BCC_DNA_Data$Metadata$Public_Patient_ID))
# Number of patients with reported genes on the xT DNA panels:
length(unique(BCC_DNA_Data$Variants$Public_Patient_ID))
# (There are 6 panels with no reported genes.)

# Add to RNA metadata, but make sure to leave values as NA for patients without
# DNA panel data.  In the RNA metadata, treat this as a patient level attribute
# (if any sample for a patient has any HR_DDR alterations, then that patient and
# any samples for that patient will be analyzed as HR_DDR_Altered_All == TRUE.)
BCC_RNA_Data$Metadata$HR_DDR_Altered_All = NA
idx = which(BCC_RNA_Data$Metadata$Public_Patient_ID %in%
              BCC_DNA_Data$Metadata$Public_Patient_ID[BCC_DNA_Data$Metadata$HR_DDR_Altered_All == TRUE])
BCC_RNA_Data$Metadata$HR_DDR_Altered_All[idx] = TRUE 
idx = which(BCC_RNA_Data$Metadata$Public_Patient_ID %in%
              BCC_DNA_Data$Metadata$Public_Patient_ID[BCC_DNA_Data$Metadata$HR_DDR_Altered_All == FALSE])
BCC_RNA_Data$Metadata$HR_DDR_Altered_All[idx] = FALSE 

length(which(BCC_RNA_Data$Metadata$HR_DDR_Altered_All == TRUE))
length(which(BCC_RNA_Data$Metadata$HR_DDR_Altered_All == FALSE))
length(which(is.na(BCC_RNA_Data$Metadata$HR_DDR_Altered_All)))

BCC_RNA_Data$Metadata$HR_DDR_Altered_BioRel = NA
idx = which(BCC_RNA_Data$Metadata$Public_Patient_ID %in%
              BCC_DNA_Data$Metadata$Public_Patient_ID[BCC_DNA_Data$Metadata$HR_DDR_Altered_BioRel == TRUE])
BCC_RNA_Data$Metadata$HR_DDR_Altered_BioRel[idx] = TRUE 
idx = which(BCC_RNA_Data$Metadata$Public_Patient_ID %in%
              BCC_DNA_Data$Metadata$Public_Patient_ID[BCC_DNA_Data$Metadata$HR_DDR_Altered_BioRel == FALSE])
BCC_RNA_Data$Metadata$HR_DDR_Altered_BioRel[idx] = FALSE 

length(which(BCC_RNA_Data$Metadata$HR_DDR_Altered_BioRel == TRUE))
length(which(BCC_RNA_Data$Metadata$HR_DDR_Altered_BioRel == FALSE))
length(which(is.na(BCC_RNA_Data$Metadata$HR_DDR_Altered_BioRel)))

#= Collapse mutation type down to a more manageable number =====================

mutation_priority = BCC_DNA_Data$Variants$Alteration_Type
unique(mutation_priority)
CheckSum = 0
FoundIdx = NULL
for (SearchTerm in c("missense", "stop", "inframe", 
                     "frameshift", "splice",
                     "Copy number loss", "Copy number gain"))
{
  idx = grep(toupper(SearchTerm), toupper(mutation_priority))
  print(unique(mutation_priority[idx]))
  if(length(idx) > 0)
  {
    mutation_priority[idx] = SearchTerm
  }
  CheckSum = CheckSum + length(idx)
  FoundIdx = c(FoundIdx, idx)
}
length(mutation_priority)
CheckSum
unique(mutation_priority)
FoundIdx = unique(FoundIdx)
unique(mutation_priority[-FoundIdx])
mutation_priority[-FoundIdx] = "other_variant"

mutation_priority[mutation_priority == "missense"] = "Missense"
mutation_priority[mutation_priority == "frameshift"] = "Frameshift"
mutation_priority[mutation_priority == "inframe"] = "Inframe_indels"
mutation_priority[mutation_priority == "Copy number gain"] = "Copy_number_gain"
mutation_priority[mutation_priority == "Copy number loss"] = "Copy_number_loss"
mutation_priority[mutation_priority == "splice"] = "Splice_region_variants"
mutation_priority[mutation_priority == "stop"] = "Stop_gains"
mutation_priority[mutation_priority == "protein_protein_contact"] = "Protein_protein_contact"
mutation_priority[mutation_priority == ""] = "UnDefined"
unique(mutation_priority)

idx = which(mutation_priority == "Protein_protein_contact")
BCC_DNA_Data$Variants$Alteration_Type[idx]
BCC_DNA_Data$Variants$Nucleic_Acid[idx]

BCC_DNA_Data$Variants$mutation_priority = mutation_priority
rm(mutation_priority)

#= Define cohorts for DNA panel data ===========================================

# Link DNA metadata to RNA metadata.
# Just one DNA specimen without an RNA specimen due to finding error in that 
# RNA specimen.  That specimen is for patient ST-00021102. We are using the RNA
# data to define some of the cohorts below, so the DNA specimen for ST-00021102
# will have to be left out.
idx = match(x=BCC_DNA_Data$Metadata$Public_Specimen_ID,
            table=BCC_RNA_Data$Metadata$Public_Specimen_ID)
BCC_DNA_Data$RnaMetadataIdx = idx

BCC_DNA_Data$LiverLungCohort = rep("Unknown", times=dim(BCC_DNA_Data$Metadata)[[1]])

BCC_DNA_Data$LiverLungCohort[(!is.na(BCC_DNA_Data$Metadata$Liver_Met_Present) & 
                             BCC_DNA_Data$Metadata$Liver_Met_Present == 1)] = "Liver"
BCC_DNA_Data$LiverLungCohort[(!is.na(BCC_DNA_Data$Metadata$Lung_Met_Present) & 
                             BCC_DNA_Data$Metadata$Lung_Met_Present == 1) &
                             !(!is.na(BCC_DNA_Data$Metadata$Liver_Met_Present) & 
                               BCC_DNA_Data$Metadata$Liver_Met_Present == 1)] = "Lung"

#BCC_DNA_Data$LiverLungCohort
#BCC_DNA_Data$Metadata$Tumor_Type
#BCC_DNA_Data$Metadata$HR_DDR_Altered_BioRel

BCC_DNA_Data$Metadata$PurIST_Subtype = BCC_RNA_Data$Metadata$PurIST_Subtype[BCC_DNA_Data$RnaMetadataIdx]
BCC_DNA_Data$Metadata$PurIST_Detailed_Subtype = BCC_RNA_Data$Metadata$PurIST_Detailed_Subtype[BCC_DNA_Data$RnaMetadataIdx]
BCC_DNA_Data$Metadata$PurIST_Score = BCC_RNA_Data$Metadata$PurIST_Score[BCC_DNA_Data$RnaMetadataIdx]

Number = floor(length(which(!is.na(BCC_DNA_Data$Metadata$PurIST_Score) &
                              BCC_DNA_Data$Metadata$Tumor_Type == "Primary")) / 4)
BCC_DNA_Data$TopBottomQuarterPrimariesPurIST = rep(NA, times=dim(BCC_DNA_Data$Metadata)[[1]])
idx = order(BCC_DNA_Data$Metadata$PurIST_Score, decreasing = TRUE)
idx = idx[BCC_DNA_Data$Metadata$Tumor_Type[idx] == "Primary"]
BCC_DNA_Data$TopBottomQuarterPrimariesPurIST[idx[1:Number]] = "TopQuarter"
idx = order(BCC_DNA_Data$Metadata$PurIST_Score, decreasing = FALSE)
idx = idx[BCC_DNA_Data$Metadata$Tumor_Type[idx] == "Primary"]
BCC_DNA_Data$TopBottomQuarterPrimariesPurIST[idx[1:Number]] = "BottomQuarter"

Number = floor(length(which(!is.na(BCC_DNA_Data$Metadata$PurIST_Score) &
                              BCC_DNA_Data$Metadata$Tumor_Type == "Met")) / 4)
BCC_DNA_Data$TopBottomQuarterMetsPurIST = rep(NA, times=dim(BCC_DNA_Data$Metadata)[[1]])
idx = order(BCC_DNA_Data$Metadata$PurIST_Score, decreasing = TRUE)
idx = idx[BCC_DNA_Data$Metadata$Tumor_Type[idx] == "Met"]
BCC_DNA_Data$TopBottomQuarterMetsPurIST[idx[1:Number]] = "TopQuarter"
idx = order(BCC_DNA_Data$Metadata$PurIST_Score, decreasing = FALSE)
idx = idx[BCC_DNA_Data$Metadata$Tumor_Type[idx] == "Met"]
BCC_DNA_Data$TopBottomQuarterMetsPurIST[idx[1:Number]] = "BottomQuarter"

idx = match(x=BCC_DNA_Data$Metadata$Public_Specimen_ID,
            table=rownames(BCC_RNA_Data$GSVA_pORGpSUB_Primaries))
BCC_DNA_Data$Metadata$GSVA_pORG_Up_55_Primaries = BCC_RNA_Data$GSVA_pORGpSUB_Primaries[idx, "pORG_Up_55"]
BCC_DNA_Data$Metadata$GSVA_pSUB_Up_51_Primaries = BCC_RNA_Data$GSVA_pORGpSUB_Primaries[idx, "pSUB_Up_51"]

idx = match(x=BCC_DNA_Data$Metadata$Public_Specimen_ID,
            table=rownames(BCC_RNA_Data$GSVA_pORGpSUB_Mets))
BCC_DNA_Data$Metadata$GSVA_pORG_Up_55_Mets = BCC_RNA_Data$GSVA_pORGpSUB_Mets[idx, "pORG_Up_55"]
BCC_DNA_Data$Metadata$GSVA_pSUB_Up_51_Mets = BCC_RNA_Data$GSVA_pORGpSUB_Mets[idx, "pSUB_Up_51"]

Number = floor(length(which(!is.na(BCC_DNA_Data$Metadata$GSVA_pORG_Up_55_Primaries))) / 4)
BCC_DNA_Data$pORG_Primaries_TopBottomQuarter = rep(NA, times=dim(BCC_DNA_Data$Metadata)[[1]])
idx = order(BCC_DNA_Data$Metadata$GSVA_pORG_Up_55_Primaries, decreasing = TRUE)
BCC_DNA_Data$pORG_Primaries_TopBottomQuarter[idx[1:Number]] = "TopQuarter"
idx = order(BCC_DNA_Data$Metadata$GSVA_pORG_Up_55_Primaries, decreasing = FALSE)
BCC_DNA_Data$pORG_Primaries_TopBottomQuarter[idx[1:Number]] = "BottomQuarter"
BCC_DNA_Data$pSUB_Primaries_TopBottomQuarter = rep(NA, times=dim(BCC_DNA_Data$Metadata)[[1]])
idx = order(BCC_DNA_Data$Metadata$GSVA_pSUB_Up_51_Primaries, decreasing = TRUE)
BCC_DNA_Data$pSUB_Primaries_TopBottomQuarter[idx[1:Number]] = "TopQuarter"
idx = order(BCC_DNA_Data$Metadata$GSVA_pSUB_Up_51_Primaries, decreasing = FALSE)
BCC_DNA_Data$pSUB_Primaries_TopBottomQuarter[idx[1:Number]] = "BottomQuarter"

Number = floor(length(which(!is.na(BCC_DNA_Data$Metadata$GSVA_pORG_Up_55_Mets))) / 4)
BCC_DNA_Data$pORG_Mets_TopBottomQuarter = rep(NA, times=dim(BCC_DNA_Data$Metadata)[[1]])
idx = order(BCC_DNA_Data$Metadata$GSVA_pORG_Up_55_Mets, decreasing = TRUE)
BCC_DNA_Data$pORG_Mets_TopBottomQuarter[idx[1:Number]] = "TopQuarter"
idx = order(BCC_DNA_Data$Metadata$GSVA_pORG_Up_55_Mets, decreasing = FALSE)
BCC_DNA_Data$pORG_Mets_TopBottomQuarter[idx[1:Number]] = "BottomQuarter"
BCC_DNA_Data$pSUB_Mets_TopBottomQuarter = rep(NA, times=dim(BCC_DNA_Data$Metadata)[[1]])
idx = order(BCC_DNA_Data$Metadata$GSVA_pSUB_Up_51_Mets, decreasing = TRUE)
BCC_DNA_Data$pSUB_Mets_TopBottomQuarter[idx[1:Number]] = "TopQuarter"
idx = order(BCC_DNA_Data$Metadata$GSVA_pSUB_Up_51_Mets, decreasing = FALSE)
BCC_DNA_Data$pSUB_Mets_TopBottomQuarter[idx[1:Number]] = "BottomQuarter"

write.table(BCC_DNA_Data$Metadata, 
            paste(AnalysisOutputPath, "DNA_Specimen_Metadata_Plus.tsv", sep=""), 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#= Generate OncoPrint matrix ===================================================

# This includes all samples with panels, even without any mutations found.
ColNames = unique(BCC_DNA_Data$Metadata$Public_Specimen_ID)
RowNames = unique(BCC_DNA_Data$Variants$Gene)
OncoPrintM = matrix(data="", nrow = length(RowNames), ncol = length(ColNames), dimnames = list(RowNames, ColNames))

PriorityIdx = which(!is.na(BCC_DNA_Data$Variants$Public_Specimen_ID))
PriorityIdx = PriorityIdx[which(BCC_DNA_Data$Variants$ImpactLevel[PriorityIdx] <= 3)]
PriorityIdx = PriorityIdx[order(BCC_DNA_Data$Variants$ImpactLevel[PriorityIdx], decreasing=TRUE)]
for (idx in PriorityIdx)
{
  if (OncoPrintM[BCC_DNA_Data$Variants$Gene[idx], BCC_DNA_Data$Variants$Public_Specimen_ID[idx]] == "")
  {
    OncoPrintM[BCC_DNA_Data$Variants$Gene[idx], BCC_DNA_Data$Variants$Public_Specimen_ID[idx]] = BCC_DNA_Data$Variants$mutation_priority[idx]
  } else
  {
    OncoPrintM[BCC_DNA_Data$Variants$Gene[idx], BCC_DNA_Data$Variants$Public_Specimen_ID[idx]] = 
      paste(OncoPrintM[BCC_DNA_Data$Variants$Gene[idx], BCC_DNA_Data$Variants$Public_Specimen_ID[idx]], BCC_DNA_Data$Variants$mutation_priority[idx], sep=";")
  }
}
NumMutated = apply(X=OncoPrintM[, ] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

#= Setup OncoPrint functions ===================================================

PlottingPath = paste(AnalysisOutputPath, "OncoPrints/", sep="")
suppressWarnings(dir.create(PlottingPath, recursive=TRUE))

col = c("Missense" = "blue", 
        "Frameshift" = "red", 
        "Inframe_indels" = "brown", 
        "Stop_gains" = "hotpink", 
        "Protein_protein_contact" = "green",
        "Splice_region_variants" = "yellow", 
        "Copy_number_gain" = "purple", 
        "Copy_number_loss" = "black",
        "Other_variant" = "orange")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h-unit(0, "pt"), 
              gp = gpar(fill = "gray", col = NA))
  },
  # big blue
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h-unit(0, "pt"), 
              gp = gpar(fill = col["Missense"], col = NA))
  },
  # big red
  Frameshift = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h-unit(0, "pt"), 
              gp = gpar(fill = col["Frameshift"], col = NA))
  },
  # big brown
  Inframe_indels = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h-unit(0, "pt"), 
              gp = gpar(fill = col["Inframe_indels"], col = NA))
  },
  # big pink
  Stop_gains = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h-unit(0, "pt"), 
              gp = gpar(fill = col["Stop_gains"], col = NA))
  },
  # small green
  Protein_protein_contact = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["Protein_protein_contact"], col = NA))
  },
  # small yellow
  Splice_region_variants = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h*0.33, 
              gp = gpar(fill = col["Splice_region_variants"], col = NA))
  },
  # small purple
  Copy_number_gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h*0.33, 
              gp = gpar(fill = col["Copy_number_gain"], col = NA))
  },
  # small black
  Copy_number_loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h*0.33, 
              gp = gpar(fill = col["Copy_number_loss"], col = NA))
  },
  # small white
  Other_variant = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h*0.33, 
              gp = gpar(fill = col["Other_variant"], col = NA))
  }
)

heatmap_legend_param = list(title = "Alternations", at = names(col), 
                            labels = gsub("_", " ", names(col)))

# See this for annotation options/ideas:
# https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html

#= All_mutations_for_legend ===========================================

ReportName = "All_mutations_for_legend"
x = unlist(strsplit(OncoPrintM, ";"))
unique(x)
xrle = rle(x[order(x)])
i = order(xrle$lengths)
cbind(xrle$lengths[i], xrle$values[i])

idx = which(BCC_DNA_Data$Variants$mutation_priority == "other_variant")
BCC_DNA_Data$Variants$Gene[idx]

SelectedGenes = which(!(rownames(OncoPrintM) %in% c(BCC_DNA_Data$Variants$Gene[idx])))

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[SelectedGenes, ],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", dim(OncoPrintM)[[2]], ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Top4th_pORG_Primaries Top10 Genes (Fig 2 - F) ===============================

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$pORG_Primaries_TopBottomQuarter == "TopQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Top4th_pORG_Primaries"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Bottom4th_pORG_Primaries Top10 Genes (Fig 2 - F) ============================

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$pORG_Primaries_TopBottomQuarter == "BottomQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Bottom4th_pORG_Primaries"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= All Samples Top10 Genes (Supp Fig 7 - A) ====================================

SelectedColumns = 1:dim(OncoPrintM)[[2]]

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "All_Samples"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=11, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")
#filetype="png")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Top4th_pORG_Mets Top10 Genes (Supp Fig 7 - B) ===============================

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Met" &
                                                    BCC_DNA_Data$pORG_Mets_TopBottomQuarter == "TopQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Top4th_pORG_Mets"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Bottom4th_pORG_Mets Top10 Genes (Supp Fig 7 - B) ============================

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Met" &
                                                    BCC_DNA_Data$pORG_Mets_TopBottomQuarter == "BottomQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Bottom4th_pORG_Mets"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Top4th_pORG_PrimariesHR_DDR_Genes (Supp Fig 7 - B) ==========================

SelectedGenes = which(rownames(OncoPrintM) %in% BCC_DNA_Data$HR_DDR_Genes)

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$pORG_Primaries_TopBottomQuarter == "TopQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Top4th_pORG_PrimariesHR_DDR_Genes"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[SelectedGenes, SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Bottom4th_pORG_PrimariesHR_DDR_Genes (Supp Fig 7 - B) =======================

SelectedGenes = which(rownames(OncoPrintM) %in% BCC_DNA_Data$HR_DDR_Genes)

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$pORG_Primaries_TopBottomQuarter == "BottomQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Bottom4th_pORG_PrimariesHR_DDR_Genes"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[SelectedGenes, SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= LiverCohortPrimaries Top10 Genes (Supp Fig 7 - C) ===========================
SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$LiverLungCohort == "Liver"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "LiverCohortPrimaries"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}
#= LungCohortPrimaries Top10 Genes (Supp Fig 7 - C) ============================
SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$LiverLungCohort == "Lung"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "LungCohortPrimaries"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}
#= LiverCohortMets Top10 Genes (Supp Fig 7 - C) ================================
SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Met" &
                                                    BCC_DNA_Data$LiverLungCohort == "Liver"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "LiverCohortMets"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}
#= LungCohortMets Top10 Genes (Supp Fig 7 - C) =================================

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Met" &
                                                    BCC_DNA_Data$LiverLungCohort == "Lung"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "LungCohortMets"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}
#= LiverCohortPrimariesHR_DDR_Genes (Supp Fig 7 - C) ===========================

SelectedGenes = which(rownames(OncoPrintM) %in% BCC_DNA_Data$HR_DDR_Genes)

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$LiverLungCohort == "Liver"])
length(SelectedColumns)

ReportName = "LiverCohortPrimariesHR_DDR_Genes"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[SelectedGenes, SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}
#= LungCohortPrimariesHR_DDR_Genes (Supp Fig 7 - C) ============================

SelectedGenes = which(rownames(OncoPrintM) %in% BCC_DNA_Data$HR_DDR_Genes)

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$LiverLungCohort == "Lung"])
length(SelectedColumns)

ReportName = "LungCohortPrimariesHR_DDR_Genes"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[SelectedGenes, SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}
#= Top4th_pSUB_Primaries Top10 Genes (Supp Fig 8 - A) ==========================

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$pSUB_Primaries_TopBottomQuarter == "TopQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Top4th_pSUB_Primaries"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Bottom4th_pSUB_Primaries Top10 Genes (Supp Fig 8 - A) =======================

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$pSUB_Primaries_TopBottomQuarter == "BottomQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Bottom4th_pSUB_Primaries"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Top4th_pSUB_Mets Top10 Genes (Supp Fig 8 - A) ===============================

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Met" &
                                                    BCC_DNA_Data$pSUB_Mets_TopBottomQuarter == "TopQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Top4th_pSUB_Mets"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Bottom4th_pSUB_Mets Top10 Genes (Supp Fig 8 - A) ============================

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Met" &
                                                    BCC_DNA_Data$pSUB_Mets_TopBottomQuarter == "BottomQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Bottom4th_pSUB_Mets"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Top4th_pSUB_PrimariesHR_DDR_Genes (Supp Fig 8 - A) ==========================

SelectedGenes = which(rownames(OncoPrintM) %in% BCC_DNA_Data$HR_DDR_Genes)

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$pSUB_Primaries_TopBottomQuarter == "TopQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Top4th_pSUB_PrimariesHR_DDR_Genes"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[SelectedGenes, SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Bottom4th_pSUB_PrimariesHR_DDR_Genes (Supp Fig 8 - A) =======================

SelectedGenes = which(rownames(OncoPrintM) %in% BCC_DNA_Data$HR_DDR_Genes)

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$pSUB_Primaries_TopBottomQuarter == "BottomQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Bottom4th_pSUB_PrimariesHR_DDR_Genes"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[SelectedGenes, SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}



#= Top4thPurIST_Primaries Top10 Genes (Supp Fig 8 - B) =========================
SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$TopBottomQuarterPrimariesPurIST == "TopQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Top4thPurIST_Primaries"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}
#= Bottom4thPurIST_Primaries Top10 Genes (Supp Fig 8 - B) ======================
SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$TopBottomQuarterPrimariesPurIST == "BottomQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Bottom4thPurIST_Primaries"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}
#= Top4thPurIST_Mets Top10 Genes (Supp Fig 8 - B) ==============================

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Met" &
                                                    BCC_DNA_Data$TopBottomQuarterMetsPurIST == "TopQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Top4thPurIST_Mets"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Bottom4thPurIST_Mets Top10 Genes (Supp Fig 8 - B) ===========================
SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Met" &
                                                    BCC_DNA_Data$TopBottomQuarterMetsPurIST == "BottomQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Bottom4thPurIST_Mets"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}
#= Top4thPurIST_PrimariesHR_DDR_Genes (Supp Fig 8 - B) =========================

SelectedGenes = which(rownames(OncoPrintM) %in% BCC_DNA_Data$HR_DDR_Genes)

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$TopBottomQuarterPrimariesPurIST == "TopQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Top4thPurIST_PrimariesHR_DDR_Genes"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[SelectedGenes, SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Bottom4thPurIST_PrimariesHR_DDR_Genes (Supp Fig 8 - B) ======================

SelectedGenes = which(rownames(OncoPrintM) %in% BCC_DNA_Data$HR_DDR_Genes)

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary" &
                                                    BCC_DNA_Data$TopBottomQuarterPrimariesPurIST == "BottomQuarter"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Bottom4thPurIST_PrimariesHR_DDR_Genes"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[SelectedGenes, SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}

#= Primaries Top10 Genes (Supp Fig 9 - F) ======================================
SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Primary"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Primaries"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}
#= Mets Top10 Genes (Supp Fig 9 - F) ===========================================

SelectedColumns = which(colnames(OncoPrintM) %in% 
                        BCC_DNA_Data$Metadata$Public_Specimen_ID[BCC_DNA_Data$Metadata$Tumor_Type == "Met"])
length(SelectedColumns)

NumMutated = apply(X=OncoPrintM[, SelectedColumns] != "", MARGIN=1, FUN=sum)
TopMutated = order(NumMutated, decreasing = TRUE)

ReportName = "Mets"

PrintToFile = GraphSetup(width=6, height=4, 
                         pointsize=12, 
                         ReportName=paste(PlottingPath, "OncoPrint_", ReportName, sep=""), 
                         filetype="cairopdf")

print(oncoPrint(OncoPrintM[TopMutated[1:10], SelectedColumns],
                get_type = function(x) strsplit(x, ";")[[1]],                 
                alter_fun = alter_fun, col = col, 
                remove_empty_columns = FALSE, remove_empty_rows = FALSE,
                column_title = paste(ReportName, " (N=", length(SelectedColumns), ")", sep=""), 
                heatmap_legend_param = heatmap_legend_param,
                column_title_gp = gpar(fontsize = 12, fontface = "bold")))

if (PrintToFile)
{
  # Stop outputting to file.
  dev.off()
}
}
#= The following produces some of the Kaplan-Meier plots used. =================
{
  PlottingPath = paste(AnalysisOutputPath, "SurvivalPlots/", sep="")
  suppressWarnings(dir.create(PlottingPath, recursive=TRUE))
  
  # BCC data -
  # We only use patients with RNA-Seq from a primary tumor resection.  We use 
  # Days_from_Diagnosis_to_FU for overall survival.  We need vital status at
  # last follow up to use for censoring.  
  idx = which(BCC_RNA_Data$Metadata$Tumor_Type == "Primary")
  BCC_KM_TestData = data.frame(Public_Patient_ID = BCC_RNA_Data$Metadata$Public_Patient_ID[idx],
                               Public_Specimen_ID = BCC_RNA_Data$Metadata$Public_Specimen_ID[idx],
                               Event = BCC_RNA_Data$Metadata$Vital_Status_at_FU[idx] == "Dead",
                               Survival = BCC_RNA_Data$Metadata$Days_from_Diagnosis_to_FU[idx],
                               ResectionSurvival = BCC_RNA_Data$Metadata$Days_from_Resection_to_FU[idx],
                               PurIST_Score = BCC_RNA_Data$Metadata$PurIST_Score[idx],
                               HR_DDR_Altered_All = BCC_RNA_Data$Metadata$HR_DDR_Altered_All[idx],
                               HR_DDR_Altered_BioRel = BCC_RNA_Data$Metadata$HR_DDR_Altered_BioRel[idx],
                               stringsAsFactors = FALSE)
  dim(BCC_KM_TestData)
  
  # Merge in pORG and pSUB scores for Primary samples.
  idx = match(x=BCC_KM_TestData$Public_Specimen_ID, rownames(BCC_RNA_Data$GSVA_pORGpSUB_Primaries))
  BCC_KM_TestData$pORG_Up_55 = BCC_RNA_Data$GSVA_pORGpSUB_Primaries[idx, "pORG_Up_55"]
  BCC_KM_TestData$pSUB_Up_51 = BCC_RNA_Data$GSVA_pORGpSUB_Primaries[idx, "pSUB_Up_51"]
  dim(BCC_KM_TestData)
  
  # Exclude any patients that are missing the survival time or event data.
  BCC_KM_TestData = BCC_KM_TestData[!is.na(BCC_KM_TestData$Event), ]
  dim(BCC_KM_TestData)
  BCC_KM_TestData = BCC_KM_TestData[!is.na(BCC_KM_TestData$Survival), ]
  dim(BCC_KM_TestData)
  
  # Exclude patients that died in the first month after resection (lived 30 days 
  # or less after resection).  These patients likely died from a complication of 
  # surgery.)
  BCC_KM_TestData = BCC_KM_TestData[is.na(BCC_KM_TestData$ResectionSurvival) | 
                                      BCC_KM_TestData$ResectionSurvival > 30, ]
  dim(BCC_KM_TestData)
  
  # Only include each patient once, there are two patients with two primary 
  # tumors.
  BCC_KM_TestData$Public_Specimen_ID[which(duplicated(BCC_KM_TestData$Public_Patient_ID))]
  BCC_KM_TestData = BCC_KM_TestData[!duplicated(BCC_KM_TestData$Public_Patient_ID), ]
  dim(BCC_KM_TestData)
  
  # Here we are selecting the best cutoff for pORG/pSUB/PurIST by maximizing the
  # Youden's Index for a ROC curve.  This is based on choosing a good vs bad
  # outcome which, in turn, is based on a "clinically relevant" survival 
  # threshold of 545 days (~ 18 months).  We are selecting the cutoff based on
  # the BCC data and using that same threshold on the public data sets for 
  # validation.
  SurvivalThreshold = 545
  predSign = -1
  DataToTest = "BCC"
  TestData = BCC_KM_TestData
  
  for (TestCol in c("pORG_Up_55", "pSUB_Up_51", "PurIST_Score"))
  {
    CensorNotice = NA
    GoodIdx = which(TestData$Survival > SurvivalThreshold)
    BadIdx = which(TestData$Survival <= SurvivalThreshold)
    # Check for patients that are alive at last follow up.  These patients 
    # should not be included in the bad outcome since they likely will live
    # long enough to be in the good outcome group.
    idx = which(TestData$Event[BadIdx] == FALSE)
    if (length(idx) > 0) BadIdx = BadIdx[-idx]
    CensorNotice = paste(DataToTest, " (", SurvivalThreshold, " Days): ", 
                         length(idx), " censored patients removed from bad",
                         " outcome group for Youden's Index optimization.",
                         sep="")
    print(CensorNotice)
    
    lab2test = c(rep(0, times=length(BadIdx)), rep(1, times=length(GoodIdx))) 
    pred2test = predSign * TestData[[TestCol]][c(BadIdx, GoodIdx)]
    ROCit_obj <- rocit(score=pred2test,class=lab2test)
    summary(ROCit_obj)
    YoudenIndex = ROCit_obj$TPR - ROCit_obj$FPR
    i = order(YoudenIndex, decreasing = TRUE)
    MaxYoudenIndexIdx = i[1]
    Cutoff = predSign * ROCit_obj$Cutoff[MaxYoudenIndexIdx]
    print(Cutoff)
    if (TestCol == "pORG_Up_55") pORG_Up_55_Cutoff = Cutoff
    if (TestCol == "pSUB_Up_51") pSUB_Up_51_Cutoff = Cutoff
    if (TestCol == "PurIST_Score") PurIST_Score_Cutoff = Cutoff
  }

  #= Load ICGC and TCGA data and calculate pORG/pSUB/PurIST scores. ============
  # Load two public PDAC data sets from ICGC and TCGA.
  ICGC = Load_ICGC_Data(paste(BaseDir, sep=""))
  TCGA = Load_TCGA_Data(paste(BaseDir, sep=""))
  
  # Calculate GSVA pORG/pSUB scores and PurIST scores for ICGC and TCGA data.
  M = ICGC$Data
  Ensembl_ID = rownames(M)
  idx = match(x=Ensembl_ID, table=BCC_RNA_Data$EnsemblAnn[, "ensembl_gene_id"])
  M = M[!is.na(idx), ]
  rownames(M) = BCC_RNA_Data$EnsemblAnn[idx[!is.na(idx)], c("hgnc_symbol")]
  # Note that we are leaving genes with blank symbols in to to preserve the 
  # "background".
  M = 2^M
  M = log2(M + 1)
  ICGC$Metadata$GSVA_pORGpSUB = ApplyGSVA(M, BCC_RNA_Data$pORGpSUB_GeneSetsList, RowCV_ThreshFactor = 0.1)
  dim(ICGC$Metadata$GSVA_pORGpSUB)
  rownames(ICGC$Metadata$GSVA_pORGpSUB) = colnames(M)
  
  TCGA$Metadata$GSVA_pORGpSUB = ApplyGSVA(log2(TCGA$Data + 1), BCC_RNA_Data$pORGpSUB_GeneSetsList, RowCV_ThreshFactor = 0.1)
  dim(TCGA$Metadata$GSVA_pORGpSUB)
  rownames(TCGA$Metadata$GSVA_pORGpSUB) = colnames(TCGA$Data)

  # Do PurIST subtyping of ICGC and TCGA samples.
  # We want gene length normalized data for PurIST since it is using the 
  # difference between pairs of genes.  We will check the correlation of ICGC
  # and TCGA data with gene length to see if it is similar to our TPM data.
  # I conclude from spearman correlations that we do NOT want to scale
  # these values for gene length (it is likely that these were 3-prime libraries).
  M = ICGC$Data
  Ensembl_ID = rownames(M)
  idx = match(x=Ensembl_ID, table=BCC_RNA_Data$EnsemblAnn[, "ensembl_gene_id"])
  M = M[!is.na(idx), ]
  rownames(M) = BCC_RNA_Data$EnsemblAnn[idx[!is.na(idx)], c("hgnc_symbol")]
  M = M[rownames(M) != "", ]
  M = 2^M
  PurIST_Result = ApplyPurIST(ScriptsPath = ScriptsPath, 
                              Data = M)
  ICGC$Metadata$PurIST_Subtype = PurIST_Result$PurIST_Subtype
  ICGC$Metadata$PurIST_Detailed_Subtype = PurIST_Result$PurIST_Detailed_Subtype
  ICGC$Metadata$PurIST_Score = PurIST_Result$PurIST_Score
  
  M = TCGA$Data
  PurIST_Result = ApplyPurIST(ScriptsPath = ScriptsPath, 
                              Data = M)
  TCGA$Metadata$PurIST_Subtype = PurIST_Result$PurIST_Subtype
  TCGA$Metadata$PurIST_Detailed_Subtype = PurIST_Result$PurIST_Detailed_Subtype
  TCGA$Metadata$PurIST_Score = PurIST_Result$PurIST_Score

  #= Set up data frames for Kaplan Meier plots. ================================
  
  # ICGC data -
  
  ICGC_KM_TestData = data.frame(Public_Patient_ID = ICGC$Metadata$ICGC_ID,
                                Event = suppressWarnings(as.numeric(ICGC$Metadata$DFS.Censor) != 1),
                                Survival = suppressWarnings(as.numeric(ICGC$Metadata$D_F_S) * 30),
                                PurIST_Score = ICGC$Metadata$PurIST_Score,
                                pORG_Up_55 = ICGC$Metadata$GSVA_pORGpSUB[, "pORG_Up_55"],
                                pSUB_Up_51 = ICGC$Metadata$GSVA_pORGpSUB[, "pSUB_Up_51"],
                                stringsAsFactors = FALSE)
  
  dim(ICGC_KM_TestData)
  
  # Exclude any patients that are missing the survival time or event data.
  ICGC_KM_TestData = ICGC_KM_TestData[!is.na(ICGC_KM_TestData$Event), ]
  dim(ICGC_KM_TestData)
  ICGC_KM_TestData = ICGC_KM_TestData[!is.na(ICGC_KM_TestData$Survival), ]
  dim(ICGC_KM_TestData)
  
  # Exclude patients that lived 30 days or less (we don't have the data for 
  # how long they lived after resection, but this should be very similar).
  #ICGC_KM_TestData = ICGC_KM_TestData[is.na(ICGC_KM_TestData$Survival) | 
  #                                    ICGC_KM_TestData$Survival > 30, ]
  #dim(ICGC_KM_TestData)
  
  # Check for duplicates
  which(duplicated(ICGC_KM_TestData$Public_Specimen_ID))
  
  # TCGA data -
  
  TCGA_KM_TestData = data.frame(Public_Patient_ID = TCGA$Metadata$Patient.ID,
                                Public_Specimen_ID = TCGA$Metadata$Sample.ID,
                                Event = TCGA$Metadata$Overall.Survival.Status == "1:DECEASED",
                                Survival = TCGA$Metadata$Overall.Survival..Months. * 30,
                                PurIST_Score = TCGA$Metadata$PurIST_Score,
                                pORG_Up_55 = TCGA$Metadata$GSVA_pORGpSUB[, "pORG_Up_55"],
                                pSUB_Up_51 = TCGA$Metadata$GSVA_pORGpSUB[, "pSUB_Up_51"],
                                stringsAsFactors = FALSE)
  dim(TCGA_KM_TestData)
  
  # Exclude any patients that are missing the survival time or event data.
  TCGA_KM_TestData = TCGA_KM_TestData[!is.na(TCGA_KM_TestData$Event), ]
  dim(TCGA_KM_TestData)
  TCGA_KM_TestData = TCGA_KM_TestData[!is.na(TCGA_KM_TestData$Survival), ]
  dim(TCGA_KM_TestData)
  
  # Exclude patients that lived 30 days or less (we don't have the data for 
  # how long they lived after resection, but this should be very similar).
  TCGA_KM_TestData = TCGA_KM_TestData[is.na(TCGA_KM_TestData$Survival) | 
                                        TCGA_KM_TestData$Survival > 30, ]
  dim(TCGA_KM_TestData)
  
  # Check for duplicates
  which(duplicated(TCGA_KM_TestData$Public_Patient_ID))
  which(duplicated(TCGA_KM_TestData$Public_Specimen_ID))
  
  for (TestCol in c("pORG_Up_55", "pSUB_Up_51", "PurIST_Score"))
  {
    if (TestCol == "pORG_Up_55") Cutoff = pORG_Up_55_Cutoff
    if (TestCol == "pSUB_Up_51") Cutoff = pSUB_Up_51_Cutoff
    if (TestCol == "PurIST_Score") Cutoff = PurIST_Score_Cutoff
    
    CensorNotice = NA
    for (DataToTest in c("BCC", "ICGC", "TCGA"))
    {
      if (DataToTest == "BCC")
      {
        TestData = BCC_KM_TestData

        for (filetype in c("cairopdf"))
        {
          width <- 6; height <- 4; pointsize <- 12
          
          PrintToFile = GraphSetup(width=width, height=height, pointsize=pointsize, 
                                   ReportName=paste(PlottingPath, "ROCit_Plot_", DataToTest, "_", TestCol, "_", SurvivalThreshold, "_Days", sep=""), 
                                   filetype=filetype)
          
          plot(ROCit_obj)
          
          if (PrintToFile) { dev.off() } # Stop outputting to file.
        }
      }
      
      if (DataToTest == "ICGC")
      {
        TestData = ICGC_KM_TestData
      }
      
      if (DataToTest == "TCGA")
      {
        TestData = TCGA_KM_TestData
      }
      
      TestData$TestState = "Undefined"
      TestData$TestState[TestData[[TestCol]] > Cutoff] = "High"
      TestData$TestState[TestData[[TestCol]] <= Cutoff] = "Low"
      length(which(TestData$TestState == "Low"))
      length(which(TestData$TestState == "High"))
      
      pval = 1
      if (length(which(TestData$TestState == "Low")) > 0 &
          length(which(TestData$TestState == "High")) > 0)
      {
        cox = coxph(Surv(Survival, Event) ~ TestState, data = TestData)
        print(summary(cox))
        pval = as.numeric(summary(cox)$sctest["pvalue"])
      }
      
      PlotName = paste(DataToTest, "_", TestCol, "_", SurvivalThreshold, "_Days", sep="")
      
      km_fit = survfit(Surv(Survival, Event) ~ TestState, data = TestData)
      
      # See: help(autoplot.survfit)
      
      for (filetype in c("cairopdf"))
      {
        width <- 6; height <- 4; pointsize <- 12
        
        PrintToFile = GraphSetup(width=width, height=height, 
                                 pointsize=pointsize, 
                                 ReportName=paste(PlottingPath, "KaplanMeier_", PlotName, sep=""), 
                                 filetype=filetype)
        
        kmPlot = autoplot(km_fit,
                          surv.size = 1,
                          conf.int = TRUE,
                          conf.int.alpha = 0.1)
        kmPlot = kmPlot + 
          ggtitle(label = paste(DataToTest, ": ", TestCol, " cutoff = ", Cutoff, " from max Youden point", sep=""),
                  subtitle = paste("pVal = ", formatC(pval, format="e", digits=4),
                                   " (", length(which(TestData$TestState == "Low")), " Low, ", 
                                   length(which(TestData$TestState == "High")), " High)", sep="")) +
          labs(x = "Timeline (Days)", y = "Probability of Survival") +
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank())
        #theme(plot.title = element_text(hjust = 0.5),
        #      plot.subtitle = element_text(hjust = 0.5))
        print(kmPlot)
        
        if (PrintToFile)
        {
          # Stop outputting to file.
          dev.off()
        }
      }
      
      Report = c("SurvivalThreshold", SurvivalThreshold, "", "")
      Report = rbind(Report, c("From BCC data:", "", "", ""))
      Report = rbind(Report, c("MaxYoudenIndex", YoudenIndex[MaxYoudenIndexIdx], "", ""))
      Report = rbind(Report, c("Cutoff", Cutoff, "", ""))
      Report = rbind(Report, c(CensorNotice, "", "", ""))
      Report = rbind(Report, c("", "", "", ""))
      Report = rbind(Report, c(DataToTest, "", "", ""))
      Report = rbind(Report, c("", "", "", ""))
      Report = rbind(Report, c(paste("coxph sctest pvalue (logrank) = ", pval, sep=""), "", "", ""))
      Report = rbind(Report, c("", "", "", ""))
      Report = rbind(Report, c(TestCol, "TestState", "Event", "Survival"))
      Report = as.matrix(Report)
      Report = rbind(Report, as.matrix(TestData[order(TestData[, TestCol]), c(TestCol, "TestState", "Event", "Survival")]))
      
      write.table(Report, 
                  paste(PlottingPath, "Report_", PlotName, ".txt", sep=""), 
                  sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    }
  }

  TestCol = "pORG_Up_55"
  DataToTest = "BCC"
  TestData = BCC_KM_TestData
  Cutoff = pORG_Up_55_Cutoff
  
  # Remove patients without xT DNA panels.
  idx = which(is.na(TestData$HR_DDR_Altered_All))
  if (length(idx) > 0) TestData = TestData[-idx, ]
  
  ### Use all alterations. ### 
  {
    # Patients with DDR Altered sample & high pORG:
    DDR_Altered_High = which(!is.na(TestData$HR_DDR_Altered_All) &
                             TestData$HR_DDR_Altered_All  & 
                             TestData[[TestCol]] > Cutoff)
    # Patients with DDR Altered sample & low pORG:
    DDR_Altered_Low = which(!is.na(TestData$HR_DDR_Altered_All) &
                            TestData$HR_DDR_Altered_All  & 
                            TestData[[TestCol]] <= Cutoff)
    # Patients with DDR Intact sample & high pORG:
    DDR_Intact_High = which(!is.na(TestData$HR_DDR_Altered_All) &
                            !TestData$HR_DDR_Altered_All  & 
                            TestData[[TestCol]] > Cutoff)
    # Patients with DDR Intact sample & low pORG:
    DDR_Intact_Low = which(!is.na(TestData$HR_DDR_Altered_All) &
                           !TestData$HR_DDR_Altered_All  & 
                           TestData[[TestCol]] <= Cutoff)

    # Missing "ST-00007307-T" and "ST-00019598-T" from altered DDR set.
    TestData$TestState = "Undefined"
    TestData$TestState[DDR_Altered_High] = "DDR_Altered_High"
    TestData$TestState[DDR_Altered_Low] = "DDR_Altered_Low" 
    TestData$TestState[DDR_Intact_High] = "DDR_Intact_High"
    TestData$TestState[DDR_Intact_Low] = "DDR_Intact_Low"
    
    length(DDR_Altered_High)
    length(DDR_Altered_Low)
    length(DDR_Intact_High)
    length(DDR_Intact_Low)
    length(which(TestData$TestState == "Undefined"))
    
    pval = 1
    if (length(which(TestData$TestState == "DDR_Altered_High")) > 0 &
        length(which(TestData$TestState == "DDR_Altered_Low")) > 0 &
        length(which(TestData$TestState == "DDR_Intact_High")) > 0 &
        length(which(TestData$TestState == "DDR_Intact_Low")) > 0)
    {
      cox = coxph(Surv(Survival, Event) ~ TestState, data = TestData)
      print(summary(cox))
      pval = as.numeric(summary(cox)$sctest["pvalue"])
    }
    
    
    PlotName = paste("HR_DDR_Altered_All_and_pORG_Up_55", sep="")
  
    km_fit = survfit(Surv(Survival, Event) ~ TestState, data = TestData)
    
    # See: help(autoplot.survfit)
    
    for (filetype in c("cairopdf"))
    {
      width <- 6; height <- 4; pointsize <- 12
      
      PrintToFile = GraphSetup(width=width, height=height, 
                               pointsize=pointsize, 
                               ReportName=paste(PlottingPath, "KaplanMeier_", PlotName, sep=""), 
                               filetype=filetype)
      
      kmPlot = autoplot(km_fit,
                        surv.size = 1,
                        conf.int = TRUE,
                        conf.int.alpha = 0.1)
      kmPlot = kmPlot + 
        ggtitle(label = paste(DataToTest, ": ", TestCol, " cutoff = ", Cutoff, " from max Youden point", sep=""),
                subtitle = paste("pVal = ", formatC(pval, format="e", digits=4),
                                 " (", length(which(TestData$TestState == "DDR_Altered_High")), " DDR_Altered_High, ", 
                                 length(which(TestData$TestState == "DDR_Altered_Low")), " DDR_Altered_Low, ", 
                                 length(which(TestData$TestState == "DDR_Intact_High")), " DDR_Intact_High, ", 
                                 length(which(TestData$TestState == "DDR_Intact_Low")), " DDR_Intact_Low)", sep="")) +
        labs(x = "Survival (Days)", y = "Survival Probability") +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())
      #theme(plot.title = element_text(hjust = 0.5),
      #      plot.subtitle = element_text(hjust = 0.5))
      print(kmPlot)
  
      if (PrintToFile)
      {
        # Stop outputting to file.
        dev.off()
      }
    }
  }
  
  ### Use Biologically Relevant alterations. ### 
  {
    # Patients with DDR Altered sample & high pORG:
    DDR_Altered_High = which(!is.na(TestData$HR_DDR_Altered_BioRel) &
                               TestData$HR_DDR_Altered_BioRel  & 
                               TestData[[TestCol]] > Cutoff)
    # Patients with DDR Altered sample & low pORG:
    DDR_Altered_Low = which(!is.na(TestData$HR_DDR_Altered_BioRel) &
                              TestData$HR_DDR_Altered_BioRel  & 
                              TestData[[TestCol]] <= Cutoff)
    # Patients with DDR Intact sample & high pORG:
    DDR_Intact_High = which(!is.na(TestData$HR_DDR_Altered_BioRel) &
                              !TestData$HR_DDR_Altered_BioRel  & 
                              TestData[[TestCol]] > Cutoff)
    # Patients with DDR Intact sample & low pORG:
    DDR_Intact_Low = which(!is.na(TestData$HR_DDR_Altered_BioRel) &
                             !TestData$HR_DDR_Altered_BioRel  & 
                             TestData[[TestCol]] <= Cutoff)
    
    # Missing "ST-00007307-T" and "ST-00019598-T" from altered DDR set.
    TestData$TestState = "Undefined"
    TestData$TestState[DDR_Altered_High] = "DDR_Altered_High"
    TestData$TestState[DDR_Altered_Low] = "DDR_Altered_Low"
    TestData$TestState[DDR_Intact_High] = "DDR_Intact_High"
    TestData$TestState[DDR_Intact_Low] = "DDR_Intact_Low"
    
    length(DDR_Altered_High)
    length(DDR_Altered_Low)
    length(DDR_Intact_High)
    length(DDR_Intact_Low)
    length(which(TestData$TestState == "Undefined"))
    
    pval = 1
    if (length(which(TestData$TestState == "DDR_Altered_High")) > 0 &
        length(which(TestData$TestState == "DDR_Altered_Low")) > 0 &
        length(which(TestData$TestState == "DDR_Intact_High")) > 0 &
        length(which(TestData$TestState == "DDR_Intact_Low")) > 0)
    {
      cox = coxph(Surv(Survival, Event) ~ TestState, data = TestData)
      print(summary(cox))
      pval = as.numeric(summary(cox)$sctest["pvalue"])
    }
    
    
    PlotName = paste("HR_DDR_Altered_BioRel_and_pORG_Up_55", sep="")
    
    km_fit = survfit(Surv(Survival, Event) ~ TestState, data = TestData)
    
    # See: help(autoplot.survfit)
    
    for (filetype in c("cairopdf"))
    {
      width <- 6; height <- 4; pointsize <- 12
      
      PrintToFile = GraphSetup(width=width, height=height, 
                               pointsize=pointsize, 
                               ReportName=paste(PlottingPath, "KaplanMeier_", PlotName, sep=""), 
                               filetype=filetype)
      
      kmPlot = autoplot(km_fit,
                        surv.size = 1,
                        conf.int = TRUE,
                        conf.int.alpha = 0.1)
      kmPlot = kmPlot + 
        ggtitle(label = paste(DataToTest, ": ", TestCol, " cutoff = ", Cutoff, " from max Youden point", sep=""),
                subtitle = paste("pVal = ", formatC(pval, format="e", digits=4),
                                 " (", length(which(TestData$TestState == "DDR_Altered_High")), " DDR_Altered_High, ", 
                                 length(which(TestData$TestState == "DDR_Altered_Low")), " DDR_Altered_Low, ", 
                                 length(which(TestData$TestState == "DDR_Intact_High")), " DDR_Intact_High, ", 
                                 length(which(TestData$TestState == "DDR_Intact_Low")), " DDR_Intact_Low)", sep="")) +
        labs(x = "Survival (Days)", y = "Survival Probability") +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())
      #theme(plot.title = element_text(hjust = 0.5),
      #      plot.subtitle = element_text(hjust = 0.5))
      print(kmPlot)
      
      if (PrintToFile)
      {
        # Stop outputting to file.
        dev.off()
      }
    }
  }
  
}
#= Generate *.gct and *.cls files for running GSEA analysis tool.  =============

# This section generates input files to run the GSEA analysis used in paper.

# We used the Broad Institute GSEA version 4.1.0 (https://www.gsea-msigdb.org)
# with halllmark gene sets: h.all.v7.2.symbols.gmt
# We ran using a CLI, but we recommend using the GUI with the following 
# parameters and the gct and cls file generated below.

# param	zip_report	false
# param	num	100
# param	scoring_scheme	weighted
# param	norm	meandiv
# param	mode	Max_probe
# param	include_only_symbols	true
# param	set_max	2000
# param	gmx	h.all.v7.2.symbols.gmt
# param	plot_top_x	20
# param	nperm	1000
# param	order	descending
# param	rpt_label	my_analysis
# param	rnd_seed	timestamp
# param	set_min	10
# param	create_svgs	false
# param	sort	real
# param	create_gcts	false
# param	help	false
# param	save_rnd_lists	false
# param	median	false
# param	metric	Signal2Noise
# param	make_sets	true
# param	rnd_type	no_balance
# param	gui	false
# param	permute	phenotype
# param	collapse	No_Collapse

# We use edgeR_TMM, may be more robust than TPMs, but is not gene length normalized.
# We collapsed our data to unique gene symbols before running GSEA.  This could
# have been handled by the GSEA software and/or MSigDB genesets with Ensembl
# IDs could have been used instead of the HUGO based sets.

M_ForGSEA = BCC_RNA_Data$edgeR_TMM
Ensembl_ID = rownames(M_ForGSEA)
idx = match(x=Ensembl_ID, table=BCC_RNA_Data$EnsemblAnn[, "ensembl_gene_id"])
AnnIdx_ForGSVA = which(!is.na(idx))
M_ForGSEA = M_ForGSEA[!is.na(idx), ]
rownames(M_ForGSEA) = BCC_RNA_Data$EnsemblAnn[idx[!is.na(idx)], c("hgnc_symbol")]
AnnIdx_ForGSVA = AnnIdx_ForGSVA[rownames(M_ForGSEA) != ""]
M_ForGSEA = M_ForGSEA[rownames(M_ForGSEA) != "", ]
AnnIdx_ForGSVA = AnnIdx_ForGSVA[!duplicated(rownames(M_ForGSEA))]
M_ForGSEA = M_ForGSEA[!duplicated(rownames(M_ForGSEA)), ]
dim(M_ForGSEA)

PathForGSEA_Files = paste(AnalysisOutputPath, "FilesForGSEA/", sep="")
suppressWarnings(dir.create(PathForGSEA_Files, recursive=TRUE))

### pORG_TopVsBottom_Pri_Up/Dn ###

SelectedSampleIDs = rownames(BCC_RNA_Data$GSVA_pORGpSUB_Primaries)
SortValue = BCC_RNA_Data$GSVA_pORGpSUB_Primaries[, "pORG_Up_55"]
Number = floor(length(SelectedSampleIDs) / 4)
idx = order(SortValue, decreasing = TRUE)
TopIDs = SelectedSampleIDs[idx[1:Number]]
idx = order(SortValue, decreasing = FALSE)
BottomIDs = SelectedSampleIDs[idx[1:Number]]

TopQuarterIdx = which(colnames(M_ForGSEA) %in% TopIDs)
BottomQuarterIdx = which(colnames(M_ForGSEA) %in% BottomIDs)

SetA = TopQuarterIdx
SetAName = "Top4thPri_pORG_Up_55"
SetB = BottomQuarterIdx
SetBName = "Bottom4thPri_pORG_Up_55"

NonZeroIdx = which(rowSums(M_ForGSEA[, c(SetA, SetB)]) > 0)
FormatInputForGSEA(Data = M_ForGSEA[NonZeroIdx, c(SetA, SetB)],
                   SetA = SetA,
                   SetAName = SetAName,
                   SetB = SetB,
                   SetBName = SetBName,
                   OutputPath = PathForGSEA_Files,
                   FileNameForGSEA = "GSEA_pORG_TopVsBottom_Pri")
  
### pSUB_TopVsBottom_Pri_Up/Dn ###

SelectedSampleIDs = rownames(BCC_RNA_Data$GSVA_pORGpSUB_Primaries)
SortValue = BCC_RNA_Data$GSVA_pORGpSUB_Primaries[, "pSUB_Up_51"]
Number = floor(length(SelectedSampleIDs) / 4)
idx = order(SortValue, decreasing = TRUE)
TopIDs = SelectedSampleIDs[idx[1:Number]]
idx = order(SortValue, decreasing = FALSE)
BottomIDs = SelectedSampleIDs[idx[1:Number]]

TopQuarterIdx = which(colnames(M_ForGSEA) %in% TopIDs)
BottomQuarterIdx = which(colnames(M_ForGSEA) %in% BottomIDs)

SetA = TopQuarterIdx
SetAName = "Top4thPri_pSUB_Up_51"
SetB = BottomQuarterIdx
SetBName = "Bottom4thPri_pSUB_Up_51"

NonZeroIdx = which(rowSums(M_ForGSEA[, c(SetA, SetB)]) > 0)
FormatInputForGSEA(Data = M_ForGSEA[NonZeroIdx, c(SetA, SetB)],
                   SetA = SetA,
                   SetAName = SetAName,
                   SetB = SetB,
                   SetBName = SetBName,
                   OutputPath = PathForGSEA_Files,
                   FileNameForGSEA = "GSEA_pSUB_TopVsBottom_Pri")

### LiverVsLungCohort_Pri_Up/Dn ###

SelectedSampleIDs = BCC_RNA_Data$Metadata$Public_Specimen_ID[BCC_RNA_Data$Metadata$Tumor_Type == "Primary"]
SetA = which(BCC_RNA_Data$Metadata$Public_Specimen_ID %in% SelectedSampleIDs &
             BCC_RNA_Data$Metadata$Liver_Met_Present == 1)
SetAName = "LiverCohortPri"
SetB = which(BCC_RNA_Data$Metadata$Public_Specimen_ID %in% SelectedSampleIDs &
             BCC_RNA_Data$Metadata$Lung_Met_Present == 1 &
             is.na(BCC_RNA_Data$Metadata$Liver_Met_Present))
SetBName = "LungNotLiverCohortPri"

NonZeroIdx = which(rowSums(M_ForGSEA[, c(SetA, SetB)]) > 0)
FormatInputForGSEA(Data = M_ForGSEA[NonZeroIdx, c(SetA, SetB)],
                   SetA = SetA,
                   SetAName = SetAName,
                   SetB = SetB,
                   SetBName = SetBName,
                   OutputPath = PathForGSEA_Files,
                   FileNameForGSEA = "GSEA_LiverVsLungCohort_Pri")

### PurIST_TopVsBottom_Pri_Up/Dn ###

SelectedSampleIDs = BCC_RNA_Data$Metadata$Public_Specimen_ID[BCC_RNA_Data$Metadata$Tumor_Type == "Primary"]
SortValue = BCC_RNA_Data$Metadata$PurIST_Score[BCC_RNA_Data$Metadata$Tumor_Type == "Primary"]
Number = floor(length(SelectedSampleIDs) / 4)
idx = order(SortValue, decreasing = TRUE)
TopIDs = SelectedSampleIDs[idx[1:Number]]
idx = order(SortValue, decreasing = FALSE)
BottomIDs = SelectedSampleIDs[idx[1:Number]]

TopQuarterIdx = which(colnames(M_ForGSEA) %in% TopIDs)
BottomQuarterIdx = which(colnames(M_ForGSEA) %in% BottomIDs)

SetA = TopQuarterIdx
SetAName = "Top4thPri_PurIST"
SetB = BottomQuarterIdx
SetBName = "Bottom4thPri_PurIST"

NonZeroIdx = which(rowSums(M_ForGSEA[, c(SetA, SetB)]) > 0)
FormatInputForGSEA(Data = M_ForGSEA[NonZeroIdx, c(SetA, SetB)],
                   SetA = SetA,
                   SetAName = SetAName,
                   SetB = SetB,
                   SetBName = SetBName,
                   OutputPath = PathForGSEA_Files,
                   FileNameForGSEA = "GSEA_PurIST_TopVsBottom_Pri")

### pORG_TopVsBottom_Met_Up/Dn ###

SelectedSampleIDs = rownames(BCC_RNA_Data$GSVA_pORGpSUB_Mets)
SortValue = BCC_RNA_Data$GSVA_pORGpSUB_Mets[, "pORG_Up_55"]
Number = floor(length(SelectedSampleIDs) / 4)
idx = order(SortValue, decreasing = TRUE)
TopIDs = SelectedSampleIDs[idx[1:Number]]
idx = order(SortValue, decreasing = FALSE)
BottomIDs = SelectedSampleIDs[idx[1:Number]]

TopQuarterIdx = which(colnames(M_ForGSEA) %in% TopIDs)
BottomQuarterIdx = which(colnames(M_ForGSEA) %in% BottomIDs)

SetA = TopQuarterIdx
SetAName = "Top4thMet_pORG_Up_55"
SetB = BottomQuarterIdx
SetBName = "Bottom4thMet_All_pORG_Up_55"

NonZeroIdx = which(rowSums(M_ForGSEA[, c(SetA, SetB)]) > 0)
FormatInputForGSEA(Data = M_ForGSEA[NonZeroIdx, c(SetA, SetB)],
                   SetA = SetA,
                   SetAName = SetAName,
                   SetB = SetB,
                   SetBName = SetBName,
                   OutputPath = PathForGSEA_Files,
                   FileNameForGSEA = "GSEA_pORG_TopVsBottom_Met")

### pSUB_TopVsBottom_Met_Up/Dn ###

SelectedSampleIDs = rownames(BCC_RNA_Data$GSVA_pORGpSUB_Mets)
SortValue = BCC_RNA_Data$GSVA_pORGpSUB_Mets[, "pSUB_Up_51"]
Number = floor(length(SelectedSampleIDs) / 4)
idx = order(SortValue, decreasing = TRUE)
TopIDs = SelectedSampleIDs[idx[1:Number]]
idx = order(SortValue, decreasing = FALSE)
BottomIDs = SelectedSampleIDs[idx[1:Number]]

TopQuarterIdx = which(colnames(M_ForGSEA) %in% TopIDs)
BottomQuarterIdx = which(colnames(M_ForGSEA) %in% BottomIDs)

SetA = TopQuarterIdx
SetAName = "Top4thMet_pSUB_Up_51"
SetB = BottomQuarterIdx
SetBName = "Bottom4thMet_pSUB_Up_51"

NonZeroIdx = which(rowSums(M_ForGSEA[, c(SetA, SetB)]) > 0)
FormatInputForGSEA(Data = M_ForGSEA[NonZeroIdx, c(SetA, SetB)],
                   SetA = SetA,
                   SetAName = SetAName,
                   SetB = SetB,
                   SetBName = SetBName,
                   OutputPath = PathForGSEA_Files,
                   FileNameForGSEA = "GSEA_pSUB_TopVsBottom_Met")

### LiverVsLungCohort_Met_Up/Dn ###

SelectedSampleIDs = BCC_RNA_Data$Metadata$Public_Specimen_ID[BCC_RNA_Data$Metadata$Tumor_Type == "Met"]
SetA = which(BCC_RNA_Data$Metadata$Public_Specimen_ID %in% SelectedSampleIDs &
               BCC_RNA_Data$Metadata$Liver_Met_Present == 1)
SetAName = "LiverCohortMet"
SetB = which(BCC_RNA_Data$Metadata$Public_Specimen_ID %in% SelectedSampleIDs &
               BCC_RNA_Data$Metadata$Lung_Met_Present == 1 &
               is.na(BCC_RNA_Data$Metadata$Liver_Met_Present))
SetBName = "LungNotLiverCohortMet"

NonZeroIdx = which(rowSums(M_ForGSEA[, c(SetA, SetB)]) > 0)
FormatInputForGSEA(Data = M_ForGSEA[NonZeroIdx, c(SetA, SetB)],
                   SetA = SetA,
                   SetAName = SetAName,
                   SetB = SetB,
                   SetBName = SetBName,
                   OutputPath = PathForGSEA_Files,
                   FileNameForGSEA = "GSEA_LiverVsLungCohort_Met")

### PurIST_TopVsBottom_Met_Up/Dn ###

SelectedSampleIDs = BCC_RNA_Data$Metadata$Public_Specimen_ID[BCC_RNA_Data$Metadata$Tumor_Type == "Met"]
SortValue = BCC_RNA_Data$Metadata$PurIST_Score[BCC_RNA_Data$Metadata$Tumor_Type == "Met"]
Number = floor(length(SelectedSampleIDs) / 4)
idx = order(SortValue, decreasing = TRUE)
TopIDs = SelectedSampleIDs[idx[1:Number]]
idx = order(SortValue, decreasing = FALSE)
BottomIDs = SelectedSampleIDs[idx[1:Number]]

TopQuarterIdx = which(colnames(M_ForGSEA) %in% TopIDs)
BottomQuarterIdx = which(colnames(M_ForGSEA) %in% BottomIDs)

SetA = TopQuarterIdx
SetAName = "Top4thMet_PurIST"
SetB = BottomQuarterIdx
SetBName = "Bottom4thMet_PurIST"

NonZeroIdx = which(rowSums(M_ForGSEA[, c(SetA, SetB)]) > 0)
FormatInputForGSEA(Data = M_ForGSEA[NonZeroIdx, c(SetA, SetB)],
                   SetA = SetA,
                   SetAName = SetAName,
                   SetB = SetB,
                   SetBName = SetBName,
                   OutputPath = PathForGSEA_Files,
                   FileNameForGSEA = "GSEA_PurIST_TopVsBottom_Met")

#= END =========================================================================