#################################################################
# Supplemental Figure 11A-B: GSVA Heatmaps
#################################################################
# Load libraries

library(dplyr)
library(pheatmap)
set.seed(1)



#################################################################
# Load metadata

subtype <- readxl::read_xlsx("Supplemental_Dataset_1.xlsx", sheet = 2, na = "NA")
recurrence <- subtype

recurrence$group <- NA
recurrence$group[is.na(recurrence$Days_from_Resection_to_FU) & is.na(recurrence$Days_from_Resection_to_Recurrence)] = "No Resection"

# "No Documented Recurrence"
recurrence$group[   !is.na(recurrence$Days_from_Resection_to_FU) &
                        is.na(recurrence$Days_from_Resection_to_Recurrence) &
                        is.na(recurrence$Lung_Met_Present) &
                        is.na(recurrence$Liver_Met_Present)] = "No Documented Recurrence"

# "Liver Metastasis" - no resection included
recurrence$group[recurrence$Liver_Met_Present==1 & !(is.na(recurrence$Liver_Met_Present))] = "Liver Cohort"

# "Lung Metastasis" - no resection included
recurrence$group[(recurrence$Lung_Met_Present==1 & !(is.na(recurrence$Lung_Met_Present))) &
                     is.na(recurrence$Liver_Met_Present)] = "Lung Cohort"

# "Other Recurrence Site"
recurrence$group[ is.na(recurrence$Lung_Met_Present) &
                      is.na(recurrence$Liver_Met_Present) &
                      !is.na(recurrence$Days_from_Resection_to_Recurrence)] = "Other Recurrence Site"



#################################################################
# GSVA HEATMAP: scores from Primaries matrix 
#################################################################

####################################
# Load GSVA results

gsva_results <- as.data.frame(readxl::read_xlsx("Supplemental_Dataset_5.xlsx", sheet = "GSVA_Hallmarks_Primaries"))
rownames(gsva_results) <- gsva_results$Public_Specimen_ID
gsva_results <- gsva_results[,-1]
gsva_results <- as.data.frame(t(gsva_results))

# Convert to numeric
for(col in colnames(gsva_results)){
    gsva_results[col] <- as.numeric(gsva_results[[col]])
    rm(col)
}

# Remove no resection samples
gsva_results <- gsva_results[colnames(gsva_results) %in% recurrence$Public_Specimen_ID[recurrence$group!="No Resection"]]



####################################
# Prepare annotation for GSVA heatmap

gsva_annot <- data.frame(group = colnames(gsva_results))
rownames(gsva_annot) <- gsva_annot$group

# Recurrence group
gsva_annot$Group <- recurrence$group[match(gsva_annot$group,recurrence$Public_Specimen_ID)]

# pORG
gsva_annot$pORG <- recurrence$pORG_Up_55_Primaries[match(gsva_annot$group,recurrence$Public_Specimen_ID)]

# Subtype
gsva_annot$Subtype <- recurrence$PurIST_Subtype[match(gsva_annot$group,recurrence$Public_Specimen_ID)]
gsva_annot$Subtype[gsva_annot$Subtype=="classical"] = "Classical"
gsva_annot$Subtype[gsva_annot$Subtype=="basal-like"] = "Basal-like"



####################################
# Significant GSEA results for GSVA heatmap

hallmarks_to_use <- c("HALLMARK_PROTEIN_SECRETION","HALLMARK_MTORC1_SIGNALING","HALLMARK_GLYCOLYSIS","HALLMARK_MITOTIC_SPINDLE","HALLMARK_G2M_CHECKPOINT","HALLMARK_MYC_TARGETS_V1","HALLMARK_E2F_TARGETS","HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_DNA_REPAIR","HALLMARK_OXIDATIVE_PHOSPHORYLATION","HALLMARK_ANDROGEN_RESPONSE","HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_MYOGENESIS")

hallmarks_to_use <- gsub("HALLMARK_","",hallmarks_to_use)
hallmarks_to_use <- gsub("_"," ",hallmarks_to_use)
hallmarks_to_use <- stringr::str_to_title(hallmarks_to_use)
hallmarks_to_use <- gsub("^E2f Targets$","E2F Targets",hallmarks_to_use)
hallmarks_to_use <- gsub("^Mtorc1 Signaling$","MTORC1 Signaling",hallmarks_to_use)
hallmarks_to_use <- gsub("^G2m Checkpoint$","G2M Checkpoint",hallmarks_to_use)
hallmarks_to_use <- gsub("^Myc Targets V1$","MYC Targets V1",hallmarks_to_use)
hallmarks_to_use <- gsub("^Myc Targets V2$","MYC Targets V2",hallmarks_to_use)
hallmarks_to_use <- gsub("^E2f Targets$","E2F Targets",hallmarks_to_use)
hallmarks_to_use <- gsub("^Dna Repair$","DNA Repair",hallmarks_to_use)



####################################
# GSVA Hallmark Heatmap - reordered by pORG and by correlation with pORG

gsva_annot_reordered <- gsva_annot
gsva_annot_reordered <- gsva_annot_reordered[order(gsva_annot_reordered$pORG, decreasing = T),]

gsva_results_reordered <- gsva_results[,gsva_annot_reordered$group]


# Color palette for heatmap
heatmap_colors = list(
    "Group" = c('Liver Cohort' = "#00A403",
                'Lung Cohort' = "#FF0000",
                'Other Recurrence Site' = "#79FF47",
                'No Documented Recurrence' = "#FFFB86"),
    "Subtype" = c('Classical' = "#CC79A7",
                  'Basal-like' = 'black'),
    "pORG" = colorRampPalette(rev(c("#D55E00","white","#0072B2")))(100)
)


rownames(gsva_results_reordered) <- gsub("HALLMARK_","",rownames(gsva_results_reordered))
rownames(gsva_results_reordered) <- gsub("_"," ",rownames(gsva_results_reordered))
rownames(gsva_results_reordered) <- stringr::str_to_title(rownames(gsva_results_reordered))
rownames(gsva_results_reordered) <- gsub("^E2f Targets$","E2F Targets",rownames(gsva_results_reordered))
rownames(gsva_results_reordered) <- gsub("^Mtorc1 Signaling$","MTORC1 Signaling",rownames(gsva_results_reordered))
rownames(gsva_results_reordered) <- gsub("^G2m Checkpoint$","G2M Checkpoint",rownames(gsva_results_reordered))
rownames(gsva_results_reordered) <- gsub("^Myc Targets V1$","MYC Targets V1",rownames(gsva_results_reordered))
rownames(gsva_results_reordered) <- gsub("^Myc Targets V2$","MYC Targets V2",rownames(gsva_results_reordered))
rownames(gsva_results_reordered) <- gsub("^E2f Targets$","E2F Targets",rownames(gsva_results_reordered))
rownames(gsva_results_reordered) <- gsub("^Dna Repair$","DNA Repair",rownames(gsva_results_reordered))



p_primaries_clustered <- pheatmap(gsva_results_reordered[rownames(gsva_results_reordered) %in% hallmarks_to_use,],
                                  cluster_rows = TRUE,
                                  cluster_cols = FALSE,
                                  show_rownames = TRUE,
                                  show_colnames = FALSE,
                                  annotation_colors = heatmap_colors,
                                  annotation_col = gsva_annot_reordered[-1],
                                  annotation_names_row = F,
                                  colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu")))(100),
                                  angle_col = 45,
                                  fontsize = 8.2,
                                  fontsize_col = 10,
                                  fontsize_row = 9,
                                  fontsize_number = 10,
                                  main = "Primary Tumor Heatmap (Clustered)\nGSVA scores from Primaries Only")



############################################################################################################
############################################################################################################
# PRINT PLOTS ##############################################################################################
############################################################################################################
############################################################################################################

pdf(width = 7, height = 3, file = "GSVAHeatmap_Clustered_PrimariesMatrix_pORG.pdf")
p_primaries_clustered
dev.off()








############################################################################################################
############################################################################################################
# pSUB PRIMARIES MATRIX ####################################################################################
############################################################################################################
############################################################################################################

####################################
# Prepare annotation for GSVA heatmap

gsva_annot <- data.frame(group = colnames(gsva_results))
rownames(gsva_annot) <- gsva_annot$group

# Recurrence group
gsva_annot$Group <- recurrence$group[match(gsva_annot$group,recurrence$Public_Specimen_ID)]

# pORG
gsva_annot$pSUB <- recurrence$pSUB_Up_51_Primaries[match(gsva_annot$group,recurrence$Public_Specimen_ID)]

# Subtype
gsva_annot$Subtype <- recurrence$PurIST_Subtype[match(gsva_annot$group,recurrence$Public_Specimen_ID)]
gsva_annot$Subtype[gsva_annot$Subtype=="classical"] = "Classical"
gsva_annot$Subtype[gsva_annot$Subtype=="basal-like"] = "Basal-like"



####################################
# GSVA Hallmark Heatmap - reordered by pSUB and by correlation with hallmarks

gsva_annot_reordered <- gsva_annot
gsva_annot_reordered <- gsva_annot_reordered[order(gsva_annot_reordered$pSUB, decreasing = T),]

gsva_results_reordered <- gsva_results[,gsva_annot_reordered$group]


heatmap_colors = list(
    "Group" = c('Liver Cohort' = "#00A403",
                'Lung Cohort' = "#FF0000",
                'Other Recurrence Site' = "#79FF47",
                'No Documented Recurrence' = "#FFFB86"),
    "Subtype" = c('Classical' = "#CC79A7",
                  'Basal-like' = 'black'),
    "pSUB" = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "PiYG")))(100)
)


rownames(gsva_results_reordered) <- gsub("HALLMARK_","",rownames(gsva_results_reordered))
rownames(gsva_results_reordered) <- gsub("_"," ",rownames(gsva_results_reordered))
rownames(gsva_results_reordered) <- stringr::str_to_title(rownames(gsva_results_reordered))


hallmarks_psub <- c("HALLMARK_GLYCOLYSIS","HALLMARK_HYPOXIA","HALLMARK_APICAL_JUNCTION","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
hallmarks_psub <- gsub("HALLMARK_","",hallmarks_psub)
hallmarks_psub <- gsub("_"," ",hallmarks_psub)
hallmarks_psub <- stringr::str_to_title(hallmarks_psub)


p_psub_primaries_clustered <- pheatmap(gsva_results_reordered[rownames(gsva_results_reordered) %in% hallmarks_psub,],
                                       cluster_rows = TRUE,
                                       cluster_cols = FALSE,
                                       show_rownames = TRUE,
                                       show_colnames = FALSE,
                                       annotation_colors = heatmap_colors,
                                       annotation_col = gsva_annot_reordered[c("Group","pSUB","Subtype")],
                                       annotation_names_row = F,
                                       colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu")))(100),
                                       angle_col = 45,
                                       fontsize = 8.2,
                                       fontsize_col = 10,
                                       fontsize_row = 9,
                                       fontsize_number = 10,
                                       main = "Primary Tumor Heatmap (Clustered)\nGSVA scores from Primaries Only")


############################################################################################################
############################################################################################################
# PRINT PLOTS ##############################################################################################
############################################################################################################
############################################################################################################

pdf(width = 7.5, height = 3.02, file = "GSVAHeatmap_Clustered_PrimariesMatrix_pSUB.pdf")
p_psub_primaries_clustered
dev.off()
