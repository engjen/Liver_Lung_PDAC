#################################################################
# Supplemental Figure 12: GO Network Plots
#################################################################
# Load libraries

library(readxl)
library(org.Hs.eg.db)
library(rtracklayer)
library(clusterProfiler)
library(enrichplot)
set.seed(1)



#################################################################
# Load GENCODE annotation

# Load GENCODE annotation (obtained from https://www.gencodegenes.org/human/release_24.html)
annotation <- import("gencode.v24.annotation.gtf.gz")

# Convert to data frame. To be used for converting regulon name to Ensembl ID
annotation <- as.data.frame(annotation@elementMetadata@listData)
annotation <- annotation[c("gene_id","gene_type","gene_name")]  # Keep only the info of interest
annotation <- unique(annotation)



#################################################################
# Regulon GO analysis for high/low pORG

# Read in regulon data
# The original regulon data can be found in Supplemental Dataset 6
# Here we focus on regulons that have the greatest abs. mean difference between high/low pORG
regulons <- read_xlsx("ViperANOVA-pORGhigh.xlsx")

# Order regulons from high-to-low MeanDiff_pORG_Up_55
regulons <- regulons[order(regulons$MeanDiff_pORG_Up_55,decreasing = T),]

# Add Ensembl ID for each regulon (for GO enrichment)
regulons$ensembl <- annotation$gene_id[match(regulons$regulons,annotation$gene_name)]
regulons$ensembl <- gsub("\\..*","",regulons$ensembl)

# Verify that the top 300 regulons are significant
identical(regulons$regulons[1:300],regulons$regulons[regulons$TopBottomQuarter_pORG_Up_55_qVal<=0.05][1:300])

# Run GO enrichment on top 300 regulons
goenrich <- enrichGO(regulons$ensembl[1:300],
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     OrgDb = org.Hs.eg.db,
                     universe = regulons$ensembl,
                     keyType = 'ENSEMBL',
                     ont = 'BP')

# Calculate term similarity
pairwise <- pairwise_termsim(goenrich)

# Plot GO network
p <- emapplot(pairwise,
              cex_label_category = 0.35,
              repel = TRUE,
              group_category = TRUE, 
              showCategory = 150,
              force = 1.5)

# Save plot to PDF
pdf(width = 16, height = 16, file = "GONetwork_Regulons_Top300MeanDiffpORGUp55_qval0.05.pdf")
p
dev.off()


# Extract GO enrichment results as table
enrichmentResults <- as.data.frame(goenrich)

# Replace Ensembl ID with regulon name
enrichmentResults$Regulons <- enrichmentResults$geneID
for(regulon in regulons$regulons){
    enrichmentResults$Regulons <- gsub(regulons$ensembl[regulons$regulons==regulon],regulon,enrichmentResults$Regulon)
}

# Reorder columns and rename to something more informative
enrichmentResults <- enrichmentResults[c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Regulons","Count")]
colnames(enrichmentResults) <- gsub("^ID$","GO_ID",colnames(enrichmentResults))
colnames(enrichmentResults) <- gsub("^Description$","GO_Term",colnames(enrichmentResults))

# Write table
write.table(enrichmentResults, file = "GOEnrichment_Regulons_Top300MeanDiffpORGUp55_qval0.05.txt", row.names = F, col.names = T, sep = '\t', quote = F)






#################################################################
# Regulon GO analysis for liver/lung cohort

# Read in regulon data (again)
# Re: original regulon data can be found in Supplemental Dataset 6
# Here we focus on regulons that have the greatest abs. mean difference between liver/lung cohort
regulons <- read_xlsx("ViperANOVA-pORGhigh.xlsx")

# Order regulons from high-to-low MeanDiff_LiverVsLungNotLiver
regulons <- regulons[order(regulons$MeanDiff_LiverVsLungNotLiver,decreasing = T),]

# Add Ensembl ID for each regulon (for GO enrichment)
regulons$ensembl <- annotation$gene_id[match(regulons$regulons,annotation$gene_name)]
regulons$ensembl <- gsub("\\..*","",regulons$ensembl)

# Check if top 300 are significant
identical(
    regulons$regulons[1:300],
    regulons$regulons[regulons$LiverVsLungNotLiver_qVal<=0.05][1:300]
)

# Top 300 are not significant. Select only the top 300 significant regulons
top <- regulons$ensembl[regulons$LiverVsLungNotLiver_qVal<=0.05][1:300]

# Run GO enrichment on top 300 regulons
goenrich <- enrichGO(top,
                     OrgDb = org.Hs.eg.db,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     universe = regulons$ensembl,
                     keyType = 'ENSEMBL',
                     ont = 'BP')

# Calculate term similarity
pairwise <- pairwise_termsim(goenrich)

# Plot GO network
p <- emapplot(pairwise,
              cex_label_category = 0.35,
              repel = TRUE,
              group_category = TRUE, 
              showCategory = 150,
              force = 1.5)

# Save plot to PDF
pdf(width = 16, height = 16, file = "GONetwork_Regulons_Top300MeanDiffLiverVsLungNotLiver_qval0.05.pdf")
p
dev.off()


# Extract GO enrichment results as table
enrichmentResults <- as.data.frame(goenrich)

# Replace EnsemblID with regulon name
enrichmentResults$Regulons <- enrichmentResults$geneID
for(regulon in regulons$regulons){
    enrichmentResults$Regulons <- gsub(regulons$ensembl[regulons$regulons==regulon],regulon,enrichmentResults$Regulon)
}

# Reorder columns and rename to something more informative
enrichmentResults <- enrichmentResults[c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","Regulons","Count")]
colnames(enrichmentResults) <- gsub("^ID$","GO_ID",colnames(enrichmentResults))
colnames(enrichmentResults) <- gsub("^Description$","GO_Term",colnames(enrichmentResults))

# Write table
write.table(enrichmentResults, file = "GOEnrichment_Regulons_Top300MeanDiffLiverVsLungNotLiver_qval0.05.txt", row.names = F, col.names = T, sep = '\t', quote = F)

