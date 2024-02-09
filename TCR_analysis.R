library(immunarch)
data(immdata) 

#load the file
file_path <-  "raw_TCR_data"
immdata_test <- repLoad(file_path)

# Compute statistics and visualise them
# Chao1 diversity measure
div_chao <- repDiversity(immdata_test$data, "chao1")

# Hill numbers
div_hill <- repDiversity(immdata_test$data, "hill")

# D50
div_d50 <- repDiversity(immdata_test$data, "d50")

# Ecological diversity measure
div_div <- repDiversity(immdata_test$data, "div")

# gini.simp
div_gini_simp <- repDiversity(immdata_test$data, "gini.simp")

#'pORG_0.2_Primary_quartiles_Site','pORG_0.2_Met_quartiles_Site',
#'Cohort_Primary_Site','Cohort_Met_Site'

s_meta <- "Cohort_Primary_Site"
p1 <- vis(div_chao, .by = s_meta, .meta = immdata_test$meta)
p2 <- vis(div_hill, .by = s_meta, .meta = immdata_test$meta)
p3 <- vis(div_d50, .by = s_meta, .meta = immdata_test$meta)
p4 <- vis(div_div, .by = s_meta, .meta = immdata_test$meta)

#save
write.csv(div_div,"results_TCR_true_diversity.csv", row.names = TRUE)
write.csv(div_chao,"results_TCR_chao_diversity.csv", row.names = TRUE)
write.csv(div_hill,"results_TCR_hill_diversity.csv", row.names = TRUE)
write.csv(div_d50,"results_TCR_d50_diversity.csv", row.names = TRUE)
write.csv(div_gini_simp,"results_TCR_ginisimp_diversity.csv", row.names = TRUE)

# repertoire overlap
imm_ov1 <- repOverlap(immdata_test$data, .method = "public", .verbose = F)
write.csv(imm_ov1,"results_TCR_public_overlap.csv", row.names = TRUE)
#imm_ov2 <- repOverlap(immdata_test$data, .method = "morisita", .verbose = F)
#write.csv(imm_ov2,"results_TCR_morisita_overlap.csv", row.names = TRUE)
imm_ov3 <- repOverlap(immdata_test$data, "jaccard", .verbose = FALSE)
write.csv(imm_ov3,"results_TCR_jaccard_overlap.csv", row.names = TRUE)
#vis(imm_ov1, "heatmap2")




