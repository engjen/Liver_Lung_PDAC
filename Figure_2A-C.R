################################################################## Figure 2A-C################################################################## Load librarieslibrary(ggplot2)library(ggsignif)set.seed(1)

################################################################## Primary Tumors - pORG
################################################################## Supplemental Dataset 1metadata <- readxl::read_xlsx("Supplemental_Dataset_1.xlsx", sheet = 1)metadata <- metadata[!(is.na(metadata$pORG_Primary)) | !(is.na(metadata$pORG_Primary_T2)),]########################################## Collect pORG & pSUB data for -T and -T2 samples

# -T samplespORG_primary <- data.frame(cohort = metadata$Cohort,                           subtype = metadata$PurIST_Subtype,                           pORG = metadata$pORG_Primary,                           pSUB = metadata$pSUB_Primary)# -T2 samplespORG_primary_T2 <- data.frame(cohort = metadata$Cohort,                              subtype = metadata$PurIST_Subtype,                              pORG = metadata$pORG_Primary_T2,                              pSUB = metadata$pSUB_Primary_T2)pORG_primary_T2 <- pORG_primary_T2[!(is.na(pORG_primary_T2$pORG)),]# Basal-like subtype -T and -T2 samplespORG_subtype_basal <- rbind(pORG_primary[pORG_primary$subtype=="basal-like",],                            pORG_primary_T2[pORG_primary_T2$subtype=="basal-like",])

# Classical subtype -T and -T2 samplespORG_subtype_classical <- rbind(pORG_primary[pORG_primary$subtype=="classical",],                                pORG_primary_T2[pORG_primary_T2$subtype=="classical",])# Collect all data togetherpORG_subtype <- rbind(pORG_subtype_basal, pORG_subtype_classical)pORG_subtype$cohort <- pORG_subtype$subtypepORG_all <- rbind(pORG_primary, pORG_primary_T2, pORG_subtype)pORG_all <- pORG_all[!(is.na(pORG_all$cohort)),]



########################################## Plot data

# Get one-way ANOVA (adjusted by BH method) p-values for plot
p.adjust(    c(        summary(aov(pORG ~ cohort, data = pORG_all[pORG_all$cohort %in% c("Liver","Lung"),] ))[[1]][["Pr(>F)"]][[1]],        summary(aov(pORG ~ cohort, data = pORG_all[pORG_all$cohort %in% c("classical","basal-like"),] ))[[1]][["Pr(>F)"]][[1]]    ),    method = "BH")


plot_data <- pORG_allplot_data$cohort <- factor(plot_data$cohort, levels = c("Liver","Lung","basal-like","classical"))pdf(width = 5.1, height = 3.125, file = "pORG_Primary.pdf")ggplot(plot_data, aes(x = cohort, y = pORG)) + geom_violin(color = "gray40") +     geom_jitter(width = 0.2, aes(color = cohort), size = 3, alpha = 0.8, show.legend = F, stroke = NA) +     scale_color_manual(values = c("Liver" = "#d55e00",                                  "Lung" = "#0072b2",                                  "basal-like" = "#000000",                                  "classical" = "#cc79a7")) +    stat_summary(fun=mean, geom="crossbar" , linewidth = 0.6) +     ylim(-0.7,1) +    theme(panel.grid.major = element_blank(),          panel.grid.minor = element_blank(),          panel.background = element_blank(),          axis.line = element_line(color = "black"),          axis.ticks = element_line(color = "black"),          axis.title.x = element_text(color = "black",size = 13),          axis.title.y = element_text(color = "black", size = 17, face = "bold"),          axis.text.x = element_text(color = "black", size = 16),          axis.text.y = element_text(color = "black", size = 17),          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +    xlab("Cohort N = 76 | PurIST Subtype N = 218") +    ylab("pORG Primary") +    geom_signif(comparisons = list(c("Liver", "Lung"),                                    c("basal-like", "classical")),                annotation = c("P=1.6e-08","P=0.38"), tip_length = 0.015, size = 1, y_position = 0.75, textsize = 5.5, vjust = -0.25)dev.off()


################################################################## Primary Tumors - pSUB
################################################################## Get one-way ANOVA (adjusted by BH method) p-values for plotp.adjust(    c(        summary(aov(pSUB ~ cohort, data = pORG_all[pORG_all$cohort %in% c("Liver","Lung"),] ))[[1]][["Pr(>F)"]][[1]],        summary(aov(pSUB ~ cohort, data = pORG_all[pORG_all$cohort %in% c("classical","basal-like"),] ))[[1]][["Pr(>F)"]][[1]]    ),    method = "BH")pdf(width = 5.1, height = 3.125, file = "pSUB_Primary.pdf")ggplot(plot_data, aes(x = cohort, y = pSUB)) + geom_violin(color = "gray40") +     geom_jitter(width = 0.2, aes(color = cohort), size = 3, alpha = 0.8, show.legend = F, stroke = NA) +     scale_color_manual(values = c("Liver" = "#d55e00",                                  "Lung" = "#0072b2",                                  "basal-like" = "#000000",                                  "classical" = "#cc79a7")) +    stat_summary(fun=mean, geom="crossbar" , linewidth = 0.6) +     ylim(-0.7,1) +    theme(panel.grid.major = element_blank(),          panel.grid.minor = element_blank(),          panel.background = element_blank(),          axis.line = element_line(color = "black"),          axis.ticks = element_line(color = "black"),          axis.title.x = element_text(color = "black",size = 13),          axis.title.y = element_text(color = "black", size = 17, face = "bold"),          axis.text.x = element_text(color = "black", size = 16),          axis.text.y = element_text(color = "black", size = 17),          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +    xlab("Cohort N = 76 | PurIST Subtype N = 218") +    ylab("pSUB Primary") +    geom_signif(comparisons = list(c("Liver", "Lung"),                                    c("basal-like", "classical")),                annotation = c("P=0.22","P=7.1e-27"), tip_length = 0.015, size = 1, y_position = 0.82, textsize = 5.5, vjust = -0.20)
dev.off()



################################################################## Metastases - pORG
#################################################################
# Supplemental Dataset 1metadata <- readxl::read_xlsx("Supplemental_Dataset_1.xlsx", sheet = 1)metadata <- metadata[!(is.na(metadata$pORG_Met)),]#########################################
# Collect pORG & pSUB data for metastases
pORG_met <- data.frame(cohort = metadata$Cohort,                       subtype = metadata$PurIST_Subtype_Met,                       pORG = metadata$pORG_Met,                       pSUB = metadata$pSUB_Met)pORG_subtype_met <- pORG_metpORG_subtype_met$cohort <- pORG_subtype_met$subtypepORG_all_met <- rbind(pORG_met, pORG_subtype_met)pORG_all_met <- pORG_all_met[!(is.na(pORG_all_met$cohort)),]



########################################## Plot data

# Get one-way ANOVA (adjusted by BH method) p-values for plotp.adjust(    c(        summary(aov(pORG ~ cohort, data = pORG_all_met[pORG_all_met$cohort %in% c("Liver","Lung"),] ))[[1]][["Pr(>F)"]][[1]],        summary(aov(pORG ~ cohort, data = pORG_all_met[pORG_all_met$cohort %in% c("classical","basal-like"),] ))[[1]][["Pr(>F)"]][[1]]    ),    method = "BH")


plot_data <- pORG_all_metplot_data$cohort <- factor(plot_data$cohort, levels = c("Liver","Lung","basal-like","classical"))

pdf(width = 5.1, height = 3.125, file = "pORG_Met.pdf")ggplot(plot_data, aes(x = cohort, y = pORG)) + geom_violin(color = "gray40") +     geom_jitter(width = 0.2, aes(color = cohort), size = 3, alpha = 0.8, show.legend = F, stroke = NA) +     scale_color_manual(values = c("Liver" = "#d55e00",                                  "Lung" = "#0072b2",                                  "basal-like" = "#000000",                                  "classical" = "#cc79a7")) +    stat_summary(fun=mean, geom="crossbar" , linewidth = 0.6) +     ylim(-0.7,1) +    theme(panel.grid.major = element_blank(),          panel.grid.minor = element_blank(),          panel.background = element_blank(),          axis.line = element_line(color = "black"),          axis.ticks = element_line(color = "black"),          axis.title.x = element_text(color = "black",size = 13),          axis.title.y = element_text(color = "black", size = 17, face = "bold"),          axis.text.x = element_text(color = "black", size = 16),          axis.text.y = element_text(color = "black", size = 17),          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +    xlab("Cohort N = 37 | PurIST Subtype N = 71") +    ylab("pORG Met") +    geom_signif(comparisons = list(c("Liver", "Lung"),                                    c("basal-like", "classical")),                annotation = c("P=0.91","P=0.39"), tip_length = 0.015, size = 1, y_position = 0.70, textsize = 5.5, vjust = -0.25)
dev.off()



################################################################## Metastases - pSUB
################################################################## Get one-way ANOVA (adjusted by BH method) p-values for plotp.adjust(    c(        summary(aov(pSUB ~ cohort, data = pORG_all_met[pORG_all_met$cohort %in% c("Liver","Lung"),] ))[[1]][["Pr(>F)"]][[1]],        summary(aov(pSUB ~ cohort, data = pORG_all_met[pORG_all_met$cohort %in% c("classical","basal-like"),] ))[[1]][["Pr(>F)"]][[1]]    ),    method = "BH")plot_data$cohort <- factor(plot_data$cohort, levels = c("Liver","Lung","basal-like","classical"))



pdf(width = 5.1, height = 3.125, file = "pSUB_Met.pdf")ggplot(plot_data, aes(x = cohort, y = pSUB)) + geom_violin(color = "gray40") +     geom_jitter(width = 0.2, aes(color = cohort), size = 3, alpha = 0.8, show.legend = F, stroke = NA) +     scale_color_manual(values = c("Liver" = "#d55e00",                                  "Lung" = "#0072b2",                                  "basal-like" = "#000000",                                  "classical" = "#cc79a7")) +    stat_summary(fun=mean, geom="crossbar" , linewidth = 0.6) +     ylim(-0.7,1) +    theme(panel.grid.major = element_blank(),          panel.grid.minor = element_blank(),          panel.background = element_blank(),          axis.line = element_line(color = "black"),          axis.ticks = element_line(color = "black"),          axis.title.x = element_text(color = "black",size = 13),          axis.title.y = element_text(color = "black", size = 17, face = "bold"),          axis.text.x = element_text(color = "black", size = 16),          axis.text.y = element_text(color = "black", size = 17),          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +    xlab("Cohort N = 37 | PurIST Subtype N = 71") +    ylab("pSUB Met") +    geom_signif(comparisons = list(c("Liver", "Lung"),                                    c("basal-like", "classical")),                annotation = c("P=0.0013","P=1.1e-08"), tip_length = 0.015, size = 1, y_position = 0.77, textsize = 5.5, vjust = -0.25)
dev.off()



################################################################## Primary vs Met - pORG
#################################################################

# Supplemental Dataset 1metadata <- readxl::read_xlsx("Supplemental_Dataset_1.xlsx", sheet = 1)metadata <- metadata[!(is.na(metadata$pORG_allPrimary)) | !(is.na(metadata$pORG_All_T2)) | !(is.na(metadata$pORG_allMet)),]



#########################################
# Collect pORG & pSUB data for primaries & metspORG_final_primary <- data.frame(Tumor = "Primary",                                 pORG = c(metadata$pORG_allPrimary[!(is.na(metadata$pORG_allPrimary))],                                          metadata$pORG_All_T2[!(is.na(metadata$pORG_All_T2))]),                                 pSUB = c(metadata$pSUB_allPrimary[!(is.na(metadata$pSUB_allPrimary))],                                          metadata$pSUB_All_T2[!(is.na(metadata$pSUB_All_T2))]))pORG_final_met <- data.frame(Tumor = "Met",                             pORG = metadata$pORG_allMet[!(is.na(metadata$pORG_allMet))],                             pSUB = metadata$pSUB_allMet[!(is.na(metadata$pSUB_allMet))])pORG_final <- rbind(pORG_final_primary,pORG_final_met)



#########################################
# Plot data

# Get one-way ANOVA p-valuesummary(aov(pORG ~ Tumor, data = pORG_final))[[1]][["Pr(>F)"]][[1]]plot_data <- pORG_finalplot_data$Tumor <- factor(plot_data$Tumor, levels = c("Primary","Met"))



pdf(width = 3.825, height = 3.7, file = "pORG_Primary_vs_Met.pdf")ggplot(plot_data, aes(x = Tumor, y = pORG)) + geom_violin(color = "gray40") +     geom_jitter(width = 0.2, aes(color = Tumor), size = 3, alpha = 0.8, show.legend = F, stroke = NA) +     scale_color_manual(values = c("Primary" = "gray60", "Met" = "black")) +    stat_summary(fun=mean, geom="crossbar" , linewidth = 0.6) +     ylim(-0.75,0.95) +    theme(panel.grid.major = element_blank(),          panel.grid.minor = element_blank(),          panel.background = element_blank(),          axis.line = element_line(color = "black"),          axis.ticks = element_line(color = "black"),          axis.title.x = element_text(color = "black",size = 14),          axis.title.y = element_text(color = "black", size = 17, face = "bold"),          axis.text.x = element_text(color = "black", size = 17),          axis.text.y = element_text(color = "black", size = 17),          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +    xlab("Specimen Type N = 289") +    ylab("pORG All Tumors") +    geom_signif(comparisons = list(c("Primary", "Met")),                annotation = c("P = 0.91"), tip_length = 0.015, size = 1, y_position = 0.70, textsize = 5.5, vjust = -0.25)dev.off()



################################################################## Primary vs Met - pSUB
#################################################################


# Get one-way ANOVA p-valuesummary(aov(pSUB ~ Tumor, data = pORG_final))[[1]][["Pr(>F)"]][[1]]plot_data <- pORG_finalplot_data$Tumor <- factor(plot_data$Tumor, levels = c("Primary","Met"))



pdf(width = 3.825, height = 3.7, file = "pSUB_Primary_vs_Met.pdf")ggplot(plot_data, aes(x = Tumor, y = pSUB)) + geom_violin(color = "gray40") +     geom_jitter(width = 0.2, aes(color = Tumor), size = 3, alpha = 0.8, show.legend = F, stroke = NA) +     scale_color_manual(values = c("Primary" = "gray60", "Met" = "black")) +    stat_summary(fun=mean, geom="crossbar" , linewidth = 0.6) +     ylim(-0.75,1) +    theme(panel.grid.major = element_blank(),          panel.grid.minor = element_blank(),          panel.background = element_blank(),          axis.line = element_line(color = "black"),          axis.ticks = element_line(color = "black"),          axis.title.x = element_text(color = "black",size = 14),          axis.title.y = element_text(color = "black", size = 17, face = "bold"),          axis.text.x = element_text(color = "black", size = 17),          axis.text.y = element_text(color = "black", size = 17),          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +    xlab("Specimen Type N = 289") +    ylab("pSUB All Tumors") +    geom_signif(comparisons = list(c("Primary", "Met")),                annotation = c("P = 0.39"), tip_length = 0.015, size = 1, y_position = 0.82, textsize = 5.5, vjust = -0.25)dev.off()