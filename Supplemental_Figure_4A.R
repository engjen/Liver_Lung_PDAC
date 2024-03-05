#################################################################
# Supplemental figure 4A
#################################################################
# Load libraries

library(ggplot2)
library(ggsignif)
set.seed(1)



#################################################################
# Primary Tumors - PurIST
#################################################################

# Supplemental Dataset 1
metadata <- readxl::read_xlsx("Supplemental_Dataset_1.xlsx", sheet = 1)

# Remove rows without PurIST scores
metadata <- metadata[!(is.na(metadata$PurIST_Primary)) | !(is.na(metadata$PurIST_Primary_T2)),]



#########################################
# Collect purIST data for -T and -T2 samples

# -T samples
pORG_primary <- data.frame(cohort = metadata$Cohort,
                           subtype = metadata$PurIST_Subtype,
                           pORG = metadata$pORG_Primary,
                           PurIST = metadata$PurIST_Primary)

# -T2 samples
pORG_primary_T2 <- data.frame(cohort = metadata$Cohort,
                              subtype = metadata$PurIST_Subtype,
                              pORG = metadata$pORG_Primary_T2,
                              PurIST = metadata$PurIST_Primary_T2)
pORG_primary_T2 <- pORG_primary_T2[!(is.na(pORG_primary_T2$PurIST)),]


# Basal-like subtype -T and -T2 samples
pORG_subtype_basal <- rbind(pORG_primary[pORG_primary$subtype=="basal-like",],
                            pORG_primary_T2[pORG_primary_T2$subtype=="basal-like",])

# Classical subtype -T and -T2 samples
pORG_subtype_classical <- rbind(pORG_primary[pORG_primary$subtype=="classical",],
                                pORG_primary_T2[pORG_primary_T2$subtype=="classical",])

# Collect all data together
pORG_subtype <- rbind(pORG_subtype_basal, pORG_subtype_classical)
pORG_subtype$cohort <- pORG_subtype$subtype
pORG_all <- rbind(pORG_primary, pORG_primary_T2, pORG_subtype)
pORG_all <- pORG_all[!(is.na(pORG_all$cohort)),]



#########################################
# Plot data

# Get one-way ANOVA (adjusted by BH method) p-values for plot
p.adjust(
    c(
        summary(aov(PurIST ~ cohort, data = pORG_all[pORG_all$cohort %in% c("Liver","Lung"),] ))[[1]][["Pr(>F)"]][[1]],
        summary(aov(PurIST ~ cohort, data = pORG_all[pORG_all$cohort %in% c("classical","basal-like"),] ))[[1]][["Pr(>F)"]][[1]]
    ),
    method = "BH")



plot_data <- pORG_all
plot_data$cohort <- factor(plot_data$cohort, levels = c("Liver","Lung","basal-like","classical"))



# Plot data
pdf(width = (5.1*1.25), height = (3.125*1.25), file = "PurIST_Primary.pdf")
ggplot(plot_data, aes(x = cohort, y = PurIST)) + geom_violin(color = "gray40") + 
    geom_jitter(width = 0.2, aes(color = cohort), size = 3, alpha = 0.8, show.legend = F, stroke = NA) + 
    scale_color_manual(values = c("Liver" = "#d55e00",
                                  "Lung" = "#0072b2",
                                  "basal-like" = "#000000",
                                  "classical" = "#cc79a7")) +
    stat_summary(fun=mean, geom="crossbar" , linewidth = 0.6) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.title.x = element_text(color = "black",size = 13),
          axis.title.y = element_text(color = "black", size = 17),
          axis.text.x = element_text(color = "black", size = 16),
          axis.text.y = element_text(color = "black", size = 17),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    xlab("Cohort N = 76 | PurIST Subtype N = 218") +
    scale_y_continuous(name = "PurIST Score (Primaries)",
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       limits = c(0, 1.15)) +
    geom_signif(comparisons = list(c("Liver", "Lung"), 
                                   c("basal-like", "classical")),
                annotation = c("P=0.17","P=1.0e-114"), tip_length = 0.015, size = 1, textsize = 5.5, vjust = 0, margin_top = 0.1)
dev.off()



#################################################################
# Metastases - PurIST
#################################################################

# Supplemental Dataset 1
metadata <- readxl::read_xlsx("Supplemental_Dataset_1.xlsx", sheet = 1)

metadata <- metadata[!(is.na(metadata$PurIST_Met)),]



#########################################
# Collect PurIST data for metastases

# Met samples
pORG_met <- data.frame(cohort = metadata$Cohort,
                       subtype = metadata$PurIST_Subtype_Met,
                       pORG = metadata$pORG_Met,
                       PurIST = metadata$PurIST_Met)

# Met samples by subtype
pORG_subtype_met <- pORG_met
pORG_subtype_met$cohort <- pORG_subtype_met$subtype

pORG_all_met <- rbind(pORG_met, pORG_subtype_met)
pORG_all_met <- pORG_all_met[!(is.na(pORG_all_met$cohort)),]



#########################################
# Plot data

# Get one-way ANOVA (adjusted by BH method) p-values for plot
p.adjust(
    c(
        summary(aov(PurIST ~ cohort, data = pORG_all_met[pORG_all_met$cohort %in% c("Liver","Lung"),] ))[[1]][["Pr(>F)"]][[1]],
        summary(aov(PurIST ~ cohort, data = pORG_all_met[pORG_all_met$cohort %in% c("classical","basal-like"),] ))[[1]][["Pr(>F)"]][[1]]
    ),
    method = "BH")



plot_data <- pORG_all_met
plot_data$cohort <- factor(plot_data$cohort, levels = c("Liver","Lung","basal-like","classical"))



# Plot data
pdf(width = (5.1*1.25), height = (3.125*1.25), file = "PurIST_Metastases.pdf")
ggplot(plot_data, aes(x = cohort, y = PurIST)) + geom_violin(color = "gray40") + 
    geom_jitter(width = 0.2, aes(color = cohort), size = 3, alpha = 0.8, show.legend = F, stroke = NA) + 
    scale_color_manual(values = c("Liver" = "#d55e00",
                                  "Lung" = "#0072b2",
                                  "basal-like" = "#000000",
                                  "classical" = "#cc79a7")) +
    stat_summary(fun=mean, geom="crossbar" , linewidth = 0.6) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.title.x = element_text(color = "black",size = 13),
          axis.title.y = element_text(color = "black", size = 17),
          axis.text.x = element_text(color = "black", size = 16),
          axis.text.y = element_text(color = "black", size = 17),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    xlab("Cohort N = 37 | PurIST Subtype N = 71") +
    scale_y_continuous(name = "PurIST Score (Metastases)",
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       limits = c(0, 1.15)) +
    geom_signif(comparisons = list(c("Liver", "Lung"), 
                                   c("basal-like", "classical")),
                annotation = c("P=0.045","P=4.9e-34"), tip_length = 0.015, size = 1, textsize = 5.5, vjust = 0, margin_top = 0.1)
dev.off()



#################################################################
# Primary vs Met - PurIST
#################################################################

# Supplemental Dataset 1
metadata <- readxl::read_xlsx("Supplemental_Dataset_1.xlsx", sheet = 1)

# Primaries (-T and -T2) and metastases
metadata <- metadata[!(is.na(metadata$PurIST_Primary)) | !(is.na(metadata$PurIST_Primary_T2)) | !(is.na(metadata$PurIST_Met)),]



#########################################
# Collect info for primaries and mets

pORG_final_primary <- data.frame(Tumor = "Primary",
                                 PurIST = c(metadata$PurIST_Primary[!(is.na(metadata$PurIST_Primary))],
                                            metadata$PurIST_Primary_T2[!(is.na(metadata$PurIST_Primary_T2))]))

pORG_final_met <- data.frame(Tumor = "Met",
                             PurIST = metadata$PurIST_Met[!(is.na(metadata$PurIST_Met))])

pORG_final <- rbind(pORG_final_primary,pORG_final_met)



#########################################
# Plot data

# Get one-way ANOVA p-value for plot
summary(aov(PurIST ~ Tumor, data = pORG_final))[[1]][["Pr(>F)"]][[1]]



plot_data <- pORG_final
plot_data$Tumor <- factor(plot_data$Tumor, levels = c("Primary","Met"))



pdf(width = 3.825, height = 3.7, file = "PurIST_Primary_vs_Met.pdf")
ggplot(plot_data, aes(x = Tumor, y = PurIST)) + geom_violin(color = "gray40") + 
    geom_jitter(width = 0.2, aes(color = Tumor), size = 3, alpha = 0.8, show.legend = F, stroke = NA) + 
    scale_color_manual(values = c("Primary" = "gray60", "Met" = "black")) +
    stat_summary(fun=mean, geom="crossbar" , linewidth = 0.6) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.title.x = element_text(color = "black",size = 14),
          axis.title.y = element_text(color = "black", size = 17),
          axis.text.x = element_text(color = "black", size = 17),
          axis.text.y = element_text(color = "black", size = 17),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    xlab("Specimen Type N = 289") +
    scale_y_continuous(name = "PurIST Score (All)",
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       limits = c(0, 1.15)) +
    geom_signif(comparisons = list(c("Primary", "Met")),
                annotation = c("P=0.40"), tip_length = 0.015, size = 1, textsize = 5.5, vjust = 0, margin_top = 0.1)
dev.off()


