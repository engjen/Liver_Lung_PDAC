#################################################################
# Figure 7D
#################################################################
# Load libraries

library(ggplot2)
library(ggsignif)
set.seed(1)



#################################################################
# Load data

# Supplemental Dataset 1
supp1 <- readxl::read_xlsx("Supplemental_Dataset_1.xlsx", sheet = 1)

# Supplemental Dataset 8
supp8 <- read.csv("Supplemental_Dataset_8.csv", header = T)


#################################################################
# Arrange data for plot

# Collect necessary data
data <- data.frame(PatientID = supp1$Public_Patient_ID,
                   TLS = supp1$TLS,
                   TLS_Met = supp1$TLS_Met)
data$PercentTumorDistinctClones <- supp8$Percent.Tumor.Distinct.Clones[match(data$PatientID,supp8$Public_Patient_ID)]
data <- data[!(is.na(data$PercentTumorDistinctClones)),]

data$group[data$TLS=="N" | data$TLS_Met=="N"] = "Tumors\nWith No TLSs"
data$group[data$TLS=="Y" | data$TLS_Met=="Y"] = "Tumors\nWith TLSs"
data <- data[!(is.na(data$group)),]



# Run t-test manually
t.test(x = data$PercentTumorDistinctClones[data$group=="Tumors\nWith No TLSs"],
       y = data$PercentTumorDistinctClones[data$group=="Tumors\nWith TLSs"],
       alternative = "two.sided")



#################################################################
# Plot data

pdf(height = 5.1, width = 5.0, file = "PctTumorDistinctClones_vs_TLS.pdf")
ggplot(data, aes(x = group, y = PercentTumorDistinctClones)) +
    geom_violin() +
    geom_jitter(width = 0.1) +
    ylab("Percentage of Tumor-Distinct Clones") +
    xlab("") +
    stat_summary(fun=median, geom="crossbar" , linewidth = 0.4, color = "blue") + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.title.x = element_text(color = "black",size = 13),
          axis.title.y = element_text(color = "black", size = 14, face = "bold"),
          axis.text.x = element_text(color = "black", size = 16),
          axis.text.y = element_text(color = "black", size = 17),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    ylim(0,105) +
    geom_signif(comparisons = list(c("Tumors\nWith No TLSs", "Tumors\nWith TLSs")),
                tip_length = c(0.03,0.4),
                size = 1,
                y_position = 97,
                textsize = 5.5,
                vjust = -0.05,
                test = "t.test")
dev.off()
