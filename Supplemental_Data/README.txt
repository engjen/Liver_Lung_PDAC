####### Supplemental Datasets #########

Supplemental_Dataset_1.csv: Patient level data including survival/outcome data, clinical covariates, metastasis site (liver, lung, other, "Cohort" column has whether patients were in liver or lung cohort [including resected and non-resected patients], "Recurrence_Sites_4" column has the metaststic recurrence sites for all resected patients), neoadjuvant treatment, primary resection, resection of lung or liver met (column name "Cohort_Met_Resection"), histology, GSVA scores for pORG, pSUB and PurIST from primary of met sample, DNA alterations in primary or met (for variants present in >9 patients)

Supplemental_Dataset_2.xlsx: Cox proportional hazards modeling results. 
Multiple sheets, labeled multi for multivariable CPH and single for single variable.
The sheet name will include S, Surf or Survival if the outcome is overall survival
The sheet name will include R, Recur or Recurrence  if the outcome is recurrence free survival
columns:
covariate: the variable modeled against survival
exp(coef):  the hazard ratio
exp(coef) lower 95%: lower 95% confidence interval of hazard ratio
exp(coef) upper 95%: upper 95% confidence interval of hazard ratio
p: P-value of the model
model: for multivariable models, indicating which covariates were grouped into one model
n: number of patients

Supplemental_Dataset_3: DESeq2 differentially expressed genes

Supplemental_Dataset_4: The DNA dataset

Supplemental_Dataset_5: GSEA hallmarks results

Supplemental_Dataset_6: VIPER results

Supplemental_Dataset_7.xlsx: multiplex immunohistochemistry data
Includes the following multiple sheets:
sample_avg_density: patient level data including the mean density of each cell type per patient, cohort, pORG and pSUB scores and high or low pORG (based on cutoff of 0.0249). 
liver_vs_lung: difference in density of cell types  (at the ROI level, not patient level) between liver and lung cohort,  FDR corrected p-values 
high_vs_low_porg: difference in density of cell types  (at the ROI level, not patient level) between high and low pORG,  FDR corrected p-values 
cell_type_gating: the gating scheme to get cell types
lymphoid_aggregates: number of lymphoid aggregates identified in each patient's slide

Supplemental_Dataset_8.csv: Tcell receptor sequencing data metrics, summarized at the patient level. 
Metrics will reflect whether they were from a blood sample or a tumor sample.
Blood_Type: 
Primary if the blood was collected before resection of a primary tumor, Met if the blood was collected after a recurrence
If blood was collected after resection and before recurrence, it was not categorized as primary or met blood
For patients without a resection, the blood was considered primary if it was collected 180 days before the patient was confirmed metastasis free 	
The blood was considered met blood if it was collected after patient had mets confirmed on imaging or within 30 days before the mets being detected on imaging

Tumor_Type: primary if the patient's primary tumor had TCR seq, Met if it was a met

Supplemental_Dataset_9.xlsx: sheet 1: Putative KRAS-specific CDR3Seqs occuring in liver and lung cohort tumor or blood
 sheet 2: Shared/clonal sequences from liver or lung **per clone data**: sum of frequncy, average frequency, # pts. present in, specific for pathogen
 sheet 3: Shared/clonal sequences from liver or lung **per patient data**: CD3R frequncy per pateint per shared clonal sequence
 sheet 4: Shared/clonal sequences from all tumors **per patient data**: CD3R frequncy per pateint per shared clonal sequence (present in 25% of patients and in top 50 clones in at least 1 patient)