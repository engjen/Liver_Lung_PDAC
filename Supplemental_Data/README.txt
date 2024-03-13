####### Supplemental Datasets #########

Supplemental_Dataset_1.xlsx: 
<<<<<<< HEAD
Sheet: Patient_Metadata - Patient level data including survival/outcome data, clinical covariates such as age, grade, stage, sex, lymph nodes positive, lymphovascular invasion and neoadjuvant treatment. "Survival" is a one hot encoding on vital status and Recurrence is a one hot encoding of follow up status for days from resection to recurrence (used for Kaplan-meier analysis). 

=======

Sheet: Patient_Metadata - Patient level data including survival/outcome data, clinical covariates such as age, grade, stage, sex, lymph nodes positive, lymphovascular invasion and neoadjuvant treatment. "Survival" is a one hot encoding on vital status and Recurrence is a one hot encoding of follow up status for days from resection to recurrence (used for Kaplan-meier analysis). 

>>>>>>> 6e72060997b914e02d61099bc9922b151473ea17
Metastasis site (i.e. liver, lung, other site).
"Cohort" column has whether patients were in liver or lung cohort (including resected and non-resected patients) 
"Recurrence_Sites_4" column has the metastatic recurrence sites for all resected patients

Also, columns indicating whether the patient had a primary tumor resection ("Resected"), resection of a lung or liver met ("Cohort_Met_Resection"), and if RNA DNA or TCR sequencing was done for the patient (columns "RNAseq_Patient", "DNAseq_Patient", "TCRseq_Tumor_Patient","TCRseq_Blood_Patient").
 
Finally, columns with results for histology analyses of primary and mets, GSVA scores for pORG, pSUB and PurIST from primary or met samples, and DNA alterations in primary or met (for variants present in >9 patients)

Sheet: RNA_Specimen_Metadata - Specimen level RNA data, including GSVA scores and specimen site. Public_Specimen_ID column is the patient ID plus a suffix indicating primary tumor "-T", a second primary tumor specimen "-T2", metastasis "-M", a second met "-M2", a fine needle aspirate from the primary "-F" and a primary from rapid autopsy "-A-T". To facilitate analysis and generating figures, this sheet also repeats selected patient level metadata and includes data derived from RNA and DNA data such as PurIST and GSVA scores and mutation status.
<<<<<<< HEAD

Sheet: DNA_Specimen_Metadata - Speciemen level DNA data: tumor mutation burden (TMB), estimated tumor cellularity (from DNA or pathologist), Microsatellite Instability, Normal Sample Source and whether or not the specimen has a homologous recombination (HR) or DNA damage repair (DDR) alteration. To facilitate analysis and generating figures, this sheet also repeats selected patient level metadata and data derived from RNA and DNA data. 
=======

Sheet: DNA_Specimen_Metadata - Speciemen level DNA data: tumor mutation burden (TMB), estimated tumor cellularity (from DNA or pathologist), Microsatellite Instability, Normal Sample Source and whether or not the specimen has a homologous recombination (HR) or DNA damage repair (DDR) alteration. To facilitate analysis and generating figures, this sheet also repeats selected patient level metadata and data derived from RNA and DNA data. 

>>>>>>> 6e72060997b914e02d61099bc9922b151473ea17
________________________________________________
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
________________________________________________
<<<<<<< HEAD
=======

>>>>>>> 6e72060997b914e02d61099bc9922b151473ea17
Supplemental_Dataset_3: RNA expression data and DESeq2 differentially expressed genes.
Sheet: pORG_Up_55 - DESeq2 statistics for pORG gene set.
Sheet: pSUB_Up_51 - DESeq2 statistics for pSUB gene set.
Sheet: pORG_FactorFromModel - DESeq2 statistics for pORG factor listing all genes passing expression filter.
Sheet: pSUB_FactorFromModel - DESeq2 statistics for pSUB factor listing all genes passing expression filter.
Sheet: Kallisto_TPM - TPM (Transcripts Per Million - normalized for sequencing depth and gene length) values for all RNA specimens (from kallisto pipeline with rows reduced to unique HUGO symbols).
Sheet: Kallisto_Counts - Gene level counts estimates for all RNA specimens (from kallisto pipeline with rows reduced to unique HUGO symbols).
Sheet: edgeR_TMM – EdgeR filtered, normalized counts using TMM (Trimmed Median of Means - more robust between sample normalization compared to TMM, but not normalized for gene length) method. (based on tximport from kallisto pipeline with rows reduced to unique HUGO symbols after filtering and normalizing with EdgeR, ).
________________________________________________
Supplemental_Dataset_4: Mutation calls from DNA xT panel dataset. 
Sheet: TempusReportedVariants – Gene variants called by Tempus using their xT panel.  Columns include “Public_Patient_ID”, “Public_Specimen_ID”, “Protein” and “Nucleic_Acid” changes, “Alteration_Type” (including LOF and GOF when known, these categories are simplified for display on OncoPrints), “Classification” (mostly Somatic with one Germline), “Significance” (Biologically relevant or VUS), “Panel_Name” (two panel versions used), “VCF_File” containing annotating the alteration, “VCF_FileIdx”, numeric index of alteration in vcf file, “CHROM” chromosome, “POS” position, “ID”, “REF” reference allele, “ALT” alternative allele, “TUMOR_REF_COUNT”, “TUMOR_ALT_COUNT”, “TUMOR_VAF” variant allele frequency for tumor, “NORMAL_REF_COUNT”, “NORMAL_ALT_COUNT”, “NORMAL_VAF” variant allele frequency for normal, “INFO”	column from vcf file, “ANN_TYPE” matching annotation type parsed out from INFO column, “ANN_IMPACT” matching annotation impact parsed out from INFO column.  Copy number gain/loss calls are included in this sheet and for these, many columns are not relevant and are left blank or set to NA.
________________________________________________
Supplemental_Dataset_5: : Gene Set Enrichment Analysis (GSEA) and Gene Set Variation Analysis (GSVA) results for HALLMARKs of cancer pathways.  Multiple Sheets. There are “TopVsBottom” quartile comparisons for pORG scores, pSUB scores, and PurIST scores in both Primary specimens and in Met specimens.  There are also comparisons between full Liver and Lung cohorts for Primary and for Met specimens.  Each comparison is split into up and down pathways and each sheet contains standard output fields from the GSEA analysis software.  The links are not active.  There are also sheets for GSVA analysis of Primary specimens, Met specimens, and all specimens together.  Each of these sheets contains relative pathway scores (for all 50 pathways) for each specimen in the respective sets.
________________________________________________
Supplemental_Dataset_6: VIPER results for primary tumors.
Sheet: VIPER for Primary Tumors - Viper Regulon scores for all Primary specimens.
________________________________________________
Supplemental_Dataset_7.xlsx: multiplex immunohistochemistry data
Includes the following multiple sheets:
sample_avg_density: patient level data including the mean density of each cell type per patient, cohort, pORG and pSUB scores and high or low pORG (based on cutoff of 0.0249). 
liver_vs_lung: difference in density of cell types  (at the ROI level, not patient level) between liver and lung cohort,  FDR corrected p-values 
high_vs_low_porg: difference in density of cell types  (at the ROI level, not patient level) between high and low pORG,  FDR corrected p-values 
cell_type_gating: the gating scheme to get cell types
lymphoid_aggregates: number of lymphoid aggregates identified in each patient's slide
________________________________________________
Supplemental_Dataset_8.csv: Tcell receptor sequencing data metrics, summarized at the patient level. 
Metrics will reflect whether they were from a blood sample or a tumor sample.
Blood_Type: 
Primary if the blood was collected before resection of a primary tumor, Met if the blood was collected after a recurrence
If blood was collected after resection and before recurrence, it was not categorized as primary or met blood
For patients without a resection, the blood was considered primary if it was collected 180 days before the patient was confirmed metastasis free 	
The blood was considered met blood if it was collected after patient had mets confirmed on imaging or within 30 days before the mets being detected on imaging
Tumor_Type: primary if the patient's primary tumor had TCR seq, Met if it was a met
________________________________________________
Supplemental_Dataset_9.xlsx: sheet 1: Putative KRAS-specific CDR3Seqs occuring in liver and lung cohort tumor or blood
 sheet 2: Shared/clonal sequences from liver or lung **per clone data**: sum of frequncy, average frequency, # pts. present in, specific for pathogen
 sheet 3: Shared/clonal sequences from liver or lung **per patient data**: CD3R frequncy per pateint per shared clonal sequence
 sheet 4: Shared/clonal sequences from all tumors **per patient data**: CD3R frequncy per pateint per shared clonal sequence (present in 25% of patients and in top 50 clones in at least 1 patient)
