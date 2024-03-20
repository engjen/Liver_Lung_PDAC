#= GetCollapseIndex() ==========================================================

GetCollapseIndex = function(X, Symbols)
{
  # Sort by gene symbol and select only the symbol with the highest expression
  # for any duplicated symbols.
  OriginalOrder = 1:dim(X)[[1]]
  idx = order(Symbols)
  X = X[idx, ]
  Symbols = Symbols[idx]
  OriginalOrder = OriginalOrder[idx]
  idx = NULL
  i = 1
  for (j in 2:length(Symbols))
  {
    if (Symbols[i] == Symbols[j])
    {
      if (sum(X[j, ]) > sum(X[i, ]))
      {
        idx = c(idx, i)
        i = j
      } else
      {
        idx = c(idx, j)
      }
    } else
    {
      i = j
    }
  }
  
  if (!is.null(idx))
  {
    OriginalOrder = OriginalOrder[-idx]
  }
  return(OriginalOrder)
}

#= LoadDnaPanelData() ============================================================

#FileName = paste(SourcePath, "Public/OriginalData/", "Supplemental_Dataset_3.xlsx", sep="")
#require(XLConnect)
#  
#wb <- loadWorkbook(FileName)
#
#RnaSeqTPM = readWorksheet(wb, 
#                          sheet = "Kallisto_TPM", 
#                          startCol = which(LETTERS=="A"), 
#                          startRow = 1, # Include header if you want column names.
#                          endCol = 0, #which(LETTERS=="A") * 26 + which(LETTERS=="S"),
#                          endRow = 0)
#dim(RnaSeqTPM)
#
#RnaSeqCounts = readWorksheet(wb, 
#                             sheet = "Kallisto_Counts", 
#                             startCol = which(LETTERS=="A"), 
#                             startRow = 1, # Include header if you want column names.
#                             endCol = 0, #which(LETTERS=="A") * 26 + which(LETTERS=="S"),
#                             endRow = 0)
#dim(RnaSeqCounts)
#
#edgeR_TMM = readWorksheet(wb, 
#                          sheet = "edgeR_TMM", 
#                          startCol = which(LETTERS=="A"), 
#                          startRow = 1, # Include header if you want column names.
#                          endCol = 0, #which(LETTERS=="A") * 26 + which(LETTERS=="S"),
#                          endRow = 0)
#dim(edgeR_TMM)
#
#edgeR_TMM_Filtered = readWorksheet(wb, 
#                                   sheet = "edgeR_TMM_Filtered", 
#                                   startCol = which(LETTERS=="A"), 
#                                   startRow = 1, # Include header if you want column names.
#                                   endCol = 0, #which(LETTERS=="A") * 26 + which(LETTERS=="S"),
#                                   endRow = 0)
#dim(edgeR_TMM_Filtered)

LoadDnaPanelData = function(SourcePath)
{
  HR_DDR_Genes = c("ARID1A", "ATM", "ATRX", "BAP1", "BARD1", "BLM", "BRCA1", 
                   "BRCA2", "BRIP1", "CHEK1", "CHEK2", "FANCA", "FANCC", "FANCD2", 
                   "FANCE", "FANCF", "FANCG", "FANCL", "MRE11", "NBN", "PALB2", 
                   "RAD50", "RAD51", "RAD51B", "WRN")
  
  # DNA_Specimen_Metadata:
  
  #DNA_Specimen_Metadata = read.table(paste(SourcePath, "OriginalData/DNA_Specimen_Metadata.tsv", sep=""),  
  #                                   sep="\t", quote="", comment.char = "", header=TRUE, stringsAsFactors=FALSE) 
  
  FileName = paste(SourcePath, "Supplemental_Data/", "Supplemental_Dataset_1.xlsx", sep="")
  require(XLConnect)
  wb <- loadWorkbook(FileName)
  
  DNA_Specimen_Metadata = readWorksheet(wb, 
                                        sheet = "DNA_Specimen_Metadata", 
                                        startCol = which(LETTERS=="A"), 
                                        startRow = 1, # Include header if you want column names.
                                        endCol = 0, #which(LETTERS=="A") * 26 + which(LETTERS=="S"),
                                        endRow = 0)
  dim(DNA_Specimen_Metadata)
  names(DNA_Specimen_Metadata)
  idx = which(names(DNA_Specimen_Metadata) %in% c("Public_Patient_ID", 
                                                  "Public_Specimen_ID", 
                                                  "Tumor_Type", 
                                                  "Specimen_Site", 
                                                  "Liver_Met_Present", 
                                                  "Lung_Met_Present", 
                                                  "Days_from_Diagnosis_to_FU", 
                                                  "Days_from_Resection_to_Recurrence", 
                                                  "Days_from_Resection_to_FU", 
                                                  "Days_from_Earliest_Recur_to_FU", 
                                                  "Vital_Status_at_FU", 
                                                  "Neoadjuvant_Treatment", 
                                                  "Percent_Tumor_Pathologist_Estimate", 
                                                  "Percent_Tumor_from_DNA", 
                                                  "Tumor_Mutation_Burden", 
                                                  "Microsatellite_Instability", 
                                                  "NormalSampleSource"))
  DNA_Specimen_Metadata = DNA_Specimen_Metadata[, idx]
  
  # TempusReportedVariants:
  FileName = paste(SourcePath, "Supplemental_Data/", "Supplemental_Dataset_4.xlsx", sep="")
  require(XLConnect)
  wb <- loadWorkbook(FileName)
  
  TempusReportedVariants = readWorksheet(wb, 
                                        sheet = "Reported_Mutations", 
                                        startCol = which(LETTERS=="A"), 
                                        startRow = 1, # Include header if you want column names.
                                        endCol = 0, #which(LETTERS=="A") * 26 + which(LETTERS=="S"),
                                        endRow = 0)
  dim(TempusReportedVariants)
  names(TempusReportedVariants)
  
  #TempusReportedVariants = read.table(paste(SourcePath, "OriginalData/TempusReportedVariants.tsv", sep=""),  
  #                                    sep="\t", quote="", comment.char = "", header=TRUE, 
  #                                    stringsAsFactors=FALSE)    

  return(list(Metadata = DNA_Specimen_Metadata,
              Variants = TempusReportedVariants,
              HR_DDR_Genes = HR_DDR_Genes))
}

#= LoadRnaSeqData() ============================================================

LoadRnaSeqData = function(SourcePath, AnalysisOutputPath)
{
  if (length(list.files(path = AnalysisOutputPath, pattern = "^BCC_RNA_Data.RData$")) == 1)
  {
    # Load saved analysis.
    load(file=paste(AnalysisOutputPath,"BCC_RNA_Data.RData", sep=""))
  } else
  {
    RnaSeqTPM = read.table(paste(SourcePath, "data/RnaSeqTPM.tsv", sep=""),  
                 sep="\t", quote="", comment.char = "", header=TRUE, 
                 row.names=1, check.names=FALSE, stringsAsFactors=FALSE)    
    RnaSeqCounts = read.table(paste(SourcePath, "data/RnaSeqCounts.tsv", sep=""),  
                 sep="\t", quote="", comment.char = "", header=TRUE, 
                 row.names=1, check.names=FALSE, stringsAsFactors=FALSE)    
    edgeR_TMM = read.table(paste(SourcePath, "data/RnaSeqDGE_All_cpms.tsv", sep=""),  
                 sep="\t", quote="", comment.char = "", header=TRUE, 
                 row.names=1, check.names=FALSE, stringsAsFactors=FALSE)    
    edgeR_TMM_Filtered = read.table(paste(SourcePath, "data/RnaSeqDGE_Filtered_cpms.tsv", sep=""),  
                 sep="\t", quote="", comment.char = "", header=TRUE, 
                 row.names=1, check.names=FALSE, stringsAsFactors=FALSE)    
    
    Result = CacheEnsemblAnnotion(RnaSeqData = RnaSeqTPM, 
                                  SourcePath = AnalysisOutputPath,
                                  Version = 111)
                                  #Version = "Original")
    EnsemblAnn = Result$EnsemblAnn
    ExtraAnn = Result$ExtraAnn
    rm(Result)
    
    # Remove ".version" numbers from Ensembl IDs if they are there.
    # It's assumed that RnaSeqTPM & RnaSeqCounts have same rows, same number and order.
    RnaSeqTPM = as.matrix(RnaSeqTPM)
    RnaSeqCounts = as.matrix(RnaSeqCounts)
    edgeR_TMM = as.matrix(edgeR_TMM)
    edgeR_TMM_Filtered = as.matrix(edgeR_TMM_Filtered)
    Target_Ensembl_ID = rownames(RnaSeqTPM)
    temp = unlist(strsplit(Target_Ensembl_ID, split=".", fixed=TRUE))
    if (length(temp) == 2*length(Target_Ensembl_ID))
    {
      idx = 2 * 1:length(Target_Ensembl_ID) - 1
      Target_Ensembl_ID = temp[idx]
    }
    rownames(RnaSeqTPM) = Target_Ensembl_ID
    rownames(RnaSeqCounts) = Target_Ensembl_ID
  
    # Match ensembl IDs from the kallisto output with the ensembl IDs
    # from MioMart annotation.
    Target_Ensembl_ID = rownames(RnaSeqTPM)
    idx = match(x=Target_Ensembl_ID, table=EnsemblAnn[, "ensembl_gene_id"])
    EnsemblAnn = EnsemblAnn[idx, ]
    dim(EnsemblAnn)
    
    ExtraAnn = ExtraAnn[!duplicated(ExtraAnn$ENSEMBL), ]
    dim(ExtraAnn)
    
    # I think these are already in order, but just to be sure...
    idx = match(x=Target_Ensembl_ID, table=ExtraAnn[, "ENSEMBL"])
    ExtraAnn = ExtraAnn[idx, ]
    dim(ExtraAnn)
    
    # Create a matrix of RnaSeqTPM data where rows are unique gene symbols and
    # any redundant symbols are resolved by selecting the variant with the highest
    # level of expression across the dataset.  Features without HUGO gene symbols
    # are excluded. This helps with pathway tools based on HUGO symbols.
    Symbols = ExtraAnn[, c("SYMBOL")]
    Symbols[is.na(Symbols)] = "NA"
    CollapseIndex = GetCollapseIndex(X=RnaSeqTPM, Symbols = Symbols)
    idx = which(Symbols[CollapseIndex] == "NA")
    if (length(idx) > 0) CollapseIndex = CollapseIndex[-idx]
    length(CollapseIndex)
    RnaSeqCollapsedTPM = RnaSeqTPM[CollapseIndex, ]
    RnaSeqCollapsedCounts = RnaSeqCounts[CollapseIndex, ]
    rownames(RnaSeqCollapsedTPM) = ExtraAnn[CollapseIndex, c("SYMBOL")]
    rownames(RnaSeqCollapsedCounts) = ExtraAnn[CollapseIndex, c("SYMBOL")]
    
    write.table(cbind(gene = rownames(RnaSeqCollapsedTPM), RnaSeqCollapsedTPM), 
                paste(AnalysisOutputPath, "RnaSeqCollapsedTPM.tsv", sep=""), 
                sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    
    write.table(cbind(gene = rownames(RnaSeqCollapsedCounts), RnaSeqCollapsedCounts), 
                paste(AnalysisOutputPath, "RnaSeqCollapsedCounts.tsv", sep=""), 
                sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    
    # RNA_Specimen_Metadata:
    
    #RNA_Specimen_Metadata = read.table(paste(SourcePath, "OriginalData/RNA_Specimen_Metadata.tsv", sep=""),  
    #             sep="\t", quote="", comment.char = "", header=TRUE, stringsAsFactors=FALSE)   
    
    FileName = paste(SourcePath, "Supplemental_Data/", "Supplemental_Dataset_1.xlsx", sep="")
    require(XLConnect)
    wb <- loadWorkbook(FileName)
    RNA_Specimen_Metadata = readWorksheet(wb, 
                                          sheet = "RNA_Specimen_Metadata", 
                                          startCol = which(LETTERS=="A"), 
                                          startRow = 1, # Include header if you want column names.
                                          endCol = 0, #which(LETTERS=="A") * 26 + which(LETTERS=="S"),
                                          endRow = 0)
    dim(RNA_Specimen_Metadata)
    names(RNA_Specimen_Metadata)
    
    idx = which(names(RNA_Specimen_Metadata) %in% c("Public_Patient_ID", 
                                                    "Public_Specimen_ID", 
                                                    "Tumor_Type", 
                                                    "Specimen_Site", 
                                                    "Percent_Tumor_from_DNA", 
                                                    "Liver_Met_Present", 
                                                    "Lung_Met_Present", 
                                                    "Days_from_Diagnosis_to_FU", 
                                                    "Days_from_Resection_to_Recurrence", 
                                                    "Days_from_Resection_to_FU", 
                                                    "Days_from_Earliest_Recur_to_FU", 
                                                    "Vital_Status_at_FU", 
                                                    "Neoadjuvant_Treatment"))
    RNA_Specimen_Metadata = RNA_Specimen_Metadata[, idx]
    
    # Make sure metadata rows match matrix columns.
    idx = match(x=colnames(RnaSeqTPM), table=RNA_Specimen_Metadata$Public_Specimen_ID)
    if (length(which(is.na(idx))) != 0)
    {
      stop("RNA_Specimen_Metadata does not match column names of RnaSeqTPM!")
    }
    if (length(idx) !=  dim(RnaSeqTPM)[[2]])
    {
      stop("RNA_Specimen_Metadata does not have the same number of specimens as RnaSeqTPM!")
    }
    RNA_Specimen_Metadata = RNA_Specimen_Metadata[idx, ]

    save(RnaSeqTPM, 
         RnaSeqCounts, 
         RnaSeqCollapsedTPM, 
         RnaSeqCollapsedCounts, 
         CollapseIndex,
         edgeR_TMM,
         edgeR_TMM_Filtered,
         EnsemblAnn,
         ExtraAnn,
         RNA_Specimen_Metadata, 
         file=paste(AnalysisOutputPath,"BCC_RNA_Data.RData", sep=""))
  }
  
  return(list(TPM = RnaSeqTPM, 
              Counts = RnaSeqCounts, 
              CollapsedTPM = RnaSeqCollapsedTPM, 
              CollapsedCounts = RnaSeqCollapsedCounts, 
              CollapseIndex = CollapseIndex,
              edgeR_TMM = edgeR_TMM,
              edgeR_TMM_Filtered = edgeR_TMM_Filtered,
              EnsemblAnn = EnsemblAnn,
              ExtraAnn = ExtraAnn,
              Metadata = RNA_Specimen_Metadata))
}

#= CalculateBCCpORGpSUBandGSVAscores() =========================================
CalculateBCCpORGpSUBandGSVAscores = function(AnalysisOutputPath, BCC_RNA_Data)
{
  if (length(list.files(path = AnalysisOutputPath, pattern = "^pORGpSUBandGSVA.RData$")) == 1)
  {
    load(file=paste(AnalysisOutputPath,"pORGpSUBandGSVA.RData", sep=""))
  } else
  {
    # Theses results are consistent between different versions of R and packages
    # tested.
    BCC_DESeq2Result = DoTwoFactorDESeq2(data = BCC_RNA_Data$CollapsedCounts, 
                                         Metadata = BCC_RNA_Data$Metadata, 
                                         FCThreshold = 1.0, # No FC filtering.
                                         TPMs = BCC_RNA_Data$CollapsedTPM,
                                         averageTPMThreshold = 0.25) # Exclude genes with very low expression across data set.
    
    For_pORG = ApplyExclude(DESeq2Results = BCC_DESeq2Result,
                            GeneAnn = NULL,
                            FCThreshold = 1.0,
                            pValThreshold = 0.01,
                            padjThreshold = 0.20,
                            UseNominalPval = FALSE,
                            GetExtraInfo = TRUE)
    
    For_pSUB = ApplyExclude(DESeq2Results = BCC_DESeq2Result,
                            GeneAnn = NULL,
                            FCThreshold = 1.0,
                            pValThreshold = 0.01,
                            padjThreshold = 0.0001,
                            UseNominalPval = FALSE,
                            GetExtraInfo = TRUE)
    
    # Write some reports that can go into supp data sets.
    write.table(cbind(gene=rownames(For_pORG[[5]]), For_pORG[[5]]), 
                paste(AnalysisOutputPath, "DESeq2_pORG_GeneList.tsv", sep=""), 
                sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    
    write.table(cbind(gene=rownames(For_pSUB[[6]]), For_pSUB[[6]]), 
                paste(AnalysisOutputPath, "DESeq2_pSUB_GeneList.tsv", sep=""), 
                sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    
    write.table(cbind(genes = rownames(BCC_DESeq2Result[[3]]), BCC_DESeq2Result[[3]]), 
                paste(AnalysisOutputPath, "DESeq2_pORG_FactorFromModel.tsv", sep=""), 
                sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    
    write.table(cbind(genes = rownames(BCC_DESeq2Result[[2]]), BCC_DESeq2Result[[2]]), 
                paste(AnalysisOutputPath, "DESeq2_pSUB_FactorFromModel.tsv", sep=""), 
                sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
    
    #= Calculate GSVA scores. ======================================================
    # Note that gsva() call has been updated.
    # Load HALLMARKS of cancer gene sets (50) and format for GSVA.
    HallmarkGeneSetsList = LoadGeneSets()
    
    pORGpSUB_GeneSets = list(pORG_Up_55 = For_pORG[[1]],
                             pSUB_Up_51 = For_pSUB[[2]])
    pORGpSUB_GeneSetsList = MakeGeneList(pORGpSUB_GeneSets)
    
    # Use DGE_cpms, may be more robust, but not gene length normalized?
    M = log2(BCC_RNA_Data$edgeR_TMM_Filtered + 1)
    Ensembl_ID = rownames(M)
    idx = match(x=Ensembl_ID, table=BCC_RNA_Data$EnsemblAnn[, "ensembl_gene_id"])
    M = M[!is.na(idx), ]
    dim(M)
    rownames(M) = BCC_RNA_Data$EnsemblAnn[idx[!is.na(idx)], c("hgnc_symbol")]
    M = M[rownames(M) != "", ]
    dim(M)

    PrimaryIdx = which(!is.na(BCC_RNA_Data$Metadata$Tumor_Type) &
                         BCC_RNA_Data$Metadata$Tumor_Type == "Primary")
    MetIdx = which(!is.na(BCC_RNA_Data$Metadata$Tumor_Type) &
                     BCC_RNA_Data$Metadata$Tumor_Type == "Met")
    AllIdx = which(!is.na(BCC_RNA_Data$Metadata$Tumor_Type) &
                     (BCC_RNA_Data$Metadata$Tumor_Type == "Primary" | BCC_RNA_Data$Metadata$Tumor_Type == "Met"))
    
    length(PrimaryIdx)
    length(MetIdx)
    length(AllIdx)
    
    GSVA_Hallmarks_All = ApplyGSVA(M[, AllIdx], HallmarkGeneSetsList, RowCV_ThreshFactor = 0.1)
    dim(GSVA_Hallmarks_All)
    rownames(GSVA_Hallmarks_All) = colnames(M[, AllIdx])
    
    GSVA_Hallmarks_Primaries = ApplyGSVA(M[, PrimaryIdx], HallmarkGeneSetsList, RowCV_ThreshFactor = 0.1)
    dim(GSVA_Hallmarks_Primaries)
    rownames(GSVA_Hallmarks_Primaries) = colnames(M[, PrimaryIdx])
    
    GSVA_Hallmarks_Mets = ApplyGSVA(M[, MetIdx], HallmarkGeneSetsList, RowCV_ThreshFactor = 0.1)
    dim(GSVA_Hallmarks_Mets)
    rownames(GSVA_Hallmarks_Mets) = colnames(M[, MetIdx])
    
    GSVA_pORGpSUB_All = ApplyGSVA(M[, AllIdx], pORGpSUB_GeneSetsList, RowCV_ThreshFactor = 0.1)
    dim(GSVA_pORGpSUB_All)
    rownames(GSVA_pORGpSUB_All) = colnames(M[, AllIdx])
    
    GSVA_pORGpSUB_Primaries = ApplyGSVA(M[, PrimaryIdx], pORGpSUB_GeneSetsList, RowCV_ThreshFactor = 0.1)
    dim(GSVA_pORGpSUB_Primaries)
    rownames(GSVA_pORGpSUB_Primaries) = colnames(M[, PrimaryIdx])
    
    GSVA_pORGpSUB_Mets = ApplyGSVA(M[, MetIdx], pORGpSUB_GeneSetsList, RowCV_ThreshFactor = 0.1)
    dim(GSVA_pORGpSUB_Mets)
    rownames(GSVA_pORGpSUB_Mets) = colnames(M[, MetIdx])
    
    save(BCC_DESeq2Result, For_pORG, For_pSUB, 
         HallmarkGeneSetsList, pORGpSUB_GeneSetsList,
         AllIdx, PrimaryIdx, MetIdx,
         GSVA_Hallmarks_All, GSVA_Hallmarks_Primaries, GSVA_Hallmarks_Mets,
         GSVA_pORGpSUB_All, GSVA_pORGpSUB_Primaries, GSVA_pORGpSUB_Mets,
         file=paste(AnalysisOutputPath,"pORGpSUBandGSVA.RData", sep=""))
  }
  
  BCC_RNA_Data$BCC_DESeq2Result = BCC_DESeq2Result
  BCC_RNA_Data$For_pORG = For_pORG
  BCC_RNA_Data$For_pSUB = For_pSUB
  BCC_RNA_Data$HallmarkGeneSetsList = HallmarkGeneSetsList
  BCC_RNA_Data$pORGpSUB_GeneSetsList = pORGpSUB_GeneSetsList
  BCC_RNA_Data$AllIdx = AllIdx
  BCC_RNA_Data$PrimaryIdx = PrimaryIdx
  BCC_RNA_Data$MetIdx = MetIdx
  BCC_RNA_Data$GSVA_Hallmarks_All = GSVA_Hallmarks_All
  BCC_RNA_Data$GSVA_Hallmarks_Primaries = GSVA_Hallmarks_Primaries
  BCC_RNA_Data$GSVA_Hallmarks_Mets = GSVA_Hallmarks_Mets
  BCC_RNA_Data$GSVA_pORGpSUB_All = GSVA_pORGpSUB_All
  BCC_RNA_Data$GSVA_pORGpSUB_Primaries = GSVA_pORGpSUB_Primaries
  BCC_RNA_Data$GSVA_pORGpSUB_Mets = GSVA_pORGpSUB_Mets
  
  return(BCC_RNA_Data)
}

#= DoTwoFactorDESeq2() =========================================================

DoTwoFactorDESeq2 = function(data, # List for tximport or matrix if not.
                             Metadata = NULL, 
                             MetStatus = NULL,
                             SubTypeStatus = NULL,
                             FCThreshold = 1.0,
                             TPMs = NULL,
                             averageTPMThreshold = 0.25, # Applied when using matrix data.
                             min.count = 10, # Applied when using tximport method (list).
                             min.total.count = 15) # Applied when using tximport method (list).
{
  if (!is.null(Metadata))
  {
    LiverCohort = which(!is.na(Metadata$Liver_Met_Present) & 
                          Metadata$Liver_Met_Present == 1)
    LungCohort = which(!is.na(Metadata$Lung_Met_Present) & 
                         Metadata$Lung_Met_Present == 1 &
                         (is.na(Metadata$Liver_Met_Present) |
                          Metadata$Liver_Met_Present == "NA"))
    
    PrimaryIdx = which(!is.na(Metadata$Tumor_Type) &
                         Metadata$Tumor_Type == "Primary")
    
    SelectedSet = PrimaryIdx
    SelectedSet = SelectedSet[which(SelectedSet %in% c(LiverCohort, LungCohort))]
    
    MetStatus = rep(NA, times=length(SelectedSet))
    MetStatus[SelectedSet %in% LiverCohort] = "Liver"
    MetStatus[SelectedSet %in% LungCohort] = "LungOnly"
    
    SubTypeStatus = Metadata$PurIST_Subtype[SelectedSet]
    SubTypeStatus[SubTypeStatus == "classical"] = "Classical"
    SubTypeStatus[SubTypeStatus == "basal-like"] = "BasalLike"
  } else
  {
    SelectedSet = which(!is.na(MetStatus) &
                          (MetStatus %in% c("Liver", "LungOnly"))) 
    MetStatus = MetStatus[SelectedSet]
    SubTypeStatus = SubTypeStatus[SelectedSet]
  }
  
  colData = data.frame(SubTypeStatus = factor(SubTypeStatus), MetStatus = factor(MetStatus))
  if (is.matrix(data))
  {
    if (!is.null(TPMs))
    {
      SelectedGeneIdx = which(rowSums(TPMs[, SelectedSet]) / length(SelectedSet) >= averageTPMThreshold)
    } else
    {
      SelectedGeneIdx = which(rowSums(data[, SelectedSet]) / length(SelectedSet) >= averageTPMThreshold)
    }
    M = as.matrix(data[SelectedGeneIdx, SelectedSet])
    M = apply(M + 0.5, MARGIN=c(1,2), FUN=floor)
    dds1 = DESeqDataSetFromMatrix(countData=M, colData=colData, design=~SubTypeStatus + MetStatus)
    featureData <- data.frame(EnsemblID=rownames(M))
  } else if (is.list(data))
  {
    TxiSelected = data
    TxiSelected$abundance = TxiSelected$abundance[, SelectedSet]
    TxiSelected$counts = TxiSelected$counts[, SelectedSet]
    TxiSelected$length = TxiSelected$length[, SelectedSet]
    cts <- TxiSelected$counts
    normMat <- TxiSelected$length
    # Obtaining per-observation scaling factors for length, adjusted to avoid
    # changing the magnitude of the counts.
    normMat <- normMat/exp(rowMeans(log(normMat)))
    normCts <- cts/normMat
    eff.lib <- calcNormFactors(normCts) * colSums(normCts)
    # Combining effective library sizes with the length factors, and calculating
    # offsets for a log-link GLM.
    normMat <- sweep(normMat, 2, eff.lib, "*")
    normMat <- log(normMat)
    # Creating a DGEList object for use in edgeR.
    y <- DGEList(cts)
    y <- scaleOffset(y, normMat)
    # filtering using the design information
    # Do not use a design here.
    design = NULL
    keep <- filterByExpr(y, design, min.count = min.count, min.total.count = min.total.count)
    SelectedGeneIdx = which(keep)
    TxiSelected$abundance = TxiSelected$abundance[SelectedGeneIdx, ]
    TxiSelected$counts = TxiSelected$counts[SelectedGeneIdx, ]
    TxiSelected$length = TxiSelected$length[SelectedGeneIdx, ]
    dds1 = DESeqDataSetFromTximport(TxiSelected, colData=colData, design=~SubTypeStatus + MetStatus)
    featureData <- data.frame(EnsemblID=rownames(TxiSelected$counts))
  }
  mcols(dds1) <- DataFrame(mcols(dds1), featureData)
  dds1$SubTypeStatus = factor(SubTypeStatus, levels=c("BasalLike", "Classical"))
  dds1$MetStatus = factor(MetStatus, levels=c("Liver", "LungOnly"))
  
  dds1 = DESeq(dds1, fitType = "parametric", betaPrior = TRUE, minReplicatesForReplace = 1000000)
  
  dd1ResultSubTypeStatus = 
    results(dds1,
            contrast=c("SubTypeStatus", "BasalLike", "Classical"), # variable, condition A, condition B => A/B
            alpha=0.1) # In manual, it says that alpha should match the adjusted p-value cutoff.
  
  dd1ResultMetStatus = 
    results(dds1,
            contrast=c("MetStatus", "Liver", "LungOnly"), # variable, condition A, condition B => A/B
            alpha=0.1) # In manual, it says that alpha should match the adjusted p-value cutoff.
  return(list(dds1, dd1ResultSubTypeStatus, dd1ResultMetStatus))
}

#= ApplyExclude ================================================================

# Apply exclude method to pORG
ApplyExclude = function(DESeq2Results,
                        GeneAnn = NULL,
                        FCThreshold = 1.0,
                        pValThreshold = 0.01,
                        padjThreshold = 0.05,
                        UseNominalPval = FALSE,
                        GetExtraInfo = FALSE)
{
  DESeq_pSUB = DESeq2Results[[2]]
  DESeq_pORG = DESeq2Results[[3]]
  
  ExcludeThreshold = DESeq_pSUB$pvalue
  ExcludeThreshold[is.na(DESeq_pSUB$pvalue)] = 1
  ExcludeThreshold = ExcludeThreshold[order(ExcludeThreshold)]
  ExcludeThreshold = ExcludeThreshold[length(ExcludeThreshold) / 2]
  
  if (UseNominalPval)
  {
    IdxA = (1:dim(DESeq_pORG)[[1]])[(!is.na(DESeq_pORG$pvalue) & DESeq_pORG$pvalue < pValThreshold) &
                                      (!is.na(DESeq_pORG$log2FoldChange) & abs(DESeq_pORG$log2FoldChange) >= log2(FCThreshold))]
    LengthPreExcludeA = length(IdxA[DESeq_pORG$log2FoldChange[IdxA] > 0])
    IdxA = IdxA[is.na(DESeq_pSUB$pvalue[IdxA]) | DESeq_pSUB$pvalue[IdxA] > ExcludeThreshold]
  } else
  {
    IdxA = (1:dim(DESeq_pORG)[[1]])[(!is.na(DESeq_pORG$padj) & DESeq_pORG$padj < padjThreshold) &
                                      (!is.na(DESeq_pORG$log2FoldChange) & abs(DESeq_pORG$log2FoldChange) >= log2(FCThreshold))]
    LengthPreExcludeA = length(IdxA[DESeq_pORG$log2FoldChange[IdxA] > 0])
    IdxA = IdxA[is.na(DESeq_pSUB$pvalue[IdxA]) | DESeq_pSUB$pvalue[IdxA] > ExcludeThreshold]
  }
  
  IdxA = IdxA[DESeq_pORG$log2FoldChange[IdxA] > 0]
  IdxA = IdxA[order(DESeq_pORG$pvalue[IdxA])]
  
  ExcludeThreshold = DESeq_pORG$pvalue
  ExcludeThreshold[is.na(DESeq_pORG$pvalue)] = 1
  ExcludeThreshold = ExcludeThreshold[order(ExcludeThreshold)]
  ExcludeThreshold = ExcludeThreshold[length(ExcludeThreshold) / 2]
  
  if (UseNominalPval)
  {
    IdxB = (1:dim(DESeq_pSUB)[[1]])[(!is.na(DESeq_pSUB$pvalue) & DESeq_pSUB$pvalue < pValThreshold) &
                                      (!is.na(DESeq_pSUB$log2FoldChange) & abs(DESeq_pSUB$log2FoldChange) >= log2(FCThreshold))]
    LengthPreExcludeB = length(IdxB[DESeq_pSUB$log2FoldChange[IdxB] > 0])
    IdxB = IdxB[(is.na(DESeq_pORG$pvalue[IdxB]) | DESeq_pORG$pvalue[IdxB] > ExcludeThreshold)]
  } else
  {
    IdxB = (1:dim(DESeq_pSUB)[[1]])[(!is.na(DESeq_pSUB$padj) & DESeq_pSUB$padj < padjThreshold) &
                                      (!is.na(DESeq_pSUB$log2FoldChange) & abs(DESeq_pSUB$log2FoldChange) >= log2(FCThreshold))]
    LengthPreExcludeB = length(IdxB[DESeq_pSUB$log2FoldChange[IdxB] > 0])
    IdxB = IdxB[(is.na(DESeq_pORG$pvalue[IdxB]) | DESeq_pORG$pvalue[IdxB] > ExcludeThreshold)]
  }
  
  IdxB = IdxB[DESeq_pSUB$log2FoldChange[IdxB] > 0]
  IdxB = IdxB[order(DESeq_pSUB$pvalue[IdxB])]
  
  if (is.null(GeneAnn))
  {
    pORG_Up_Genes = rownames(DESeq_pORG)[IdxA]
    
    pSUB_Up_Genes = rownames(DESeq_pSUB)[IdxB]
  } else
  {
    Ensembl_ID = rownames(DESeq_pORG)[IdxA]
    idx = match(x=Ensembl_ID, table=GeneAnn[, "ensembl_gene_id"])
    IdxA = IdxA[!is.na(idx)]
    pORG_Up_Genes = GeneAnn[idx[!is.na(idx)], c("hgnc_symbol")]
    IdxA = IdxA[pORG_Up_Genes != ""]
    pORG_Up_Genes = pORG_Up_Genes[pORG_Up_Genes != ""]
    IdxA = IdxA[!duplicated(pORG_Up_Genes)]
    pORG_Up_Genes = pORG_Up_Genes[!duplicated(pORG_Up_Genes)]
    
    Ensembl_ID = rownames(DESeq_pSUB)[IdxB]
    idx = match(x=Ensembl_ID, table=GeneAnn[, "ensembl_gene_id"])
    IdxB = IdxB[!is.na(idx)]
    pSUB_Up_Genes = GeneAnn[idx[!is.na(idx)], c("hgnc_symbol")]
    IdxB = IdxB[pSUB_Up_Genes != ""]
    pSUB_Up_Genes = pSUB_Up_Genes[pSUB_Up_Genes != ""]
    IdxB = IdxB[!duplicated(pSUB_Up_Genes)]
    pSUB_Up_Genes = pSUB_Up_Genes[!duplicated(pSUB_Up_Genes)]
  }
  print(paste("length(pORG_Up_Genes) = ", length(pORG_Up_Genes), sep=""))
  print(paste("length(pSUB_Up_Genes) = ", length(pSUB_Up_Genes), sep=""))
  if (GetExtraInfo)
  {
    return(list(pORG_Up_Genes, pSUB_Up_Genes, LengthPreExcludeA, LengthPreExcludeB,
                DESeq_pORG[IdxA, ], DESeq_pSUB[IdxB, ]))
  } else
  {
    return(list(pORG_Up_Genes, pSUB_Up_Genes, LengthPreExcludeA, LengthPreExcludeB))
  }
}

#= GraphSetup() ================================================================

GraphSetup <- function(width=6, height=6, pointsize=12, 
                       ReportName="test", filetype="none", quality=75)
{
  # All pair wise combinations.
  PrintToFile <- FALSE
  if (filetype == "png")
  {
    PrintToFile <- TRUE
    file <- paste(ReportName, ".png", sep="")
    # Write plots as image files.
    png(filename=file, width=96*width, height=96*height, pointsize=pointsize, 
        bg="white")
  } 
  
  if (filetype == "wmf")
  {
    PrintToFile <- TRUE
    file <- paste(ReportName, ".wmf", sep="")
    # Write plots as windows meta files.
    win.metafile(filename=file, width=width, height=height, pointsize=pointsize)
  } 
  
  if (filetype == "postscript")
  {
    PrintToFile <- TRUE
    file <- paste(ReportName, ".eps", sep="")
    # Write plots as postscript files.
    postscript(file=file, width=width, height=height, pointsize=pointsize,
               onefile=FALSE, horizontal=FALSE, paper="special", family="Times")
  } 
  
  if (filetype == "pdf")
  {
    PrintToFile = TRUE
    file = paste(ReportName, ".pdf", sep="")
    # Write plots as pdf files.
    pdf(file=file, width=width, height=height, pointsize=pointsize,
        onefile=FALSE, paper="special", family="Times")
  } 
  
  if (filetype == "cairopdf")
  {
    PrintToFile = TRUE
    file = paste(ReportName, ".pdf", sep="")
    # Write plots as pdf files.
    CairoPDF(file = paste(ReportName, ".pdf", sep=""),
             width=width, height=height, onefile=TRUE, family = "Arial",
             title="R Graphics Output", fonts=NULL, version="1.5",
             paper="special", encoding, bg="white", fg, pointsize=pointsize, pagecentre)
  } 
  
  if (filetype == "cairojpeg")
  {
    PrintToFile = TRUE
    # Write plots as jpeg files.
    CairoJPEG(file = paste(ReportName, ".jpeg", sep=""),
              width=96*width, height=96*height, pointsize=pointsize,
              quality = quality, bg = "white", res = NA)
  } 
  
  if (filetype == "cairosvg")
  {
    PrintToFile = TRUE
    # Write plots as SVG files.
    CairoSVG(file = paste(ReportName, ".svg", sep=""),
             width=width, height=height, pointsize=pointsize,
             onefile=TRUE, bg="white")
  } 
  
  if (filetype == "svg")
  {
    #install.packages("RSvgDevice")
    library("RSvgDevice")
    
    PrintToFile = TRUE
    file = paste(ReportName, ".svg", sep="")
    # Write plots as pdf files.
    devSVG(file=file, width=width, height=height, onefile=FALSE)
  } 
  
  return(PrintToFile)
}

#= ApplyPurIST() ===============================================================

ApplyPurIST = function(ScriptsPath, Data)
{
  ### Calculate PurIST scores based on example that Hannah set up. ###
  
  # Author: Hannah Holly
  # Contact: holly@ohsu.edu
  # Based on instructions provided by Naim Rashid
  # PurIST github repo: https://github.com/naimurashid/PurIST
  # PATENT PENDING; PAPER ACCEPTED TO AACR.
  # Do not distribute without citing yeh, rashid, and moffit's new paper.
  # Do not use for profit.
  
  # load classifier objects
  load(paste(ScriptsPath, "PurIST-master/fitteds_public_2019-02-12.Rdata", sep=""))
  
  # Extract classifier object of interest
  classifier = classifs[[1]]
  
  # source the helper functions
  source(paste(ScriptsPath, "PurIST-master/functions.R", sep=""))
  
  # apply classifier to our samples
  predictions = apply_classifier(data = Data, classifier = classifier)
  return(data.frame(PurIST_Subtype = as.character(predictions$Subtype),
                    PurIST_Detailed_Subtype = as.character(predictions$Subtype_graded),
                    PurIST_Score = predictions$Pred_prob_basal,
                    stringsAsFactors = FALSE))
}

#= MakeGeneList() ==============================================================

MakeGeneList = function(GeneSets)
{
  gs_list = list(NULL)
  gs_list_idx = 0
  
  SetNames = names(GeneSets)
  
  for (idx in 1:length(SetNames))
  {
    gs = GeneSet(GeneSets[[SetNames[idx]]])
    setName(gs) = paste(SetNames[idx], sep="")
    geneIdType(gs) = SymbolIdentifier()
    collectionType(gs) = BroadCollection(category="c2")
    gs_list_idx = gs_list_idx + 1
    gs_list[[gs_list_idx]] = gs
  }
  
  return(GeneSetCollection(gs_list))
}

#= LoadGeneSets() ==============================================================

LoadGeneSets = function()
{
  # Get human collection of gene sets.
  m_df_all = msigdbr(species = "Homo sapiens")
  
  # Grep out just the HALLMARKS.
  SelectedIdx = (1:length(m_df_all$gs_name))[grep("^HALLMARK_", m_df_all$gs_name)]
  
  # Now convert these gene sets into a GSVA compatible GeneSetCollection object.
  m_df_rle = rle(m_df_all$gs_name[SelectedIdx])
  idx_start = 1
  gs_list = list(NULL)
  for (i in 1:length(m_df_rle$values))
  {
    idx = idx_start:(idx_start + m_df_rle$lengths[i] - 1)
    idx_start = idx_start + m_df_rle$lengths[i]
    gs = GeneSet(unique(m_df_all$human_gene_symbol[SelectedIdx[idx]]))
    setName(gs) = m_df_rle$values[i]
    geneIdType(gs) = SymbolIdentifier()
    collectionType(gs) = BroadCollection(category="h")
    gs_list[[i]] = gs
  }
  
  gs_collection = GeneSetCollection(gs_list)
  return(gs_collection)
}

#= ApplyGSVA() =================================================================

ApplyGSVA = function(M, gs_collection, RowCV_ThreshFactor = 0.5)
{
  # Coefficient of Variation (CV).
  RowMax = apply(M, MARGIN = 1, FUN = max)
  RowMean = apply(M, MARGIN = 1, FUN = mean)
  RowSd = apply(M, MARGIN = 1, FUN = sd)
  RowCV = RowMean
  RowCV[RowMean != 0] = RowSd[RowMean != 0] / RowMean[RowMean != 0]
  
  RowCV_Thresh = RowCV_ThreshFactor * median(RowCV)
  length(RowCV[RowCV >= RowCV_Thresh])
  
  M_GSVA = gsva(gsvaParam(exprData = M[RowCV >= RowCV_Thresh, ],
                          geneSets = gs_collection,
                          minSize = 3,
                          maxSize = 5000,
                          kcdf="Gaussian",
                          tau = 1,
                          maxDiff = TRUE,
                          absRanking = FALSE), 
                verbose=TRUE)
  
  #M_GSVA = gsva(M[RowCV >= RowCV_Thresh, ], gs_collection,
  #              min.sz=3, max.sz=5000, verbose=TRUE)
  
  return(t(M_GSVA))
}

#= Load_ICGC_Data() ============================================================
Load_ICGC_Data = function(SourcePath)
{
  # Load 96 ICGC PDAC samples.
  ICGC_Data = read.table(paste(SourcePath, "PublicDataUsed/", "nature16965-s2-PDAC_normalised_exp.tsv", sep=""),  
                         sep="\t", quote="", comment.char = "", header=TRUE, row.names=1, stringsAsFactors=FALSE)
  dim(ICGC_Data)
  names(ICGC_Data)
  ICGC_Data <- as.matrix(ICGC_Data)
  
  ICGCMetadata = read.table(paste(SourcePath, "data/", "ICGC_APGI_Survival_Metadata.tsv", sep=""),  
                            sep="\t", quote="", comment.char = "", header=TRUE, stringsAsFactors=FALSE)
  dim(ICGCMetadata)
  names(ICGCMetadata)
  
  idx = match(x=colnames(ICGC_Data), table=ICGCMetadata$ICGC_ID)
  ICGC_Data = ICGC_Data[, !is.na(idx)]
  ICGCMetadata = ICGCMetadata[idx[!is.na(idx)], ]
  
  return(list(Data = ICGC_Data,
              Metadata = ICGCMetadata))
}

#= Load_TCGA_Data() ============================================================

Load_TCGA_Data = function(SourcePath)
{
  # Get TCGA PDAC data from cBioPortal download (paad_tcga_pan_can_atlas_2018).
  Data = 
    read.table(paste(SourcePath, "PublicDataUsed/", 
                     "data_RNA_Seq_v2_expression_median.txt", 
                     sep=""),
               sep="\t", quote="", comment.char = "", header=TRUE)
  
  dim(Data)
  head(colnames(Data))
  
  ### TCGA samples to use from Johnathan. ###
  
  TCGA_ListToUse = 
    read.table(paste(SourcePath, "data/", 
                     "TCGA_ListToUse.tsv", 
                     sep=""),
               sep="\t", quote="", comment.char = "", header=TRUE)
  
  dim(TCGA_ListToUse)
  head(colnames(TCGA_ListToUse))
  
  # There are duplicate gene symbols as expected, but there are also duplicate 
  # Entrez IDs.  So I will collaplse Entrez IDs, then gene symbols.
  Symbols = as.character(Data$Entrez_Gene_Id)
  Symbols[is.na(Symbols)] = "NA"
  M = as.matrix(Data[, -c(1:2)])
  CollapseIdx = GetCollapseIndex(M, Symbols)
  length(CollapseIdx)
  idx = c(grep("^NA$", Symbols[CollapseIdx]), grep("^$", Symbols[CollapseIdx]))
  if (length(idx) > 0) { CollapseIdx = CollapseIdx[-idx] }
  length(CollapseIdx)
  Symbols = as.character(Data$Hugo_Symbol[CollapseIdx])
  Symbols[is.na(Symbols)] = "NA"
  M = as.matrix(Data[CollapseIdx, -c(1:2)])
  CollapseIdx = GetCollapseIndex(M, Symbols)
  length(CollapseIdx)
  idx = c(grep("^NA$", Symbols[CollapseIdx]), grep("^$", Symbols[CollapseIdx]))
  if (length(idx) > 0) { CollapseIdx = CollapseIdx[-idx] }
  length(CollapseIdx)
  M = M[CollapseIdx, ]
  rownames(M) = Symbols[CollapseIdx]
  
  RnaSeqMedianData = M
  
  Metadata = 
    read.table(paste(SourcePath, "PublicDataUsed/", 
                     "paad_tcga_pan_can_atlas_2018_clinical_data.tsv", 
                     sep=""),
               sep="\t", quote="", comment.char = "", header=TRUE)
  
  dim(Metadata)
  head(colnames(Metadata))
  
  idx = which(Metadata$Patient.ID %in% TCGA_ListToUse$bcr_patient_barcode)
  Metadata = Metadata[idx, ]
  
  # Metadata does not line up with RnaSeqData, so match it now.
  ModifiedSampleNames = gsub("-", ".", Metadata$Sample.ID)
  idx = match(x=colnames(RnaSeqMedianData), table=ModifiedSampleNames)
  colnames(RnaSeqMedianData)[is.na(idx)]
  
  RnaSeqMedianData = RnaSeqMedianData[, which(!is.na(idx))]
  Metadata = Metadata[idx[!is.na(idx)], ]
  rownames(Metadata) = ModifiedSampleNames[idx[!is.na(idx)]]
  
  dim(RnaSeqMedianData)
  dim(Metadata)
  
  return(list(Data = RnaSeqMedianData,
              Metadata = Metadata))
}

#= tximport_Example() ==========================================================

tximport_Example = function(PathToH5Files)
{
  # Load h5 files from kallisto pipeline to use for TMM normalization.
  library(tximport)
  library(rhdf5)
  library(edgeR)
  
  # Import from h5 files which should be faster.
  Samples = list.dirs(PathToH5Files, full.names = FALSE, recursive = FALSE)
  # The following is only used if ordering samples to match RnaSeqCounts.
  #idx = match(x=colnames(RnaSeqCounts), table=Samples)
  #Samples = Samples[idx]
  h5Files = paste(PathToH5Files, Samples, "/abundance.h5", sep="")
  
  # Use cool tricks from Hannah's code to split by "|" and build data frame from
  # first two elements of each row which are the ensembl transcript ID and 
  # ensembl gene ID (many to 1 mapping that defines transcripts for each gene).
  #h5ls(file=h5Files[1])
  X = h5read(file=h5Files[1], name="/aux/ids")
  X = data.frame(t(sapply(strsplit(as.character(X),'\\|'),'['))) 
  tx2geneMap = X[1:2]
  names(tx2geneMap) = c('transcript', 'gene')
  
  RnaSeq_tximport = tximport(files = h5Files, 
                             type = "kallisto", 
                             tx2gene = tx2geneMap, txOut = FALSE,
                             ignoreAfterBar=TRUE)
  
  colnames(RnaSeq_tximport$abundance) = Samples
  colnames(RnaSeq_tximport$counts) = Samples
  colnames(RnaSeq_tximport$length) = Samples
  
  Ensembl_ID = rownames(RnaSeq_tximport$counts)
  # Remove ".version" numbers from Ensembl IDs if they are there.
  idx = 2 * 1:length(Ensembl_ID) - 1
  Ensembl_ID = unlist(strsplit(Ensembl_ID, split=".", fixed=TRUE))[idx]
  
  rownames(RnaSeq_tximport$abundance) = Ensembl_ID
  rownames(RnaSeq_tximport$counts) = Ensembl_ID
  rownames(RnaSeq_tximport$length) = Ensembl_ID
  
  # Use edgeR to calculate cpms from this txi object.
  # This is the TMM method.
  cts <- RnaSeq_tximport$counts
  normMat <- RnaSeq_tximport$length
  
  # Obtaining per-observation scaling factors for length, adjusted to avoid
  # changing the magnitude of the counts.
  normMat <- normMat/exp(rowMeans(log(normMat)))
  normCts <- cts/normMat
  
  # Computing effective library sizes from scaled counts, to account for
  # composition biases between samples.
  eff.lib <- calcNormFactors(normCts) * colSums(normCts)
  
  # Combining effective library sizes with the length factors, and calculating
  # offsets for a log-link GLM.
  normMat <- sweep(normMat, 2, eff.lib, "*")
  normMat <- log(normMat)
  
  # Creating a DGEList object for use in edgeR.
  RnaSeqDGEList <- DGEList(cts)
  RnaSeqDGEList <- scaleOffset(RnaSeqDGEList, normMat)
  dim(RnaSeqDGEList$counts)
  
  # y is now ready for estimate dispersion functions see edgeR User's Guide
  # For creating a matrix of CPMs within edgeR, the following code chunk can be used:
  edgeR_TMM <- edgeR::cpm(RnaSeqDGEList, offset = RnaSeqDGEList$offset, log = FALSE)
  dim(edgeR_TMM)
  
  # filtering using the design information
  # Do not use a design here.
  design = NULL
  keep <- filterByExpr(RnaSeqDGEList, design)
  length(which(keep))
  RnaSeqDGEList <- RnaSeqDGEList[keep, ]
  # y is now ready for estimate dispersion functions see edgeR User's Guide
  # For creating a matrix of CPMs within edgeR, the following code chunk can be used:
  edgeR_TMM_Filtered <- edgeR::cpm(RnaSeqDGEList, offset = RnaSeqDGEList$offset, log = FALSE)
  dim(edgeR_TMM_Filtered)
  
  return(list(edgeR_TMM = edgeR_TMM,
              edgeR_TMM_Filtered = edgeR_TMM_Filtered))
}

#= Viper_ANOVA() ===============================================================

Viper_ANOVA = function(BCC_RNA_Data)
{
  # Use viper all primary matrix and test significant regulons between top/bottom
  # 1/4 by pORG_Up_55, pSUB_Up_51, and PurIST score.
  
  ViperMetaData = BCC_RNA_Data$Metadata
  
  ViperData = read.table(paste(SourcePath, "OriginalData/ViperData.tsv", sep=""),  
                         sep="\t", quote="", comment.char = "", header=TRUE, 
                         row.names=1, check.names=FALSE, stringsAsFactors=FALSE)    
  
  idx = match(x=colnames(ViperData), table=ViperMetaData$Public_Specimen_ID)
  # Are there any samples missing from metadata?
  colnames(ViperData)[is.na(idx)]
  ViperMetaData = ViperMetaData[idx, ]
  
  idx = match(x=colnames(ViperData), table=rownames(BCC_RNA_Data$GSVA_pORGpSUB_Primaries))
  # Are there any samples missing from this metadata?
  colnames(ViperData)[is.na(idx)]
  ViperGSVA_pORGpSUB_Primaries = BCC_RNA_Data$GSVA_pORGpSUB_Primaries[idx, ]
  
  Number = floor(dim(ViperData)[[2]] / 4) 
  idx = order(ViperGSVA_pORGpSUB_Primaries[, "pORG_Up_55"], decreasing=TRUE)
  TopQuarterIdx = idx[1:Number]
  idx = order(ViperGSVA_pORGpSUB_Primaries[, "pORG_Up_55"], decreasing=FALSE)
  BottomQuarterIdx = idx[1:Number]
  
  Values = ViperData[, c(TopQuarterIdx, BottomQuarterIdx)]
  Factors = as.factor(c(rep("N", times=length(TopQuarterIdx)), rep("D", times=length(BottomQuarterIdx))))
  
  fm = aov(t(Values) ~ Factors, data=data.frame(t(Values),Factors))
  AOV_pVal = unlist(summary(fm))[1:dim(Values)[[1]] * 10 - 1]
  AOV_pVal[is.na(AOV_pVal)] = 1
  AOV_qVal = fdrtool(AOV_pVal, statistic="pvalue")$qval
  print(length(AOV_pVal[AOV_pVal < 0.01]))
  print(length(AOV_qVal[AOV_qVal < 0.01]))
  
  Report = cbind(regulons = rownames(ViperData), 
                 TopBottomQuarter_pORG_Up_55_pVal = AOV_pVal, 
                 TopBottomQuarter_pORG_Up_55_qVal = AOV_qVal)
  
  idx = order(ViperGSVA_pORGpSUB_Primaries[, "pSUB_Up_51"], decreasing=TRUE)
  TopQuarterIdx = idx[1:Number]
  idx = order(ViperGSVA_pORGpSUB_Primaries[, "pSUB_Up_51"], decreasing=FALSE)
  BottomQuarterIdx = idx[1:Number]
  
  Values = ViperData[, c(TopQuarterIdx, BottomQuarterIdx)]
  Factors = as.factor(c(rep("N", times=length(TopQuarterIdx)), rep("D", times=length(BottomQuarterIdx))))
  
  fm = aov(t(Values) ~ Factors, data=data.frame(t(Values),Factors))
  AOV_pVal = unlist(summary(fm))[1:dim(Values)[[1]] * 10 - 1]
  AOV_pVal[is.na(AOV_pVal)] = 1
  AOV_qVal = fdrtool(AOV_pVal, statistic="pvalue")$qval
  print(length(AOV_pVal[AOV_pVal < 0.01]))
  print(length(AOV_qVal[AOV_qVal < 0.01]))
  
  Report = cbind(Report, 
                 TopBottomQuarter_pSUB_Up_51_pVal = AOV_pVal, 
                 TopBottomQuarter_pSUB_Up_51_qVal = AOV_qVal)
  Number = floor(dim(ViperData)[[2]] / 4) 
  idx = order(ViperGSVA_pORGpSUB_Primaries[, "pORG_Up_55"], decreasing=TRUE)
  TopQuarterIdx = idx[1:Number]
  idx = order(ViperGSVA_pORGpSUB_Primaries[, "pORG_Up_55"], decreasing=FALSE)
  BottomQuarterIdx = idx[1:Number]
  
  MeanDiff = apply(ViperData[, TopQuarterIdx], MARGIN=1, FUN=mean) - 
    apply(ViperData[, BottomQuarterIdx], MARGIN=1, FUN=mean)
  
  Values = ViperData[, c(TopQuarterIdx, BottomQuarterIdx)]
  Factors = as.factor(c(rep("N", times=length(TopQuarterIdx)), rep("D", times=length(BottomQuarterIdx))))
  
  fm = aov(t(Values) ~ Factors, data=data.frame(t(Values),Factors))
  AOV_pVal = unlist(summary(fm))[1:dim(Values)[[1]] * 10 - 1]
  AOV_pVal[is.na(AOV_pVal)] = 1
  AOV_qVal = fdrtool(AOV_pVal, statistic="pvalue")$qval
  print(length(AOV_pVal[AOV_pVal < 0.01]))
  print(length(AOV_qVal[AOV_qVal < 0.01]))
  
  Report = cbind(regulons = rownames(ViperData),
                 MeanDiff_pORG_Up_55 = MeanDiff,
                 TopBottomQuarter_pORG_Up_55_pVal = AOV_pVal, 
                 TopBottomQuarter_pORG_Up_55_qVal = AOV_qVal)
  
  
  idx = order(ViperGSVA_pORGpSUB_Primaries[, "pSUB_Up_51"], decreasing=TRUE)
  TopQuarterIdx = idx[1:Number]
  idx = order(ViperGSVA_pORGpSUB_Primaries[, "pSUB_Up_51"], decreasing=FALSE)
  BottomQuarterIdx = idx[1:Number]
  
  MeanDiff = apply(ViperData[, TopQuarterIdx], MARGIN=1, FUN=mean) - 
    apply(ViperData[, BottomQuarterIdx], MARGIN=1, FUN=mean)
  
  Values = ViperData[, c(TopQuarterIdx, BottomQuarterIdx)]
  Factors = as.factor(c(rep("N", times=length(TopQuarterIdx)), rep("D", times=length(BottomQuarterIdx))))
  
  fm = aov(t(Values) ~ Factors, data=data.frame(t(Values),Factors))
  AOV_pVal = unlist(summary(fm))[1:dim(Values)[[1]] * 10 - 1]
  AOV_pVal[is.na(AOV_pVal)] = 1
  AOV_qVal = fdrtool(AOV_pVal, statistic="pvalue")$qval
  print(length(AOV_pVal[AOV_pVal < 0.01]))
  print(length(AOV_qVal[AOV_qVal < 0.01]))
  
  Report = cbind(Report, 
                 MeanDiff_pSUB_Up_51 = MeanDiff,
                 TopBottomQuarter_pSUB_Up_51_pVal = AOV_pVal, 
                 TopBottomQuarter_pSUB_Up_51_qVal = AOV_qVal)
  
  ViperLiverCohort = which(!is.na(ViperMetadata$Liver_Met_Present) & 
                           ViperMetadata$Liver_Met_Present == 1)
  ViperLungCohort = which(!is.na(ViperMetadata$Lung_Met_Present) & 
                          ViperMetadata$Lung_Met_Present == 1 &
                          (is.na(ViperMetadata$Liver_Met_Present) |
                           ViperMetadata$Liver_Met_Present == "NA"))

  MeanDiff = apply(ViperData[, ViperLiverCohort], MARGIN=1, FUN=mean) - 
    apply(ViperData[, ViperLungCohort], MARGIN=1, FUN=mean)
  
  Values = ViperData[, c(ViperLiverCohort, ViperLungCohort)]
  Factors = as.factor(c(rep("N", times=length(ViperLiverCohort)), rep("D", times=length(ViperLungCohort))))
  
  fm = aov(t(Values) ~ Factors, data=data.frame(t(Values),Factors))
  AOV_pVal = unlist(summary(fm))[1:dim(Values)[[1]] * 10 - 1]
  AOV_pVal[is.na(AOV_pVal)] = 1
  AOV_qVal = fdrtool(AOV_pVal, statistic="pvalue")$qval
  print(length(AOV_pVal[AOV_pVal < 0.01]))
  print(length(AOV_qVal[AOV_qVal < 0.01]))
  
  Report = cbind(Report, 
                 MeanDiff_LiverVsLungNotLiver = MeanDiff,
                 LiverVsLungNotLiver_pVal = AOV_pVal, 
                 LiverVsLungNotLiver_qVal = AOV_qVal)
  
  idx = order(ViperMetaData[, "PurIST_Score"], decreasing=TRUE)
  TopQuarterIdx = idx[1:Number]
  idx = order(ViperMetaData[, "PurIST_Score"], decreasing=FALSE)
  BottomQuarterIdx = idx[1:Number]
  
  MeanDiff = apply(ViperData[, TopQuarterIdx], MARGIN=1, FUN=mean) - 
    apply(ViperData[, BottomQuarterIdx], MARGIN=1, FUN=mean)
  
  Values = ViperData[, c(TopQuarterIdx, BottomQuarterIdx)]
  Factors = as.factor(c(rep("N", times=length(TopQuarterIdx)), rep("D", times=length(BottomQuarterIdx))))
  
  fm = aov(t(Values) ~ Factors, data=data.frame(t(Values),Factors))
  AOV_pVal = unlist(summary(fm))[1:dim(Values)[[1]] * 10 - 1]
  AOV_pVal[is.na(AOV_pVal)] = 1
  AOV_qVal = fdrtool(AOV_pVal, statistic="pvalue")$qval
  print(length(AOV_pVal[AOV_pVal < 0.01]))
  print(length(AOV_qVal[AOV_qVal < 0.01]))
  
  Report = cbind(Report, 
                 MeanDiff_PurIST_Score = MeanDiff,
                 TopBottomQuarter_PurIST_Score_pVal = AOV_pVal, 
                 TopBottomQuarter_PurIST_Score_qVal = AOV_qVal)
  
  # Write some reports that can go into supp data sets.
  write.table(Report, 
              paste(AnalysisOutputPath, "Viper_ANOVA.tsv", sep=""), 
              sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  return(Report)
}

# Example to generate viper regulon statistics used in Fig...
# TestResult = Viper_ANOVA(BCC_RNA_Data)

#= FormatInputForGSEA() ========================================================

FormatInputForGSEA = function(Data,
                              SetA,
                              SetAName,
                              SetB,
                              SetBName,
                              OutputPath = "./",
                              FileNameForGSEA = "GSEA_Example")
{
  # *.gct file.
  TmpFileName = paste(OutputPath, FileNameForGSEA, ".gct", sep="")
  
  LineData = NULL
  OneLine = paste("#1.2", sep="\t")
  LineData = c(LineData, OneLine)
  OneLine = paste(dim(Data)[[1]], dim(Data)[[2]], sep="\t")
  LineData = c(LineData, OneLine)
  
  OneLine = paste("NAME", "Description", sep="\t")
  for(SampleName in colnames(Data))
  {
    OneLine = paste(OneLine, SampleName, sep="\t")
  }
  LineData = c(LineData, OneLine)
  
  Output = cbind(rownames(Data),
                 rownames(Data),
                 Data)
  
  LineData = c(LineData, apply(Output, MARGIN=1, FUN=paste, collapse="\t"))
  writeLines(LineData, TmpFileName)
  
  # *.cls file.
  TmpFileName = paste(OutputPath, FileNameForGSEA, ".cls", sep="")
  
  LineData = NULL
  OneLine = paste(length(c(SetA, SetB)), 2, 1, sep="\t")
  LineData = c(LineData, OneLine)
  OneLine = paste("#", SetAName, SetBName, sep="\t")
  LineData = c(LineData, OneLine)
  OneLine = paste(paste(rep(SetAName, times=length(SetA)), collapse="\t"),
                  paste(rep(SetBName, times=length(SetB)), collapse="\t"), sep="\t")
  LineData = c(LineData, OneLine)
  writeLines(LineData, TmpFileName)
}

#= CacheEnsemblAnnotion() ======================================================

CacheEnsemblAnnotion = function(RnaSeqData,
                                SourcePath,
                                Species = "human",
                                Version = "Current")
{
  CachePath = paste(SourcePath, sep="")
  CacheFile = paste("BioMartAnn_", Species, "_", Version, ".RData", sep="")
  
  Target_Ensembl_ID = rownames(RnaSeqData)
  
  # Remove ".version" numbers from Ensembl IDs if they are there.
  temp = unlist(strsplit(Target_Ensembl_ID, split=".", fixed=TRUE))
  if (length(temp) == 2*length(Target_Ensembl_ID))
  {
    idx = 2 * 1:length(Target_Ensembl_ID) - 1
    Target_Ensembl_ID = temp[idx]
  }
  
  if (length(list.files(CachePath, 
                        pattern=paste("^", CacheFile, "$", sep=""))) == 1) 
  {
    load(paste(CachePath, CacheFile, sep=""))
    
    # Match the ensembl IDs RnaSeqData with the ensembl IDs
    # of the annotation just loaded.
    idx = match(x=Target_Ensembl_ID, table=EnsemblAnn[, "ensembl_gene_id"])
    EnsemblAnn = EnsemblAnn[idx, ]
    
  } else
  {
    if (Version == "Current")
    {
      # Default, current version.
      ensembl = useEnsembl(biomart="Ensembl", dataset="hsapiens_gene_ensembl")
    } else
    {
      # Use listEnsemblArchives() for list of versions.
      ensembl = useEnsembl(biomart="Ensembl", dataset="hsapiens_gene_ensembl", version=Version)
    }
    
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    require("BiocManager", quietly = TRUE)
    
    if (!require("org.Hs.eg.db", quietly = TRUE))
      BiocManager::install("org.Hs.eg.db")    # See listAttributes(ensembl)
    require("org.Hs.eg.db", quietly = TRUE)
    
    if (!require("biomaRt", quietly = TRUE))
      BiocManager::install("biomaRt")    # See listAttributes(ensembl)
    require("biomaRt", quietly = TRUE)
    
    EnsemblAnnDB = getBM(attributes=c('entrezgene_id',
                                      'ensembl_gene_id', 
                                      'hgnc_symbol', # 'hgnc_symbol', 
                                      # 'mgi_symbol'
                                      'gene_biotype',
                                      'chromosome_name', 
                                      'start_position', 
                                      'end_position'), 
                         mart = ensembl)
    
    # This is needed to get nice gene descriptions...
    # Some examples: 
    # https://rdrr.io/bioc/GenomicFeatures/man/select-methods.html
    #columns(org.Mm.eg.db)
    # Note that this gives a 1:many mapping.
    ExtraAnnDB = select(org.Hs.eg.db, Target_Ensembl_ID, 
                        columns=c("ENSEMBL", "SYMBOL", "GENENAME"), 
                        keytype="ENSEMBL")
    
    save(EnsemblAnnDB, ExtraAnnDB, 
         file=paste(CachePath, CacheFile, sep=""))
  }
  
  # Match the cleaned ensembl IDs from RnaSeqData with the ensembl IDs
  # of the annotation just loaded.
  idx = match(x=Target_Ensembl_ID, table=EnsemblAnnDB[, "ensembl_gene_id"])
  EnsemblAnn = EnsemblAnnDB[idx, ]
  
  ExtraAnn = ExtraAnnDB[!duplicated(ExtraAnnDB$ENSEMBL), ]
  dim(ExtraAnn)
  idx = match(x=Target_Ensembl_ID, table=ExtraAnn[, "ENSEMBL"])
  ExtraAnn = ExtraAnn[idx, ]
  
  return(list(Ensembl_ID = Target_Ensembl_ID,
              EnsemblAnn = EnsemblAnn,
              ExtraAnn = ExtraAnn))
}

#= END =========================================================================