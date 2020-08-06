if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
getGDCInfo()

# 1. Search data in GDC
query = GDCquery(project="TCGA-LIHC",
                 data.category="Transcriptome Profiling",
                 data.type="Gene Expression Quantification",
                 workflow.type="HTSeq - Counts")
# 2. Download from GDC repository
GDCdownload(query)

# 3. Make R object from the downloaded data
data = GDCprepare(query)
# 4. Extract Gene expression matrix
library(SummarizedExperiment)
eset = assay(data)
# 5. Save the matrix as .csv format
write.csv(eset, file="GE.csv")
