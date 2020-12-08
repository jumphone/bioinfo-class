
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html
# Install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")


#####################################
# Install ComplexHeatmap and dependencies...
install.packages('circlize')
install.packages('seriation')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
######################################


######################################
# Install GSEA

# Please check https://www.gsea-msigdb.org/gsea/index.jsp
