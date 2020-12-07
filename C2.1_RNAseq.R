##########################
# Set working directory
setwd('C:/Users/user/Desktop')


###########################
# Load Data
ORIG_MAT=read.table('COUNT_MAT.txt',row.names=1,header=TRUE,sep='\t')


###########################
# Check header
head(ORIG_MAT)


###########################
# Check size 
dim(ORIG_MAT)


###########################
GENE_LENGTH=ORIG_MAT[,1]
MAT_COUNT=ORIG_MAT[,2:ncol(ORIG_MAT)]


##############################
# Remove non-expressed genes
R_SUM=apply(MAT_COUNT,1,sum)
GENE_LENGTH=GENE_LENGTH[which(R_SUM>0)]
MAT_COUNT=MAT_COUNT[which(R_SUM>0),]
dim(MAT_COUNT)


##################################
# Normalization (归一化)
#####
.normRPM <- function(x){
    x.sum=sum(x)
    if(x.sum==0){
        y=x
    }else{
        y=x / x.sum *1000000
        }
    return(y)
    }

.normRPKM <- function(x, L){
    x_sum=sum(x)
    #print(x_sum)
    if(x_sum==0){
        y=x
    }else{
        y=x  / L * 1000 / x_sum * 1000000
        }
    return(y)
    }
#####

MAT_RPM=apply(MAT_COUNT, 2, .normRPM)
MAT_RPKM=apply(MAT_COUNT, 2, .normRPKM, GENE_LENGTH)

sd(MAT_RPM[,1])
sd(MAT_RPKM[,1])
sd(MAT_RPM[1,])
sd(MAT_RPKM[1,])
########################


##################################
# Log

#####
MAT_RPKM_1=MAT_RPKM+1
MAT_RPKM_LOG=apply(MAT_RPKM_1, 2, log, 10)

par(mfrow=c(2,2))
plot(density(MAT_RPKM[,1]))
plot(density(MAT_RPKM_LOG[,1]))
plot(density(MAT_RPKM[1,]))
plot(density(MAT_RPKM_LOG[1,]))
par(mfrow=c(1,1))
########################



#######################################
# Standardization (标准化)

####
MAT_RPKM_STD=t(apply(MAT_RPKM, 1, scale))
MAT_RPKM_STD[1:5,1:5]
colnames(MAT_RPKM_STD)=colnames(MAT_RPKM)
rownames(MAT_RPKM_STD)=rownames(MAT_RPKM)

MAT_RPKM_LOG_STD=t(apply(MAT_RPKM_LOG, 1, scale))
MAT_RPKM_LOG_STD[1:5,1:5]
colnames(MAT_RPKM_LOG_STD)=colnames(MAT_RPKM)
rownames(MAT_RPKM_LOG_STD)=rownames(MAT_RPKM)

par(mfrow=c(2,2))
plot(density(MAT_RPKM_STD[,1]))
plot(density(MAT_RPKM_LOG_STD[,1]))
plot(density(MAT_RPKM_STD[1,]))
plot(density(MAT_RPKM_LOG_STD[1,]))
par(mfrow=c(1,1))



############################################
# PCA
pca=prcomp(t(MAT_RPKM_LOG_STD),center=FALSE,scale=FALSE)

plot(pca$x[,1],pca$x[,2],xlim=c(-140,160),ylim=c(-100,100))
set.seed(123)
text(label=rownames(pca$x),x=pca$x[,1], y=pca$x[,2], pos=sample(c(1,2,3,4),nrow(pca$x),replace = TRUE))


###########################
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html
# Install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")


#################################
# DESeq2
library(DESeq2)
mydata=MAT_COUNT[,1:4]
type<-factor(c(rep("D0",2),rep("D7",2)),levels = c("D0","D7"))
database<-round(as.matrix(mydata))
coldata<-data.frame(row.names = colnames(database),type)
dds<-DESeqDataSetFromMatrix(database,coldata,design = ~type)
dds<-DESeq(dds)
res = results(dds, contrast=c("type", "D0", "D7"))
res=as.matrix(res)

head(res)

mydata_rpkm=MAT_RPKM[,1:4]
colnames(mydata_rpkm)=paste0('RPKM_',colnames(mydata_rpkm))

###########################
all_res=cbind(mydata, mydata_rpkm, res)

head(all_res)

write.table(all_res,file='DESeq2_result.txt',sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)


###########################
LOG2FC=all_res$log2FoldChange
PADJ=all_res$padj



########################
# Draw Figure
plot(LOG2FC, -log(PADJ,10), pch=16, col='grey50', xlim=c(-15,15), ylim=c(0,300))
#######################
FOLD_CUT=5
PADG_CUT=0.01
#######################
abline(v=log(FOLD_CUT,2),lty=3,lwd=1)
abline(v=-log(FOLD_CUT,2),lty=3,lwd=1)
abline(h=-log(PADG_CUT,10),lty=3,lwd=1)
############################
select_point=which( abs(LOG2FC)>log(FOLD_CUT,2) & PADJ < PADG_CUT)
points(LOG2FC[select_point], -log(PADJ[select_point],10), pch=16,col='red')
#####################
text(x=12, y=200,labels='D0',cex=1.5)
text(x=-12, y=200,labels='D7',cex=1.5)



########################
pdf('f1_Volcano.pdf',width=6,height=6)
# Draw Figure
plot(LOG2FC, -log(PADJ,10), pch=16, col='grey50', xlim=c(-15,15), ylim=c(0,300))
#######################
FOLD_CUT=5
PADG_CUT=0.01
#######################
abline(v=log(FOLD_CUT,2),lty=3,lwd=1)
abline(v=-log(FOLD_CUT,2),lty=3,lwd=1)
abline(h=-log(PADG_CUT,10),lty=3,lwd=1)
############################
select_point=which( abs(LOG2FC)>log(FOLD_CUT,2) & PADJ < PADG_CUT)
points(LOG2FC[select_point], -log(PADJ[select_point],10), pch=16,col='red')
#####################
text(x=12, y=200,labels='D0',cex=1.5)
text(x=-12, y=200,labels='D7',cex=1.5)
dev.off()


#####################################
# Output result
D0_OUT=all_res[which( LOG2FC> log(FOLD_CUT,2) & PADJ < PADG_CUT),]
D7_OUT=all_res[which( LOG2FC< -log(FOLD_CUT,2) & PADJ < PADG_CUT),]


write.table(D0_OUT,file='D0_DESeq2_result.txt',sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)
write.table(D7_OUT,file='D7_DESeq2_result.txt',sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)



#####################################
# Install ComplexHeatmap ...
install packages('circlize')
install packages('seriation')

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
######################################


library('ComplexHeatmap')
library('circlize')
library('seriation')


HM_MAT=all_res[which( abs(LOG2FC) > log(FOLD_CUT,2) & PADJ < PADG_CUT),]
HM_MAT=HM_MAT[order(HM_MAT$log2FoldChange),]

USED_HM_MAT=HM_MAT[,5:8]
USED_HM_MAT_LOG=log(USED_HM_MAT+1,10)
USED_HM_MAT_LOG_STD=t(apply(USED_HM_MAT_LOG,1,scale))
colnames(USED_HM_MAT_LOG_STD)=colnames(USED_HM_MAT_LOG)



mat=USED_HM_MAT_LOG_STD
mat=as.matrix(mat)
colnames(mat)=c('D0R1','D0R2','D7R1','D7R2')
color_fun =colorRamp2(c(-1,-0.5,0,0.5,1 ), c('royalblue3','white','white','white','indianred3'))

SHOW_GENE= c('CGA','KRT23','LEFTY1','POU5F1')
SHOW_INDEX=which(rownames(mat) %in% SHOW_GENE)
SHOW_NAME=rownames(mat)[SHOW_INDEX]
ha = rowAnnotation(foo = anno_mark(at = SHOW_INDEX ,  labels =SHOW_NAME))

Heatmap(mat, row_title='', name="V",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE ,
         right_annotation = ha
	      )


pdf('f2_Heatmap.pdf',width=4,height=6)
Heatmap(mat, row_title='', name="V",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE ,
         right_annotation = ha
	      )
dev.off()

########################################



mat=USED_HM_MAT_LOG
mat=as.matrix(mat)
colnames(mat)=c('D0R1','D0R2','D7R1','D7R2')
color_fun =colorRamp2(c(0,1,3 ), c('royalblue3','yellow1','yellow3'))

SHOW_GENE= c('CGA','KRT23','LEFTY1','POU5F1')
SHOW_INDEX=which(rownames(mat) %in% SHOW_GENE)
SHOW_NAME=rownames(mat)[SHOW_INDEX]
ha = rowAnnotation(foo = anno_mark(at = SHOW_INDEX ,  labels =SHOW_NAME))

Heatmap(mat, row_title='', name="V",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=FALSE,
	      col=color_fun, border = TRUE ,
         right_annotation = ha
	      )





##############################################################################
# GSEA

.getGSEAinput <- function( DATA, TAG, PATH ){
    DATA=as.matrix(DATA)
    TAG=TAG
    VAR=apply(DATA,1,var)
    DATA=DATA[which(VAR>0),]
    ########################
    EXP.FILE=paste0(PATH,'.EXP.txt')
    colnames(DATA)=paste0(TAG,'_',colnames(DATA))
    OUT=cbind(toupper(rownames(DATA)),rep('NO',nrow(DATA)),DATA)
    colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
    write.table(OUT, EXP.FILE, sep='\t',quote=F,row.names=F,col.names=T)
    
    #########################
    PT.FILE=paste0(PATH,'.PT.cls')
    PT=t(as.character(TAG))
    cat(paste0(length(TAG),' 2 1'),file=PT.FILE,sep="\n") 
    cat(paste(c('#',unique(TAG)),collapse=' '),file=PT.FILE,sep='\n',append=TRUE)
    cat(PT,file=PT.FILE,sep=' ',append=TRUE)
    ##########################  
    }


.getGSEAinput(DATA=MAT_RPKM[,1:4], TAG=c('D0','D0','D7','D7'), PATH='./GSEA' )

# https://www.gsea-msigdb.org/gsea/index.jsp












