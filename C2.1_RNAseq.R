##########################
# Set working directory
setwd('C:/Users/user/Desktop')


###########################
# Load Data
COUNT_MAT=read.table('COUNT_MAT.txt',row.names=1,header=TRUE,sep='\t')


###########################
# Check header
head(COUNT_MAT)


###########################
# Check size 
dim(COUNT_MAT)


###########################
GENE_LENGTH=COUNT_MAT[,1]
MAT_COUNT=COUNT_MAT[,2:ncol(COUNT_MAT)]


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







