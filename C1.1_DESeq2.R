##########################
# Set working directory
setwd('C:/Users/user/Desktop')

###########################
# Load Data
COUNT_MAT=read.table('ALL.SE.MAT.txt',row.names=1,header=TRUE,sep='\t')

###########################
# Check header
head(COUNT_MAT)

###########################
# Check size 
dim(COUNT_MAT)

###########################
# DESeq2
library(DESeq2)
mydata=COUNT_MAT
type<-factor(c(rep("D0",2),rep("D7",2)),levels = c("D0","D7"))
database<-round(as.matrix(mydata))
coldata<-data.frame(row.names = colnames(database),type)
dds<-DESeqDataSetFromMatrix(database,coldata,design = ~type)
dds<-DESeq(dds)
res = results(dds, contrast=c("type", "D0", "D7"))
res=as.matrix(res)

###########################
# Get RPM
.norm <- function(x){
    y=x/sum(x)*1000000
    return(y)
    }

NORM_COUNT_MAT=apply(COUNT_MAT, 2,.norm)
colnames(NORM_COUNT_MAT)=paste0('RPM_',colnames(COUNT_MAT))

###########################
all_res=cbind(COUNT_MAT, NORM_COUNT_MAT, res)
write.table(all_res,file='DESeq2_result.txt',sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)


###########################
LOG2FC=all_res$log2FoldChange
PVALUE=all_res$pvalue



########################
# Draw Figure
plot(LOG2FC, -log(PVALUE,10), pch=16, col='grey50', xlim=c(-1,1), ylim=c(0,3))
#######################
abline(v=log(1.2,2),lty=3,lwd=1)
abline(v=-log(1.2,2),lty=3,lwd=1)
abline(h=-log(0.05,10),lty=3,lwd=1)
############################
select_point=which( abs(LOG2FC)>log(1.2,2) & PVALUE < 0.05)
points(LOG2FC[select_point], -log(PVALUE[select_point],10), pch=16,col='red')
#####################
text(x=0.5, y=0.5,labels='D0',cex=1.5)
text(x=-0.5, y=0.5,labels='D7',cex=1.5)



########################
pdf('f1_Volcano.pdf',width=6,height=6)
# Draw Figure
plot(LOG2FC, -log(PVALUE,10), pch=16, col='grey50', xlim=c(-1,1), ylim=c(0,3))
#######################
abline(v=log(1.2,2),lty=3,lwd=1)
abline(v=-log(1.2,2),lty=3,lwd=1)
abline(h=-log(0.05,10),lty=3,lwd=1)
############################
select_point=which( abs(LOG2FC)>log(1.2,2) & PVALUE < 0.05)
points(LOG2FC[select_point], -log(PVALUE[select_point],10), pch=16,col='red')
#####################
text(x=0.5, y=0.5,labels='D0',cex=1.5)
text(x=-0.5, y=0.5,labels='D7',cex=1.5)
dev.off()



D0_OUT=rownames(all_res)[which( LOG2FC> log(1.2,2) & PVALUE < 0.05)]
D7_OUT=rownames(all_res)[which( LOG2FC< -log(1.2,2) & PVALUE < 0.05)]

D0_OUT=strsplit(D0_OUT, '\\.')
D7_OUT=strsplit(D7_OUT, '\\.')

.list2mat <- function(LIST){
    MAT=c()
    i=1
    while(i<=length(LIST)){
        MAT=cbind(MAT,LIST[[i]])
        i=i+1}
    MAT=t(MAT)
    return(MAT)
    }


D0_BED=.list2mat(D0_OUT)
D7_BED=.list2mat(D7_OUT)

write.table(D0_BED,file='DESeq2_D0.bed',sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(D7_BED,file='DESeq2_D7.bed',sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)





