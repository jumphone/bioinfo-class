setwd('/home/database/data/RNAseq_H9d0tod45/hisat2')
data=read.table('MAT.txt',header=T,row.names=1,sep='\t')
pc.gene=read.table('/home/database/annotation/hg19/Homo_sapiens.GRCh37.75.chr.pc.gene.SYM.bed',sep='\t',header=TRUE,row.names=NULL)
pc.gene=as.character(pc.gene[,4])

data=data[which(row.names(data) %in% pc.gene),]
#data=data[,2:ncol(data)]

used.data=data[,1:7]
colnames(used.data)=c('LENGTH','D0_R1','D0_R2','D7_R1','D7_R2','D14_R1','D14_R2')


write.table(used.data,file='./class_use/COUNT_MAT.txt',col.names=T,row.names=T,sep='\t',quote=FALSE)
