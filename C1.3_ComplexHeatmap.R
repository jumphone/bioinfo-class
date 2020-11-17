'''
colnames(MAT)=paste0('C',1:ncol(MAT))
write.table(MAT,'DEMO_MAT.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
CM1=DevStage
CM1[which(DevStage=='Adult')]='Ctrl'
CM1[which(DevStage=='Fetal')]='Case'
CM1[which(DevStage=='NA')]='Other'
CM2=ST_EXP
COL.META=cbind(SCORE, CM1, CM2)
rownames(COL.META)=colnames(MAT)
write.table(COL.META,'DEMO_COL.META.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
ROW.META=cbind(as.numeric(rownames(MAT)),189:1)
colnames(ROW.META)=c('RM1','RM2')
rownames(ROW.META)=rownames(MAT)
write.table(ROW.META,'DEMO_ROW.META.txt',sep='\t',quote=FALSE,row.names=TRUE,col.names=TRUE)
'''





