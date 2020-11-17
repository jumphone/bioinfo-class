#https://jokergoo.github.io/ComplexHeatmap-reference/book/

setwd('  ')


library('ComplexHeatmap')
library('circlize')
library('seriation')


MAT=read.table('DEMO_MAT.txt',sep='\t',row.names=1,header=TRUE)
ROW.MATE=read.table('DEMO_ROW.META.txt',sep='\t',row.names=1,header=TRUE)
COL.MATE=read.table('DEMO_COL.META.txt',sep='\t',row.names=1,header=TRUE)


MAT[1:5,1:5]
head(ROW.MATE)
head(COL.MATE)



ORDER.COL=order(COL.MATE$SCORE)
ORDER.ROW=order(ROW.MATE$RM1)



############################################
ha_top = HeatmapAnnotation(  
     SCORE =anno_lines( COL.MATE$SCORE[ORDER.COL] )
     )
#############################################




#############################################
color_fun_1=colorRamp2(c(0,0.01,0.025,0.03,0.05 ), c('royalblue3','royalblue3','white','indianred3','indianred3'))
ha_bot = HeatmapAnnotation(  
     CM1 = COL.MATE$CM1[ORDER.COL],
     CM2 = COL.MATE$CM2[ORDER.COL], 
     col = list(CM1=c('Ctrl'='royalblue3','Case'='indianred3','Other'='white'),
                CM2=color_fun_1 )
     )
##############################################


##############################################
color_fun_2=colorRamp2(c(1,95,189 ), c('royalblue3','white','indianred3'))
ha_left = rowAnnotation( RM1 =ROW.MATE$RM1,
                        col = list(RM1=color_fun_2)
                         )
##############################################


################################################
o.mat=as.matrix(MAT[ORDER.ROW,ORDER.COL])
color_fun_3 =colorRamp2(c(-1,-0.5,0,0.5,1 ), c('royalblue3','white','white','white','indianred3'))
Heatmap(o.mat,row_title='',name="C",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=FALSE, show_row_names=FALSE,
	      col=color_fun_3, border = TRUE,
        top_annotation = ha_top,
	      bottom_annotation = ha_bot,
        left_annotation = ha_left
	      )
###########################################


################################################
o.mat=as.matrix(MAT[ORDER.ROW,ORDER.COL])
color_fun_3 =colorRamp2(c(-1,-0.5,0,0.5,1 ), c('royalblue3','white','white','white','indianred3'))
ht = Heatmap(o.mat,row_title='',name="C",
        cluster_columns=TRUE,  cluster_rows=TRUE,
	      show_column_dend = TRUE, show_row_dend = TRUE, 
	      show_column_names=FALSE, show_row_names=FALSE,
	      col=color_fun_3, border = TRUE,
        top_annotation = ha_top,
	      bottom_annotation = ha_bot,
        
	      )

ht = draw(ht)
row_order(ht)
ht.dend=as.dendrogram(ht@ht_list$C@row_dend_list[[1]])
plot(ht.dend)
###########################################















