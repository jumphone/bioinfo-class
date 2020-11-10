bedtools intersect -a ./H3K27ac_ChIP_H9/ALL.bed -b DESeq2_D7.bed -wa -wb | cut -f 1,2,3,7 | bedtools sort -i - | uniq - > DESeq2_D7_SE_EN.bed
bedtools intersect -a ./H3K27ac_ChIP_H9/ALL.bed -b DESeq2_D0.bed -wa -wb | cut -f 1,2,3,7 | bedtools sort -i - | uniq - > DESeq2_D0_SE_EN.bed



findMotifsGenome.pl DESeq2_D7.bed hg19 DESeq2_D7.bed.motif 
findMotifsGenome.pl DESeq2_D0.bed hg19 DESeq2_D0.bed.motif 
findMotifsGenome.pl DESeq2_D7_SE_EN.bed hg19 DDESeq2_D7_SE_EN.bed.motif 
findMotifsGenome.pl DESeq2_D0_SE_EN.bed hg19 DESeq2_D0_SE_EN.bed.motif 



annotatePeaks.pl DESeq2_D7.bed hg19 > DESeq2_D7.bed.anno.txt
annotatePeaks.pl DESeq2_D0.bed hg19 > DESeq2_D0.bed.anno.txt
annotatePeaks.pl DESeq2_D7_SE_EN.bed hg19 > DESeq2_D7_SE_EN.bed.anno.txt
annotatePeaks.pl DESeq2_D0_SE_EN.bed hg19 > DESeq2_D0_SE_EN.bed.anno.txt




