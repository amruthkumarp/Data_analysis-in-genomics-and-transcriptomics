# Load packages

library(SGSeq)
library(pheatmap)

# Define sample info for 4 BAM files

file_bam<-c("Con_A1_sorted.bam","Con_A2_sorted.bam","Con_B1_sorted.bam","Con_B2_sorted.bam")
sample_name<-c("mt1","mt2","WT1","WT2")
v<-data.frame(sample_name,file_bam)
Baminfo=getBamInfo(v)
#Import Annotation file and converting the formats
x<-importTranscripts(file = "R64_genomic.gtf")

#Convert	a	TxDb	object	or	a	GRangesList	of	exons	grouped	by	transcripts	to	a	TxFeatures	object.

TxFeat<-convertToTxFeatures(x)
TxFeat
type(TxFeat)
head(txName(TxFeat))
head(geneName(TxFeat))

sgf_ucsc <- convertToSGFeatures(TxFeat)
head(sgf_ucsc)
#analyzeFeatures   
sgfc_ucsc <- analyzeFeatures(Baminfo, features = TxFeat)
sgfc_ucsc



colData(sgfc_ucsc)
rowRanges(sgfc_ucsc)
(counts(sgfc_ucsc))
head(FPKM(sgfc_ucsc))

#df <- plotFeatures(sgfc_ucsc, geneID = 30)
sgfc_pred <- analyzeFeatures(Baminfo)
#annotate features
A2=annotate(sgfc_pred,TxFeat)
head(rowRanges(A2))
df <- plotFeatures(A2, geneID = 10,color_novel="red")
