GilbertCounts = read.table(file = "20150112_ctxdta_crispri_mergedcountstable.txt", sep = "\t", header = TRUE)
GilbertCounts = GilbertCounts[-which(rowSums(GilbertCounts[ ,6:9]) == 0), ]
counts = GilbertCounts[ ,grep("treat", colnames(GilbertCounts))]
colData = data.frame(replicate = rep(c("JEV", "MH"), times = 2), treatment = rep(c("treated", "untreated"), each = 2))
colData$treatment = factor(colData$treatment, levels = c("untreated", "treated"))
rownames(colData) = colnames(counts)
counts.DESeq = DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~treatment)
counts.DESeq = DESeq2::DESeq(counts.DESeq)
counts.DESeq = DESeq2::results(counts.DESeq)
log2fc = counts.DESeq$log2FoldChange
log2fc.negCtrl = log2fc[which(GilbertCounts$gene == "negative_control")]
log2fc.geneTargeting = log2fc[which(GilbertCounts$gene != "negative_control")]
geneIds = GilbertCounts$gene[which(GilbertCounts$gene != "negative_control")]
geneIds = factor(geneIds, levels = unique(geneIds))
library(CRISPhieRmix)
GilbertTreatedVsUntreatedCRISPhieRmix = CRISPhieRmix(x = log2fc.geneTargeting, geneIds = geneIds, negCtrl = log2fc.negCtrl, PLOT = TRUE, BIMODAL = TRUE)

GilbertTreatedVsUntreatedCRISPhieRmix.GeneScores = data.frame(gene = GilbertTreatedVsUntreatedCRISPhieRmix$genes, 
                                                              positivePosteriorProb = GilbertTreatedVsUntreatedCRISPhieRmix$posGenePosteriors, 
                                                              negativePosteriorProb = GilbertTreatedVsUntreatedCRISPhieRmix$negGenePosteriors,
                                                              FDR = sapply(GilbertTreatedVsUntreatedCRISPhieRmix$locfdr, function(x) 
                                                                mean(GilbertTreatedVsUntreatedCRISPhieRmix$locfdr[which(GilbertTreatedVsUntreatedCRISPhieRmix$locfdr <= x)])))
write.table(GilbertTreatedVsUntreatedCRISPhieRmix.GeneScores, 
            file = "GilbertTreatedVsUntreatedCRISPhieRmixGeneScores.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
