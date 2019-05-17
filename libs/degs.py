"""
    deferential analysis
"""
import os
import subprocess


class DESeqHandler:
    def __init__(self, config):
        self.config = config

    def model_script(self):
        diretory = self.config.get('REPORT')
        rscript = f"""library(DESeq2)
        library(airway)
        library(pheatmap2)
        library(RColorBrewer)
        # load and prepare count data
        directory <- '{directory}'
        sampleFiles <- grep('treated', list.files(directory), value=True)
        sampleCondition <- sub(...)
        sampleTable <- data.frame(sampleName = sampleFiles,
                                  fileName = sampleFiles,
                                  condition = sampleCondition)
        ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                               directory = diretory,
                                               design =~ condition)
        # pre-filtering
        keep <- rowSums(counts(dds)) >= 10
        dds <- dds[keep,]
        # note on factor levels
        dds$condition <- relevel(dds$condition, ref={ref})
        # differential analysis
        dds <- DESeq(dds)
        res <- result(dds)
        resSig <- subset(res, padj < 0.05)
        # save the result
        tx_count <- counts(dds, normalized=TRUE)
        output <- merge(res, tx_count, by=0)
        write.csv(output, filename='degs.csv')
        # plots
        plotMA(res, ylim=c(-2, 2))
        rld <- rlog(dds, blind=FALSE)
        rldSig <- rld[rownames(resSig),]
        pheadmap(assay(rldSig))

        sampleDistsf <- dist(t(assay(rld)))
        sampleDistMatrix <- as.matrix(sampleDists)
        rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep='-')
        clonames(sampleDistMatrix) <- NULL
        colors <- colorRampPalette( rev(brewer.pal(9, 'Blues')))(255)
        pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 col=colors)
        plotPCA(rld, intgroup=c('condition', 'type'))
        """
        with open('deseq_manager.R', 'w') as f:
            f.write(rscript)
        return script

    def process(self, script):
        try:
            subprocess.run([
                Rscript,
                script
            ])
        except subprocess.CalledProcessError:
            pass
