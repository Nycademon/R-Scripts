### Load programs
library("cummeRbund")
# Need plyr to merge gene expression and annotation tables, since gene expression
# (diffData) tables don't have unique IDs.
library("plyr")

### Read Cufflinks data
# The gtf file can be via a symlink, genome is required with the gtfFile option,
# 'rebuild = T' is sometimes necessary to incorprate changes.
cuff <- readCufflinks(gtfFile = "merged.gtf", genome = "CriGri_1.0")
cuff

### Various plots describing the data set
pdf("dispersionPlot.pdf", width=10, height=10)
print(dispersionPlot(genes(cuff)))
dev.off()
pdf("genesSCV.pdf", width=10, height=10)
print(fpkmSCVPlot(genes(cuff)))
dev.off()
pdf("isoformsSCV.pdf", width=10, height=10)
print(fpkmSCVPlot(isoforms(cuff)))
dev.off()
pdf("csDensity.pdf", width=10, height=10)
print(csDensity(genes(cuff)))
dev.off()
pdf("csDensity_reps.pdf", width=10, height=10)
print(csDensity(genes(cuff),replicates=T))
dev.off()
pdf("csBoxplot.pdf", width=10, height=10)
print(csBoxplot(genes(cuff)))
dev.off()
pdf("csBoxplot_reps.pdf", width=10, height=10)
print(csBoxplot(genes(cuff),replicates=T))
dev.off()
pdf("csScatterMatrix.pdf", width=10, height=10)
print(csScatterMatrix(genes(cuff)))
dev.off()
pdf("csDendro_reps.pdf", width=10, height=10)
print(csDendro(genes(cuff),replicates=T))
dev.off()
pdf("csDendro.pdf", width=10, height=10)
print(csDendro(genes(cuff)))
dev.off()
pdf("csVolcanoMatrix.pdf", width=10, height=10)
print(csVolcanoMatrix(genes(cuff)))
dev.off()
pdf("PCA.pdf", width=10,height=10)
print(PCAplot(genes(cuff),"PC1","PC2"))
dev.off()
pdf("PCAreps.pdf", width=10,height=10)
print(PCAplot(genes(cuff),"PC1","PC2",replicates=T))
dev.off()
pdf("csDistHeat.pdf", width=10,height=10)
print(csDistHeat(genes(cuff)))
dev.off()
pdf("csDistHeatreps.pdf", width=10,height=10)
print(csDistHeat(genes(cuff),replicates=T))
dev.off()

### Gene specific output
jeans <- c("Adsl","Adss","Atic","Gart","Pkm")
l <- length(jeans)
for (i in 1:l){
    myGeneID <- as.character(jeans[i])
    myGene<-getGene(cuff,myGeneID)
    pdf(paste0(myGeneID,"_exp.pdf"), width=10,height=10)
    print(expressionPlot(myGene))
    dev.off()
    pdf(paste0(myGeneID,"_exp_reps.pdf"), width=10,height=10)
    print(expressionPlot(myGene,replicates=T))
    dev.off()
    pdf(paste0(myGeneID,"_exp_iso_reps.pdf"), width=10,height=10)
    print(expressionPlot(isoforms(myGene),replicates=T))
    dev.off()
    pdf(paste0(myGeneID,"_exp_cds_reps.pdf"), width=10,height=10)
    print(expressionPlot(CDS(myGene),replicates=T))
    dev.off()
    # Need this prior to makeGeneRegionTrack (default expects "regular" chromosome names).
    options(ucscChromosomeNames=FALSE)
    genetrack<-makeGeneRegionTrack(myGene)
    pdf(paste0(myGeneID,"_genetrack.pdf"), width=10,height=10)
    print(plotTracks(genetrack))
    dev.off()
}

### Gene set specific output

source(file="GeneSets.R")

# List of Gene Sets
GeneSet <- list(EifIds,RpIds,MrpIds,PurPyrIds,PsomeIds,tRNASynIds,ATPaseIds,NADHDeIds,TFacIds,
                CompI_Ids,CompII_Ids,CompIII_Ids,CompIV_Ids,CompV_Ids)

# Vector of Gene Set Names
GeneSetName <- c("EifIds","RpIds","MrpIds","PurPyrIds","PsomeIds","tRNASynIds","ATPaseIds","NADHDeIds","TFacIds",
                    "CompI_Ids","CompII_Ids","CompIII_Ids","CompIV_Ids","CompV_Ids")

y<-length(GeneSet)

for(j in 1:y){
    ## OK, THE BELOW WORKS, but could be better. You don't need ATPaseIds, etc.
    # NOTE: GeneIds is a "List of 1". Solved using unlist.
    # GeneIds<-GeneSet[j]
    GeneIds<-unlist(GeneSet[j])
    head(GeneIds)
    Genes<-getGenes(cuff,GeneIds)
    Genes
    Genes_diff<-diffData(Genes)
    Genes_featureNames<-featureNames(Genes)
    head(Genes_featureNames)
    # Need to have the same column name in both tables to do a join, so rename tracking_id to gene_id.
    names(Genes_featureNames)<-sub('tracking_id','gene_id',names(Genes_featureNames))
    # Then combine the tables using gene_id for the join
    Genes_merge <- join(Genes_diff,Genes_featureNames,by='gene_id',type='left',match='all')
    write.table(Genes_merge, file=paste0(GeneSetName[j],"_AdeC_gene_exp.txt"), sep="\t", col.names=T)
    pdf(paste0(GeneSetName[j],"_cluster3.pdf"), width=10, height=10)
    # k=6 seems to cause problems sometimes.
    print(csClusterPlot(csCluster(Genes,k=3)))
    dev.off()
    postscript(paste0(GeneSetName[j],"_heatmap_bw.ps"), paper="letter", horizontal = F)
    print(csHeatmap(Genes,heatscale = c(low='white',high='black'),clustering = 'none'))
    dev.off()
    postscript(paste0(GeneSetName[j],"_heatmap.ps"), paper="letter", horizontal = F)
    print(csHeatmap(Genes,heatscale = c(low='red',mid='white',high='blue'),clustering = 'none'))
    dev.off()
    postscript(paste0(GeneSetName[j],"_heatmap_wb.ps"), paper="letter", horizontal = F)
    print(csHeatmap(Genes,heatscale = c(low='white',high='blue'),clustering = 'none'))
    dev.off()
}
