# Read in table (from norm_volumes.txt, from prot2D)
PCAread<-read.table("norm_volumes.txt", header=T, row.names=1, sep="\t")
# Transpose the matrix, make data.frame
PCAtab<-data.frame(t(PCAread))
# Check
head(PCAtab)
# Do PCA
PCAtab.pca=prcomp(~., data=PCAtab, center=T, scale=T)
#Get PCA results summary
summary(PCAtab.pca)
#Plot the variances associated with each PCA axis
pdf("variances.pdf", width=10, height=10)
plot(PCAtab.pca)
dev.off()
#Plot the scores associated with PC1 and PC2 (by default), the cex option makes the "points" smaller. The ablines puts horizontal and vertical lines through the plot origin. The text command replaces the circles with text labels.
pdf("PCA_scatter.pdf", width=10, height=10)
plot(PCAtab.pca$x, type="p", cex=1.0, col="red")
abline(h=0,col=8)
abline(v=0,col=8)
text(PCAtab.pca$x,labels=rownames(PCAtab.pca$x),cex=0.5)
dev.off()
#Produce a SVG file from the plot, that you can manipulate with Illustrator.
svg("scatterplot.svg", width=10, height=10)
plot(PCAtab.pca$x,type="p",cex=0.5,col="red")
dev.off()

#Produce a 3D plot, with PC1, PC2, and PC3
pdf("PCA_3D.pdf", width=10, height=10)
library("scatterplot3d", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
scatterplot3d(PCAtab.pca$x[,1],PCAtab.pca$x[,2],PCAtab.pca$x[,3], main="3D PCA Plot")
dev.off()

#Plot the "samples" (mice) which indicate how the sample variation is aligned with variation in principal component axes.
pdf("PCA_Samples.pdf", width=10, height=10)
plot(PCAtab.pca$rotation)
dev.off()

#Can get row names with str()
str(PCAtab)

#Column labels
g=c("6CD","6CD","6CD","6CD","6CD","6CD","6CD","6CD","6CT","6CT","6CT","6CT","6CT","6CT","6CT","6CT","6HD","6HD","6HD","6HD","6HD","6HD","6HD","6HD","6HT","6HT","6HT","6HT","6HT","6HT","6HT","6HT","12CD","12CD","12CD","12CD","12CD","12CD","12CD","12CD","12CD","12CD","12CD","12CD","12CD","12CD","12CT","12CT","12CT","12CT","12CT","12CT","12CT","12CT","12CT","12CT","12CT","12CT","12CT","12CT","12HD","12HD","12HD","12HD","12HD","12HD","12HD","12HD","12HD","12HD","12HD","12HD","12HD","12HD","12HD","12HT","12HT","12HT","12HT","12HT","12HT","12HT","12HT","12HT","12HT","12HT","12HT","12HT","12HT","12HT")
identify(PCAtab$rotation, labels=g)

#biplot() overlays loading plot and scores plot in one figure (cex resizes the text to a reasonable level)
biplot(PCAtab.pca,cex=0.6)

#To create plots or biplots using PC axes other than default . . .
plot(PCAtab.pca$x[,1],PCAtab.pca$x[,3],xlab="PC1",ylab="PC3")
biplot(PCAtab.pca,choices=c(1,3))


# NOTES:
    