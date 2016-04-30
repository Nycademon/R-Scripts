# Load volume data into R
file1 <- "duncan_norm_vols.txt"
normvols <- read.table(file1, header=T, row.names=1, sep="\t")
# Get number of columns
cols <- ncol(normvols)-3
# Get column names
n <- colnames(normvols)
# Create matrix
# a <- c("comparison","Df","Sum Sq","Mean Sq","F value","Pr(>F)","significant")
# m <- data.frame()
# Start for loop
for(Z in 1:cols){
    i
    y <- anova(lm(normvols[,n[Z]]~age*somy*region,data=normvols))
    y
#    write.table(y,file="results.txt",append = TRUE)
}
# Export table
# file2 <- "3wayanova_duncan_norm_vols.txt"
# write.table(m, file = file2, quote=F, col.names= F, row.names=F, sep = "\t")
