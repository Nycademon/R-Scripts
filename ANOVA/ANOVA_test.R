## Load volume data into R
# Get data file name
file <- "values.txt"
# Create data frame
peakvols <- read.table(file, header=T, row.names=1, sep="\t")
# Get number of columns
cols <- ncol(peakvols)-3
# Get column names
n <- colnames(peakvols)
# Create results table, add master columns
b <- c("Df","Sum_Sq","Mean_Sq","F_value","Pr(>F)","Sample")
m <- matrix(b,ncol=6)
write.table(m,file="results.txt",quote=F,col.names=F,sep="\t")
# Start for loop
for(i in 1:cols){
    # Calculate ANOVAs
    y <- anova(lm(peakvols[,n[i]]~age*somy*region,data=peakvols))
    # Populate column 6 with peak name
    y[,6] <-n[i]
    # append results to output file
    write.table(y,file="results.txt",quote=F,col.names=F,sep="\t",append = TRUE)
}
