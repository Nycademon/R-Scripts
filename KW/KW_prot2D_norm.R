# Load volume data into R
file1 <- "prot2d_norm_volumes_all.txt"
normvols <- read.table(file1, header=T, row.names=1, sep="\t")
# Get number of columns
cols <- ncol(normvols)-3
# Get column names
n <- colnames(normvols)
# Create four column matrix
a <- c("Spot","comparison","chi-squared","df","p-value")
m <- matrix(a,ncol=5)
# Start for loop
for(i in 1:cols){
    x <- kruskal.test(normvols[,n[i]]~age,data=normvols)
    y <- kruskal.test(normvols[,n[i]]~region,data=normvols)
    z <- kruskal.test(normvols[,n[i]]~somy,data=normvols)
    d <- as.character(normvols[,n[i]])
    p <- c(as.character(c(n[i],"age")),as.numeric(as.vector(c(x[1],x[2],x[3]))))
    q <- c(as.character(c(n[i],"region")),as.numeric(as.vector(c(y[1],y[2],y[3]))))
    r <- c(as.character(c(n[i],"somy")),as.numeric(as.vector(c(z[1],z[2],z[3]))))
    m <- rbind(m,p,q,r)
}
# Export table
file2 <- "kw_prot2D_norm_vols.txt"
write.table(m, file = file2, quote=F, col.names= F, row.names=F, sep = "\t")
