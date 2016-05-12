# Load coin package
library("coin")
# Get enzyme activity data
sample <- read.table("samples.txt", header=T, sep="\t")
sample2 <- read.table("samples2.txt", header=T, sep="\t")

## Run wilcox.test
# null hypothesis of equality
sample_wc <- wilcox_test(Act ~ Group, data=sample, distribution="exact",conf.int=TRUE)
sample2_wc <- wilcox_test(Act ~ Group, data=sample2, distribution="exact",conf.int=TRUE)

# alternative = "less"
sample_wc_less <- wilcox_test(Act ~ Group, data=sample, distribution="exact",conf.int=TRUE, alternative = "less")
sample2_wc_less <- wilcox_test(Act ~ Group, data=sample2, distribution="exact",conf.int=TRUE, alternative = "less")

# alternative = "greater"
sample_wc_greater <- wilcox_test(Act ~ Group, data=sample, distribution="exact",conf.int=TRUE, alternative = "greater")
sample2_wc_greater <- wilcox_test(Act ~ Group, data=sample2, distribution="exact",conf.int=TRUE, alternative = "greater")

## Get statistics
# null hypothesis
statistic(sample_wc, "linear")
statistic(sample2_wc, "linear")

# "less"
statistic(sample_wc_less, "linear")
statistic(sample2_wc_less, "linear")

# "greater"
statistic(sample_wc_greater, "linear")
statistic(sample2_wc_greater, "linear")

# Write out results
colz<-c("samples","samples2")
pval<-c(pvalue(sample_wc),pvalue(sample2_wc))
write.table(pval,file="wc_pvalues.txt",quote=F,row.names=colz,col.names=F,sep="\t")
colz_less<-c("samples_less","samples2_less")
pval_less<-c(pvalue(sample_wc_less),pvalue(sample2_wc_less))
write.table(pval_less,file="wc_pvalues_less.txt",quote=F,row.names=colz_less,col.names=F,sep="\t")
colz_greater<-c("samples_greater","samples2_greater")
pval_greater<-c(pvalue(sample_wc_greater),pvalue(sample2_wc_greater))
write.table(pval_greater,file="wc_pvalues_greater.txt",quote=F,row.names=colz_greater,col.names=F,sep="\t")

