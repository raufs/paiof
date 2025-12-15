library(ggplot2)

dat <- read.table("Boxplot_Inputs.txt", header=T, sep='\t')
# comparison      ahba_genes      aai

dat.filt <- dat[dat$ahba_genes == 'True',]

pdf("Percent_Identity_Boxplots_Paiof.pdf", height=5, width=15)
ggplot(dat, aes(x=comparison, y=aai)) + geom_boxplot(color='black') + geom_point(data=dat.filt, aes(x=comparison, y=aai), size=7, color='red', shape=4)+ coord_flip() + theme_classic() + xlab("") + ylab("")
dev.off()
