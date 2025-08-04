library(ggplot2)

dat <- read.table("Plotting_Input.txt", header=T, sep='\t')
pdf("AAI_Plot_Simple_vs_TwoMulti.pdf", height=5, width=10)
ggplot(dat, aes(x=AAI, y=Shared_Prop, color=Comparator_Group)) + geom_point(alpha=0.75) + theme_bw() + facet_wrap(~Focal, nrow=2) + scale_color_manual(values=c('#9fc4cc', '#41668a'))
dev.off()
