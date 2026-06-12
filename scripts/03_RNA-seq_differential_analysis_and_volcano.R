# Differential expression analysis and volcano plot
# Dataset: GSE126848

library(DESeq2)
library(ggplot2)

# 1. Differential analysis
sel=samples.info$CONDITION %in% c("DISEASE","CONTROL") #select the column CONDITION 

data=Gdata
data$counts =data$counts[,sel]
data$abundance=data$abundance[,sel]
data$length=data$length[,sel]

samples.info$CONDITION <- factor(
	samples.info$CONDITION, 
	levels = c(
		"DISEASE", 
		"CONTROL"
	)
) #transform characters in factors

# create a new dds and keep only gene with count more than 10
dds = DESeqDataSetFromTximport(
	data, 
	samples.info[sel,],
	~CONDITION
)
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

dds_d = DESeq(dds)
resultsNames(dds_d) # check the name
res=data.frame(
	results(
		dds_d, 
		name="CONDITION_CONTROL_vs_DISEASE", 
		lfcThreshold=1
		)
	)

res<-res[!is.na(res$padj),]
res<-cbind(
	res, 
	annot[match(rownames(res),annot[,1]),c(3,4,5)]
)

write.table(res, file="results/res",sep="\t")

# select only values with p-value<0,05 and |fold change| =>1 (doubled or halved)
res_sel<-res[res$padj<0.05 & abs(res$log2FoldChange)>=1,] 

# save the differentually expressed genes table
write.table(res_sel, file="results/res_sel",sep="\t")

# 2. Shrinkage
res_shrinked<-lfcShrink(
	dds_d, 
	coef="CONDITION_CONTROL_vs_DISEASE", 
	type="apeglm" 
)
res_shrinked<-res_shrinked[!is.na(res_shrinked$padj),]

## annotation
res_shrinked<-cbind(
	res_shrinked, 
	annot[match(rownames(res_shrinked), annot[,1]),c(3,4,5)]
)

# filter with padj and log2FC
res_shrinked_sel<-res_shrinked[res_shrinked$padj<0.05 & abs(res_shrinked$log2FoldChange)>=1,]

write.table(res_shrinked_sel,file="results/res_shrinked_sel",sep="\t")

# 3. Volcano plot
df_clean<-df[complete.cases(res), ]
df_clean<-df_clean[df_clean$padj > 0,]
df_clean$diffexpressed <- "NO"
df_clean$diffexpressed[
	df_clean$log2FoldChange > log2_fc_thr & 
	df_clean$padj < p_value_thr
] <- "UP"

df_clean$diffexpressed[
	df_clean$log2FoldChange < -log2_fc_thr & 
	df_clean$padj < p_value_thr
] <- "DOWN"

df_clean<-rbind(
	df_clean[df_clean$diffexpressed == "NO",],
	df_clean[df_clean$diffexpressed == "UP",],
	df_clean[df_clean$diffexpressed == "DOWN",]
)

color_values <- c(
	"DOWN" = "#00AFBB",
	"NO" = "grey",
	"UP" = "#bb0c00"
)

color_labels<-c(color_labels,"NO")


xlimit <- ceiling(max(c(
	abs(max(df_clean$log2FoldChange)),
	abs(min(df_clean$log2FoldChange))
)) / 2.5) * 2.5

ylimit <- -log10(min(df_clean$padj))

p<-ggplot2::ggplot(
	data = df_clean, 
	ggplot2::aes(
		x = log2FoldChange,
		y = -log10(padj),
		col = diffexpressed,
	)
) +

ggplot2::geom_vline(xintercept = c(-log2_fc_thr, log2_fc_thr), col = "gray", linetype = 'dashed') + #border up and down
 ggplot2::geom_hline(yintercept = -log10(p_value_thr), col = "gray", linetype = 'dashed')+ #significatività 

ggplot2::geom_point(size = 0.5) + #desing a gene like a 0.5 size dot

ggplot2::scale_color_manual( #color the dots
	values = color_values, 
	labels = color_labels
) +
ggplot2::coord_cartesian(ylim = c(0, ylimit), xlim = c(-xlimit, xlimit)) + #limits of the axes
ggplot2::labs( #label
	color = 'Expression', #legend_title
	x = expression("log"[2]*"FC"),
	y = expression("-log"[10]*"p-value adj.")
) +

ggplot2::ggtitle("GSE126848") #assign the title to the plot

print(p)

ggsave(
	filename = "plots/GSE126848_volcano_plot.png",
	plot = p,
	dpi = 300,
	width = 1800, 
	height = 2400, 
	units = "px"
)

ggsave(
	filename = "plots/GSE126848_volcano_plot.tiff",
	plot = p,
	dpi = 300,
	width = 1800, 
	height = 2400, 
	units = "px",
	compression = "lzw"
)

ggsave(
	filename = "plots/GSE126848_volcano_plot.pdf",
	plot = p,
	width = 1800, 
	height = 2400, 
	units = "px"
)

up_genes<-df_clean$Symbol[df_clean$diffexpressed=="UP"]
down_genes<-df_clean$Symbol[df_clean$diffexpressed=="DOWN"]

