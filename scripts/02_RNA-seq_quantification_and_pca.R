# Import Salmon quantification and perform PCA
# Dataset: GSE126848

library(tximport)
library(biomaRt)
library(DESeq2)
library(ggplot2)
library(cowplot)


# 1. Import Salmon quantification files

files=paste0(dir()[2:58],"/salmon_results/quant.sf")  
Tdata=tximport(
	files,
	type="salmon",
	txOut=T,
	varReduce=T
)

# modify rownames in Tdata to make the match with annot
newrownames=gsub("\\.\\d+","",rownames(Tdata$counts))
rownames(Tdata$counts)=newrownames
rownames(Tdata$abundance)=newrownames
rownames(Tdata$length)=newrownames
rownames(Tdata$variance)=newrownames

# 2. Retrieve biomart annotation
mart <- biomaRt::useDataset(
	dataset = "hsapiens_gene_ensembl",
	mart = useMart(
		"ENSEMBL_MART_ENSEMBL",
		host = "https://www.ensembl.org"
	)
)

annot <- biomaRt::getBM(
	attributes=c(
		"ensembl_gene_id", 
		"ensembl_transcript_id", 
		"external_gene_name", 
		"entrezgene_id", 
		"description", 
		"strand"
	), 
	mart=mart
) 

names(annot)=c(
	"Gene",
	"Transcript",
	"Symbol",
	"Entrez",
	"Description",
	"Strand"
) 

annot=annot[match(rownames(Tdata$counts),annot$Transcript),] #match between counts and trascripts

# 3. Summarize transcript level data to gene-level data

Gdata=summarizeToGene(
	Tdata,
	annot[,c(2,1)]
) #trascript to gene

#4. Import sample metadata

samples.info=read.delim(
	"metadata/SRR_info", 
	sep="\t", 
	header=F
) #import informations about samples (GSM, GSE, SRR, condition)

# label columns and row in samples.info
colnames(samples.info)<-c(
	"GSM",
	"GSE",
	"CONDITION", 
	"SRR"
)
rownames(samples.info) <- samples.info$SRR

# label columnns in Gdata
colnames(Gdata$counts)<-samples.info$SRR 
colnames(Gdata$abundance)<-samples.info$SRR
colnames(Gdata$length)<-samples.info$SRR

# 5. PCA

dds = DESeqDataSetFromTximport(
	Gdata, 
	samples.info, 
	~CONDITION
)

# applying a variance stabilizing transformation (VST) to normalize the raw data
vst_data <- vst(dds, blind = TRUE)
vst_mat <- assay(vst_data)
v <- apply(vst_mat, 1, var)
vst_mat <- vst_mat[v > 0, ]

# transpose the matrix
pca <- prcomp(t(vst_mat),scale=TRUE)

# create parameter for the plot
percentVar <- pca$sdev^2 / sum(pca$sdev^2)

match_indexes = match(colnames(Gdata),rownames(samples.info))

xlab = paste0("PC1: ", round(percentVar[1]*100, 1), "% variance")
ylab = paste0("PC2: ", round(percentVar[2]*100, 1), "% variance")

pca_df <- data.frame(
	PC1 = pca$x[,1],
	PC2 = pca$x[,2]
)
pca_df<-pca_df[complete.cases(pca_df), ]

x_min <- min(pca_df[,1])
x_max <- max(pca_df[,1])
y_min <- min(pca_df[,2])
y_max <- max(pca_df[,2])
max <- abs(round(max(x_max,y_max) + 10))
min <- abs(round(min(x_min,y_min) - 10))
range=c(-min,max)

png("plot/pca.png", width = 500, height = 500) #save the plot
plot.new() #open
par(bg = 'white')
condition_colors<-c(
	"DISEASE"="blue",
	"HEALTHY"="red"
)

cols<-condition_colors[samples.info$CONDITION]
plot(
	pca_df[,1], pca_df[,2],
        xlim = range,
        ylim = range,
        xlab = xlab,
        ylab = ylab,
	main = "GSE126848",
        pch = 19,
	col=cols
)

legend(
	"bottomleft",
	legend=names(condition_colors),
	col=condition_colors,
	pch=19,
	title="Condition"
)
dev.off()

