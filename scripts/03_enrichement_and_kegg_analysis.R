# Enrichment and KEGG analysis
# Dataset: GSE126848

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(KEGGREST)
library(VennDiagram)

# 1. Prepare ranked genes list for GSEA
res <- results(dds_d, contrast = c("CONDITION","DISEASE","CONTROL"))
res_all <- as.data.frame(res)

res_all$Symbol <- annot$Symbol[match(rownames(res_all), annot$Gene)]
res_all <- res_all[!is.na(res_all$Symbol) & !is.na(res_all$log2FoldChange),]

res_gene_unique <- aggregate(
	res_all$log2FoldChange, 
	by = list(res_all$Symbol), 
	FUN = max, 
	na.rm = TRUE
)

colnames(res_gene_unique) <- c("Symbol","log2FoldChange")

res_gene_filtered <- res_gene_unique[abs(res_gene_unique$log2FoldChange) > 1,]

write.table(res_gene_filtered, file="results/transcriptomic_symbol_logFC.txt", sep="\t", quote=F, col.names=T, row.names=F)

gene_list <- res_gene_unique$log2FoldChange
names(gene_list) <- res_gene_unique$Symbol
gene_list <- gene_list[!is.na(names(gene_list)) & names(gene_list) != ""]
gene_list <- sort(gene_list, decreasing=TRUE)

# 2. GSEA

organism="org.Hs.eg.db"
library(organism, character.only = TRUE)

gse <- gseGO(
	geneList=gene_list, 
	ont="BP", 
	keyType="SYMBOL", 
	nPerm=1000, 
	minGSSize =3, 
	maxGSSize=800, 
	pvalueCutoff=0.05, 
	verbose=TRUE, 
	OrgDb = organism, 
	pAdjustMethod="none"
)

gse_MF <- gseGO(
	geneList=gene_list, 
	ont="MF", 
	keyType="SYMBOL", 
	nPerm=1000, 
	minGSSize =3, 
	maxGSSize=800, 
	pvalueCutoff=0.05, 
	verbose=TRUE, 
	OrgDb = organism, 
	pAdjustMethod="none"
)

png("plots/GSEA.png", width = 1200, height = 800)
gseaplot(
	gse, 
	by="all", 
	title=gse$Description[1], 
	geneSetID=1
)
dev.off()

png("plots/ridgeplot.png", width = 1200, height = 800)
ridgeplot(gse, showCategory = 15) +
 labs(x = "enrichment distribution") +
 theme(
	axis.text.y = element_text(size = 10), 
	axis.text.x = element_text(size = 10), 
	axis.title.x = element_text(size = 14))
dev.off()

png("plots/dotplot_BP.png", width=3000, height=3000, res=300)
dotplot(
	gse, 
	showCategory=15, 
	orderBy="NES", 
	font.size=10,
	title = "Biological Process Enriched in DEGs"
)
dev.off()

png("plots/dotplot_MF.png", width=3000, height=3000, res=300)
dotplot(
	gse_MF, 
	showCategory=15, 
	orderBy="NES", 
	font.size=10
)
dev.off()

# 3. ORA
gene_sets <- clusterProfiler::read.gmt("metadata/c5.go.bp.v2024.1.Hs.symbols.gmt")
res_gene_rank <- res_all[order(res$padj),]

# split genes by direction
up_genes <- rownames(res_gene_rank)[res_gene_rank$log2FoldChange > 0]
down_genes <- rownames(res_gene_rank)[res_gene_rank$log2FoldChange < 0]

# run enrichment
enriched_up <- clusterProfiler::enricher(
	gene = up_genes,
	TERM2GENE = gene_sets
)

enriched_down <- clusterProfiler::enricher(
	gene = down_genes,
	TERM2GENE = gene_sets
)

enriched_all <- clusterProfiler::enricher(
	gene = res_gene_rank$Symbol,
	TERM2GENE = gene_sets
)

write.csv(
	as.data.frame(enriched_all),
	"results/ORA_all_NASH_vs_HEALTHY.csv",
	row.names = FALSE
)

write.csv(
	as.data.frame(enriched_up),
	"results/ORA_up_NASH_vs_HEALTHY.csv",
	row.names = FALSE
)
write.csv(
	as.data.frame(enriched_down),
	"results/ORA_down_NASH_vs_HEALTHY.csv",
	row.names = FALSE
)

png("plots/dotplot_enrich.png", width=1200, height=1200, res=150)
dotplot(enriched_all)
dev.off()

png("plots/barplot_up.png", width=1200, height=1200, res=150)
barplot(enriched_up, title = "GO Biological Process Enriched in Upregulated DEGs")
dev.off()

png("plots/barplot_down.png", width=1200, height=1200, res=150)
barplot(enriched_down, title = "GO Biological Process Enriched in Downregulated DEGs")
dev.off()

# . KEGG
Human.GRCh38.p13.annot <- read.table(
	"metadata/Human.GRCh38.p13.annot.tsv", 
	header=T, 
	sep="\t", 
	fill= T, 
	comment.char ="#"
)

# remove NA or invalid values
res_gene_rank <- res_gene_rank[
	res_gene_rank$Symbol != "" & 
	!is.na(res_gene_rank$Symbol), 
]

# merge to obtain an object with symbol and entrez of each gene
res_gene_entrez <- merge(
	res_gene_rank, 
	Human.GRCh38.p13.annot[, c("Symbol", "GeneID")], 
	by.x ="Symbol", 
	by.y = "Symbol", 
	all.x =T
)

path <- keggGet("hsa04932") #get the pathway
nafld_genes<- path[[1]]$GENE #extract the genes

nafld_entrez <- nafld_genes[seq(1, length(nafld_genes), 2)] #extract the entrez
nafld_entrez <- sub(";.*", "", nafld_entrez) #remove ;

## overlap to see if the gene I obtained are involved in the pathology
gene_list <- unique(na.omit(res_gene_entrez$GeneID)) #extract the gene obtained from transcriptomic
overlap <- intersect(
	gene_list, 
	nafld_entrez
)

length(overlap)

png("plots/venn_kegg.png", width = 800, height = 800)
venn.plot <- venn.diagram(
	x = list(
		Transcriptomic = gene_list, 
		"NAFLD pathway" = nafld_entrez
	),
	filename = NULL, 
	scaled=FALSE, # to have circle of the same dimension
	fill = c("#FF9900", "#3276FF"), 
	alpha = 0.5, 
	cex = 2, 
	cat.cex = 2, 
	cat.pos = c(-20, 20), 
	cat.dist = c(0.05, 0.05), 
	main = "Overlap Transcriptomic vs KEGG NAFLD", 
	main.cex=2
)
grid.newpage()
grid.draw(venn.plot)
dev.off()

overlap_symbols <- Human.GRCh38.p13.annot[Human.GRCh38.p13.annot$GeneID %in% overlap, c("GeneID", "Symbol", "Description")]

overlap_symbols <- overlap_symbols[order(overlap_symbols$Symbol), ]

write.csv(overlap_symbols, "new/NASH_overlap_genes.csv", row.names = FALSE, quote=F)







