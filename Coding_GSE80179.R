library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(illuminaHumanv4.db)
library(AnnotationDbi)
library(umap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(KEGGREST)

gset <- getGEO("GSE80179", GSEMatrix=TRUE)[[1]]

ex <- exprs(gset)
pd <- pData(gset)

summary(ex[1:5,1:5])

ex <- log2(ex - min(ex) + 1)
quantile(ex)

group <- factor(pd$`category:ch1`)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

contrast_matrix <- makeContrasts(Control - RSV, levels=design)

fit <- lmFit(ex, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, number=Inf)

sum(results$adj.P.Val < 0.05)

summary(ex)

topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  number = Inf,)
head(topTableResults)

probe_ids <- rownames(topTableResults)


gene_annotation <- AnnotationDbi::select(illuminaHumanv4.db, keys = probe_ids, columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID")


topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE)


head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])


#Boxplot

group_colors <- as.numeric(group)

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2(ex-min(ex)+1))")

legend(
  "topright",
  legend = levels(group),
  fill = unique(group_colors),
  cex = 0.8)

#Density plot

expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(group, each = nrow(ex)))

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen (Control)",
    x = "Expression Value (log2(ex-min(ex)+1))",
    y = "Density")


#UMAP

umap_input <- t(ex)

umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = group)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot: RSV vs Control",
    x = "UMAP 1",
    y = "UMAP 2")


#Volcano plot 

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL)

volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 0.1 & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -0.1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG RSV")


#Heatmap 

topTableResults <- topTableResults[order(topTableResults$adj.P.Val)]


top50 <- head(topTableResults,50)

mat_heatmap <- ex[top50$PROBEID,]

gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,      
  top50$SYMBOL)

rownames(mat_heatmap) <- gene_label

mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,]

gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

annotation_col <- data.frame(
  Group = group)

rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale = "row",            
  annotation_col = annotation_col,
  show_colnames = FALSE,         
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes")

#GO KEGG

gene_list <- topTableResults[["SYMBOL"]]
# Ambil gen signifikan saja
sig_genes <- subset(topTableResults,
                    adj.P.Val < 0.01 & abs(logFC) > 2)
length(sig_genes)
genes <- sig_genes$SYMBOL

# Convert ke ENTREZ
gene_df <- bitr(
  genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

entrez_ids <- gene_df$ENTREZID

length(entrez_ids)
gene_df <- bitr(
  gene_list,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db)

head(gene_df)

gene_list <- gene_df$ENTREZID

ego_bp <- enrichGO(
  gene          = gene_list,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

head(ego_bp)

##KEGG
ekegg <- enrichKEGG(
  gene         = gene_list,
  organism     = "hsa",   # human
  pvalueCutoff = 0.05
)

head(ekegg)

ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db)

dotplot(ego_bp, showCategory=15)
dotplot(ekegg, showCategory=15)
barplot(ego_bp, showCategory=10)

length(entrez_ids)

entrez_ids <- gene_df$ENTREZID

kegg_ids <- paste0("hsa:", entrez_ids)

ekegg <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa")

dotplot(ekegg)

#Menyimpan hasil

write.csv(topTableResults, "Hasil_GSE80179_DEG.csv")
write.csv(top50, "output_GSE80179_top50.csv")
write.csv(gene_annotation, "output_GSE80179_volcano data.csv")
message("Analisis selesai. File hasil telah disimpan.")
