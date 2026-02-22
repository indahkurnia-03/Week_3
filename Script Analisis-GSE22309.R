#Title: 	Expression data from human skeletal muscle
#Dataset: GSE22309
#Platform: Microarray (Affymetrix Human Genome U95A Array)
# Tujuan: Mengidentifikasi Differential Expression Gene (DEG)

#1. Install BiocManager (manajer paket Bioconductor) 
#IF adalah struktur logika : “jika kondisi terpenuhi, lakukan aksi”

if (!require("BiocManager", quietly = TRUE))  {
  install.packages("BiocManager") 
}

# 2. Install paket Bioconductor (GEOquery & limma) 
#GEOquery: mengambil data dari database GEO 
#limma: analisis statistik ekspresi gen 

BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE) 

#Install annotation package sesuai platform
#GPL91: Affymetrix Human Genome U95A
BiocManager::install("hgu95a.db", ask = FALSE, update = FALSE)

#3. Install paket CRAN untuk visualisasi dan manipulasi data 
#phetmap: heatmap ekspresi gen 
#ggplot2: grafik (volcano plot)
#dplyr: manipulasi tabel data 

install.packages(c("pheatmap", "ggplot2", "dplyr"))

#4. Memanggil library 
#library() digunakan agar fungsi di dalam package bisa digunakan 
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu95a.db)
library(AnnotationDbi)
library(umap)

#PART C. PENGAMBILAN DATA DARI GEO 

Data GEO tidak didownload langsung karena mengalami error. setdata GEO diuduh 
diunduh manual kemudian dimasukkan dalam R

gse <- getGEO(filename = "GSE22309_series_matrix.txt.gz")

#PART D. PRE-PROCESSING DATA EKSPRESI 

ex <- exprs(gse)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

#LogTransform

LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

#IF statement:
#Jika LogTransform = TRUE, maka lakukan log2
if (LogTransform) {
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

#PART E. DEFINISI KELOMPOK SAMPEL

#pData(): metadata sampel
group_info <- pData(gse)[["agent:ch1"]]

#make.names():
groups <- make.names(group_info)

#factor():
gse$group <- factor(groups)

#levels():
nama_grup <- levels(gse$group)
print(nama_grup)

#PART F. DESIGN MATRIX (KERANGKA STATISTIK) 

#model.matrix():
design <- model.matrix(~0 + gse$group)

#colnames():
colnames(design) <- levels(gse$group)

#Menentukan perbandingan biologis

grup_insulin    <- nama_grup[1]
grup_untreated  <- nama_grup[2]

contrast_formula <- paste(grup_insulin, "-", grup_untreated)
print(paste("Kontras yang dianalisis:", contrast_formula))

#PART G. ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)

#lmFit():
#Membangun model linear untuk setiap gen
fit <- lmFit(ex, design)

#makeContrasts(): mendefinisikan perbandingan antar grup
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

#contrasts.fit(): menerapkan kontras ke model
fit2 <- contrasts.fit(fit, contrast_matrix)

#eBayes():
#Empirical Bayes untuk menstabilkan estimasi varians
fit2 <- eBayes(fit2)

#topTable():
#Mengambil hasil akhir DEG
#adjust = "fdr" -> koreksi multiple testing
#p.value = 0.01  -> gen sangat signifikan
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)

head(topTableResults)


#PART H. ANOTASI NAMA GEN 

#Mengambil ID probe dari hasil DEG
probe_ids <- rownames(topTableResults)

#Mapping probe -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(
  hgu95a.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

#Gabungkan dengan hasil limma
topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

#Cek hasil anotasi
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])
#PART I.1 BOXPLOT DISTRIBUSI NILAI EKSPRESI

group_colors <- as.numeric(gse$group)

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gse$group),
  fill = unique(group_colors),
  cex = 0.8
)
#PART I.2 DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT) 

#Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gse$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )

#PART I.3 UMAP (VISUALISASI DIMENSI RENDAH)

umap_input <- t(ex)

#Jalankan UMAP
umap_result <- umap(umap_input)

#Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gse$group
)

#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

#Klasifikasi status gen

volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 0 & volcano_data$adj.P.Val < 0.01] <- 
  "UP"
volcano_data$status[volcano_data$logFC < 0 & volcano_data$adj.P.Val < 0.01] <- 
  "DOWN"

#Visualisasi

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(0, 0), linetype = "dashed") + # Sesuaikan garis vline
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG Insulin vs Untreated")

#PART J.2 VISUALISASI HEATMAP

#Pilih 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)

#Ambil matriks ekspresi untuk gen terpilih
mat_heatmap <- ex[top50$PROBEID, ]

#Gunakan Gene Symbol (fallback ke Probe ID)
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,      # jika SYMBOL kosong → probe ID
  top50$SYMBOL        # jika ada → gene symbol
)

rownames(mat_heatmap) <- gene_label

#Pembersihan data (WAJIB agar tidak error hclust)
#Hapus baris dengan NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

#Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

#Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

#Anotasi kolom (kelompok sampel)
annotation_col <- data.frame(
  Group = gse$group
)

rownames(annotation_col) <- colnames(mat_heatmap)

#Visualisasi heatmap 
pheatmap(
  mat_heatmap,
  scale = "row",                 # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE,         # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

#PART K. MENYIMPAN HASIL 

# write.csv(): menyimpan hasil analisis ke file CSV
write.csv(topTableResults, "Hasil_GSE22309csv")

message("Analisis selesai. File hasil telah disimpan.")

