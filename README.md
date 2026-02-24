# Case Study 3: Comprehensive Differential Gene Expression and Functional Enrichment Analysis in Breast cancer
## Background

Studi kasus ini bertujuan untuk melakukan analisis Differentially Expressed Genes (DEGs) pada dataset transcriptomics kanker payudara menggunakan pendekatan bioinformatika berbasis R. Dataset yang digunakan adalah GSE15852 yang diperoleh dari Gene Expression Omnibus (GEO), yang membandingkan jaringan kanker payudara dengan jaringan normal.

Kanker payudara berkembang akibat akumulasi perubahan molekuler yang memengaruhi proliferasi sel, regulasi apoptosis, dan reprogramming metabolik. Analisis transcriptomics memungkinkan identifikasi gen-gen kunci yang terlibat dalam patogenesis tersebut.

## Project Objective

Proyek ini bertujuan untuk:
- Mengidentifikasi gen yang mengalami upregulation dan downregulation pada kanker payudara
- Memvisualisasikan 50 DEGs teratas dalam bentuk heatmap
- Melakukan analisis enrichment (GO dan KEGG) untuk memahami implikasi biologis
- Menginterpretasikan perubahan molekuler yang membedakan jaringan kanker dan normal
- Mengimplementasikan pipeline analisis DEG secara reproducible menggunakan R

## Methods
Dataset
- Dataset: GSE15852
- Sumber: GEO
- Tipe data: Microarray

Analytical Environment

Analisis dilakukan menggunakan R dengan package:
- GEOquery
- limma
- dplyr
- ggplot2
- pheatmap
- umap
- clusterProfiler
- org.Hs.eg.db

Differential Expression Analysis
- Model linear: limma
- Kriteria signifikan:
  1. Adjusted p-value (FDR) < 0.05
  2. |log2FC| > 1
- Gen dengan log2FC positif → upregulated
- Gen dengan log2FC negatif → downregulated

Visualizations
- Boxplot & density plot → evaluasi distribusi data
- UMAP → pemisahan sampel
- Volcano plot → hubungan magnitude vs signifikansi
- Heatmap → 50 DEGs teratas

Functional Enrichment
- Gene Ontology (Biological Process)
- KEGG pathway analysis
- Signifikansi enrichment: adjusted p-value < 0.05

## Key Results
### Differentially Expressed Genes

Analisis menghasilkan sejumlah gen signifikan dengan perbedaan ekspresi kuat antara jaringan kanker dan normal.

Gen yang mengalami upregulation signifikan antara lain:
- KRT19
- CD24
- EPCAM
- TACSTD2
- KRT18

Gen-gen ini berkaitan dengan proliferasi sel epitel, adhesi sel, dan karakteristik tumor.

Gen yang mengalami downregulation signifikan antara lain:
- RBP4
- ACACB
- LPL
- ADIPOQ
- PPP1R1A

Sebagian besar gen ini terkait metabolisme lipid dan fungsi jaringan adiposa.

Volcano plot menunjukkan distribusi sistematis dengan banyak gen melewati ambang signifikansi dan magnitude perubahan.

### Heatmap: 50 Top DEGs

Heatmap menunjukkan pemisahan klaster yang jelas antara sampel kanker dan normal.
- Gen upregulated mendominasi ekspresi pada klaster tumor
- Gen metabolik dan adiposa cenderung tereduksi pada jaringan kanker

Pola ini menunjukkan sinyal biologis yang kuat dan konsisten antar sampel.

### Functional Enrichment Analysis
**Gene Ontology (GO)**

Enrichment menunjukkan keterlibatan gen dalam:
- Regulasi metabolisme lipid
- Proliferasi sel
- Respons hormonal
- Organisasi struktur sel

**KEGG Pathway**

Pathway signifikan meliputi:
- AMPK signaling pathway
- PPAR signaling pathway
- ECM–receptor interaction
- Adipocytokine signaling pathway

Hasil ini mendukung adanya metabolic reprogramming serta perubahan interaksi sel dengan lingkungan mikro tumor.

## Biological Interpretation

Perubahan ekspresi gen pada kanker payudara mencerminkan:
- Aktivasi jalur proliferatif dan epitelial
- Supresi jalur metabolisme lipid normal
- Disrupsi homeostasis jaringan
- Adaptasi metabolik untuk mendukung pertumbuhan tumor

Secara keseluruhan, hasil ini menunjukkan adanya reprogramming molekuler sistemik yang konsisten dengan karakteristik sel kanker.

## Skills Developed

Melalui studi ini, dikembangkan kemampuan:
- Implementasi analisis DEG berbasis limma
- Visualisasi data transcriptomics secara komprehensif
- Interpretasi volcano plot dan heatmap
- Functional enrichment analysis (GO & KEGG)
- Penyusunan pipeline analisis reproducible berbasis R
