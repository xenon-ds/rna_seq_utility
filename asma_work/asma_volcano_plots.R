library(ggrepel)

annotation_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/cross_sectional/annotation.csv", stringsAsFactors = FALSE)
head(annotation_df)

only_in_dr_negative_degs_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports/dr_negative_degs_which_never_came_in_total.csv", stringsAsFactors = FALSE)
head(only_in_dr_negative_degs_df)

dr_0_vs_2_degs_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/unannotated/dr_negative_0_vs_2hr_degs_df.csv", stringsAsFactors = FALSE )
dr_0_vs_24_degs_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/unannotated/dr_negative_0_vs_24hr_degs_df.csv", stringsAsFactors = FALSE )
dr_0_vs_96_degs_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/unannotated/dr_negative_0_vs_96hr_degs_df.csv", stringsAsFactors = FALSE )

dr_0_vs_2_degs_df = dr_0_vs_2_degs_df %>% dplyr::mutate(clean_ensembl_id = getCleanEnsembleIds(X)) %>%  dplyr::inner_join(annotation_df, by = c("clean_ensembl_id" = "ensembl_gene_id"))
dr_0_vs_24_degs_df = dr_0_vs_24_degs_df %>% dplyr::mutate(clean_ensembl_id = getCleanEnsembleIds(X)) %>%  dplyr::inner_join(annotation_df, by = c("clean_ensembl_id" = "ensembl_gene_id"))
dr_0_vs_96_degs_df = dr_0_vs_96_degs_df %>% dplyr::mutate(clean_ensembl_id = getCleanEnsembleIds(X)) %>%  dplyr::inner_join(annotation_df, by = c("clean_ensembl_id" = "ensembl_gene_id"))

get_volcano_plot(dr_0_vs_2_degs_df)
get_volcano_plot(dr_0_vs_24_degs_df)
get_volcano_plot(dr_0_vs_96_degs_df)

get_volcano_plot(dr_0_vs_2_degs_df) + geom_label_repel(data = dr_0_vs_2_degs_df %>% dplyr::filter(hgnc_symbol %in% (only_in_dr_negative_degs_df %>% dplyr::filter(time == 2) %>% dplyr::select(gene_symbol)%>% pull())) , aes( x = log2FoldChange, y = -log10(padj), label = hgnc_symbol))
get_volcano_plot(dr_0_vs_24_degs_df) + geom_label_repel(data = dr_0_vs_24_degs_df %>% dplyr::filter(hgnc_symbol %in% (only_in_dr_negative_degs_df %>% dplyr::filter(time == 24) %>% dplyr::select(gene_symbol)%>% pull())) , aes( x = log2FoldChange, y = -log10(padj), label = hgnc_symbol))
get_volcano_plot(dr_0_vs_96_degs_df) + geom_label_repel(data = dr_0_vs_96_degs_df %>% dplyr::filter(hgnc_symbol %in% (only_in_dr_negative_degs_df %>% dplyr::filter(time == 96) %>% dplyr::select(gene_symbol)%>% pull())) , aes( x = log2FoldChange, y = -log10(padj), label = hgnc_symbol))
