library(DESeq2)

####################################################################################################################################
#analysis with vaccination status as a contributing factor
#the analysis excludes eptb , CD154 negative , and CD154 poitive pooled samples
####################################################################################################################################

########################################
#loading ful counts matrix file
########################################
temp_combined_counts_for_deseq = read.delim("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/counts/combined_counts.csv",sep = ",",stringsAsFactors = FALSE, row.names = "X")
head(temp_combined_counts_for_deseq)
############################################
#changing the columnnames of the counts matrix
############################################
names(temp_combined_counts_for_deseq) = gsub("X.media.cidr.d7416dce.cdf6.43ed.9df7.978f8a9438d8.databackup_2019_07_11.work.SO_9292.bam_files.",x =  names(temp_combined_counts_for_deseq), replacement = "")

############################################
#loading metadata only for group 1,2,3,4
############################################
metadata_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/metadata/vaccination_status_metadata.csv", stringsAsFactors = FALSE)
vaccination_status_metadata = metadata_df %>% dplyr::filter(Group %in% 1:4)

############################################
#loading counts matrix for bam files belonging to group 1,2,3,4
############################################
counts_matrix = temp_combined_counts_for_deseq[, colnames(temp_combined_counts_for_deseq) %in% vaccination_status_metadata$bam_files ]

############################################
#creating vector corresponding to the sample bam files. grouping info will be used in the design matrix
############################################
grouping = factor(vaccination_status_metadata$vaccination_status)
names(grouping) = vaccination_status_metadata$bam_files
grouping = grouping[colnames(counts_matrix)]

############################################
#creating group name vector corresponding to the sample bam files. will be used for visualization
############################################
group_names = vaccination_status_metadata$Group
names(group_names) = vaccination_status_metadata$bam_files
group_names = group_names[colnames(counts_matrix)]

############################################
#creating design matrix
############################################
design = model.matrix(~vaccination_status, data = vaccination_status_metadata)

dds = DESeqDataSetFromMatrix(counts_matrix, colData = vaccination_status_metadata, design = formula(~ vaccination_status))
#deseq_res = DESeq(dds, test = "LRT", reduced = ~ cell_type + vaccination_status)
deseq_res = DESeq(dds)
resultsNames(deseq_res)
DESeq2::plotDispEsts(deseq_res)
res = results(deseq_res, alpha = 0.05)
mcols(res)
DESeq2::plotMA(res)
res = lfcShrink(deseq_res, coef = "vaccination_status_vaccinated_vs_not_vaccinated", res = res)
degs_df = as.data.frame(res)
degs_df = annotate_degs_df_offline(degs_df)
write.csv(degs_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/deseq_vaccinated_vs_unvaccinated_without_eptb_and_cd154_negative_and_cd154_pooled_samples/unfiltered_vaccinated_vs_unvaccinated.csv", row.names = FALSE)

degs_df = degs_df %>% dplyr::filter(padj<0.05) %>% dplyr::arrange(desc(log2FoldChange))
nrow(degs_df)

write.csv(degs_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/deseq_vaccinated_vs_unvaccinated_without_eptb_and_cd154_negative_and_cd154_pooled_samples/vaccinated_vs_unvaccinated.csv", row.names = FALSE)
get_degs_with_pathways(degs_csv_path ="/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/deseq_vaccinated_vs_unvaccinated_without_eptb_and_cd154_negative_and_cd154_pooled_samples/vaccinated_vs_unvaccinated.csv", excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/excel_reports/vaccinated_vs_unvaccinated", excel_report_filename = "vaccinated_vs_unvaccinated.xlsx" )
