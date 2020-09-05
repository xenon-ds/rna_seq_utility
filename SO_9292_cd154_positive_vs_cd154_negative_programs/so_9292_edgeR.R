library(edgeR)
library(dplyr)

########################################
#loading counts matrix file
########################################
temp_combined_counts_for_edgeR = read.delim("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/counts/combined_counts.csv",sep = ",",stringsAsFactors = FALSE, row.names = "X")
head(temp_combined_counts_for_edgeR)

########################################
#loading metadata file
########################################
metadata_edgeR = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/metadata/vaccination_status_metadata.csv", stringsAsFactors = FALSE)

############################################
#changing the columnnames of the counts matrix
############################################  
names(temp_combined_counts_for_edgeR) = gsub("X.media.cidr.d7416dce.cdf6.43ed.9df7.978f8a9438d8.databackup_2019_07_11.work.SO_9292.bam_files.",x =  names(temp_combined_counts_for_edgeR), replacement = "")


############################################
#creating values groups and cell types for the current column names of the counts matrix 
############################################
groups_edgeR = sapply(names(temp_combined_counts_for_edgeR), function(x){
  return(metadata_edgeR %>% dplyr::filter(bam_files == x) %>% pull(Group))
})

cell_type_edgeR = sapply(names(temp_combined_counts_for_edgeR), function(x){
  return(metadata_edgeR %>% dplyr::filter(bam_files == x) %>% pull(cell_type))
})

####################################################################################################################################
#analysis with group as a contributing factor
####################################################################################################################################
edgeR_object = DGEList(temp_combined_counts_for_edgeR, group = groups_edgeR)
edgeR_object = calcNormFactors(edgeR_object)
plotMDS.DGEList(edgeR_object, labels = groups_edgeR)
############################################
#creating design matrix
############################################
design = model.matrix(~0+groups_edgeR)
############################################
#visualizing the dispersion within the data
############################################
edgeR_object = estimateDisp(edgeR_object, design = design, robust = TRUE)
plotBCV(edgeR_object)
############################################
# fitting the model
############################################
fit = glmQLFit(edgeR_object, design)
############################################
# checking the parameters in the model
############################################
colnames(fit)

####################################################################################################################################
# conducting significance analysis tests
####################################################################################################################################
############################################
# getting contrast for group 1 cd154+ vs cd154-
############################################
qlf <- glmQLFTest(fit, contrast = makeContrasts(groups_edgeR1-groups_edgeRGroup1, levels = design))
res = as.data.frame(topTags(qlf, nrow(qlf$table)))
head(res)
#nrow(as.data.frame(res)[res$FDR<0.05,])
res = annotate_degs_df_offline(res, filter_by_padj = T, degs_coming_from = "edgeR", arranged = FALSE)
res_group_1_cd154_positive_negative = res
head(res_group_1_cd154_positive_negative)
write.csv(res_group_1_cd154_positive_negative, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/edgeR_degs/groupwise_cd154_positive_vs_cd154_negative/group_1_cd154_positive_vs_cd154_negative.csv", row.names = FALSE)

############################################
# getting contrast for group 2 cd154+ vs cd154-
############################################
qlf <- glmQLFTest(fit, contrast = makeContrasts(groups_edgeR2-groups_edgeRGroup2, levels = design))
res = as.data.frame(topTags(qlf, nrow(qlf$table)))
head(res)
#nrow(as.data.frame(res)[res$FDR<0.05,])
res = annotate_degs_df_offline(res, filter_by_padj = T, degs_coming_from = "edgeR", arranged = FALSE)
res_group_2_cd154_positive_negative = res
head(res_group_2_cd154_positive_negative)
write.csv(res_group_2_cd154_positive_negative, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/edgeR_degs/groupwise_cd154_positive_vs_cd154_negative/group_2_cd154_positive_vs_cd154_negative.csv", row.names = FALSE)

############################################
# getting contrast for group 3 cd154+ vs cd154-
############################################
qlf <- glmQLFTest(fit, contrast = makeContrasts(groups_edgeR3-groups_edgeRGroup3, levels = design))
res = as.data.frame(topTags(qlf, nrow(qlf$table)))
head(res)
#nrow(as.data.frame(res)[res$FDR<0.05,])
res = annotate_degs_df_offline(res, filter_by_padj = T, degs_coming_from = "edgeR", arranged = FALSE)
res_group_3_cd154_positive_negative = res
head(res_group_3_cd154_positive_negative)
write.csv(res_group_3_cd154_positive_negative, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/edgeR_degs/groupwise_cd154_positive_vs_cd154_negative/group_3_cd154_positive_vs_cd154_negative.csv", row.names = FALSE)

############################################
# getting contrast for group 4 cd154+ vs cd154-
############################################
qlf <- glmQLFTest(fit, contrast = makeContrasts(groups_edgeR4-groups_edgeRGroup4, levels = design))
res = as.data.frame(topTags(qlf, nrow(qlf$table)))
head(res)
#nrow(as.data.frame(res)[res$FDR<0.05,])
res = annotate_degs_df_offline(res, filter_by_padj = T, degs_coming_from = "edgeR", arranged = FALSE)
res_group_4_cd154_positive_negative = res
head(res_group_4_cd154_positive_negative)
write.csv(res_group_4_cd154_positive_negative, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/edgeR_degs/groupwise_cd154_positive_vs_cd154_negative/group_4_cd154_positive_vs_cd154_negative.csv", row.names = FALSE)

############################################
# getting contrast for group 5 cd154+ vs cd154-
############################################
qlf <- glmQLFTest(fit, contrast = makeContrasts(groups_edgeREPTB-groups_edgeRGroup5, levels = design))
res = as.data.frame(topTags(qlf, nrow(qlf$table)))
head(res)
#nrow(as.data.frame(res)[res$FDR<0.05,])
res = annotate_degs_df_offline(res, filter_by_padj = T, degs_coming_from = "edgeR", arranged = FALSE)
res_group_5_cd154_positive_negative = res
head(res_group_5_cd154_positive_negative)
write.csv(res_group_5_cd154_positive_negative, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/edgeR_degs/groupwise_cd154_positive_vs_cd154_negative/group_5_cd154_positive_vs_cd154_negative.csv", row.names = FALSE)


res_group_1_cd154_positive_negative %>% dplyr::inner_join(res_group_2_cd154_positive_negative, by = c("ensembl_id" = "ensembl_id"))%>% dplyr::inner_join(res_group_3_cd154_positive_negative, by = c("ensembl_id" = "ensembl_id"))%>% dplyr::inner_join(res_group_4_cd154_positive_negative, by = c("ensembl_id" = "ensembl_id"))%>% dplyr::inner_join(res_group_5_cd154_positive_negative, by = c("ensembl_id" = "ensembl_id"))

nrow(res_group_1_cd154_positive_negative)
nrow(res_group_2_cd154_positive_negative)
nrow(res_group_3_cd154_positive_negative)
nrow(res_group_4_cd154_positive_negative)
nrow(res_group_5_cd154_positive_negative)

####################################################################################################################################
#analysis with cell_type as a contributing factor
#the analysis is conducted with cd154 positive and negative 
####################################################################################################################################

cell_type_metadata = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/metadata/vaccination_status_metadata.csv", stringsAsFactors = FALSE)
not_eptb_not_cd154_positive_pooled_selected_samples = sapply(names(temp_combined_counts_for_edgeR), function(x){
  #print(so_9292_metadata_df %>% dplyr::filter(bam_files == x))
  group_name = as.character(cell_type_metadata %>% dplyr::filter(bam_files == x) %>% dplyr::select(Group) %>% dplyr::pull())
  not_eptb_not_cd154_positive_pooled_groups = c(1:4, paste("Group", 1:4,sep = ""))
  if(group_name %in% not_eptb_not_cd154_positive_pooled_groups){
    #print("selected")
    return(TRUE)
  }
  #print("NOT selected")
  print(metadata_edgeR %>% dplyr::filter(bam_files == x))
  return(FALSE)
})

not_eptb_not_cd154_positive_pooled_selected_samples_metadata_df = cell_type_metadata %>% dplyr::filter(bam_files %in% names(temp_combined_counts_for_edgeR[,not_eptb_not_cd154_positive_pooled_selected_samples])) %>% dplyr::select(bam_files,Group, cell_type) %>% dplyr::mutate(custom_group = ifelse(Group %in% 1:4, paste(Group,cell_type, sep = "."), cell_type))

rownames(not_eptb_not_cd154_positive_pooled_selected_samples_metadata_df) = not_eptb_not_cd154_positive_pooled_selected_samples_metadata_df$bam_files

not_eptb_not_cd154_positive_pooled_counts_matrix = temp_combined_counts_for_edgeR[,not_eptb_not_cd154_positive_pooled_selected_samples]
not_eptb_not_cd154_positive_pooled_selected_samples_metadata_df$custom_group = relevel(factor(not_eptb_not_cd154_positive_pooled_selected_samples_metadata_df$custom_group), ref = "CD154_negative")


design = model.matrix(~custom_group, data = not_eptb_not_cd154_positive_pooled_selected_samples_metadata_df)

edgeR_object = DGEList(not_eptb_not_cd154_positive_pooled_counts_matrix, group = not_eptb_not_cd154_positive_pooled_selected_samples_metadata_df[names(not_eptb_not_cd154_positive_pooled_counts_matrix),]$custom_group)
edgeR_object = calcNormFactors(edgeR_object)
plotMDS.DGEList(edgeR_object, labels = not_eptb_not_cd154_positive_pooled_selected_samples_metadata_df[names(not_eptb_not_cd154_positive_pooled_counts_matrix),]$Group)

edgeR_object = estimateDisp(edgeR_object, design = design, robust = TRUE)
plotBCV(edgeR_object)

fit = glmQLFit(edgeR_object, design)
colnames(fit)

for(i in 2:5)
{
  qlf <- glmQLFTest(fit, coef=i)
  res = as.data.frame(topTags(qlf, n = nrow(qlf$table)))
  head(res)
  #nrow(as.data.frame(res)[res$FDR<0.05,])
  res = annotate_degs_df_offline(res, filter_by_padj = T, degs_coming_from = "edgeR", arranged = TRUE)
  print(nrow(res))
  write.csv(res, paste("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/edgeR_degs/degs_with_combined_cd154_negative_samples/cd154_positive_vs_cd154_negative_not_eptb_not_pooled/group_",i-1,"_cd154_positive_vs_combined_cd154_negative.csv"), row.names = FALSE)
}

####################################################################################################################################
#analysis with vaccination status as a contributing factor
#the analysis excludes eptb , CD154 negative , and CD154 poitive pooled samples
####################################################################################################################################

############################################
#loading metadata only for group 1,2,3,4
############################################
metadata_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/metadata/vaccination_status_metadata.csv", stringsAsFactors = FALSE)
vaccination_status_metadata = metadata_df %>% dplyr::filter(Group %in% 1:4)

############################################
#loading counts matrix for bam files belonging to group 1,2,3,4
############################################
counts_matrix = temp_combined_counts_for_edgeR[, colnames(temp_combined_counts_for_edgeR) %in% vaccination_status_metadata$bam_files ]

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

############################################
#initialising edger objects
############################################
edgeR_object = DGEList(counts_matrix, group = grouping)
edgeR_object = calcNormFactors(edgeR_object)
plotMDS.DGEList(edgeR_object, labels = group_names)

############################################
#estimating dispersion
############################################
edgeR_object = estimateDisp(edgeR_object, design = design, robust = TRUE)
plotBCV(edgeR_object)

############################################
#fitting model
############################################
fit = glmQLFit(edgeR_object, design)
colnames(fit)

############################################
#performing significance test
############################################
qlf <- glmQLFTest(fit, coef=2)

############################################
#annotating results
############################################
res = as.data.frame(topTags(qlf, n = nrow(qlf$table)))
head(res)
#nrow(as.data.frame(res)[res$FDR<0.05,])
res = annotate_degs_df_offline(res, filter_by_padj = T, degs_coming_from = "edgeR", arranged = TRUE)
print(nrow(res))
write.csv(res, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/edgeR_degs/degs_with_combined_cd154_negative_samples/vaccinated_vs_unvaccinated_withcd154_positive_not_eptb_notpooled/vaccinated_vs_unvaccinated.csv", row.names = FALSE)

################################################################################################################################
#generating degs with pathways excel reports for cd154+ and negative files 
################################################################################################################################
degs_file_paths = list.files("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/edgeR_degs/degs_with_combined_cd154_negative_samples/cd154_positive_vs_cd154_negative_not_eptb_not_pooled", pattern = "[.]csv", full.names = T)
degs_file_names = list.files("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/edgeR_degs/degs_with_combined_cd154_negative_samples/cd154_positive_vs_cd154_negative_not_eptb_not_pooled/", pattern = "[.]csv")
excel_report_names = gsub(x = degs_file_names, pattern = "csv", replacement = "xlsx")

for(i in 1:length(degs_file_paths))
{
  get_degs_with_pathways(degs_csv_path = degs_file_paths[i], excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/excel_reports/edgeR/cd154_positive_vs_cd154_negative_not_eptb_not_pooled/", excel_report_filename = excel_report_names[i])
}

################################################################################################################################
#generating degs with pathways excel reports for group wise cd154+ and negative files 
#note that this comparison is 5 vs 1 sample
################################################################################################################################
degs_file_paths = list.files("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/edgeR_degs/groupwise_cd154_positive_vs_cd154_negative/", pattern = "[.]csv", full.names = T)
degs_file_names = list.files("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/degs_csv/edgeR_degs/groupwise_cd154_positive_vs_cd154_negative/", pattern = "[.]csv")
excel_report_names = gsub(x = degs_file_names, pattern = "csv", replacement = "xlsx")

for(i in 1:length(degs_file_paths))
{
  get_degs_with_pathways(degs_csv_path = degs_file_paths[i], excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/SO_9292/excel_reports/edgeR/group_wise_cd154_positive_vs_cd154_negative/", excel_report_filename = excel_report_names[i])
}