library(dplyr)

#################################
# loading cross sectional data
#################################
cross_0hr_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/cross_sectional/counts_total_0hr_vs_dr_negative_0hr.csv", stringsAsFactors = FALSE)
cross_2hr_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/cross_sectional/counts_total_2hr_vs_dr_negative_2hr.csv", stringsAsFactors = FALSE)
cross_24hr_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/cross_sectional/counts_total_24hr_vs_dr_negative_24hr.csv", stringsAsFactors = FALSE)
cross_96hr_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/cross_sectional/counts_total_96hr_vs_dr_negative_96hr.csv", stringsAsFactors = FALSE)

#################################
# adding time column to asmas data
#################################
cross_0hr_df = cross_0hr_df %>% dplyr::mutate(time = 0)
cross_2hr_df = cross_2hr_df %>% dplyr::mutate(time = 2)
cross_24hr_df = cross_24hr_df %>% dplyr::mutate(time = 24)
cross_96hr_df = cross_96hr_df %>% dplyr::mutate(time = 96)

head(cross_0hr_df)
head(cross_2hr_df)
head(cross_24hr_df)
head(cross_96hr_df)

#for_viewing_counts_df  = rbind.data.frame(cross_0hr_df, cross_2hr_df, cross_24hr_df, cross_96hr_df)
#################################
# modifying column names 
#################################

fix_asma_data_columnm_names = function(counts_df, time = ""){
  names(counts_df) = gsub("X.media.cidr.d7416dce.cdf6.43ed.9df7.978f8a9438d8.databackup_2019_07_11.work.asma_whole_blood_rnaseq.asma_durbar_analysis_raw_reads_26_jun_2020.sam_files.sam_","",names(counts_df))
  names(counts_df) = gsub("_trimmed_Aligned.out.sam","",names(counts_df))
  return(counts_df)
}

cross_0hr_df = fix_asma_data_columnm_names(cross_0hr_df)
cross_2hr_df = fix_asma_data_columnm_names(cross_2hr_df)
cross_24hr_df = fix_asma_data_columnm_names(cross_24hr_df)
cross_96hr_df = fix_asma_data_columnm_names(cross_96hr_df)

cross_0hr_df = cross_0hr_df %>% tidyr::gather(key = "sample_id",value = "counts", - X,-time)
cross_2hr_df = cross_2hr_df %>% tidyr::gather(key = "sample_id",value = "counts", - X,-time)
cross_24hr_df = cross_24hr_df %>% tidyr::gather(key = "sample_id",value = "counts", - X,-time)
cross_96hr_df = cross_96hr_df %>% tidyr::gather(key = "sample_id",value = "counts", - X,-time)

combined_cross_data = rbind.data.frame(cross_0hr_df, cross_2hr_df, cross_24hr_df, cross_96hr_df)
for_viewing_counts_df  = rbind.data.frame(cross_0hr_df, cross_2hr_df, cross_24hr_df, cross_96hr_df)
combined_cross_data = combined_cross_data %>% dplyr::mutate(clean_ensembl_id = getCleanEnsembleIds(X))

##############################
# getting annottaion from biomart 
##############################
# ensembl_ids = unique(combined_cross_data$clean_ensembl_id)
# length(ensembl_ids)
# biomart_description_df = getAnnotationFromBiomart(ids = ensembl_ids, filter = "ensembl_gene_id", columns = c("ensembl_gene_id", "hgnc_symbol", "description"))
# write.csv(biomart_description_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/cross_sectional/annotation.csv", row.names = FALSE)
biomart_description_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/cross_sectional/annotation.csv", stringsAsFactors = FALSE)
head(combined_cross_data) %>% dplyr::inner_join(biomart_description_df, by = c("clean_ensembl_id" = "ensembl_gene_id"))
combined_cross_data = combined_cross_data %>% dplyr::inner_join(biomart_description_df, by = c("clean_ensembl_id" = "ensembl_gene_id"))

head(for_viewing_counts_df) %>% dplyr::mutate(clean_ensembl_id = getCleanEnsembleIds(X)) %>% dplyr::inner_join(biomart_description_df, by = c("clean_ensembl_id" = "ensembl_gene_id"))

##############################
# importing run info
##############################
run_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/metadata/run_df.csv", stringsAsFactors = FALSE)
run_df = run_df %>% dplyr::select(Run, Library.Name, Sample.Name, sample_type)
head(run_df)
names(run_df)
head(combined_cross_data) %>% dplyr::inner_join(run_df, by = c("sample_id" = "Run"))

combined_cross_data = combined_cross_data %>% dplyr::inner_join(run_df, by = c("sample_id" = "Run"))

write.csv(combined_cross_data, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_data_for_plotting/counts_data_for_plotting.csv", row.names = FALSE)
gene_symbol = "CCL4"

combined_cross_data = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_data_for_plotting/counts_data_for_plotting.csv", stringsAsFactors = FALSE)
combined_cross_data %>% dplyr::filter(hgnc_symbol == gene_symbol)%>% dplyr::mutate(sample = factor(sample_type, levels = c("total","dr_negative"))) %>% ggplot(aes(x = factor(time), y = counts, color = Sample.Name)) + geom_point(alpha = 0.4) +geom_line(aes(group = Sample.Name))+ scale_y_log10()+ facet_wrap(~sample) + ggtitle(gene_symbol) + theme_bw()

gene_symbols = c("IL2","IFNG","IL17A","IL22","CCL3L3","CCL3","CCL4","CD274","CD46","TRAF1","TRAF3","FASLG","TNFRSF8")
for(gene_symbol in sort(gene_symbols))
{
  print(gene_symbol)
  print(combined_cross_data %>% dplyr::filter(hgnc_symbol == gene_symbol)%>% dplyr::mutate(sample = factor(sample_type, levels = c("total","dr_negative"))) %>% ggplot(aes(x = factor(time), y = counts, color = Sample.Name)) + geom_point(alpha = 0.4) +geom_line(aes(group = Sample.Name))+ scale_y_log10()+ facet_wrap(~sample) + ggtitle(gene_symbol) + theme_bw())
}

head(for_viewing_counts_df)
generate_counts_report = function(gene_symbols, excel_file_path)
{
  loaded_packages = (.packages())
  if("xlsx" %in% loaded_packages)
  {
    detach("package:xlsx", unload = TRUE)
  }
  if(!("openxlsx" %in% loaded_packages))
  {
    library("openxlsx")
  }
  if(!("dplyr" %in% loaded_packages))
  {
    library("dplyr")
  }
  
  if(!("ggplot2" %in% loaded_packages))
  {
    library("ggplot2")
  }
  combined_cross_data = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_data_for_plotting/counts_data_for_plotting.csv", stringsAsFactors = FALSE)
  wb = createWorkbook()
  for(gene_symbol in sort(gene_symbols))
  {
    print(gene_symbol)
    gene_counts_df = combined_cross_data %>% dplyr::filter(hgnc_symbol == gene_symbol) %>% dplyr::mutate(sample_group = paste(Sample.Name, sample_type, sep = "_")) %>% dplyr::select(hgnc_symbol, description, counts, time, sample_group,) %>% tidyr::spread(sample_group, counts) %>% dplyr::select(hgnc_symbol, description, time, SC168_total,SC169_total,SC170_total,SC176_total,SC189_total,SC168_dr_negative,SC169_dr_negative,SC170_dr_negative,SC176_dr_negative,SC189_dr_negative)
    if(nrow(gene_counts_df)!=0)
    {
      plot = combined_cross_data %>% dplyr::filter(hgnc_symbol == gene_symbol)%>% dplyr::mutate(sample = factor(sample_type, levels = c("total","dr_negative"))) %>% ggplot(aes(x = factor(time), y = counts, color = Sample.Name)) + geom_point(alpha = 0.4) +geom_line(aes(group = Sample.Name))+ scale_y_log10()+ facet_wrap(~sample) + ggtitle(gene_symbol) + theme_bw()
      print(plot)
      plot_path = paste("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/plots/",gene_symbol,".jpg", sep = "")
      ggsave(filename = plot_path,dpi = 400)
      addWorksheet(wb, sheetName = gene_symbol)
      insertImage(wb = wb, sheet = gene_symbol, file = plot_path, startCol = 5, startRow = 10,width=25,height = 15, units = "cm" ,dpi = 400)
      #********************************
      #adding sheet to workbook
      #********************************
      writeData(wb=wb,sheet=gene_symbol,x=gene_counts_df)
    }
  }
  saveWorkbook(wb,file = excel_file_path,overwrite = TRUE)
  
  #*****************************************
  # unloading the openxlsx package
  #*****************************************
  print("unloading the openxlsx package")
  detach("package:openxlsx", unload = TRUE)
  
  gc()
}
generate_counts_report(gene_symbols = gene_symbols, excel_file_path = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/plots_excel_reports/asma_selected_gene_counts.xlsx")
asma_selected_df = combined_cross_data %>% dplyr::filter(hgnc_symbol %in% gene_symbols) %>% dplyr::mutate(sample_group = paste(Sample.Name, sample_type, sep = "_")) %>% dplyr::select(hgnc_symbol, description, counts, time, sample_group,) %>% tidyr::spread(sample_group, counts) %>% dplyr::select(hgnc_symbol, description, time, SC168_total,SC169_total,SC170_total,SC176_total,SC189_total,SC168_dr_negative,SC169_dr_negative,SC170_dr_negative,SC176_dr_negative,SC189_dr_negative)
write.csv(asma_selected_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_data_for_plotting/asma_selected_gene_counts.csv", row.names = FALSE)
