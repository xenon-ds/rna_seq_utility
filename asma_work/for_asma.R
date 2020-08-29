library(stringr)
library(Rsubread)


differentialExpressionAnalysis=function(X,Y,gtfFilePath,countsFilePath,differentialExpressionResultsFilePath)
{
  library(Rsubread)
  #***********************************************************
  #counting raw reads
  #***********************************************************
  rawReadCounts=featureCounts(files=c(X,Y),GTF.featureType = "exon", annot.ext = gtfFilePath,isGTFAnnotationFile = T,nthreads = 20)
  
  #***********************************************************
  #saving raw reads into a file
  #***********************************************************
  df=as.data.frame(rawReadCounts$counts)
  print(paste("counts data : ",countsFilePath))
  write.csv(df,file = countsFilePath,row.names = T,col.names = T)
  colnames(df)
  
  
  #***********************************************************
  #conducting differential expression analysis
  #***********************************************************
  colData = DataFrame(condition=factor(c("ctrl","ctrl", "ctrl","ctrl", "ctrl","treat", "treat", "treat", "treat", "treat")))
  dds = DESeqDataSetFromMatrix(rawReadCounts$counts, colData, formula(~ condition))
  
  # run DEseq
  res= DESeq(dds)
  
  #plotMA(res)
  
  #diffrentiallyExpressed 
  deGenes <- results(res,alpha=0.05)
  
  # order by BH adjusted p-value
  deGenesOrdered <- deGenes[order(deGenes$padj),]
  #deGenesOrdered=deGenesOrdered[!is.na(deGenesOrdered$padj),]
  #deGenesOrdered <- deGenesOrdered[deGenesOrdered$padj<0.05,]
  print(paste("diffex data data : ",differentialExpressionResultsFilePath))
  write.csv(deGenesOrdered,file=differentialExpressionResultsFilePath,row.names = T,col.names = T)
  # top of ordered matrix
  #head(deGenesOrdered)
}

generate_counts_file=function(X,Y,gtfFilePath = "/home/cidr/Documents/work/test_analysis/homoSapien/forHisat/gencode.v31.annotation.gtf",countsFilePath)
{
  library(Rsubread)
  #***********************************************************
  #counting raw reads
  #***********************************************************
  rawReadCounts=featureCounts(files=c(X,Y),GTF.featureType = "exon", annot.ext = gtfFilePath,isGTFAnnotationFile = T,nthreads = 20)
  
  #***********************************************************
  #saving raw reads into a file
  #***********************************************************
  df=as.data.frame(rawReadCounts$counts)
  print(paste("counts data : ",countsFilePath))
  write.csv(df,file = countsFilePath,row.names = T,col.names = T)
}


##########################################################################################################################
#loading metadata
##########################################################################################################################
run_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/metadata/SraRunTable.txt", stringsAsFactors = FALSE)
head(run_df)

time_point = sapply(run_df$Library.Name, function(x){
  if(str_detect(x, "0hr"))
  {
    return(0)
  }
  if(str_detect(x, "2hr"))
  {
    return(2)
  }
  if(str_detect(x, "24hr"))
  {
    return(24)
  }
  if(str_detect(x, "4D"))
  {
    return(96)
  }
})

sample_type = sapply(run_df$Library.Name, function(x){
  if(str_detect(x, "DR"))
  {
    return("dr_negative")
  }
  if(str_detect(x, "Total"))
  {
    return("total")
  }
})

run_df = cbind.data.frame(run_df, time_point, sample_type)
head(run_df)


##########################################################################################################################
#
##########################################################################################################################
sam_files = list.files(path = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/sam_files", full.names = TRUE)
file_path = sapply(run_df$Run, function(x){
  return(sam_files[str_detect(sam_files, x)])
})
run_df = cbind.data.frame(run_df, file_path)

#run_df %>% dplyr::filter(time_point == 96, sample_type == "dr_negative" )
write.csv(run_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/metadata/run_df.csv", row.names = FALSE)

##########################################################################################
#counting
##########################################################################################
#************************
#initialising the ids
#************************
total_0hr_run_ids = as.vector(run_df %>% dplyr::filter(time_point == 0, sample_type == "total") %>% dplyr::select(Library.Name) %>% dplyr::pull())
dr_0hr_run_ids = as.vector(run_df %>% dplyr::filter(time_point == 0, sample_type == "dr_negative") %>% dplyr::select(Library.Name) %>% dplyr::pull())

total_2hr_run_ids = as.vector(run_df %>% dplyr::filter(time_point == 2, sample_type == "total") %>% dplyr::select(Library.Name) %>% dplyr::pull())
dr_2hr_run_ids = as.vector(run_df %>% dplyr::filter(time_point == 2, sample_type == "dr_negative")%>% dplyr::select(Library.Name) %>% dplyr::pull())

total_24hr_run_ids = as.vector(run_df %>% dplyr::filter(time_point == 24, sample_type == "total") %>% dplyr::select(Library.Name) %>% dplyr::pull())
dr_24hr_run_ids = as.vector(run_df %>% dplyr::filter(time_point == 24, sample_type == "dr_negative")%>% dplyr::select(Library.Name) %>% dplyr::pull())

total_96hr_run_ids = as.vector(run_df %>% dplyr::filter(time_point == 96, sample_type == "total") %>% dplyr::select(Library.Name) %>% dplyr::pull())
dr_96hr_run_ids = as.vector(run_df %>% dplyr::filter(time_point == 96, sample_type == "dr_negative") %>% dplyr::select(Library.Name) %>% dplyr::pull())

gtfFilePath="/home/cidr/Documents/work/test_analysis/homoSapien/forHisat/gencode.v31.annotation.gtf"
counts_dir= "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/"
degs_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/"


###############################################################
#cross sectional counting
###############################################################
dr_negative_0hr_files = as.vector(run_df %>% dplyr::filter(time_point == 0, sample_type == "dr_negative") %>% dplyr::select(file_path) %>% dplyr::pull())
total_0hr_files = as.vector(run_df %>% dplyr::filter(time_point == 0, sample_type == "total") %>% dplyr::select(file_path) %>% dplyr::pull())

countsFilePath = getPath(filename = "counts_total_0hr_vs_dr_negative_0hr.csv", directory = counts_dir)
degsFilePath = getPath(filename = "degs_total_0hr_vs_dr_negative_0hr.csv", directory = degs_dir)

generate_counts_file(X = total_0hr_files, Y = dr_negative_0hr_files, gtfFilePath = gtfFilePath, countsFilePath = countsFilePath)

dr_negative_2hr_files = as.vector(run_df %>% dplyr::filter(time_point == 2, sample_type == "dr_negative") %>% dplyr::select(file_path) %>% dplyr::pull())
total_2hr_files = as.vector(run_df %>% dplyr::filter(time_point == 2, sample_type == "total") %>% dplyr::select(file_path) %>% dplyr::pull())
countsFilePath_2hr = getPath(filename = "counts_total_2hr_vs_dr_negative_2hr.csv", directory = counts_dir)
generate_counts_file(X = total_2hr_files, Y = dr_negative_2hr_files, gtfFilePath = gtfFilePath, countsFilePath = countsFilePath_2hr)


dr_negative_24hr_files = as.vector(run_df %>% dplyr::filter(time_point == 24, sample_type == "dr_negative") %>% dplyr::select(file_path) %>% dplyr::pull())
total_24hr_files = as.vector(run_df %>% dplyr::filter(time_point == 24, sample_type == "total") %>% dplyr::select(file_path) %>% dplyr::pull())
countsFilePath_24hr = getPath(filename = "counts_total_2hr_vs_dr_negative_24hr.csv", directory = counts_dir)
generate_counts_file(X = total_24hr_files, Y = dr_negative_24hr_files, gtfFilePath = gtfFilePath, countsFilePath = countsFilePath_24hr)

dr_negative_96hr_files = as.vector(run_df %>% dplyr::filter(time_point == 96, sample_type == "dr_negative") %>% dplyr::select(file_path) %>% dplyr::pull())
total_96hr_files = as.vector(run_df %>% dplyr::filter(time_point == 96, sample_type == "total") %>% dplyr::select(file_path) %>% dplyr::pull())
countsFilePath_96hr = getPath(filename = "counts_total_2hr_vs_dr_negative_96hr.csv", directory = counts_dir)
generate_counts_file(X = total_96hr_files, Y = dr_negative_96hr_files, gtfFilePath = gtfFilePath, countsFilePath = countsFilePath_96hr)


###############################################################
#total longitudinal counting
###############################################################
countsFilePath_total_0_vs_2hr = getPath(filename = "total_0hr_vs_total_2hr.csv", directory = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/longitudinal/")
generate_counts_file(X = total_0hr_files, Y = total_2hr_files, gtfFilePath = gtfFilePath, countsFilePath = countsFilePath_total_0_vs_2hr)

countsFilePath_total_0_vs_24hr = getPath(filename = "total_0hr_vs_total_24hr.csv", directory = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/longitudinal/")
generate_counts_file(X = total_0hr_files, Y = total_24hr_files, gtfFilePath = gtfFilePath, countsFilePath = countsFilePath_total_0_vs_24hr)

countsFilePath_total_0_vs_96hr = getPath(filename = "total_0hr_vs_total_96hr.csv", directory = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/longitudinal/")
generate_counts_file(X = total_0hr_files, Y = total_96hr_files, gtfFilePath = gtfFilePath, countsFilePath = countsFilePath_total_0_vs_96hr)

###############################################################
#dr_negative longitudinal counting
###############################################################
countsFilePath_dr_negative_0_vs_2hr = getPath(filename = "dr_negative_0hr_vs_dr_negative_2hr.csv", directory = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/longitudinal/")
generate_counts_file(X = dr_negative_0hr_files, Y = dr_negative_2hr_files, gtfFilePath = gtfFilePath, countsFilePath = countsFilePath_dr_negative_0_vs_2hr)

countsFilePath_dr_negative_0_vs_24hr = getPath(filename = "dr_negative_0hr_vs_dr_negative_24hr.csv", directory = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/longitudinal/")
generate_counts_file(X = dr_negative_0hr_files, Y = dr_negative_24hr_files, gtfFilePath = gtfFilePath, countsFilePath = countsFilePath_dr_negative_0_vs_24hr)

countsFilePath_dr_negative_0_vs_96hr = getPath(filename = "dr_negative_0hr_vs_dr_negative_96hr.csv", directory = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/longitudinal/")
generate_counts_file(X = dr_negative_0hr_files, Y = dr_negative_96hr_files, gtfFilePath = gtfFilePath, countsFilePath = countsFilePath_dr_negative_0_vs_96hr)


###############################################################
#generating longitudinal venn between for total 
###############################################################

get_mega_report(degs1_filepath = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/annotated/annotated_total_0_vs_2hr_degs_df.csv", degs2_filepath = '/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/annotated/annotated_total_0_vs_24hr_degs_df.csv', set1_name = "total_0_vs_2", set2_name = "total_0_vs_24", excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports")

get_mega_report(degs1_filepath = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/annotated/annotated_total_0_vs_2hr_degs_df.csv", degs2_filepath = '/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/annotated/annotated_total_0_vs_96hr_degs_df.csv', set1_name = "total_0_vs_2", set2_name = "total_0_vs_96", excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports")

get_mega_report(degs1_filepath = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/annotated/annotated_total_0_vs_24hr_degs_df.csv", degs2_filepath = '/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/annotated/annotated_total_0_vs_96hr_degs_df.csv', set1_name = "total_0_vs_24", set2_name = "total_0_vs_96", excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports")

###############################################################
#generating longitudinal venn between for dr negative 
###############################################################

get_mega_report(degs1_filepath = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/annotated/annotated_dr_negative_0_vs_2hr_degs_df.csv", degs2_filepath = '/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/annotated/annotated_dr_negative_0_vs_24hr_degs_df.csv', set1_name = "dr_0_vs_2", set2_name = "dr_0_vs_24", excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports")

get_mega_report(degs1_filepath = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/annotated/annotated_dr_negative_0_vs_2hr_degs_df.csv", degs2_filepath = '/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/annotated/annotated_dr_negative_0_vs_96hr_degs_df.csv', set1_name = "dr_0_vs_2", set2_name = "dr_0_vs_96", excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports")

get_mega_report(degs1_filepath = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/annotated/annotated_dr_negative_0_vs_24hr_degs_df.csv", degs2_filepath = '/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/annotated/annotated_dr_negative_0_vs_96hr_degs_df.csv', set1_name = "dr_0_vs_24", set2_name = "dr_0_vs_96", excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports")

###############################################################
#generating longitudinal venn between total and dr_negative
###############################################################
get_mega_report(degs1_filepath = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/annotated/annotated_total_0_vs_2hr_degs_df.csv", degs2_filepath = '/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/annotated/annotated_dr_negative_0_vs_2hr_degs_df.csv', set1_name = "total_0_vs_2", set2_name = "dr_neg_0_vs_2", excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports")
get_mega_report(degs1_filepath = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/annotated/annotated_total_0_vs_24hr_degs_df.csv", degs2_filepath = '/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/annotated/annotated_dr_negative_0_vs_24hr_degs_df.csv', set1_name = "total_0_vs_24", set2_name = "dr_neg_0_vs_24", excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports")
get_mega_report(degs1_filepath = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/annotated/annotated_total_0_vs_96hr_degs_df.csv", degs2_filepath = '/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/annotated/annotated_dr_negative_0_vs_96hr_degs_df.csv', set1_name = "total_0_vs_96", set2_name = "dr_neg_0_vs_96", excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports")


only_in_total_0_vs_2 = read_excel("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports/mega_excelReport_venn_total_0_vs_2_dr_neg_0_vs_2.xlsx", sheet = "total_0_vs_2_degs")
only_in_dr_neg_0_vs_2 = read_excel("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports/mega_excelReport_venn_total_0_vs_2_dr_neg_0_vs_2.xlsx", sheet = "only_dr_neg_0_vs_2_degs")

only_in_total_0_vs_24 = read_excel("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports/mega_excelReport_venn_total_0_vs_24_dr_neg_0_vs_24.xlsx", sheet = "total_0_vs_24_degs")
only_in_dr_neg_0_vs_24 = read_excel("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports/mega_excelReport_venn_total_0_vs_24_dr_neg_0_vs_24.xlsx", sheet = "only_dr_neg_0_vs_24_degs")

only_in_total_0_vs_96 = read_excel("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports/mega_excelReport_venn_total_0_vs_96_dr_neg_0_vs_96.xlsx", sheet = "total_0_vs_96_degs")
only_in_dr_neg_0_vs_96 = read_excel("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports/mega_excelReport_venn_total_0_vs_96_dr_neg_0_vs_96.xlsx", sheet = "only_dr_neg_0_vs_96_degs")


only_in_total_0_vs_2 = only_in_total_0_vs_2 %>% dplyr::mutate(time = 2)
only_in_total_0_vs_24 = only_in_total_0_vs_24 %>% dplyr::mutate(time = 24)
only_in_total_0_vs_96 = only_in_total_0_vs_96 %>% dplyr::mutate(time = 96)

only_in_dr_neg_0_vs_2 = only_in_dr_neg_0_vs_2 %>% dplyr::mutate(time = 2)
only_in_dr_neg_0_vs_24 = only_in_dr_neg_0_vs_24 %>% dplyr::mutate(time = 24)
only_in_dr_neg_0_vs_96 = only_in_dr_neg_0_vs_96 %>% dplyr::mutate(time = 96)

only_in_total = rbind.data.frame(only_in_total_0_vs_2, only_in_total_0_vs_24, only_in_total_0_vs_96)
only_in_dr_neg = rbind.data.frame(only_in_dr_neg_0_vs_2, only_in_dr_neg_0_vs_24, only_in_dr_neg_0_vs_96)

head(only_in_total)
head(only_in_dr_neg)

only_in_dr_neg %>% dplyr::filter(!(gene_symbol %in% only_in_total$gene_symbol)) %>% dplyr::arrange(desc(log2FoldChange))
nrow(only_in_dr_neg %>% dplyr::filter(!(gene_symbol %in% only_in_total$gene_symbol)) %>% dplyr::arrange(desc(log2FoldChange)))
write.csv(only_in_dr_neg %>% dplyr::filter(!(gene_symbol %in% only_in_total$gene_symbol)) %>% dplyr::arrange(desc(log2FoldChange)) %>% dplyr::select(gene_symbol, description, log2FoldChange, padj, time), "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports/dr_negative_degs_which_never_came_in_total.csv", row.names = FALSE)

only_in_dr_neg %>% dplyr::filter((gene_symbol %in% only_in_total$gene_symbol))

interested_df = only_in_dr_neg %>% dplyr::filter(!(gene_symbol %in% only_in_total$gene_symbol)) %>% dplyr::arrange(desc(log2FoldChange))
get_kegg_pathways_from_gene_symbols(interested_df$gene_symbol)
getReactomePathways(interested_df$gene_symbol)
write.csv(get_kegg_pathways_from_gene_symbols(interested_df$gene_symbol), "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports/kegg_pathways_from_dr_negative_degs_which_never_came_in_total.csv")
write.csv(getReactomePathways(interested_df$gene_symbol), "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/mega_reports/reactome_pathways_from_dr_negative_degs_which_never_came_in_total.csv")





#*************************************************
#comparison of rna-seq results with 
#*************************************************
library(readxl)
readxl::excel_sheets("/home/cidr/Documents/work/programs/avlab/asma_work/ppat.1007289.s014.xlsx")

genotypic_dr_negative_0_vs_2_df = read_excel("/home/cidr/Documents/work/programs/avlab/asma_work/ppat.1007289.s014.xlsx", sheet ="HLA-DR- Teff (0 h vs 2 h)")
genotypic_dr_negative_0_vs_24_df = read_excel("/home/cidr/Documents/work/programs/avlab/asma_work/ppat.1007289.s014.xlsx", sheet ="HLA-DR- Teff (0 h vs 24 h)")
genotypic_dr_negative_0_vs_96_df = read_excel("/home/cidr/Documents/work/programs/avlab/asma_work/ppat.1007289.s014.xlsx", sheet ="HLA-DR- Teff (0 h vs 96 h)")

genotypic_total_0_vs_2_df = read_excel("/home/cidr/Documents/work/programs/avlab/asma_work/ppat.1007289.s014.xlsx", sheet ="Total Teff (0 h vs 2 h)")
genotypic_total_0_vs_24_df = read_excel("/home/cidr/Documents/work/programs/avlab/asma_work/ppat.1007289.s014.xlsx", sheet ="Total Teff (0 'vs' 24 h)")
genotypic_total_0_vs_96_df = read_excel("/home/cidr/Documents/work/programs/avlab/asma_work/ppat.1007289.s014.xlsx", sheet ="Total Teff (0 'vs' 96 h)")

head(genotypic_dr_negative_0_vs_2_df)
head(genotypic_dr_negative_0_vs_24_df)
head(genotypic_dr_negative_0_vs_96_df)

head(genotypic_total_0_vs_2_df)
head(genotypic_total_0_vs_24_df)
head(genotypic_total_0_vs_96_df)

nrow(genotypic_dr_negative_0_vs_2_df)
nrow(genotypic_dr_negative_0_vs_24_df)
nrow(genotypic_dr_negative_0_vs_96_df)

nrow(genotypic_total_0_vs_2_df)
nrow(genotypic_total_0_vs_24_df)
nrow(genotypic_total_0_vs_96_df)

genotypic_total_0_vs_2_df %>% dplyr::filter(`Q value` < 0.05)
genotypic_total_0_vs_24_df %>% dplyr::filter(`Q value` < 0.05)
genotypic_total_0_vs_96_df %>% dplyr::filter(`Q value` < 0.05)

genotypic_dr_negative_0_vs_2_df %>% dplyr::filter(`Q value` < 0.05)
genotypic_dr_negative_0_vs_24_df %>% dplyr::filter(`Q_value` < 0.05)
genotypic_dr_negative_0_vs_96_df %>% dplyr::filter(`Q value` < 0.05)

total_0_vs_2_df = read_excel("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/total_0_vs_2hr_degs.xlsx", "all_degs" ) %>% dplyr::filter(padj<0.05)
total_0_vs_24_df = read_excel("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/total_0_vs_24hr_degs.xlsx", "all_degs" ) %>% dplyr::filter(padj<0.05)
total_0_vs_96_df = read_excel("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/total_0_vs_96hr_degs.xlsx", "all_degs" ) %>% dplyr::filter(padj<0.05)

nrow(total_0_vs_2_df)
nrow(total_0_vs_24_df)
nrow(total_0_vs_96_df)

nrow(genotypic_total_0_vs_2_df %>% dplyr::filter(`Q value` < 0.05))
nrow(genotypic_dr_negative_0_vs_2_df %>% dplyr::filter(`Q value` < 0.05))
nrow(genotypic_dr_negative_0_vs_24_df %>% dplyr::filter(`Q_value` < 0.05))
nrow(genotypic_dr_negative_0_vs_96_df %>% dplyr::filter(`Q value` < 0.05))


dr_negative_0_vs_2_df = read_excel("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/dr_negative_0_vs_2hr.xlsx", "all_degs" ) %>% dplyr::filter(padj<0.05)
dr_negative_0_vs_24_df = read_excel("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/dr_negative_0_vs_24hr.xlsx", "all_degs" ) %>% dplyr::filter(padj<0.05)
dr_negative_0_vs_96_df = read_excel("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/dr_negative_0_vs_96hr.xlsx", "all_degs" ) %>% dplyr::filter(padj<0.05)

nrow(total_0_vs_2_df %>% dplyr::inner_join(genotypic_total_0_vs_2_df, by = c("gene_symbol"="Gene Name")))
nrow(total_0_vs_24_df %>% dplyr::inner_join(genotypic_total_0_vs_24_df, by = c("gene_symbol"="Gene Name")))

nrow(total_0_vs_2_df %>% dplyr::inner_join(genotypic_total_0_vs_2_df %>% dplyr::filter(`Q value` < 0.05), by = c("gene_symbol"="Gene Name")))
nrow(total_0_vs_24_df %>% dplyr::inner_join(genotypic_total_0_vs_24_df %>% dplyr::filter(`Q value` < 0.05), by = c("gene_symbol"="Gene Name")))
nrow(total_0_vs_96_df %>% dplyr::inner_join(genotypic_total_0_vs_96_df %>% dplyr::filter(`Q value` < 0.05), by = c("gene_symbol"="Gene Name")))

nrow(dr_negative_0_vs_2_df %>% dplyr::inner_join(genotypic_dr_negative_0_vs_2_df %>% dplyr::filter(`Q value` < 0.05), by = c("gene_symbol"="Gene Name")))
nrow(dr_negative_0_vs_24_df %>% dplyr::inner_join(genotypic_dr_negative_0_vs_24_df %>% dplyr::filter(`Q_value` < 0.05), by = c("gene_symbol"="Gene Name")))
nrow(dr_negative_0_vs_96_df %>% dplyr::inner_join(genotypic_dr_negative_0_vs_96_df %>% dplyr::filter(`Q value` < 0.05), by = c("gene_symbol"="Gene Name")))
