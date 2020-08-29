library(ggrepel)

get_volcano_plot = function(degs_df, title= "")
{
  library(ggplot2)
  library(dplyr)
  #library(ggrepel)
  return(degs_df %>% dplyr::filter(!is.na(padj)) %>% ggplot( aes(x = log2FoldChange, -log10(padj), color = factor(padj>0.05))) + geom_point(alpha = 0.3) + geom_hline(yintercept = -log10(0.05)) +ggtitle(title))
}
get_annotated_deg_df = function(deg_file_path,degs_df=NA, filter_by_padj = TRUE)
{
  library(dplyr)
  library(stringr)
  if(nrow(degs_df)==0)
  {
    return(degs_df)
  }
  if(is.na(degs_df))
  {
    if(str_ends(deg_file_path,"csv"))
    {
      degs_df=read.csv(file = deg_file_path,header = T,stringsAsFactors = F)
    }
    if(str_ends(deg_file_path,"xlsx") | str_ends(deg_file_path,"xls"))
    {
      degs_df=read_excel(path =deg_file_path,sheet = 1)
    }
  }
  if(filter_by_padj){
    degs_df=filterByPadjValue(degs_df,cutoff = 0.05)
  }
  print(head(degs_df))
  # columns=colnames(degs_df)
  # columns[1]="ensembl_gene_id"
  # colnames(degs_df)=columns
  # print(head(degs_df))
  #names(degs_df)
  #rownames(degs_df)=degs_df$ensembl_gene_id
  my_ensemble_ids=rownames(degs_df)
  clean_ensemble_ids=getCleanEnsembleIds(my_ensemble_ids)
  degs_df$clean_ensemble_id=clean_ensemble_ids
  degs_df$ensemble_id=my_ensemble_ids
  #length(clean_ensemble_ids)
  #length(unique(clean_ensemble_ids))
  
  #**************************************************
  # fetching annotation from biomart
  #**************************************************
  columns=c("ensembl_gene_id","entrezgene_id","hgnc_symbol","description","gene_biotype")
  annotation_df=getAnnotationFromBiomart(ids = clean_ensemble_ids,filter = "ensembl_gene_id",columns = columns)
  nrow(annotation_df)
  
  
  allowed_biotypes=c("IG_C_gene","IG_D_gene","IG_J_gene","IG_LV_gene","IG_V_gene","TR_C_gene","TR_J_gene","TR_V_gene","TR_D_gene ","protein_coding")
  
  annotation_df=annotation_df %>% dplyr::filter(gene_biotype %in% allowed_biotypes)
  #sum(duplicated(annotation_df$ensembl_gene_id))
  
  
  annotation_summary=as.data.frame(table(annotation_df$ensembl_gene_id))
  #names(annotation_summary)
  #annotation_summary %>% dplyr::filter(Freq>1)
  
  # duplicated_genes=annotation_df[duplicated(annotation_df$ensembl_gene_id),1]
  # annotation_df[annotation_df$ensembl_gene_id %in% duplicated_genes,]
  # degs_df[degs_df$clean_ensemble_id %in% duplicated_genes,]
  
  
  #little experimet
  #degs_df[degs_df$clean_ensemble_id=="ENSG00000143226",] %>% inner_join(annotation_df,c("clean_ensemble_id"="ensembl_gene_id"))
  
  #nrow(degs_df %>% inner_join(annotation_df,c("clean_ensemble_id"="ensembl_gene_id")))
  print(head(degs_df))
  print(head(annotation_df))
  print(class(degs_df))
  print(class(annotation_df))
  degs_df=degs_df %>% dplyr::inner_join(annotation_df,by=c("clean_ensemble_id"="ensembl_gene_id"))
  names(degs_df)
  print(head(degs_df))
  degs_df=degs_df %>% dplyr::select(hgnc_symbol,description,log2FoldChange,padj,pvalue,ensemble_id,clean_ensemble_id,entrezgene_id) %>% dplyr::rename(gene_symbol=hgnc_symbol) %>% arrange(padj)
  
  return(degs_df)
}

get_excel_degs_with_pathways=function(my_degs_df,excel_report_output_dir,excel_report_filename)
{
  library(dplyr)
  library(xlsx)
  #my_degs_df=read.csv(file=degs_csv_path,stringsAsFactors = F,header = T)
  my_degs_df=my_degs_df %>% dplyr::arrange(desc(log2FoldChange))
  my_degs_df=my_degs_df %>% group_by(gene_symbol) %>% dplyr::filter(row_number()==1) %>% ungroup()
  
  #print(head(my_degs_df))
  up_degs_df=my_degs_df %>% dplyr::filter(log2FoldChange >= 1) %>% dplyr::arrange(desc(log2FoldChange))
  down_degs_df=my_degs_df %>% dplyr::filter(log2FoldChange <= -1) %>% dplyr::arrange(log2FoldChange)
  #getting the details dataframe
  up_down_details_df=get_up_down_details_df(my_degs_df)
  
  #preparing the gene list as required for the pathway enrichment analysis
  for_gsea_df=my_degs_df
  for_gsea_df=for_gsea_df %>% dplyr::arrange(log2FoldChange)
  #for_gsea_annotation_df=getAnnotationFromBiomart(ids=for_gsea_df$gene_symbol,filter = "hgnc_symbol",columns = c("hgnc_symbol","entrezgene_id"))
  #for_gsea_df=for_gsea_df %>% inner_join(for_gsea_annotation_df,c("gene_symbol"="hgnc_symbol"))
  for_gsea_gene_list=for_gsea_df$log2FoldChange
  names(for_gsea_gene_list)=for_gsea_df$entrezgene_id
  up_gene_list=for_gsea_gene_list[for_gsea_gene_list>=1]
  down_gene_list=for_gsea_gene_list[for_gsea_gene_list<=-1]
  
  #adding results to workbook
  wb=createWorkbook()
  
  sheet=createSheet(wb,sheetName="details")
  addDataFrame(x=as.data.frame(up_down_details_df),sheet=sheet,row.names = F)
  
  sheet=createSheet(wb,sheetName = "all_degs")
  addDataFrame(x=as.data.frame(my_degs_df),sheet=sheet,row.names = F)
  try(doGSEA(geneList = for_gsea_gene_list,wb=wb,sheetNamePrefix = "pthwys_all_genes_"),silent=T)
  
  sheet=createSheet(wb,sheetName = "up_degs")
  addDataFrame(x=as.data.frame(up_degs_df),sheet=sheet,row.names = F)
  try(doGSEA(geneList = up_gene_list,wb=wb,sheetNamePrefix = "pthwys_up_gene_"),silent=T)
  
  sheet=createSheet(wb,sheetName = "down_degs")
  addDataFrame(x=as.data.frame(down_degs_df),sheet=sheet,row.names = F)
  try(doGSEA(geneList = down_gene_list,wb=wb,sheetNamePrefix = "pthwys_down_gene_"),silent=T)
  
  excel_report_filepath=getPath(filename = excel_report_filename,directory = excel_report_output_dir)
  saveWorkbook(wb,file=excel_report_filepath)
  gc()
}
get_differential_expression_results_df=function(counts_df, filter_by_padj = TRUE)
{
  library(stringr)
  library(DESeq2)
  
  
  
  print(paste("nrow(counts_df[rowSums(counts_df)<10,])",nrow(counts_df[rowSums(counts_df)<10,])))
  print(paste("nrow(counts_df[rowSums(counts_df)<6,])",nrow(counts_df[rowSums(counts_df)<6,])))
  print(paste("nrow(counts_df[rowSums(counts_df)<3,])",nrow(counts_df[rowSums(counts_df)<2,])))
  
  
  #counts_df=counts_df[rowSums(counts_df)>rowsum_cutoff,]
  
  columns=colnames(counts_df)
  
  #preparing the details dataframe
  sam_file_details=as.data.frame(cbind(columns[1:(ncol(counts_df)/2)],columns[((ncol(counts_df)/2)+1): ncol(counts_df)]))
  colnames(sam_file_details)=c("control","treatment")
  
  
  #***********************************************************
  #conducting differential expression analysis
  #***********************************************************
  colData = DataFrame(condition=factor(c("ctrl","ctrl", "ctrl","ctrl", "ctrl","treat", "treat", "treat","treat", "treat")))
  dds = DESeqDataSetFromMatrix(counts_df, colData, formula(~ condition))
  
  # run DEseq
  res= DESeq(dds, minReplicatesForReplace = 3)
  
  #plotMA(res)
  
  #diffrentiallyExpressed 
  deGenes = results(res,alpha=0.05)
  
  # order by BH adjusted p-value
  deGenesOrdered = deGenes[order(deGenes$padj),]
  if(filter_by_padj){
    deGenesOrdered=filterByPadjValue(df=deGenesOrdered,cutoff = 0.05)
  }
  return(as.data.frame(deGenesOrdered))
}


##############################################################################################################################
#differential expression analysis
##############################################################################################################################
run_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/metadata/run_df.csv", stringsAsFactors = FALSE)
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

#************************************
#total longitudinal preparation of counts files 
#************************************
counts_total_0_vs_2hr_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/longitudinal/total_0hr_vs_total_2hr.csv", stringsAsFactors = FALSE)
names(counts_total_0_vs_2hr_df)
counts_total_0_vs_24hr_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/longitudinal/total_0hr_vs_total_24hr.csv", stringsAsFactors = FALSE)
names(counts_total_0_vs_24hr_df)
counts_total_0_vs_96hr_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/longitudinal/total_0hr_vs_total_96hr.csv", stringsAsFactors = FALSE)
names(counts_total_0_vs_96hr_df)


rownames(counts_total_0_vs_2hr_df) = counts_total_0_vs_2hr_df[,1]
counts_total_0_vs_2hr_df = counts_total_0_vs_2hr_df[,-1]
head(counts_total_0_vs_2hr_df)
colnames(counts_total_0_vs_2hr_df) = c(total_0hr_run_ids, total_2hr_run_ids)
counts_total_0_vs_2hr_df = median_based_count_matrix_filtering(counts_total_0_vs_2hr_df)
head(counts_total_0_vs_2hr_df)

rownames(counts_total_0_vs_24hr_df) = counts_total_0_vs_24hr_df[,1]
counts_total_0_vs_24hr_df = counts_total_0_vs_24hr_df[,-1]
head(counts_total_0_vs_24hr_df)
colnames(counts_total_0_vs_24hr_df) = c(total_0hr_run_ids, total_24hr_run_ids)
counts_total_0_vs_24hr_df = median_based_count_matrix_filtering(counts_total_0_vs_24hr_df)
head(counts_total_0_vs_24hr_df)

rownames(counts_total_0_vs_96hr_df) = counts_total_0_vs_96hr_df[,1]
counts_total_0_vs_96hr_df = counts_total_0_vs_96hr_df[,-1]
head(counts_total_0_vs_96hr_df)
colnames(counts_total_0_vs_96hr_df) = c(total_0hr_run_ids, total_96hr_run_ids)
counts_total_0_vs_96hr_df = median_based_count_matrix_filtering(counts_total_0_vs_96hr_df)
head(counts_total_0_vs_96hr_df)

#************************************
#dr_negative longitudinal preparation of counts files 
#************************************
counts_dr_negative_0_vs_2hr_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/longitudinal/dr_negative_0hr_vs_dr_negative_2hr.csv", stringsAsFactors = FALSE)
names(counts_dr_negative_0_vs_2hr_df)
counts_dr_negative_0_vs_24hr_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/longitudinal/dr_negative_0hr_vs_dr_negative_24hr.csv", stringsAsFactors = FALSE)
names(counts_dr_negative_0_vs_24hr_df)
counts_dr_negative_0_vs_96hr_df = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/longitudinal/dr_negative_0hr_vs_dr_negative_96hr.csv", stringsAsFactors = FALSE)
names(counts_dr_negative_0_vs_96hr_df)


rownames(counts_dr_negative_0_vs_2hr_df) = counts_dr_negative_0_vs_2hr_df[,1]
counts_dr_negative_0_vs_2hr_df = counts_dr_negative_0_vs_2hr_df[,-1]
head(counts_dr_negative_0_vs_2hr_df)
colnames(counts_dr_negative_0_vs_2hr_df) = c(dr_0hr_run_ids, dr_2hr_run_ids)
counts_dr_negative_0_vs_2hr_df = median_based_count_matrix_filtering(counts_dr_negative_0_vs_2hr_df)
head(counts_dr_negative_0_vs_2hr_df)

rownames(counts_dr_negative_0_vs_24hr_df) = counts_dr_negative_0_vs_24hr_df[,1]
counts_dr_negative_0_vs_24hr_df = counts_dr_negative_0_vs_24hr_df[,-1]
head(counts_dr_negative_0_vs_24hr_df)
colnames(counts_dr_negative_0_vs_24hr_df) = c(dr_0hr_run_ids, dr_24hr_run_ids)
counts_dr_negative_0_vs_24hr_df = median_based_count_matrix_filtering(counts_dr_negative_0_vs_24hr_df)
head(counts_dr_negative_0_vs_24hr_df)

rownames(counts_dr_negative_0_vs_96hr_df) = counts_dr_negative_0_vs_96hr_df[,1]
counts_dr_negative_0_vs_96hr_df = counts_dr_negative_0_vs_96hr_df[,-1]
head(counts_dr_negative_0_vs_96hr_df)
colnames(counts_dr_negative_0_vs_96hr_df) = c(dr_0hr_run_ids, dr_96hr_run_ids)
counts_dr_negative_0_vs_96hr_df = median_based_count_matrix_filtering(counts_dr_negative_0_vs_96hr_df)
head(counts_dr_negative_0_vs_96hr_df)

#************************************************************************
#cross sectional analysis 
#************************************************************************

counts_df_0hr = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/counts_total_0hr_vs_dr_negative_0hr.csv", stringsAsFactors = FALSE)

counts_df_2hr = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/counts_total_2hr_vs_dr_negative_2hr.csv", stringsAsFactors = FALSE)

counts_df_24hr = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/counts_total_24hr_vs_dr_negative_24hr.csv", stringsAsFactors = FALSE)

counts_df_96hr = read.csv("/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/counts_csv/counts_total_96hr_vs_dr_negative_96hr.csv", stringsAsFactors = FALSE)


#*********************************
#count based filtering of cross sectional data
#*********************************
names(counts_df_0hr)
rownames(counts_df_0hr) = counts_df_0hr[,1]
counts_df_0hr = counts_df_0hr[,-1]
colnames(counts_df_0hr) = c(total_0hr_run_ids, dr_0hr_run_ids)
head(counts_df_0hr)
counts_df_0hr = median_based_count_matrix_filtering(counts_df_0hr)
head(counts_df_0hr)


names(counts_df_2hr)
rownames(counts_df_2hr) = counts_df_2hr[,1]
counts_df_2hr = counts_df_2hr[,-1]
colnames(counts_df_2hr) = c(total_2hr_run_ids, dr_2hr_run_ids)
head(counts_df_2hr)
counts_df_2hr = median_based_count_matrix_filtering(counts_df_2hr)
head(counts_df_2hr)

names(counts_df_24hr)
rownames(counts_df_24hr) = counts_df_24hr[,1]
counts_df_24hr = counts_df_24hr[,-1]
colnames(counts_df_24hr) = c(total_24hr_run_ids, dr_24hr_run_ids)
head(counts_df_24hr)
counts_df_24hr = median_based_count_matrix_filtering(counts_df_24hr)
head(counts_df_24hr)

head(counts_df_96hr)
rownames(counts_df_96hr) = counts_df_96hr[,1]
counts_df_96hr = counts_df_96hr[,-1]
colnames(counts_df_96hr) = c(total_96hr_run_ids, dr_96hr_run_ids)
head(counts_df_96hr)
counts_df_96hr = median_based_count_matrix_filtering(counts_df_96hr)
head(counts_df_96hr)

############################################################
#PCA Analysis 
############################################################
#******************************************************
#pca total 0 vs 2 
#******************************************************
colData = DataFrame(condition=factor(c("Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative")))
dds = DESeqDataSetFromMatrix(counts_total_0_vs_2hr_df, colData, formula(~ condition))
vsd = vst(dds, blind = FALSE)
head(assay(vsd))

print(plotPCA(vsd, returnData = F))
rv = plotPCA(vsd, returnData = TRUE)

ggplot(data = rv, aes(x= PC1, y = PC2, color = group, label = name)) + geom_point() + expand_limits(x = c(-30,70)) + geom_label_repel() + theme_classic()

#******************************************************
#pca total 0 vs 24
#******************************************************
colData = DataFrame(condition=factor(c("Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative")))
dds = DESeqDataSetFromMatrix(counts_total_0_vs_2hr_df, colData, formula(~ condition))
vsd = vst(dds, blind = FALSE)
head(assay(vsd))

print(plotPCA(vsd, returnData = F))
rv = plotPCA(vsd, returnData = TRUE)

ggplot(data = rv, aes(x= PC1, y = PC2, color = group, label = name)) + geom_point() + expand_limits(x = c(-30,70)) + geom_label_repel() + theme_classic()

colData = DataFrame(condition=factor(c("Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative")))
dds = DESeqDataSetFromMatrix(counts_dr_negative_0_vs_96hr_df, colData, formula(~ condition))
vsd = vst(dds, blind = FALSE)
head(assay(vsd))

print(plotPCA(vsd, returnData = F))
rv = plotPCA(vsd, returnData = TRUE)

ggplot(data = rv, aes(x= PC1, y = PC2, color = group, label = name)) + geom_point() + expand_limits(x = c(-30,70)) + geom_label_repel() + theme_classic()

#******************************************************
#pca total 0 vs 96
#******************************************************
colData = DataFrame(condition=factor(c("Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative")))
dds = DESeqDataSetFromMatrix(counts_total_0_vs_2hr_df, colData, formula(~ condition))
vsd = vst(dds, blind = FALSE)
head(assay(vsd))

print(plotPCA(vsd, returnData = F))
rv = plotPCA(vsd, returnData = TRUE)

ggplot(data = rv, aes(x= PC1, y = PC2, color = group, label = name)) + geom_point() + expand_limits(x = c(-30,70)) + geom_label_repel() + theme_classic()

colData = DataFrame(condition=factor(c("Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative")))
dds = DESeqDataSetFromMatrix(counts_dr_negative_0_vs_96hr_df, colData, formula(~ condition))
vsd = vst(dds, blind = FALSE)
head(assay(vsd))

print(plotPCA(vsd, returnData = F))
rv = plotPCA(vsd, returnData = TRUE)

ggplot(data = rv, aes(x= PC1, y = PC2, color = group, label = name)) + geom_point() + expand_limits(x = c(-30,70)) + geom_label_repel() + theme_classic()


#******************************************************
#pca dr negative 0 vs 2
#******************************************************
colData = DataFrame(condition=factor(c("Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative")))
dds = DESeqDataSetFromMatrix(counts_dr_negative_0_vs_2hr_df, colData, formula(~ condition))
vsd = vst(dds, blind = FALSE)
head(assay(vsd))

print(plotPCA(vsd, returnData = F))
rv = plotPCA(vsd, returnData = TRUE)

ggplot(data = rv, aes(x= PC1, y = PC2, color = group, label = name)) + geom_point() + expand_limits(x = c(-30,70)) + geom_label_repel() + theme_classic()

colData = DataFrame(condition=factor(c("Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative")))
dds = DESeqDataSetFromMatrix(counts_dr_negative_0_vs_96hr_df, colData, formula(~ condition))
vsd = vst(dds, blind = FALSE)
head(assay(vsd))

print(plotPCA(vsd, returnData = F))
rv = plotPCA(vsd, returnData = TRUE)

ggplot(data = rv, aes(x= PC1, y = PC2, color = group, label = name)) + geom_point() + expand_limits(x = c(-30,70)) + geom_label_repel() + theme_classic()

#******************************************************
#pca dr negative 0 vs 24
#******************************************************
colData = DataFrame(condition=factor(c("Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative")))
dds = DESeqDataSetFromMatrix(counts_dr_negative_0_vs_24hr_df, colData, formula(~ condition))
vsd = vst(dds, blind = FALSE)
head(assay(vsd))

print(plotPCA(vsd, returnData = F))
rv = plotPCA(vsd, returnData = TRUE)

ggplot(data = rv, aes(x= PC1, y = PC2, color = group, label = name)) + geom_point() + expand_limits(x = c(-30,70)) + geom_label_repel() + theme_classic()

colData = DataFrame(condition=factor(c("Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative")))
dds = DESeqDataSetFromMatrix(counts_dr_negative_0_vs_96hr_df, colData, formula(~ condition))
vsd = vst(dds, blind = FALSE)
head(assay(vsd))

print(plotPCA(vsd, returnData = F))
rv = plotPCA(vsd, returnData = TRUE)

ggplot(data = rv, aes(x= PC1, y = PC2, color = group, label = name)) + geom_point() + expand_limits(x = c(-30,70)) + geom_label_repel() + theme_classic()

#******************************************************
#pca dr negative 0 vs 96
#******************************************************
colData = DataFrame(condition=factor(c("Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative")))
dds = DESeqDataSetFromMatrix(counts_dr_negative_0_vs_96hr_df, colData, formula(~ condition))
vsd = vst(dds, blind = FALSE)
head(assay(vsd))

print(plotPCA(vsd, returnData = F))
rv = plotPCA(vsd, returnData = TRUE)

ggplot(data = rv, aes(x= PC1, y = PC2, color = group, label = name)) + geom_point() + expand_limits(x = c(-30,70)) + geom_label_repel() + theme_classic()

colData = DataFrame(condition=factor(c("Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative")))
dds = DESeqDataSetFromMatrix(counts_dr_negative_0_vs_96hr_df, colData, formula(~ condition))
vsd = vst(dds, blind = FALSE)
head(assay(vsd))

print(plotPCA(vsd, returnData = F))
rv = plotPCA(vsd, returnData = TRUE)

ggplot(data = rv, aes(x= PC1, y = PC2, color = group, label = name)) + geom_point() + expand_limits(x = c(-30,70)) + geom_label_repel() + theme_classic()


#******************************************************
#pca combined samples
#******************************************************
combined_asma_df = counts_df_0hr %>% dplyr::mutate(X = rownames(counts_df_0hr)) %>% dplyr::inner_join(counts_df_2hr %>% dplyr::mutate(X = rownames(counts_df_2hr)), by = c("X"="X") ) %>% dplyr::inner_join(counts_df_24hr %>% dplyr::mutate(X = rownames(counts_df_24hr)), by = c("X"="X") ) %>% dplyr::inner_join(counts_df_96hr %>% dplyr::mutate(X = rownames(counts_df_96hr)), by = c("X"="X") ) 

rownames(combined_asma_df) = combined_asma_df[,c("X")]
combined_asma_df$X= NULL

colData = DataFrame(condition=factor(c("Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative","Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative","Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative","Total","Total", "Total","Total", "Total","DR_negative", "DR_negative", "DR_negative", "DR_negative", "DR_negative")))
dds = DESeqDataSetFromMatrix(combined_asma_df, colData, formula(~ condition))
vsd = vst(dds, blind = FALSE)
head(assay(vsd))

print(plotPCA(vsd, returnData = F))
rv = plotPCA(vsd, returnData = TRUE)

ggplot(data = rv, aes(x= PC1, y = PC2, color = group, label = name)) + geom_point() + expand_limits(x = c(-30,70)) + geom_label_repel() + theme_classic()


################
#deseq2 cross sectional analysis
################
degs_df_0hr = get_differential_expression_results_df(counts_df_0hr, filter_by_padj = FALSE)
degs_df_2hr = get_differential_expression_results_df(counts_df_2hr, filter_by_padj = FALSE)
degs_df_24hr = get_differential_expression_results_df(counts_df_24hr, filter_by_padj = FALSE)
degs_df_96hr = get_differential_expression_results_df(counts_df_96hr, filter_by_padj = FALSE)

annotated_degs_df_0hr = get_annotated_deg_df(degs_df = degs_df_0hr,filter_by_padj = FALSE )
annotated_degs_df_2hr = get_annotated_deg_df(degs_df = degs_df_2hr,filter_by_padj = FALSE )
annotated_degs_df_24hr = get_annotated_deg_df(degs_df = degs_df_24hr,filter_by_padj = FALSE )
annotated_degs_df_96hr = get_annotated_deg_df(degs_df = degs_df_96hr,filter_by_padj = FALSE )

get_excel_degs_with_pathways(annotate_degs_df_0hr, excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/", excel_report_filename = "cross_0hr.xlsx")
get_excel_degs_with_pathways(annotate_degs_df_2hr, excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/", excel_report_filename = "cross_2hr.xlsx")
get_excel_degs_with_pathways(annotate_degs_df_24hr, excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/", excel_report_filename = "cross_24hr.xlsx")
get_excel_degs_with_pathways(annotate_degs_df_96hr, excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/", excel_report_filename = "cross_96hr.xlsx")
################
#deseq2 longitudinal for Total
################
total_0_vs_2hr_degs_df = get_differential_expression_results_df(counts_total_0_vs_2hr_df, filter_by_padj = FALSE)
total_0_vs_24hr_degs_df = get_differential_expression_results_df(counts_total_0_vs_24hr_df, filter_by_padj = FALSE )
total_0_vs_96hr_degs_df = get_differential_expression_results_df(counts_total_0_vs_96hr_df, filter_by_padj = FALSE)

write.csv(total_0_vs_2hr_degs_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/unannotated/total_0_vs_2hr_degs_df.csv")
write.csv(total_0_vs_24hr_degs_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/unannotated/total_0_vs_24hr_degs_df.csv")
write.csv(total_0_vs_96hr_degs_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/unannotated/total_0_vs_96hr_degs_df.csv")

dr_negative_0_vs_2hr_degs_df = get_differential_expression_results_df(counts_dr_negative_0_vs_2hr_df, filter_by_padj = FALSE)
dr_negative_0_vs_24hr_degs_df = get_differential_expression_results_df(counts_dr_negative_0_vs_24hr_df, filter_by_padj = FALSE )
dr_negative_0_vs_96hr_degs_df = get_differential_expression_results_df(counts_dr_negative_0_vs_96hr_df, filter_by_padj = FALSE)

write.csv(dr_negative_0_vs_2hr_degs_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/unannotated/dr_negative_0_vs_2hr_degs_df.csv")
write.csv(dr_negative_0_vs_24hr_degs_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/unannotated/dr_negative_0_vs_24hr_degs_df.csv")
write.csv(dr_negative_0_vs_96hr_degs_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/unannotated/dr_negative_0_vs_96hr_degs_df.csv")



annotated_total_0_vs_2hr_degs_df = get_annotated_deg_df(degs_df = total_0_vs_2hr_degs_df)
annotated_total_0_vs_24hr_degs_df = get_annotated_deg_df(degs_df = total_0_vs_24hr_degs_df)
annotated_total_0_vs_96hr_degs_df = get_annotated_deg_df(degs_df = total_0_vs_96hr_degs_df)


write.csv(annotated_total_0_vs_2hr_degs_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/annotated/annotated_total_0_vs_2hr_degs_df.csv", row.names = FALSE)
write.csv(annotated_total_0_vs_24hr_degs_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/annotated/annotated_total_0_vs_24hr_degs_df.csv", row.names = FALSE)
write.csv(annotated_total_0_vs_96hr_degs_df, "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/total_longitudinal/annotated/annotated_total_0_vs_96hr_degs_df.csv", row.names = FALSE)

annotated_dr_negative_0_vs_2hr_degs_df = get_annotated_deg_df(degs_df = dr_negative_0_vs_2hr_degs_df)
annotated_dr_negative_0_vs_24hr_degs_df = get_annotated_deg_df(degs_df = dr_negative_0_vs_24hr_degs_df)
annotated_dr_negative_0_vs_96hr_degs_df = get_annotated_deg_df(degs_df = dr_negative_0_vs_96hr_degs_df)

write.csv(annotated_dr_negative_0_vs_2hr_degs_df,"/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/annotated/annotated_dr_negative_0_vs_2hr_degs_df.csv", row.names = FALSE)
write.csv(annotated_dr_negative_0_vs_24hr_degs_df,"/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/annotated/annotated_dr_negative_0_vs_24hr_degs_df.csv", row.names = FALSE)
write.csv(annotated_dr_negative_0_vs_96hr_degs_df,"/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/degs_csv/dr_negative_longitudinal/annotated/annotated_dr_negative_0_vs_96hr_degs_df.csv", row.names = FALSE)

get_excel_degs_with_pathways(my_degs_df = annotated_total_0_vs_2hr_degs_df, excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/", excel_report_filename = "total_0_vs_2hr_degs.xlsx")
get_excel_degs_with_pathways(my_degs_df = annotated_total_0_vs_24hr_degs_df, excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/", excel_report_filename = "total_0_vs_24hr_degs.xlsx")
get_excel_degs_with_pathways(my_degs_df = annotated_total_0_vs_96hr_degs_df, excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/", excel_report_filename = "total_0_vs_96hr_degs.xlsx")

get_excel_degs_with_pathways(my_degs_df = annotated_dr_negative_0_vs_2hr_degs_df, excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/", excel_report_filename = "dr_negative_0_vs_2hr.xlsx")
get_excel_degs_with_pathways(my_degs_df = annotated_dr_negative_0_vs_24hr_degs_df, excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/", excel_report_filename = "dr_negative_0_vs_24hr.xlsx")
get_excel_degs_with_pathways(my_degs_df = annotated_dr_negative_0_vs_96hr_degs_df, excel_report_output_dir = "/media/cidr/d7416dce-cdf6-43ed-9df7-978f8a9438d8/databackup_2019_07_11/work/asma_whole_blood_rnaseq/asma_durbar_analysis_raw_reads_26_jun_2020/excel_reports/", excel_report_filename = "dr_negative_0_vs_96hr.xlsx")

