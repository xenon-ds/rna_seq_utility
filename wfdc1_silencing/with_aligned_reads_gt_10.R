library(stringr)


#scr 0 vs 2
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/longitudinal_igra_negative/counts_igra_negative_scr_0hr_vs_igra_negative_scr_2hr.csv",header = TRUE,stringsAsFactors = FALSE)


#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

columns=colnames(counts_df)
columns[str_detect(columns,pattern = "SCO149")]="SCO149"
columns[str_detect(columns,pattern = "SCO68")]="SCO68"
columns[str_detect(columns,pattern = "SCO71")]="SCO71"
columns[1:3]=paste("scr_0hr_",columns[1:3],sep="")
columns[4:6]=paste("scr_2hr_",columns[4:6],sep="")
colnames(counts_df)=columns

head(counts_df)
counts_df = median_based_count_matrix_filtering(counts_df, expected_median = 10)
degs_df=get_differential_expression_results_df(counts_df, filter_by_padj = FALSE)
get_volcano_plot(degs_df, title = "IGRA -ve SCR 0 vs 2 hr")
degs_df=get_annotated_deg_df(degs_df = degs_df, filter_by_padj = FALSE)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_2hr.csv",row.names = FALSE)


#scr 0 vs 24
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/longitudinal_igra_negative/counts_igra_negative_scr_0hr_vs_igra_negative_scr_24hr.csv",header = TRUE,stringsAsFactors = FALSE)


#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

columns=colnames(counts_df)
columns[str_detect(columns,pattern = "SCO149")]="SCO149"
columns[str_detect(columns,pattern = "SCO68")]="SCO68"
columns[str_detect(columns,pattern = "SCO71")]="SCO71"
columns[1:3]=paste("scr_0hr_",columns[1:3],sep="")
columns[4:6]=paste("scr_24hr_",columns[4:6],sep="")
colnames(counts_df)=columns

head(counts_df)
counts_df = median_based_count_matrix_filtering(counts_df, expected_median = 10)
degs_df=get_differential_expression_results_df(counts_df, filter_by_padj = FALSE)
get_volcano_plot(degs_df, title = "IGRA -ve SCR 0 vs 24")
degs_df=get_annotated_deg_df(degs_df = degs_df, filter_by_padj = FALSE)

write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_24hr.csv",row.names = FALSE)


########################################
#modifying longitudinal counts for b3 
########################################
#*******
# b3 0 vs 2
#*******
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/longitudinal_igra_negative/counts_igra_negative_b3_0hr_vs_igra_negative_b3_2hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

columns=colnames(counts_df)
columns[str_detect(columns,pattern = "SCO149")]="SCO149"
columns[str_detect(columns,pattern = "SCO68")]="SCO68"
columns[str_detect(columns,pattern = "SCO71")]="SCO71"
columns[1:3]=paste("b3_0hr_",columns[1:3],sep="")
columns[4:6]=paste("b3_2hr_",columns[4:6],sep="")
colnames(counts_df)=columns

counts_df$b3_0hr_SCO71=round((counts_df$b3_0hr_SCO149 + counts_df$b3_0hr_SCO68)/2)
counts_df$b3_2hr_SCO71=round((counts_df$b3_2hr_SCO149 + counts_df$b3_2hr_SCO68)/2)
counts_df = median_based_count_matrix_filtering(counts_df, expected_median = 10)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df, filter_by_padj = FALSE)
get_volcano_plot(degs_df, title = "IGRA -ve B3 0 vs 2 hr")
degs_df=get_annotated_deg_df(degs_df = degs_df, filter_by_padj = FALSE)

write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_2hr.csv",row.names = FALSE)

#*******
# b3 0 vs 24
#*******
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/longitudinal_igra_negative/counts_igra_negative_b3_0hr_vs_igra_negative_b3_24hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

columns=colnames(counts_df)
columns[str_detect(columns,pattern = "SCO149")]="SCO149"
columns[str_detect(columns,pattern = "SCO68")]="SCO68"
columns[str_detect(columns,pattern = "SCO71")]="SCO71"
columns[1:3]=paste("b3_0hr_",columns[1:3],sep="")
columns[4:6]=paste("b3_24hr_",columns[4:6],sep="")
colnames(counts_df)=columns


counts_df$b3_0hr_SCO71=round((counts_df$b3_0hr_SCO149 + counts_df$b3_0hr_SCO68)/2)
counts_df$b3_24hr_SCO71=round((counts_df$b3_24hr_SCO149 + counts_df$b3_24hr_SCO68)/2)
counts_df = median_based_count_matrix_filtering(counts_df, expected_median = 10)



#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df, filter_by_padj = FALSE)
get_volcano_plot(degs_df, title = "IGRA -ve B3 0 vs 24 hr")
degs_df=get_annotated_deg_df(degs_df = degs_df, filter_by_padj = FALSE)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_24hr.csv",row.names = FALSE)

#write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/with_outlier_diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_24hr.csv",row.names = FALSE)

#scr 0 hr vs b3 0 hr
########################################
#modifying cross sectional counts for b3 
########################################

#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/cross_sectional_igra_negative/counts_igra_negative_scr_0hr_vs_igra_negative_b3_0hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

columns=colnames(counts_df)
columns[str_detect(columns,pattern = "SCO149")]="SCO149"
columns[str_detect(columns,pattern = "SCO68")]="SCO68"
columns[str_detect(columns,pattern = "SCO71")]="SCO71"
columns[1:3]=paste("scr_0hr_",columns[1:3],sep="")
columns[4:6]=paste("b3_0hr_",columns[4:6],sep="")
colnames(counts_df)=columns


counts_df$b3_0hr_SCO71=round((counts_df$b3_0hr_SCO149 + counts_df$b3_0hr_SCO68)/2)
counts_df = median_based_count_matrix_filtering(counts_df, expected_median = 10)
write.csv(counts_df, "/home/cidr/Documents/work/programs/avlab/modified_counts/igra_negative_cross_0hr.csv")

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df, filter_by_padj = FALSE)
get_volcano_plot(degs_df, title = "IGRA -ve SCR 0 vs B3 0 hr")
degs_df=get_annotated_deg_df(degs_df = degs_df, filter_by_padj = FALSE)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_b3_0hr.csv",row.names = FALSE)

#scr 2 hr vs b3 2 hr
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/cross_sectional_igra_negative/counts_igra_negative_scr_2hr_vs_igra_negative_b3_2hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

columns=colnames(counts_df)
columns[str_detect(columns,pattern = "SCO149")]="SCO149"
columns[str_detect(columns,pattern = "SCO68")]="SCO68"
columns[str_detect(columns,pattern = "SCO71")]="SCO71"
columns[1:3]=paste("scr_2hr_",columns[1:3],sep="")
columns[4:6]=paste("b3_2hr_",columns[4:6],sep="")
colnames(counts_df)=columns

counts_df$b3_2hr_SCO71=round((counts_df$b3_2hr_SCO149 + counts_df$b3_2hr_SCO68)/2)
counts_df = median_based_count_matrix_filtering(counts_df, expected_median = 10)
write.csv(counts_df, "/home/cidr/Documents/work/programs/avlab/modified_counts/igra_negative_cross_2hr.csv")


#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df, filter_by_padj = FALSE)
get_volcano_plot(degs_df, title = "IGRA -ve SCR 2 vs B3 2 hr")
degs_df=get_annotated_deg_df(degs_df = degs_df, filter_by_padj = FALSE)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_2hr_vs_igra_negative_b3_2hr.csv",row.names = FALSE)

#scr 24 vs b3 24
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/cross_sectional_igra_negative/counts_igra_negative_scr_24hr_vs_igra_negative_b3_24hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

columns=colnames(counts_df)
columns[str_detect(columns,pattern = "SCO149")]="SCO149"
columns[str_detect(columns,pattern = "SCO68")]="SCO68"
columns[str_detect(columns,pattern = "SCO71")]="SCO71"
columns[1:3]=paste("scr_24hr_",columns[1:3],sep="")
columns[4:6]=paste("b3_24hr_",columns[4:6],sep="")
colnames(counts_df)=columns

counts_df$b3_24hr_SCO71=round((counts_df$b3_24hr_SCO149 + counts_df$b3_24hr_SCO68)/2)
counts_df = median_based_count_matrix_filtering(counts_df, expected_median = 10)
write.csv(counts_df, "/home/cidr/Documents/work/programs/avlab/modified_counts/igra_negative_cross_24hr.csv")


#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df, filter_by_padj = FALSE)
get_volcano_plot(degs_df, title = "IGRA -ve SCR 24 vs B3 24 hr")
degs_df=get_annotated_deg_df(degs_df = degs_df, filter_by_padj = FALSE)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_24hr_vs_igra_negative_b3_24hr.csv",row.names = FALSE)


####################################
#generate labelled volcano plots
####################################
library(ggrepel)
annotated_cross_sectional_0hr_df = read.csv("/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_b3_0hr.csv", stringsAsFactors = FALSE)
head(annotated_cross_sectional_0hr_df)

annotated_cross_sectional_2hr_df = read.csv("/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_2hr_vs_igra_negative_b3_2hr.csv", stringsAsFactors = FALSE)
head(annotated_cross_sectional_2hr_df)

annotated_cross_sectional_24hr_df = read.csv("/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_24hr_vs_igra_negative_b3_24hr.csv", stringsAsFactors = FALSE)
head(annotated_cross_sectional_24hr_df)

get_volcano_plot(annotated_cross_sectional_0hr_df, title = "0 hr") + theme_bw() + geom_label_repel(data = annotated_cross_sectional_0hr_df %>% dplyr::filter(-log10(padj)>=2),aes( x = log2FoldChange,y = -log10(padj), label = gene_symbol) , label.size = 0, size = 3 )

get_volcano_plot(annotated_cross_sectional_2hr_df, title = "2 hr") + theme_bw() + geom_label_repel(data = annotated_cross_sectional_2hr_df %>% dplyr::filter(-log10(padj)>=2),aes( x = log2FoldChange,y = -log10(padj), label = gene_symbol) , label.size = 0, size = 3 )

get_volcano_plot(annotated_cross_sectional_24hr_df, title = "24 hr") + theme_bw() + geom_label_repel(data = annotated_cross_sectional_24hr_df %>% dplyr::filter(padj<0.05),aes( x = log2FoldChange,y = -log10(padj), label = gene_symbol) , label.size = 0, size = 3 )

####################################
#generate mega reports
####################################
get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_2hr.csv", degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_24hr.csv",set1_name = "scr_0_vs_2", set2_name = "scr_0_vs_24", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/mega_reports/", pathwaysType ="gsea")

get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_2hr.csv", degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_24hr.csv",set1_name = "b3_0_vs_2", set2_name = "b3_0_vs_24", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/mega_reports/", pathwaysType ="gsea")

get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_2hr.csv", degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_2hr.csv",set1_name = "scr_0_v_2", set2_name = "b3_0_v_2", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/mega_reports/", pathwaysType ="gsea")

get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_24hr.csv", degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_24hr.csv",set1_name = "scr_0_v_24", set2_name = "b3_0_v_24", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/mega_reports/", pathwaysType ="gsea")

#sanity check
get_mega_report(degs2_filepath = "/home/cidr/Documents/work/fresh/standard_degs/standard_diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_2hr.csv", degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_2hr.csv",set1_name = "scr_0_vs_2", set2_name = "prv_scr_0_v_2", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/mega_reports/", pathwaysType ="gsea")

get_mega_report(degs2_filepath = "/home/cidr/Documents/work/fresh/standard_degs/standard_diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_24hr.csv", degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_24hr.csv",set1_name = "scr_0_vs_24", set2_name = "prv_scr_0_v_24", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/mega_reports/", pathwaysType ="gsea")

get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_2hr.csv",degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/degs_csv_files/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_2hr.csv",set1_name = "b3_0_v_2", set2_name = "prv_b3_0_v_2", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/mega_reports/", pathwaysType ="gsea")

get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_24hr.csv",degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/degs_csv_files/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_24hr.csv",set1_name = "b3_0_v_24", set2_name = "prv_b3_0_v_24", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/mega_reports/", pathwaysType ="gsea")

#cross sectional_sanity_
get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_b3_0hr.csv",degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/degs_csv_files/diffEx_igra_negative_scr_0hr_vs_igra_negative_b3_0hr.csv",set1_name = "crs_0", set2_name = "prv_crs_0", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/mega_reports/", pathwaysType ="gsea")

get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_2hr_vs_igra_negative_b3_2hr.csv",degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/degs_csv_files/diffEx_igra_negative_scr_2hr_vs_igra_negative_b3_2hr.csv",set1_name = "crs_2", set2_name = "prv_crs_2", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/mega_reports/", pathwaysType ="gsea")

get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_24hr_vs_igra_negative_b3_24hr.csv",degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/degs_csv_files/diffEx_igra_negative_scr_24hr_vs_igra_negative_b3_24hr.csv",set1_name = "crs_24", set2_name = "prv_crs_24", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/mega_reports/", pathwaysType ="gsea")



#get excel report from the differential expression
#excel reports of scr longitudinal files
degs_filepath="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_2hr.csv"
get_degs_with_pathways(degs_csv_path = degs_filepath,excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/excel_reports/longitudinal/scr",excel_report_filename = "degs_and_pathways_igra_negative_scr_0hr_vs_igra_negative_scr_2hr.xlsx")

degs_filepath="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_24hr.csv"
get_degs_with_pathways(degs_csv_path = degs_filepath,excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/excel_reports/longitudinal/scr",excel_report_filename = "degs_and_pathways_igra_negative_scr_0hr_vs_igra_negative_scr_24hr.xlsx")

#excel reports of b3 longitudinal files 
degs_filepath="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_2hr.csv"
get_degs_with_pathways(degs_csv_path = degs_filepath,excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/excel_reports/longitudinal/b3",excel_report_filename = "degs_and_pathways_igra_negative_b3_0hr_vs_igra_negative_b3_2hr.xlsx")

degs_filepath="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_24hr.csv"
get_degs_with_pathways(degs_csv_path = degs_filepath,excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/excel_reports/longitudinal/b3",excel_report_filename = "degs_and_pathways_igra_negative_b3_0hr_vs_igra_negative_b3_24hr.xlsx")

#excel reports of cross sectional data
degs_filepath="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_b3_0hr.csv"
get_degs_with_pathways(degs_csv_path = degs_filepath,excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/excel_reports/cross_sectional",excel_report_filename = "degs_and_pathways_igra_negative_scr_0hr_vs_igra_negative_b3_0hr.xlsx")

degs_filepath="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_2hr_vs_igra_negative_b3_2hr.csv"
get_degs_with_pathways(degs_csv_path = degs_filepath,excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/excel_reports/cross_sectional",excel_report_filename = "degs_and_pathways_igra_negative_scr_2hr_vs_igra_negative_b3_2hr.xlsx")

degs_filepath="/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_24hr_vs_igra_negative_b3_24hr.csv"
get_degs_with_pathways(degs_csv_path = degs_filepath,excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/excel_reports/cross_sectional",excel_report_filename = "degs_and_pathways_igra_negative_scr_24hr_vs_igra_negative_b3_24hr.xlsx")


#get_degs_with_pathways("/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/with_outlier_diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_24hr.csv", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/excel_reports/longitudinal/b3", excel_report_filename = "with_outlier_degs_and_pathways_igra_negative_b3_0hr_vs_igra_negative_b3_24hr.xlsx")
##############################################################################
#tb reports
##############################################################################
#*******
# tb scr 0 vs 2
#*******
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_tb_scr_0hr_vs_tb_scr_2hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

counts_df = median_based_count_matrix_filtering(counts_df)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df)
degs_df=get_annotated_deg_df(degs_df = degs_df)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_0hr_vs_tb_scr_2hr.csv",row.names = FALSE)
get_degs_with_pathways(degs_csv_path = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_0hr_vs_tb_scr_2hr.csv",excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/excel_reports",excel_report_filename = "degs_and_pathways_tb_scr_0hr_vs_tb_scr_2hr.csv.xlsx")


#tb scr 0 vs 24
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_tb_scr_0hr_vs_tb_scr_24hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

counts_df = median_based_count_matrix_filtering(counts_df)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df)
degs_df=get_annotated_deg_df(degs_df = degs_df)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_0hr_vs_tb_scr_24hr.csv",row.names = FALSE)
get_degs_with_pathways(degs_csv_path = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_0hr_vs_tb_scr_24hr.csv",excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/excel_reports",excel_report_filename = "degs_and_pathways_tb_scr_0hr_vs_tb_scr_24hr.csv.xlsx")

#*******
# tb b3 0 vs 2
#*******
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_tb_b3_0hr_vs_tb_b3_2hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

counts_df = median_based_count_matrix_filtering(counts_df)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df)
degs_df=get_annotated_deg_df(degs_df = degs_df)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_b3_0hr_vs_tb_b3_2hr.csv",row.names = FALSE)
get_degs_with_pathways(degs_csv_path = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_b3_0hr_vs_tb_b3_2hr.csv",excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/excel_reports",excel_report_filename = "degs_and_pathways_tb_b3_0hr_vs_tb_b3_2hr.csv.xlsx")

#*******
# tb b3 0 vs 24
#*******
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_tb_b3_0hr_vs_tb_b3_24hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

counts_df = median_based_count_matrix_filtering(counts_df)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df)
degs_df=get_annotated_deg_df(degs_df = degs_df)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_b3_0hr_vs_tb_b3_24hr.csv",row.names = FALSE)
get_degs_with_pathways(degs_csv_path = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_b3_0hr_vs_tb_b3_24hr.csv",excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/excel_reports",excel_report_filename = "degs_and_pathways_tb_b3_0hr_vs_tb_b3_24hr.csv.xlsx")
#cross section analysis

#*******
# tb scr 0 vs b3 0
#*******
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_tb_scr_0hr_vs_tb_b3_0hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

counts_df = median_based_count_matrix_filtering(counts_df)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df)
degs_df=get_annotated_deg_df(degs_df = degs_df)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_0hr_vs_tb_b3_0hr.csv",row.names = FALSE)
get_degs_with_pathways(degs_csv_path = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_0hr_vs_tb_b3_0hr.csv",excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/excel_reports",excel_report_filename = "degs_and_pathways_tb_scr_0hr_vs_tb_b3_0hr.xlsx")
#*******
# tb scr 2 vs b3 2
#*******
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_tb_scr_2hr_vs_tb_b3_2hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

counts_df = median_based_count_matrix_filtering(counts_df)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df)
degs_df=get_annotated_deg_df(degs_df = degs_df)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_2hr_vs_tb_b3_2hr.csv",row.names = FALSE)

#*******
# tb scr 24 vs b3 24
#*******
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_tb_scr_24hr_vs_tb_b3_24hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

counts_df = median_based_count_matrix_filtering(counts_df)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df)
degs_df=get_annotated_deg_df(degs_df = degs_df)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_24hr_vs_tb_b3_24hr.csv",row.names = FALSE)
 

#cross sectional analysis with igra negative 
#*******
# igra -ve scr 0 vs tb 0
#*******
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_igra_negative_scr_0hr_vs_tb_scr_0hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

counts_df = median_based_count_matrix_filtering(counts_df)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df)
degs_df=get_annotated_deg_df(degs_df = degs_df)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_igra_negative_scr_0hr_vs_tb_scr_0hr.csv",row.names = FALSE)
get_degs_with_pathways(degs_csv_path = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_igra_negative_scr_0hr_vs_tb_scr_0hr.csv",excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/excel_reports",excel_report_filename = "degs_and_pathways_igra_negative_scr_0hr_vs_tb_scr_0hr.xlsx")
#*******
# igra -ve scr 2 vs tb 2
#*******
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_igra_negative_scr_2hr_vs_tb_scr_2hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

counts_df = median_based_count_matrix_filtering(counts_df)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df)
degs_df=get_annotated_deg_df(degs_df = degs_df)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_igra_negative_scr_2hr_vs_tb_scr_2hr.csv",row.names = FALSE)
get_degs_with_pathways(degs_csv_path = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_igra_negative_scr_2hr_vs_tb_scr_2hr.csv",excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/excel_reports",excel_report_filename = "degs_and_pathways_igra_negative_scr_2hr_vs_tb_scr_2hr.xlsx")

# #*******
# # igra -ve scr 24 vs tb 24
# #*******
# #load the counts file
# counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_igra_negative_scr_2hr_vs_tb_scr_2hr.csv",header = TRUE,stringsAsFactors = FALSE)
# 
# #modify the counts data 
# rownames(counts_df)=counts_df[,1]
# counts_df=counts_df[,-1]
# 
# counts_df = median_based_count_matrix_filtering(counts_df)
# 
# #get differential expression standardised 
# degs_df=get_differential_expression_results_df(counts_df)
# degs_df=get_annotated_deg_df(degs_df = degs_df)
# write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_igra_negative_scr_2hr_vs_tb_scr_2hr.csv",row.names = FALSE)

#*******
# igra -ve scr 24 vs tb 24
#*******
#load the counts file
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_igra_negative_scr_24hr_vs_tb_scr_24hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]

counts_df = median_based_count_matrix_filtering(counts_df)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df)
degs_df=get_annotated_deg_df(degs_df = degs_df)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_igra_negative_scr_24hr_vs_tb_scr_24hr.csv",row.names = FALSE)
get_degs_with_pathways(degs_csv_path = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_igra_negative_scr_24hr_vs_tb_scr_24hr.csv",excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/excel_reports",excel_report_filename = "degs_and_pathways_igra_negative_scr_24hr_vs_tb_scr_24hr.xlsx")
############################################################################################
#igra -ve b3 0 hr vs tb b3 0 hr
############################################################################################

#*************************
#load 0 hr b3 counts file 
#*************************
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_igra_negative_b3_0hr_vs_tb_b3_0hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]
counts_df = median_based_count_matrix_filtering(counts_df)
columns=colnames(counts_df)
columns[str_detect(columns,pattern = "SCO149")]="SCO149"
columns[str_detect(columns,pattern = "SCO68")]="SCO68"
columns[str_detect(columns,pattern = "SCO71")]="SCO71"
colnames(counts_df)=columns
counts_df$SCO71=round((counts_df$SCO149 + counts_df$SCO68)/2)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df)
degs_df=get_annotated_deg_df(degs_df = degs_df)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_igra_negative_b3_0hr_vs_tb_b3_0hr.csv",row.names = FALSE)
get_degs_with_pathways(degs_csv_path = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_igra_negative_b3_0hr_vs_tb_b3_0hr.csv",excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/excel_reports",excel_report_filename = "degs_and_pathways_igra_negative_b3_0hr_vs_tb_b3_0hr.xlsx")
#*************************
#load 2 hr b3 counts file 
#*************************
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_igra_negative_b3_2hr_vs_tb_b3_2hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]
counts_df = median_based_count_matrix_filtering(counts_df)
columns=colnames(counts_df)
columns[str_detect(columns,pattern = "SCO149")]="SCO149"
columns[str_detect(columns,pattern = "SCO68")]="SCO68"
columns[str_detect(columns,pattern = "SCO71")]="SCO71"
colnames(counts_df)=columns
counts_df$SCO71=round((counts_df$SCO149 + counts_df$SCO68)/2)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df)
degs_df=get_annotated_deg_df(degs_df = degs_df)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_igra_negative_b3_2hr_vs_tb_b3_2hr.csv",row.names = FALSE)
get_degs_with_pathways(degs_csv_path = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_igra_negative_b3_2hr_vs_tb_b3_2hr.csv",excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/excel_reports",excel_report_filename = "degs_and_pathways_igra_negative_b3_2hr_vs_tb_b3_2hr.xlsx")
#*************************
#load 24 hr b3 counts file 
#*************************
counts_df=read.csv(file="/home/cidr/Documents/work/fresh/automatedCounts/counts_igra_negative_b3_24hr_vs_tb_b3_24hr.csv",header = TRUE,stringsAsFactors = FALSE)

#modify the counts data 
rownames(counts_df)=counts_df[,1]
counts_df=counts_df[,-1]
counts_df = median_based_count_matrix_filtering(counts_df)
columns=colnames(counts_df)
columns[str_detect(columns,pattern = "SCO149")]="SCO149"
columns[str_detect(columns,pattern = "SCO68")]="SCO68"
columns[str_detect(columns,pattern = "SCO71")]="SCO71"
colnames(counts_df)=columns
counts_df$SCO71=round((counts_df$SCO149 + counts_df$SCO68)/2)

#get differential expression standardised 
degs_df=get_differential_expression_results_df(counts_df)
degs_df=get_annotated_deg_df(degs_df = degs_df)
write.csv(degs_df,file="/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_igra_negative_b3_24hr_vs_tb_b3_24hr.csv",row.names = FALSE)
get_degs_with_pathways(degs_csv_path = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_igra_negative_b3_24hr_vs_tb_b3_24hr.csv",excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/excel_reports",excel_report_filename = "degs_and_pathways_igra_negative_b3_24hr_vs_tb_b3_24hr.xlsx")


#generating mega reports
#comparing tb scr longitudinal data
get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_0hr_vs_tb_scr_2hr.csv", degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_0hr_vs_tb_scr_24hr.csv",set1_name = "scr_0_vs_2", set2_name = "scr_0_vs_24", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/mega_reports", pathwaysType ="gsea")

#comparing tb b3 longitudinal data
get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_b3_0hr_vs_tb_b3_2hr.csv", degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_b3_0hr_vs_tb_b3_24hr.csv",set1_name = "b3_0_vs_2", set2_name = "b3_0_vs_24", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/mega_reports", pathwaysType ="gsea")

#comparing longitudinal data between scr and b3
#comparing tb scr and tb b3 longitudinal data
get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_0hr_vs_tb_scr_2hr.csv", degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_b3_0hr_vs_tb_b3_2hr.csv",set1_name = "scr_0_vs_2", set2_name = "b3_0_vs_2", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/mega_reports", pathwaysType ="gsea")

#comparing tb scr and tb b3 longitudinal data
get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_0hr_vs_tb_scr_24hr.csv", degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_b3_0hr_vs_tb_b3_24hr.csv",set1_name = "scr_0_vs_24", set2_name = "b3_0_vs_24", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/mega_reports", pathwaysType ="gsea")

#comparing igra negative longitudinal data and tb longitudinal data
#comparing scr longitudinal
#comparing tb scr and tb b3 longitudinal data
#2hr 
get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_2hr.csv", degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_0hr_vs_tb_scr_2hr.csv",set1_name = "ign_scr_0_v_2", set2_name = "tb_scr_0_v_2", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/mega_reports", pathwaysType ="gsea")

#24hr
get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_scr_0hr_vs_igra_negative_scr_24hr.csv", degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_scr_0hr_vs_tb_scr_24hr.csv",set1_name = "ign_scr_0_v_24", set2_name = "tb_scr_0_v_24", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/mega_reports", pathwaysType ="gsea")

#comparing with igra negative b3 data to tb b3 data
#2hr 
get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_2hr.csv", degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_b3_0hr_vs_tb_b3_2hr.csv",set1_name = "ign_b3_0_v_2", set2_name = "tb_b3_0_v_2", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/mega_reports", pathwaysType ="gsea")

#24hr
get_mega_report(degs1_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/with_total_trascripts_gt_10/degs_csv/diffEx_igra_negative_b3_0hr_vs_igra_negative_b3_24hr.csv", degs2_filepath = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/degs_csv_files/diffEx_tb_b3_0hr_vs_tb_b3_24hr.csv",set1_name = "ign_b3_0_v_24", set2_name = "tb_b3_0_v_24", excel_report_output_dir = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/mega_reports", pathwaysType ="gsea")

#analysing igra_negative cross section with TB data
#searching for the genes whicha re also present in the signature

ig_tb_cross_scr_0hr = read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/raw_data_degs_csv_files/cross_sectional_with_igra_negative/cross_section_with__scr/diffEx_igra_negative_scr_0hr_vs_tb_scr_0hr.csv", stringsAsFactors = FALSE)
ig_tb_cross_scr_2hr = read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/raw_data_degs_csv_files/cross_sectional_with_igra_negative/cross_section_with__scr/diffEx_igra_negative_scr_2hr_vs_tb_scr_2hr.csv", stringsAsFactors = FALSE)
ig_tb_cross_scr_24hr = read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/raw_data_degs_csv_files/cross_sectional_with_igra_negative/cross_section_with__scr/diffEx_igra_negative_scr_24hr_vs_tb_scr_24hr.csv", stringsAsFactors = FALSE)

ig_tb_cross_b3_0hr = read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/raw_data_degs_csv_files/cross_sectional_with_igra_negative/cross_section_with_b3/diffEx_igra_negative_b3_0hr_vs_tb_b3_0hr.csv", stringsAsFactors = FALSE)
ig_tb_cross_b3_2hr = read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/raw_data_degs_csv_files/cross_sectional_with_igra_negative/cross_section_with_b3/diffEx_igra_negative_b3_2hr_vs_tb_b3_2hr.csv", stringsAsFactors = FALSE)
ig_tb_cross_b3_24hr = read.csv(file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/raw_data_degs_csv_files/cross_sectional_with_igra_negative/cross_section_with_b3/diffEx_igra_negative_b3_24hr_vs_tb_b3_24hr.csv", stringsAsFactors = FALSE)


found_sig_genes_scr_0hr_df = ig_tb_cross_scr_0hr %>% dplyr::filter(gene_symbol %in% selected_genes)
found_sig_genes_scr_2hr_df = ig_tb_cross_scr_2hr %>% dplyr::filter(gene_symbol %in% selected_genes)
found_sig_genes_scr_24hr_df = ig_tb_cross_scr_24hr %>% dplyr::filter(gene_symbol %in% selected_genes)

found_sig_genes_b3_0hr_df = ig_tb_cross_b3_0hr %>% dplyr::filter(gene_symbol %in% selected_genes)
found_sig_genes_b3_2hr_df = ig_tb_cross_b3_2hr %>% dplyr::filter(gene_symbol %in% selected_genes)
found_sig_genes_b3_24hr_df = ig_tb_cross_b3_24hr %>% dplyr::filter(gene_symbol %in% selected_genes)


write.csv(as.data.frame(found_sig_genes_scr_0hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_393_genes_ig_vs_tb_degs_scr_0hr.csv",row.names = FALSE)

write.csv(as.data.frame(found_sig_genes_scr_2hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_393_genes_ig_vs_tb_degs_scr_2hr.csv",row.names = FALSE)

write.csv(as.data.frame(found_sig_genes_scr_24hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_393_genes_ig_vs_tb_degs_scr_24hr.csv",row.names = FALSE)


write.csv(as.data.frame(found_sig_genes_b3_0hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_393_genes_ig_vs_tb_degs_b3_0hr.csv",row.names = FALSE)
write.csv(as.data.frame(found_sig_genes_b3_2hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_393_genes_ig_vs_tb_degs_b3_2hr.csv",row.names = FALSE)
write.csv(as.data.frame(found_sig_genes_b3_24hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_393_genes_ig_vs_tb_degs_b3_24hr.csv",row.names = FALSE)

#16 gene signature
selected_genes = c("ANKRD22","APOL1","BATF2","ETV7","FCGR1A","FCGR1B","GBP1","GBP2","GBP4","GBP5","SCARF1","SEPT4","SERPING1","STAT1","TAP1","TRAFD1")
selected_genes = sort(selected_genes)

write.csv(as.data.frame(found_sig_genes_scr_0hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_16_genes_ig_vs_tb_degs_scr_0hr.csv",row.names = FALSE)

write.csv(as.data.frame(found_sig_genes_scr_2hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_16_genes_ig_vs_tb_degs_scr_2hr.csv",row.names = FALSE)

write.csv(as.data.frame(found_sig_genes_scr_24hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_16_genes_ig_vs_tb_degs_scr_24hr.csv",row.names = FALSE)


write.csv(as.data.frame(found_sig_genes_b3_0hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_16_genes_ig_vs_tb_degs_b3_0hr.csv",row.names = FALSE)
write.csv(as.data.frame(found_sig_genes_b3_2hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_16_genes_ig_vs_tb_degs_b3_2hr.csv",row.names = FALSE)
write.csv(as.data.frame(found_sig_genes_b3_24hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_16_genes_ig_vs_tb_degs_b3_24hr.csv",row.names = FALSE)

selected_genes = c("AGMAT","ANKRD13A","ATG3","ATP1B3","ATP6V0E1","ATP6V0E1","BRSK1","C20orf24","C20orf24","C4orf34","C9orf109","C9orf127","CALCOCO2","CAPS","CASP4","CCR2","CCR7","CD74","CHI3L2","CLEC1A","DBI","DEPDC5","DHRS9","DMXL2","EFR3A","EIF4E3","EPB41L3","FCRL3","FER1L3","GLRX","GSTK1","HLA-F","KIAA1632","KPNB1","LACTB","LOC149448","LOC284701","LOC391811","LOC440348","LOC644086","LRRC37A4","LYSMD2","MCL1","MRPL44","MS4A6A","MTHFD2","NAT1","NFKBIB","NPC2","P2RY5","PDE7A","PIGU","POLB","POMP","PPP1R3D","PRCP","PSMA4","PSMA4","PSMB10","PSMB3","PTPRE","RAC1","RBCK1","RBMS1","RNASEL","SEC23B","SLC16A6","SLC7A6","THEM2","TLR7","TMC6","TMEM51","TNFRSF14","TRIB2","TYROBP","USP15","USP47","VCPIP1","WDR33","WSB2","ZMYND15")
selected_genes = sort(selected_genes)

write.csv(as.data.frame(found_sig_genes_scr_0hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_86_genes_ig_vs_tb_degs_scr_0hr.csv",row.names = FALSE)

write.csv(as.data.frame(found_sig_genes_scr_2hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_86_genes_ig_vs_tb_degs_scr_2hr.csv",row.names = FALSE)

write.csv(as.data.frame(found_sig_genes_scr_24hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_86_genes_ig_vs_tb_degs_scr_24hr.csv",row.names = FALSE)


write.csv(as.data.frame(found_sig_genes_b3_0hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_86_genes_ig_vs_tb_degs_b3_0hr.csv",row.names = FALSE)
write.csv(as.data.frame(found_sig_genes_b3_2hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_86_genes_ig_vs_tb_degs_b3_2hr.csv",row.names = FALSE)
write.csv(as.data.frame(found_sig_genes_b3_24hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_86_genes_ig_vs_tb_degs_b3_24hr.csv",row.names = FALSE)

selected_genes = c("FCGR1A","HK3","RAB13","RBBP8","IFI44L","TIMM10","BCL6","SMARCD3","CYP4F3","SLPI")
selected_genes = sort(selected_genes)
write.csv(as.data.frame(found_sig_genes_scr_0hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_10_genes_ig_vs_tb_degs_scr_0hr.csv",row.names = FALSE)

write.csv(as.data.frame(found_sig_genes_scr_2hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_10_genes_ig_vs_tb_degs_scr_2hr.csv",row.names = FALSE)

write.csv(as.data.frame(found_sig_genes_scr_24hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_10_genes_ig_vs_tb_degs_scr_24hr.csv",row.names = FALSE)


write.csv(as.data.frame(found_sig_genes_b3_0hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_10_genes_ig_vs_tb_degs_b3_0hr.csv",row.names = FALSE)
write.csv(as.data.frame(found_sig_genes_b3_2hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_10_genes_ig_vs_tb_degs_b3_2hr.csv",row.names = FALSE)
write.csv(as.data.frame(found_sig_genes_b3_24hr_df), file = "/home/cidr/Documents/work/fresh/SCO71_removed/reports_tb_2019_10_23/consolidated_tb_reports_10_feb_2020/found_sign_10_genes_ig_vs_tb_degs_b3_24hr.csv",row.names = FALSE)

#t cell specific markers search 
#from scr data
ig_tb_cross_scr_0hr %>% dplyr::filter(gene_symbol %in% ig_tb_cross_scr_2hr$gene_symbol[ig_tb_cross_scr_2hr$gene_symbol %in% ig_tb_cross_scr_24hr$gene_symbol])

#from b3 data
ig_tb_cross_b3_0hr %>% dplyr::filter(gene_symbol %in% ig_tb_cross_b3_2hr$gene_symbol[ig_tb_cross_b3_2hr$gene_symbol %in% ig_tb_cross_b3_24hr$gene_symbol])
