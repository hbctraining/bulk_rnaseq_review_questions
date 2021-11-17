Running the complete RNA-seq workflow 
#Using bcbio-nextgen to run the QC analysis and Salmon  
#move data to /n/scratch3
$ cp /n/groups/hbctraining/kidney_fibrosis_rnaseq.tar.gz /n/scratch3/ecommons_id/
#create csv form of metadata according to original metadata 
#create configuration template 
$ vim smoc2-template.yaml 
#Apply the template to all samples 
$bcbio_nextgen.py –w template smoc2-template.yaml smoc2_project.csv *.fq 
#copy and modify sbatch job to run bcbio-nextgen on all samples 
$vim submit_bcbio.sbatch 
 
#Run RNA-seq analysis using automated workflow 
#create indice for  mm10 transcriptome 
$ salmon index -t Mus_musculus.GRCm38.cdna.all.fa -i transcripts_index --type quasi -k 31 
#adjust the scripts of rnaseq_analysis_on_input_file.sh and rnaseq_analysis_on_allfiles_for_slurm.sh 
#run MultiQC $ multiqc -n multiqc_report_rnaseq \ > /n/scratch2/xh67/homeworK/*zip \ > /n/scratch2/xh67/homeworK/results/STAR/*Log.final.out \ > /n/scratch2/xh67/homeworK/results/qualimap/* \ > /n/scratch2/xh67/homeworK/results/salmon/*salmon 
 
#DESeq Analysis 
#import data from salmon, create DESeq object, get normalized counts 
#list all directories containing data 
samples <- list.files(path = "./data/salmon", full.name = T, pattern = "\\.salmon$") 
samples <- paste0(samples, "/quant") 
#Obtain a vector of all filenames including the path 
files <- file.path(samples, "quant.sf") 
#have names for each element 
names(files) <- str_replace(samples, "./data/salmon/", "") %>%  
  str_replace(".salmon/quant", "") 
files 
#Look at the tx2gene table for GrCm38 
grcm38_tx2gene 
 
#import gene level salmon data through tximport 
txi <- tximport(files, type="salmon", tx2gene=grcm38_tx2gene, countsFromAbundance = 
"lengthScaledTPM") 
 
#Create metadata 
sampletype <- factor(c(rep("Smoc2_overexpression", 7), rep("Wildtype", 7))) 
meta <- data.frame(sampletype, row.names = colnames(txi$counts)) 
#add treatment info to metadata 
treatment <- c(rep("none", 3), rep("7day_UUO", 4), rep("none", 3), rep("7day_UUO", 4)) 
meta <- cbind(meta, treatment) 
 
#check row names of metadata equal to column names of the raw counts data 
all(colnames(txi$counts) == rownames(meta)) 
#create DESeq2Dataset object, group out treatment effect 
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~treatment+sampletype) 
 
View(counts(dds)) 
#mediam of ratio normalization 
dds <- estimateSizeFactors(dds) 
sizeFactors(dds) 
normalized_counts <- counts(dds, normalized=TRUE) 
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names = NA) 
 
#Run analysis 
dds <-DESeq(dds) 
 
#Define constrsts, extract results table, shrink the log2 fold changes 
contrast_oe <- c("sampletype", "Smoc2_overexpression", "Wildtype") 
res_tableOE_unshrunken <- results(dds, contrast = contrast_oe, alpha=0.05) 
res_tableOE <- lfcShrink(dds, contrast = contrast_oe, res=res_tableOE_unshrunken) 
plotMA(res_tableOE, ylim=c(-2,2)) 
 > summary(res_tableOE, alpha=0.05)  
out of 31947 with nonzero total read count 
adjusted p-value < 0.05 LFC > 0 (up)       : 287, 0.9% 
LFC < 0 (down)     : 159, 0.5% 
outliers [1]       : 135, 0.42% low counts [2]     : 10857, 34% 
(mean count < 3) 
  
#Set thresholds 
padj.cutoff <- 0.05 
 
#extract significant differentially expressed genes 
res_tableOE_tb <- res_tableOE %>%  
  data.frame() %>%  
  rownames_to_column(var="gene") %>%  
  as_tibble() 
sig_Smoc2 <- res_tableOE_tb %>%  
  filter(padj < padj.cutoff) 
 
write.csv(sig_Smoc2, file = "DE_genes_from_Smoc2_overexpression_vs_Wildtype.csv") 
 
#PLOT TOP 20 GENES 
 
##Sleuth Analysis 
 
#Obtain a vector of all filenames including the path 
sf_dirs <- file.path(samples, "abundance.h5") 
#have names for each element 
names(sf_dirs) <- str_replace(samples, "./data/salmon/", "") %>%  
  str_replace(".salmon/quant", "") 
sf_dirs 
all(names(sf_dirs) == rownames(meta)) 
meta$sample <- rownames(meta) 
meta 
meta$path <- sf_dirs 
meta 
design <- ~treatment + sampletype 
 
t2g <- grcm38 %>%  
  select(symbol, ensgene) %>%  
  dplyr :: inner_join(grcm38_tx2gene) %>%  
  rename(target_id = enstxp, ens_gene = ensgene, ext_gene = symbol) 
#remove duplicated transcript IDs 
t2g <- t2g[which(!(duplicated(t2g$target_id))),] 
#Create sleuth object for analysis 
#relevel basal level 
meta$sampletype <- relevel(meta$sampletype, ref = "Wildtype") 
meta$treatment <- relevel(meta$treatment, ref = "none") 
so <- sleuth_prep(meta, full_model = design, target_mapping = t2g, read_bootstrap_tpm = T, 
extra_bootstrap_summary = T, transformation_function = function(x) log2(x+0.5)) 
so <- sleuth_fit(so) 
models(so) 
 
Smoc2 <- sleuth_wt(so, which_beta = 'sampletypeSmoc2_overexpression') 
sleuth_results_Smoc2 <- sleuth_results(Smoc2, test = 'sampletypeSmoc2_overexpression', 
show_all = T)  
 
plot_pca(Smoc2, color_by = 'sampletype', text_labels = T) 
 
#PCA doesn’t show a good cluster trend according to sampletype 
 
#identify the significant different transcripts 
sig_transcripts <- sleuth_results_Smoc2 %>% filter(qval < 0.05) 
 
write.csv(sig_transcripts, file = "Smoc2_sig_transcripts_sleuth.csv") 
plot_transcript_heatmap(Smoc2, transcripts = sig_transcripts$target_id) 
 
