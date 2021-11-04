# Running the complete RNA-seq workflow 

We have downloaded the raw FASTQ files from the SRA for the sequencing data used in the paper: Silencing SMOC2 ameliorates kidney fibrosis by inhibiting fibroblast to myofibroblast transformation (Links to an external site.). The paper explores kidney fibrosis in wildtype and SMOC2-overexpressing mice. 

a. Copy the `/n/groups/hbctraining/kidney_fibrosis_rnaseq.tar.gz` folder to your own `/n/scratch2` directory. Look inside the directory, you should find the following:

* data folder with the raw fastq files
* meta folder with a metadata file containing information about each of the samples
* reference_transcriptome folder containing the transcriptome FASTA

b. Setup a file system structure your project (i.e. create subdirectories and additional directories where you feel is necessary). 

c. Using the workflow and submission scripts we generated in class, parallelize the RNA-Seq analysis of the all files in this dataset. For each FASTQ file you will need to perform the following:

* Run FastQC
* Generate abundance estimates with Salmon (be sure to include the parameter to generate bootstraps)
* Evaluate the MultiQC report
**HINT: You will need to create a mouse index for Salmon. 

d. Use the Salmon quant files to perform differential gene expression analysis with DESeq2. Report the DE genes identified for SMOC2 over-expressing samples relative to wild type using an FDR threshold of 0.05, by uploading a csv file containing the results table output by DESeq2,but only for the significant genes.

Also generate a heatmap of the top 20 differentially expressed genes and include an annotation bar assigning samples to each group.

e. Use the Salmon output to run Sleuth and identify differentially expressed transcripts. Use Sleuth functions to create a PCA plot. Save and upload the image. How well do samples cluster?

Report the DE transcripts identified for SMOC2 over-expressing samples relative to wild type using an FDR threshold of 0.05, by uploading a csv file containing the results table output by Sleuth,but only for the significant transcripts.

Use Sleuth functions to plot a heatmap of the top 20 differentially expressed transcripts. Upload the image.

 
