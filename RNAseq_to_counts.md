# Abundance estimation review

## Exploring the STAR parameters when creating an index

Using the STAR manual, investigate the parameters of the STAR command to answer the questions following each STAR command script.

1. Below is the script we ran in class to create the genome index for alignment with STAR. We have added line numbers to the script to answer the questions that follow (you will need to refer to the online manual for STAR to answer these questions):

  ```
  (1)  #!/bin/bash
  
  (2)  #SBATCH -p priority # partition name
  (3)  #SBATCH -t 0-2:00 # days-hours:minutes runlimit after which job will be killed.
  (4)  #SBATCH -c 6 # number of cores requested
  (5)  #SBATCH --job-name STAR_index         # Job name
  (6)  #SBATCH -o %j.out       # File to which standard out will be written
  (7)  #SBATCH -e %j.err       # File to which standard err will be written
  
  (8)  STAR --runThreadN 6 \
  (9)  --runMode genomeGenerate \
  (10) --genomeDir my_genome_index \
  (11) --genomeFastaFiles chr1.fa \
  (12) --sjdbGTFfile chr1-hg19_genes.gtf \
  (13) --sjdbOverhang 99
  ```

  a. Which line number(s) in the script would you change if you were aligning your reads to the entire genome (instead of chr1)?

  b. Provide the modified lines you would use in the script to use reference files from `/n/groups/shared_databases` on O2.

  c. Could you run this indexing/alignment in STAR if you had no GTF annotation file? If so, which line(s) would change, and how would the line(s) change? 

  d. Which line(s) in the script would you change if the length of your reads were 150nt? How would the line(s) change? 



## Exploring the Salmon quasi-alignment parameters

2. Below is the script containing the Salmon command to obtain abundance estimates for our sequencing reads mapped to the genome. We have added line numbers to the script to answer the questions that follow (you may need to refer to the [online manual for Salmon](https://salmon.readthedocs.io/en/latest/) to answer these questions):
  
  ```
  (1)  #!/bin/bash
  (2)  #SBATCH -p priority       # Partition to submit
  (3)  #SBATCH -c 6                  # Number of cores
  (4)  #SBATCH  -t 0-1:30               # Runtime in D-HH:MM (or use minutes)
  (5)  #SBATCH  —job-name  salmon_mov10         # Job name
  (6)  #SBATCH -o %j.out       # File to which standard out will be written
  (7)  #SBATCH -e %j.err       # File to which standard err will be written
  
  (8) salmon quant -i /n/groups/hbctraining/ngs-data-analysis-longcourse/rnaseq/salmon.ensembl38.idx \
  (9)  -l A \
  (10) -r ~/rnaseq/raw_data/Mov10_oe_1.subset.fq \
  (11) -o Mov10_oe_1.subset.salmon \
  (12) -p 6 \
  (13) —writeMappings=salmon.out \
  (14) —seqBias \
  (15) —useVBOpt
  ```

  a. Which line(s) in the script would you change if you wanted to change the number of cores you use to use only 3 cores? How would the line(s) change? 

  b. Would any line(s) in the script change if you used an unstranded library preparation method for sequencing? If so, how would the line(s) change? 

  c. Would any line(s) in the script change if you had paired-end data? If so, which line(s) would change, and how would the line(s) change?  

  d. After performing QC you realized that you have GC bias present in your data. Is there anything you would change to address this when mapping and quantifying with Salmon?

## Positional parameters for scripting

3. Use vim to create a script by copying and pasting the contents of [this page](https://steve-parker.org/sh/eg/var3.sh.txt). Call it `pos_param_test.sh`. Run the new script as follows and report back on what the contents of the variable `$0`, `$1`, `$2` and `$@` are in each case.

  a. `sh pos_param_test.sh ngs_course 2018 "4" 6`

  b. `sh pos_param_test.sh this is easy to understand`

  c. `sh pos_param_test.sh 3 x 5 = 15 , 5 x 5 = 25`


4. Open up vim to create a shell script called `run_salmon_single sample.sh`. Add a shebang line to the top of your script. Copy and paste the Salmon command below into this script and assume it is running on a file called `Mov10_oe_1.subset.fq`. 

  ```
  #!/bin/bash
  
  salmon quant -i /n/groups/hbctraining/ngs-data-analysis-longcourse/rnaseq/salmon.ensembl38.idx \
   -l A \
   -r ~/rnaseq/raw_data/Mov10_oe_1.subset.fq \
   -o Mov10_oe_1.subset.salmon \
   --writeMappings=salmon.out \
   --seqBias \
   --useVBOpt 
  ```

Now make the following modifications:

* Remove the argument that writes the mappings to a SAM file
* Before the Salmon command add a line to create a variable called `fq` and assign it the value of `$1` (the filename that we expect the user to provide)
* Change your command in the appropriate place(s) to now accept the value from this positional parameter as input to Salmon (which is now stored in the `fq` variable). _**Hint:** be sure to reference the variable with a `$`_
* Use the `fq` variable to create a new variable called `base` that stores a prefix for the output. Change the Salmon command to utilize the `base` variable for naming out output directory.

