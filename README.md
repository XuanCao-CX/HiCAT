HiCAT
=====
Tool for Hi-C data analysis.

![image](https://user-images.githubusercontent.com/57889560/113291933-54762000-9326-11eb-9100-c730be12e39d.png)

Software Requirements
---

Bowtie: the latest version of BWA should be installed from http://bowtie-bio.sourceforge.net
Shell: awk
python 2.7


Usage
---
         -1              The input R1 fastq or fastq.gz files. Files seperated by comma.
         -2              The input R2 fastq or fastq.gz files. Files seperated by comma.
                         The order of files can paired with R1 files.
         -o              The output dir.'
         --label         The lable for each paired R1_R2 fastq files, seperated by comma.
         -g              The bowtie index.
         -p              Cpus. Default:1.
         -n              The mismatch number for bowtie mapping. Default=2.
         -5              Trim bases number of 5 prime.
         -3              Trim bases number of 3 prime.
         --re            The resites for cut unmapped fastq files. It is the cut site from 5'-> 3'.
                         Default, DPNII:GATC.'
        --genomic_dis   The length for remove genomic DNA.Default:500.
        --start_search_REsite       The starting step is search REsite.If set, you should set old_file_dir.
        --start_REsite_mapping      The starting step is search RE reads mapping.If set, you should set old_file_dir.
        --old_file_dir              If not start at the full length map step,set this to get former files.
        --keep_media_fq_files       If set, keep the median sam and fastq files. Default delete them.

Run
---
	HiCAT.py -1  A.R1.fastq,B.R1.fastq  -2  A.R2.fastq,B.R2.fastq  -o  out  --label A,B  -g  bowtie_index -p 10 -n 2 --re GATC



