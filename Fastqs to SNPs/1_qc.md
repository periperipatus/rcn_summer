# Quality Control On Reads (QC)

This tutorial was conceived by Chris Balakrishnan and further modified by Peri Bolton.

# Table of Contents

* [Files Preparation](#files-preparation)
* [FastQC](#fastqc)
* [Quality Trimming](#quality-trimming)
* [Compare with results from resequencing data](#compare-with-results-from-resequencing-data]
* [QC resources](#qc-resources)

# Learning Objectives
* Understand structure of .fastq files
* Conduct subsampling of .fastq files
* Interpret quality control reports on .fastq files. 
* Understand and implement .fastq quality trimming
* Understand paired-end reads and implement appropriate commands


# Files Preparation

Move into your named directory on the server, then create a new directory called ```MAVI_reads``` and move into it.

Before we get onto the quality control, why don't you have a look at the fastq file. These are the raw reads provided by the sequencing facility
```zcat``` lets you view zipped files, which are compressed and thus unreadable

```bash
zcat ../../Bioinformatics_Workshop/mavi_reads/7_ATTACTCG-ACGTCCTG_L00M_R1_001.fastq.gz | head
```

**Question 1:** what does the ```../../``` mean?

Because we don't want to take up too much space on the server, and I want to leave time for you all to run your assigned projects, we are going to cut it to the first 400k lines of the file. 
Please note that you would NOT trim your real data like this. 

```bash
zcat ../../Bioinformatics_Workshop/mavi_reads/7_ATTACTCG-ACGTCCTG_L00M_R1_001.fastq.gz | head -n 400000 > ./7_MAVI_SH_JB1_F.fastq
```

**Question 2:** How many reads does this represent?


Now refer to the file ```MAVI_sample_key.xlsx```. This will help to explain the logic of the file name we just made, also explained below:

7 = identification of individual manakin tissue sample in a dataset that originally included multiple species
MAVI = Manacus vitellinus, the species. Golden-collared Manakin
SH = tissue type. SH is scapulohumeralis caudalis, and is involved the fast twitch movements that we know and love. 
PEC = pectoralis, involved in flight
JB1 = identification of individual manakin.
F = file contains forward reads (R1). R is reverse (R2).


**On your own:** Now copy across all the files in ```~/Bioinformatics_Workshop/mavi_reads/``` using the same process described above. *Remember* to keep track of Forward and Reverse reads. 


<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>


# FastQC

We need to run a program for quality checking of our raw reads, the standard is [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

```bash
~/programs/fastqc  7_MAVI_SH_JB1_F.fastq
```

Open a new terminal tab, or a new terminal window.
Then, download the .html file to your your folder for the tutorial.

```bash
scp -P 1200 <USERNAME>@<SERVER-IP>:~/<YOUR DIR>/MAVI_read_mapping/7_MAVI_SH_F_fastqc.html ./
```
Check the slack for your login and server IP address. 

double click on the html file in your desktop




**A Note For Windows Users:** Because your Linux WSL represents everything into a single unified file heirarchy it's a bit different to access your desktop through the commandline
```/mnt/c/Users/<YOUR-USERNAME>/Desktop``` is usually where it is. Ubuntu is also shortening the Username so if you can't find it right away just navigate to the `/mnt/c/Users/` directory and look for what looks right for your computer.

**Question 3:** What type of quality encoding scheme is used for this file? 

**Question 4:** What is the length of the reads?

**Question 5:** Does is pass all QC? If no, what fails?


**On your own:** Run fastQC on all the fastq files in the directory. Think about how you would do this without entering an unique command for every single file. 

** Question 6:** What command did you use to run fastQC on all the files?

**On your own:** The program [MultiQC](https://multiqc.info/) invoked using ``multiqc`` can take multiple FastQC reports and summarise them into one report. Use this on the server to analyse all your reports and summarise the results. 


<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>


# Quality Trimming

There are many programs that perform quality trimming. Including Trimmomatic and AdapterRemoval. We are using trim_galore, which is based off cutadapt and made by the same people as FastQC. 
trim_galore has a handy feature that it includes a database of the adapter sequences. 

Going back to the terminal that is connected to the server, we are going to trim out low quality sequences. 

```bash
trim_galore --paired --length 40 7_MAVI_SH_JB1_F.fastq 7_MAVI_SH_JB1_R.fastq
```

**Question 7:** Why are there two files here? Why do they need to be processed together?

**Question 8:** How many bases were quality trimmed? Look at the output from trim_galore on the commandline screen.

**Question 9:** How many times was the adaptor found? 

Then conduct a fastqc analysis on one of the trimmed files. 

**Question 10:** How does it look now? any improvement in terms of pass/fail?

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>


## Repeat on the remaining samples

If you haven't already, copy across the other ```.fastq.gz``` files in the folder ```Bioinformatics_Workshop/mavi_reads/``` in the way described in the previous section. Then run quality trimming.  

**On your own:** Use your loop writing knowledge to write a shell script that will loop over all the pairs of fastq files and trim them. 

**Question 11:** Share your code for this with the group. 

# Compare with results from resequencing data

In this ```~/Bioinformatics_Workshop/zf_reads/``` directory there are results from paired end Illumina data on wild zebra finches from Singhal et al. (2015) Stable recombination hotspots in birds. Science 350:928-932. 
All data for this paper are deposited here, including raw files and alignments https://www.ebi.ac.uk/ena/browser/view/PRJEB10586

To directly download these data I downloaded a .tsv report and copied the links to the forward and reverse files like so:
```bash
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/009/ERR1013179/ERR1013179_2.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/009/ERR1013179/ERR1013179_1.fastq.gz
``` 
DO NOT RUN the above code: these files are large and take forever to download (4hrs each).


Now, make a subsample of that file in your directory.
```bash
seqtk sample -s42 ERR1013179_1.fastq.gz 0.1 > <YOUR DIR>/ERR1013179_1_subs.fastq
seqtk sample -s42 ERR1013179_2.fastq.gz 0.1 > <YOUR DIR>/ERR1013179_2_subs.fastq
``` 
**Question 12:** What is the difference between the subsampling procedure using ```seqtk``` here and the one using ```head``` above?

**On your own:** Run ```fastqc``` on the subsampled files in your directory and download the .html file to your desktop to view. 

**Question 13:** Observe any differences in the fastqc results between the resequencing (Zebra finch) and the RNA-seq data (Manacus).


# QC resources

* BioInfo-Core has a great general discussion on assessing read sequence quality using FastQC [here](http://bioinfo-core.org/index.php/9th_Discussion-28_October_2010)
* HBCtraining has a great workshop series on RNAseq data including interpreting FastQC results [here] (https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html)
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/) documentation is also always a good place to start.
* Sheng et al (2017) Multi-perspective quality control of Illumina RNA sequencing data analysis. Briefings in Functional Genomics
* MacManes et al (2014) On the optimal trimming of high-throughput mRNA sequence data. Frontiers in Genetics https://doi.org/10.3389/fgene.2014.00013
* Fabbro et al (2013) An extensive evaluation of read trimming effects on Illumina NGS data analysis. PLoS One https://doi.org/10.1371/journal.pone.0085024
* Bush et al (2020) Read trimming has minimal effect on bacterial SNP-calling accuracy. Microbial Genomics 6: https://doi.org/10.1099/mgen.0.000434