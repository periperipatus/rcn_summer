# Read Mapping

This tutorial was conceived by Chris Balakrishnan and Peri Bolton.

# Table of Contents

* [Mapping resequencing data](#mapping-resequencing-data)
* [RNAseq mapping with STAR](#RNAseq-mapping-with-star)

# Learning Objectives
* Understand different types of indexing
* Understand the concept and when to use splice-aware alignment 
* Understand information encoding in .sam and .bam files 
* View information in .sam and .bam files
* Conduct steps necessary to align fastq files to a reference genome
* Implement commands to sort and index the bam file required for SNP calling and viewing. 

# Mapping resequencing data

In this tutorial  we are using the software [bwa](https://github.com/lh3/bwa), which stands for Burrows-Wheeler Aligner. Recall the Burrows-Wheeler index methodology from the lecture on mapping.
Another popular aligner for this application is ```Bowtie2```

Latest ZF genome from Vertebrate Genomes Project, downloaded from [NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/565/GCF_003957565.2_bTaeGut1.4.pri/) into this directory: ```/home/ngsclass/Bioinformatics_Workshop/zf_genome/```

```bash
nohup bwa index GCF_003957565.2_bTaeGut1.4.pri_genomic.fna &
``` 
No need to conduct this step. It takes a long time to create the index for bwa. Just go to the directory and look at the extra files created by the indexing process. They all have the same file name stem, but a different extension. 

Now we are going to map some subsampled reads that I made. The subsample we are using today is smaller than the subsamples we made yesterday because my test run on the 10% subsample was VERY slow. 
(I used the method for the manakins, by taking the first 40k lines). 

Before you start, make sure you've got a new folder in your directory called alignments
```bash
mkdir ~/<YOURDIR>/alignments/
```

```bash
REF=/home/ngsclass/Bioinformatics_Workshop/zf_genome/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna
FASTQS=/home/ngsclass/Bioinformatics_Workshop/zf_reads/subs_10k
bwa mem $REF $FASTQS/ERR1013179_1_subs_val_1.fq $FASTQS/ERR1013179_2_subs_val_2.fq > <YOURDIR>/alignments/ERR1013179.sam & 

less -S ERR1013179.sam
```

**Question 1.1:** Use your Unix loop writing skills to align all the subsampled fastq files in that directory. Share with the group. 

Recall our lecture that went through the SAM file specifications. 
Now we want to convert that SAM file into a BAM file which is a binary format that saves a lot of space. 

```bash
samtools view -b -S -o ERR1013179.bam ERR1013179.sam
```

**Question 1.2:** What do the -b and -S arguments mean?

Now, let's look at the BAM... 

```bash
less -S ERR1013179.bam
```

Errr? 
bamfiles are binary and thus not human readable.
what if you want to look at a bam file? you can convert it back to sam format, or you can do this:


```bash
samtools view ERR1013179.bam | less -S
```

Samtools includes a bunch of tools that can be used to interrogate and filter SAM and BAM files
For example the `-f` argument enables you to interrogate what information is in the ['flag'](https://broadinstitute.github.io/picard/explain-flags.html).
For example the flag 0x4 (or 4) indicates unmapped reads. 

```bash
samtools view -c -f 4 ERR1013179.bam
```

**Question 1.3:** How many reads are unmapped? What does the `-c` argument do?

```bash
samtools view -c -q 42 ERR1013179.bam
```

**Question 1.4:** How many reads are aligned with a quality of 42 (max quality score from bwa)?

**Question 1.5:** Consider how you might use these functions to filter out low quality alignments.

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

## Indexing and Sorting the BAM

For some applications, such as SNP calling, sam files have to be sorted by genomic position

```bash
samtools sort -o ERR1013179.bam.sorted ERR1013179.bam
```

If you have a library that was made with PCR, then its a good idea to mark PCR duplicates (the exact same fragment amplified multiply), so they can be easily ignored by downstream processes. Here, we will use the `java` program `PicardTools`
VALIDATION_STRINGENCY is set to SILENT as PicardTools documentation suggests it improves performance when you have variable length read data.
```bash
java -Xmx8G -jar ~/programs/picard.jar MarkDuplicates \
INPUT=ERR1013179.bam.sorted \
OUTPUT=ERR1013179_sorted_marked.bam \
METRICS_FILE=ERR1013179_metrics.txt \
ASSUME_SORT_ORDER=coordinate \
VALIDATION_STRINGENCY=SILENT
```


Like genomes, BAM and SAM files can be indexed also, and is required by many downstream applications. 
This command will create another file with the same name with a `.bai` appended to the end. It is important that these file names have the same prefix (everything before .bai) otherwise programs that require an indexed BAM wont know where to look.

```bash
samtools index ERR1013179_sorted_marked.bam
```

Sorted and Indexed files are required for viewing the alignment in IGV. The sorted and indexed alignment can also be viewed in the CLI using the following command:


```bash
REF=/home/ngsclass/Bioinformatics_Workshop/zf_genome/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna

samtools tview -p NC_054767.1:477814 ERR1013179_sorted_marked.bam $REF
```
"." indicates matches on the forward strand, and "," indicates a match on the reverse strand. Use your arrows to navigate. 
Underlined reads mean secondary mapping or an unpaired read... Colours represent some kind of base level quality score, but exactly what the breakdown is I don't know. Welcome to the world of bad bioinformatics documentation.

What happens if you just view the alignment without the positional information?


To view the alignment in [IGV](http://software.broadinstitute.org/software/igv/alignmentdata) you will need to first create an `.fai` index of the genome in samtools. 

I already did this with the following command: 

```bash
samtools faidx GCF_003957565.2_bTaeGut1.4.pri_genomic.fna
```

Then you will need to download the reference `fasta` file, its samtools index, and the sorted bam file and its index to your desktop, load it into IGV and explore. 


## Alignment QC

Now, let's use `qualimap` to have a look at our alignment. 

Let's make a new directory for our quality reports and then run the program. 

```bash
mkdir qc

qualimap bamqc \
-outdir qc/ERR1013179 \
-bam ERR1013179_sorted_marked.bam \
--java-mem-size=8G
```

**Question 1.6:** What do the \ do in the code?

Download the all the results to your desktop to view. This time you need to download the whole folder.

```bash
scp -r -P 1200 ngsclass@<IP.ADRESS>:~/<YOURDIR>/alignments/qc/ERR1013179/ .
```
Open up the results in your browser. 

**Question 1.7:** What do you notice about the alignments?


<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

# RNAseq mapping with STAR

In yesterday's tutorial we also conducted QC on RNAseq data. If you would like to compare the results from the whole-genome resequencing to RNAseq do this part of the tutorial. 

If you are working on species with a genome, you can do splice aware mapping.
Some commonly used splice-aware aligners are [STAR](https://github.com/alexdobin/STAR) and [HISAT2](http://daehwankimlab.github.io/hisat2/about/).
Here we are going to use STAR (Spliced Transcripts Alignment to a Reference) to map our *Manacus* reads to the reference genome (Yes, I lied above, we have a genome assembly). 


Just like ```bwa```, the ```STAR``` aligner needs to make a genome index so that it can efficiently access the genome. 
However, STAR is a little slow at doing the index. So genome index has already been made for you here ```~/Bioinformatics_Workshop/mavi_genome```

Here's the code used to arrive at this point. 

```bash
cd ~/Bioinformatics_Workshop/
mkdir mavi_genome/
cd mavi_genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/715/985/GCF_001715985.3_ASM171598v3/GCF_001715985.3_ASM171598v3_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/715/985/GCF_001715985.3_ASM171598v3/GCF_001715985.3_ASM171598v3_genomic.gtf.gz

gunzip *.gz #uncompress files because STAR doesn't like compressed files. 


STAR --runThreadN 50 \
--runMode genomeGenerate \
--genomeDir mavi_index \
--genomeFastaFiles GCF_001715985.3_ASM171598v3_genomic.fna \
--sjdbGTFfile GCF_001715985.3_ASM171598v3_genomic.gtf \
--sjdbOverhang 99  \
--genomeSAindexNbases 13 

#this took about 20 min
```


**Question 2.1:** What does the ```--sjdbOverhang``` argument do?



Now, you can run one of the trimmed fastqs from before. 

```bash

STAR --genomeDir ~/Bioinformatics_Workshop/mavi_genome/mavi_index/ \
--runThreadN 6 \
--readFilesIn <YOURDIR>/<SUBDIR>/7_MAVI_SH_JB1_F_val_1.fq 7_MAVI_SH_JB1_R_val_2.fq \
--outFileNamePrefix ../alignments/7_MAVI_SH_JB1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes Standard \
--quantMode GeneCounts &
```

**Question 2.2:** What do the different arguments do in this? refer to the STAR manual


Now, let's move into our ```alignments/``` folder and look at what we've produced. There should be a number of files there including those ending in ```Aligned.sortedByCoord.out.bam```, ```Log.final.out``` and ```PerGene.out.tab```

First we can look at the bam file. This is our reads aligned. Because they are compressed in BGZF format you need to use a special program to make it human readable. 


```bash
samtools view 7_MAVI_SH_JB1Aligned.sortedByCoord.out.bam | head
```
This command shows the first few lines of the alignment section of the file. if you want to view the header you need to use the ```-h``` flag.

**Question 2.3:** What do the first 5 columns of the alignment tell you?

Now, let's open the ```Log.final.out``` file. 

**Question 2.4:** What is the percent of unmapped, multi-mapped and uniquely mapped reads? What could affect these numbers?

**Question 2.5:** What is the Number of splices? What does this mean?

Now let's have a look at the contents of  ```ReadsPerGene.out.tab```

```bash
less -S 7_MAVI_SH_JB1ReadsPerGene.out.tab
```


**Question 2.6:** What does this file tell us? What do the three columns mean?

Don't worry that there are a lot of zeros - don't forget you're using a truncated fastq file.  

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

## Alignment QC

Now, let's use `qualimap` to have a look at our alignment. 

Let's make a new directory for our quality reports and then run the program. 

```bash
mkdir qc

qualimap rnaseq \
-outdir qc/7_MAVI_SH_JB1_F_quali \
-a proportional \
-bam 7_MAVI_SH_JB1Aligned.sortedByCoord.out.bam \
-p strand-specific-reverse \
-gtf ~/Bioinformatics_Workshop/mavi_genome/GCF_001715985.3_ASM171598v3_genomic.gtf \
--java-mem-size=8G
```
Download the all the results to your desktop to view. This time you need to download the whole folder.

```bash
scp -r -P 1200 ngsclass@<IP.ADRESS>:~/<YOURDIR>/alignments/qc/7_MAVI_SH_JB1_F_quali/ .
```

**Question 2.7:** what does the `-r` flag do?

**Question 2.8:** Think about the quality of the mapping. Do we have good mapping rates? Why/why not? Do we have good representation in exonic regions?

