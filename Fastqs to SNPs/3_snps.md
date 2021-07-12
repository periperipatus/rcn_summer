# SNP Calling and Filtering

# Pre-reading

* Prof Rachel White (2021) BIO 336 Week 3: DNA sequencing coverage & depth (length 10:49) https://www.youtube.com/watch?v=sKlfGZnSn5g

# Learning Objectives for tutorial and lecture
* Implement basic GATK HaplotypeCaller pipeline
* Understand genotype likelihoods
* Understand the information encoded in a VCF file. 
* Understand how sequencing depth influences distribution of variants
* Understand ReadGroups
* Understand data quality and filtering procedures 
* Implement basic filtering principles

# Table of Contents
* [Tutorial](#Tutorial)
	* [Pre-Processing](#pre-processing)
	* [Mapping](#mapping)
	* [Prepping Alignments](#prepping-alignments)
	* [SNP calling](#snp-calling)
	* [Filtering](#Filtering)
		
* [Useful Resources](#useful-resources)

# Tutorial

For today's tutorial we will start with some bam files that I have processed for input into GATK snp calling following their [germline best practices pipeline](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

You do NOT need to repeat the steps until the end of the section on SNP calling, and the [FILTERING](#Filtering) section. HOWEVER, there ARE questions in Mapping and Prepping Alignments/ReadGroups throughout the tutorial so please review the code. 


## Pre-Processing
1. subset.sh
```bash
cd /home/ngsclass/Bioinformatics_Workshop/zf_reads/
FASTQ=*.fastq.gz
for f in $FASTQ
sample=$( echo $f | sed 's/.fastq.gz//')
seqtk sample -s42 $f 0.1 > subs/${sample}_subs.fq
done
```
2. trimming.sh
```bash
cd subs/
FORWARD=*_1_subs.fq

for file1 in $FORWARD
do
        file2=$( echo $file1 | sed 's/_1/_2/') ;
        trim_galore --paired --length 40  $file1  $file2 ;
done
```

## Mapping
3. mapping.sh
```bash
READS=*val_1.fq
REF=/home/ngsclass/Bioinformatics_Workshop/zf_genome/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna
OUTDIR=/home/ngsclass/Bolton/ZF/alignments
for file1 in $READS
        do
                file2=$( echo $file1 | sed 's/1_subs_val_1/2_subs_val_2/')
                sample=$( echo $file1 | sed 's/_val_1.fq//')
                echo "bwa mem $REF $file1 $file2 | samtools view -b -S -o $OUTDIR/$sample.bam &" >> $sample.mapping.sh
                nohup sh $sample.mapping.sh
        done
```

This script was written so that each mapping command can be run simultaneously. 
**Question 1:** Have a look at this script. How has this been run so that each file is run roughly contemporaneously? What are the drawbacks of running things using this approach? 
**Question 2:** What does the pipe to `samtools` do?

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

## Prepping Alignments


```cd /home/ngsclass/Bolton/ZF/alignments/```

Here, we follow mostly the same steps in the previous mapping tutorial, with one modification that marks Read Groups, a requirement for GATK. 


### 1. sort bams
```bash
BAMS=*.bam 
for b in $BAMS
	do 
		sample=$( echo $b | sed 's/.bam//') ;
		samtools sort $b -o $sample.sorted.bam ;
	done 
```

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

### 2. Read Groups 

When you get reads back from the sequencing facility, they are often run over multiple lanes, and this might be reflected in a single sample broken up over multiple files.

GATK requires some readgroup information for SNP calling, and can use it to differentiate batches and libraries in the SNP calling process. This can be particularly handy if you are including an individual that was prepped over multiple libraries, flowcells or technologies. 

Here, we are going to assume that a single sample represents a single library, run on a single lane. See the GATK page about ReadGroups for (more info)(https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)

But essentially, we have fields as follows:

* `LB` what library did this group of reads come from
* `PL` what sequencing machine did the reads come from (there are set values that are allowed). 
* `PU` What combination of Flowcell and Lane? 
* `SM` Sample ID - this will be used in the VCF file to identify your individual
* `ID` Not included in this tutorial, but can be used to make a master readgroup ID that is by default the Flowcell + Lane. However can be expanded to include the minimum unit of any suspected batch effects (e.g. libraries).


```bash
BAMS=*.sorted.bam

PU="120524_SN553_0254_BD0V62ACXX"

for b in $BAMS
do
sample=$( echo $b | sed 's/.marked.sorted.bam//')
LB=LIB.$sample
java -Xmx4g -jar ~/programs/picard.jar AddOrReplaceReadGroups \
I=$b \
O=$sample.sorted.rg.bam \
RGLB=$LB \
RGPL=ILLUMINA \
RGPU=$PU \
RGSM=$sample 
done
```

Now go to an example .bam file downloaded from EBI. 
```
INDIR=~/Bioinformatics_Workshop/zf_reads/
samtools view -H $INDIR/101.recal.bam | grep '^@RG'
```
**Question 3:** What does `samtools view -H` do?
**Question 4:** How many readgroups are there in the original ZF file?

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>


### 3. MarkDuplicates

We are running it AFTER Read Group assignment. It will mark duplicates according to the RGLB tag set in ```AddOrReplaceReadGroups```. If you have a bam that contains multiple libraries then you need to do it in this order. HOWEVER, this will not be true for most cases. 

```
BAMS=*.sorted.rg.bam
for b in $BAMS
do
	sample=$( echo $b | sed 's/1_subs.sorted.bam//')
	java -Xmx8G -jar ~/programs/picard.jar MarkDuplicates \
	INPUT=$b \
	OUTPUT=$sample.sorted.rg.marked.bam \
	METRICS_FILE=$sample.metrics.txt \
	ASSUME_SORT_ORDER=coordinate \
	VALIDATION_STRINGENCY=SILENT
done
```

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>


## SNP calling

Before we do anything we need to make a sequence dictionary - another kind of index for Genome Analysis Toolkit (GATK) to access the reference info

```bash
REFDIR=/home/ngsclass/Bioinformatics_Workshop/zf_genome
java -jar ~/programs/picard.jar CreateSequenceDictionary R=GCF_003957565.2_bTaeGut1.4.pri_genomic.fna O=GCF_003957565.2_bTaeGut1.4.pri_genomic.dict &
```

The next steps were actually run in a shell script. See the shell script associated with this folder. 


To start we need to activate the gatk conda environment. This will load all the necessary dependencies to run gatk.

```conda activate gatk``` 

To run this within a shell script you need to activate conda environment differently ```source /home/ngsclass/miniconda3/bin/activate gatk```

We are going to skip Base quality realignment (BQSR), because we don't have good genomic resources and so we're going to skip it. However, if you have high depth data for your own project you can use that and see the tutorials [here](http://protocols.faircloth-lab.org/en/latest/protocols-computer/analysis/analysis-gatk-parallel.html), and [here](https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html#bqsr).


Versions of GATK <3.7 will need to do an indel realignment step prior to using HaplotypeCaller. See this [link to explain how this works](https://qcb.ucla.edu/wp-content/uploads/sites/14/2016/03/GATKwr12-3-IndelRealignment.pdf). 
The indel realignment step is now part of the HaplotypeCaller in newer versions of gatk. 

```bash
BAMS=*.sorted.rg.marked.bam ;
OUTDIR=~/Bolton/ZF/alignments/tmp_gvcf ;
REF=/home/ngsclass/Bioinformatics_Workshop/zf_genome/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna ; 
for b in $BAMS
do
sample=$( echo $b | sed 's/.sorted.rg.marked.bam//')
gatk --java-options "-Xmx8g" HaplotypeCaller  \
   -R $REF \
   -I $b \
   -O $OUTDIR/$sample.g.vcf.gz \
   -ERC GVCF
done 
```

Despite the small size of these files this took about a day to run.

Then we need to make a database of our individual vcfs, but first we need to make a file that will contain a list of names and file locations of the g.vcfs.

```bash
VCFS=tmp_gvcf/*.vcf.gz
for v in $VCFS; do
    filename=$(basename $v);
    name="${filename%%.*}";
    echo -e "$name\t$v" >> gvcfs.sample_map;
done
```

GATK offers two ways to combine the individual GVCFs into one for genotype calling. The first/older method is `CombineGVCFs` or using `GenomicsDBImport`. I had trouble getting DBImport to work on the server, and `CombineGVCFs` is acceptable for use in small projects like this one. 

```bash
REF=/home/ngsclass/Bioinformatics_Workshop/zf_genome/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna ;
INDIR=/home/ngsclass/Bolton/ZF/alignments/tmp_gvcf
 gatk --java-options "-Xmx24g" CombineGVCFs \
   -R $REF \
   --variant $INDIR/ERR1013161_1_subs.g.vcf.gz  \
   --variant $INDIR/ERR1013162_1_subs.g.vcf.gz \
   --variant $INDIR/ERR1013163_1_subs.g.vcf.gz \
   --variant $INDIR/ERR1013164_1_subs.g.vcf.gz \
   --variant $INDIR/ERR1013169_1_subs.g.vcf.gz \
   --variant $INDIR/ERR1013179_1_subs.g.vcf.gz \
   -O ../vcf/zfs.g.vcf.gz
```

Then we need to assign genotypes to the individuals.

```bash
gatk --java-options "-Xmx64g" GenotypeGVCFs \
      -R $REF \
      -V $INDIR/zfs.g.vcf.gz \
      -O /home/ngsclass/Bolton/ZF/vcf/zf_snps_raw.vcf.gz 
```


<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

### Taking stock

Now, let's have a look at the raw data. I have moved the vcf file to `~/Bioinformatics_workshop/vcfs` for you to use. 

In there there is another file from a whole genome resequencing project from the McRae Lab `bluebird_subset.vcf.gz`. 
Please copy these to your directory for the Filtering Tutorial.

```
bcftools stats <INVCF.vcf.gz> <OUTFILE>
```

**Question 6:** Use the above command for both VCF files. How many SNPs are there? How many Indels are there? How many multi-allelic sites?

**Question 7:**  Look at both the files in the directory. What fields are different in INFO and FORMAT? What program was used to call the genotypes for the bluebird data?

# Filtering

See the .Rmd  and the .html in this folder for the filtering exercise. You can access the .html file from the OneDrive folder and view it using google chrome. 

<br/>
<div align="right">
    <b><a href="#table-of-contents">^ back to TOC</a></b>
</div>
<br/>

# Useful resources

## GATK
https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html
https://evodify.com/gatk-in-non-model-organism/
https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
http://protocols.faircloth-lab.org/en/latest/protocols-computer/analysis/analysis-gatk-parallel.html
https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format
## Read Groups 
https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
https://eriqande.github.io/eca-bioinf-handbook/alignment-of-sequence-data-to-a-reference-genome-and-associated-steps.html
## Filtering
https://speciationgenomics.github.io/filtering_vcfs/