# Sequence extraction and alignment for comparative genomics

# Schedule

* Part 1: genome browsers and gene extraction (Thursday June 10th 1pm EDT)
* Part 2: Multiple alignment (Friday June 11th 2pm EDT) 
* Part 3: Debrief and present your work (Tuesday June 15th 11am EDT)

# Pre-viewing / Revision

* Transcription and Translation and coding vs template strand
	* MIT OpenCourseWare 9:34 video (https://www.youtube.com/watch?v=x_vlxGFrZLY)
	* Nikolays Genetics Lessons 11:01 video (https://www.youtube.com/watch?v=gAm1ASjAMf8)
* Sequence Alignment
	* How Local Alignment Works - Smith-Waterman Algorithm (2020) 7:55 (https://www.youtube.com/watch?v=lu9ScxSejSE)
	* How BLAST works - a series (2018) each video is short (<10 min). at minimum watch 1,2,6,7,9,10 (https://www.youtube.com/watch?v=8A-msg23u0w) 
	* NCBI BLAST FAQs (2020) 3:39 (https://www.youtube.com/watch?v=dzRq-5BrGD4)
	* How Global Alignment Works (2018) Doborah Thurtle Schmidt 14:27 (https://www.youtube.com/watch?v=LhpGz5--isw)
* Phylogenetic tree thinking (2015) Matthew Taylor 12:42 https://www.youtube.com/watch?v=5AElW2omY6k

# Discussion Prompts (Tuesday)

* Did you have any trouble finding sequences for your species? Were they all complete sequences?
* Present your original alignment, and describe what sequence alignment method you used.
* Present your final alignment and describe your quality control procedure.
* What are the next steps from your alignments today in terms of your research specifically?


# Required Software
All commandline software required for this tutorial is installed on the ngsclass login of the class server. 
Please install one of the following GUIs for alignment viewing. 
* MEGA (free, all platforms https://www.megasoftware.net/) 


**If** you want to run this on your personal computer, the necessary Unix commandline tools are
* [EMBOSS](http://emboss.open-bio.org/) `apt-get install emboss`
* [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html)
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [samtools](http://www.htslib.org/)

The necessary R packages are:
```R
install.packages("BiocManager")
BiocManager::install("GenomicRanges")
```

# Useful Resources & Further Reading
* Homology types: https://useast.ensembl.org/info/genome/compara/homology_types.html
* NLM (2016) Webinar: Retrieving Exon and coding sequences 28 minutes. (https://www.youtube.com/watch?v=Xac7anudTD0)
	This webinar includes how to use NCBI commandline tools to do Entrez searches NOT covered in this tutorial.
* Selecting things in the genome data viewer (https://www.youtube.com/watch?v=EVNABwOCt_s)
* Nichio et al (2017) New Tools in Orthology Analysis: A Brief Review of Promising Perspectives. Frontiers in Genetics 8: 165
* Vallender (2009) Bioinformatic approaches to identifying orthologs and assessing evolutionary relationships. Methods 49: 50-55
