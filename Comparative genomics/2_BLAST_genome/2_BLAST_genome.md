# BLAST a genome

## Learning objectives
* Conduct a BLAST search on a genome assembly using the webservice and the commandline interface
* Understand .bed format
* Use file-format documentation
* Extract sequences from a genome using coordinates.


# Tutorial

Now, if you have a genome that is totally unannotated, or you suspect there are unannotated sequences in there you can look for them using BLAST of just that genome.

Now you can use your translated sequences to conduct a tblastn search of the of interest. 
NCBI offers a commandline version of the of BLAST software, which enables novel blast databases and incorporation into pipelines. 
You can search your sequences against a reference genome using the commandline, or using the web function (as described in the lecture)

The `makeblastdb` takes forever, there is no need to do this component today. You can follow the webservice protocol, or you can use the second line of the below code, noting that the genome and its blastdb is stored in `~/Bioinformatics_Workshop/gardern_warbler_genome/`
```bash
nohup makeblastdb -in GCA_014839755.1_bSylBor1.pri_genomic.fna -parse_seqids -hash_index -dbtype nucl &
nohup tblastn -evalue 1e-20 -query  -db GCA_014839755.1_bSylBor1.pri_genomic -out SYBO_results -outfmt 4 -num_alignments 200000 -num_descriptions 200000 -num_threads 20
```
Note that we are only taking results with e-values below 1e-20. 

I also additionally decided to filter results based on a minimum exon length derived from the Zonotrichia albicollis gene model for OPN3 (minimum length 210 nt). 

There are multiple ways to extract out the regions of interest. You can download the individual fasta files from the web results OR you can parse your results file.

If you decide to do this from the text results then you need to:

1. Read the results text file some how. e.g. visually or in R. 
Note that the standard tblastn results headers are:
```queryacc.ver,subjectacc.ver,%identity,alignmentlength,mismatches,gapopens,q.start,q.end,s.start,s.end,evalue,bitscore,%positives,queryframes,subjectframes```

2. Find the coordinates of your best hits, and exclude overlapping results. 
I used the R code in the file in this github directory called `parsing_BLAST_OPN3.R` to extract the results and make a .bed file. This file is derived from the [Otter Genome Project](https://github.com/LohmuellerLab/OtterGenomeProject/tree/master/MolecularEvolution/olfactionAnalyses/FunctionalOlfactoryReceptorGene_Analysis/ferret_ORGs/mfur/step_1_c_CleanedByEvalue_ABScript)
This can be run remotely in Unix like so:

```bash
nohup R CMD BATCH parsing_BLAST_OPN3.R
```

For OPN3 results, the output was very simple and a bed file could be made manually by viewing s.start and s.end on the results and making a file that meets [.bed specifications](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)


3. Use the coordinates from the bed file to extract your sequences into a fasta file. 

```bash
REFERENCE=GCA_014839755.1_bSylBor1.pri_genomic.fna
spp="SYBO"
gene="OPN3"
bedfile=${spp}_${gene}_Hits.bed
bedtools getfasta -s -name -bed $bedfile -fi $REFERENCE -fo ${spp}_${gene}_sequences.fasta
```
Now, make sure your multi-fasta file is in the right order (largest numbers per chrom/scaff first if on `-` strand)  make it into a single entry and single line fasta file (text editor is ok). 



