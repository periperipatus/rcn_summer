# Genome browsers and sequence extraction.

## Learning Objectives
* Understand gene models in genome browsers
* Extract sequences from NCBI
* Extract CDS from NCBI transcripts.
* Conduct BLAST searches and understand the results. 
* Translate sequences 
* Understand homology terms

## Preamble

In this section you will be extracting some genes from NCBI using Entrez and BLAST. For today, you will only extract sequences from the following species:

* _Manacus vitellinus_
* _Pipra filicauda_
* _Chiroxiphia lanceolata_
* _Corapipo altera_
* _Lepidothrix coronata_
* _Neopelma chrysocephalum_
* _Empidonax trailii_
* _Taeniopygia guttata_
* _Zonotrichia albicollis_
* _Ficedula albicollis_
* _Corvus monoduloides_
* _Calypte anna_
* _Aquila chrysaetos_
* _Aythya fuliga_
* _Gallus gallus_

Here are the list of genes you will work on.

* TAP1 - Elsie
* AR - Camilo
* CIITA - Emily
* TH (Tyrosine Hydroxylase) - Diego



## Basic Instructions

More detailed instructions are included in the interactive lecture component of this module.

1. Conduct a nucleotide Entrez search of your gene (e.g. OPN3)
2. Select the CDS sequence of your gene and copy that to a fasta file in your text editor.
3. Use that sequence to BLAST the nucleotide database, but then refine the organism field to Aves.
4. View your sequences in the multiple sequence alignment viewer to get an idea of how they align to your query.
Then, for each species above: 
5. Extract the __CDS__ for the longest transcript of that gene or the transcript variant that has the best coverage for your query sequence (as viewed in MSA viewer). Do this by clicking "CDS" on the record and then clicking FASTA in the bar that appears at the bottom corresponding to the highlighted sequence. 
6. Enter your sequences into a text file, and modify the header in a way that preserves the most important information and makes sense to you.
	e.g `>GAGA OPN3 TranscriptX1 XM_426139.6:114-1301 Gallus gallus`
7. Translate those sequences using MEGA, or the following command with EMBOSS' `transeq` function on the `ngsclass` server.

```bash
transeq <infile.fa> <outfile.fa>
```

# Pulling a sequence from a genome annotation

## Learning Objectives
* Understand the different files in the FTP directory in a genome assembly
* Search for specific genes in genome annotations

# Tutorial
Sometimes you will have access to genomes before they are accessioned into Genbank. 
In some cases, the gene sequences will not be available in the nucleotide database. But you can pull this information from the annotations. 

For example, as of 8th June 2021, the Amazon Umbrellabird has an annotated genome but it's not integrated into the databases. 
We are going to extract the example sequence (OPN3) from the genome. 

Download the annotation files to your directory on `ngsclass`. 
`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/396/775/GCA_013396775.1_ASM1339677v1/GCA_013396775.1_ASM1339677v1_cds_from_genomic.fna.gz`

Take a peak at the `_genomic.fna.gz` and `_genomic.gff.gz` files as well. 

We can search for our gene using

```bash
grep 'OPN3' GCA_013396775.1_ASM1339677v1_cds_from_genomic.fna
```
It's there! 
If you use `less -S` on the file, you can see that the sequence information is over multiple lines. 
There are multiple ways to extract the sequence below the header. An easy way is to convert the multiple lines of sequence to a single line using `FASTX-toolkit`

```bash
fasta_formatter -i GCA_013396775.1_ASM1339677v1_cds_from_genomic.fna -o GCA_013396775.1_ASM1339677v1_cds_from_genomic_1line.fna -w 0
grep -i -A 1 "Opn3" GCA_013396775.1_ASM1339677v1_cds_from_genomic_1line.fna > Cephalopterus_OPN3.fna
```
Now you can look at the file and copy and paste the sequence into your multi-fasta file. 
