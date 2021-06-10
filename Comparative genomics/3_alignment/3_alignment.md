# Aligning Sequences

### Learning Objectives

* Understand and implement codon-aware alignment
* Awareness of multiple programs and sequence alignment algorithms
* Quality control sequence alignments considering reading frame and sequencing or assembly errors. 

# Tutorial
 
In order to conduct analyses of dN/dS and infer selection, we need to ensure that our codon structure is preserved when we conduct a multiple alignment of our genes.
This is the online component of a mostly lecture based tutorial. Please refer to lecture recording for more information. 

There are multiple tools that can conduct codon-aware alignment, some of these include:
* MEGA (can also edit alignments)
* [translatorX](http://translatorx.co.uk/)
* PRANK
* [PAL2NAL](http://www.bork.embl.de/pal2nal/)
* [Guidance2](http://guidance.tau.ac.il/) program provides codon-aware alignment, as well as automated quality scores for sequences, sites and individual residues.

Some of these programs are unique algorithms for alignment (e.g. PRANK), others allow you to produce alignments using common alignment methods e.g. MAFFT, MUSCLE, ClustalW, PRANK, T-coffee.
We like MAFFT and MUSCLE. MAFFT tends to be better for sequences with lower homology, or sequences with long gaps, and can handle ~30k sequences. 
We will use MEGA today because it also allows alignment editing. You can also use MAFFT in a program like [GUIDANCE](http://guidance.tau.ac.il/) or [translatorX](http://translatorx.co.uk/) and import the .fasta file into MEGA (or another sequence editor) for viewing and editing. 

From today's tutorial you will:

1. Create a codon-aware alignment on your sequences downloaded from parts 1 and 2 of this tutorial.
2. Perform quality control on that alignment through visual inspection (also feel free to use Guidance2 if you want a more quantitative approach). 
	Make sure you take notes on the positions within the alignment that you are deleting!
3. Create another codon aware alignment containing the sequences you extracted from the unannotated genome. 

## Quality Control Checklist
 
* Ensure a shared start codon (where possible)
* Trim to the nearest codon on 3' end
* Check it's in the same reading frame as when you downloaded it! (shared start codon) 
* USER DISCRETION: Delete or mask sequences that are indels, poorly aligned or  highly divergent regions if they do not appear to be biologically real
* USER DISCRETION: It may be necessary to delete a species from an alignment if the sequence cannot be aligned. 


# Discussion Prompts (Tuesday)
When presenting your findings, consider these:

* Did you have any trouble finding sequences for your species? Were they all complete sequences?
* Present your original alignment, and describe what sequence alignment method you used.
* Present your final alignment and describe your quality control procedure.
* Describe the results from the BLAST experiment on the unannotated genome. How well did it align? What sort of trimming did you do?
* What are the next steps from your alignments today in terms of your research specifically?


# Other alignment editing software

I think the GUI for MEGA is slow, but it is for all the platforms... Here are some other alignment editing software options
* eBioX (Mac) - If you can get your hands on this it is SIMPLE but very nice to use
* BioEdit (PC) - Old, Ugly, but effective. 
* Geneious - Free Trial can do all sorts of editing and WAY more. Check if your institution has a site licesnse.
* JalView - I have never used this
