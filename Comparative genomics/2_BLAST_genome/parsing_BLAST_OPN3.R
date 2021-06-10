setwd("C:/Users/perif/OneDrive - East Carolina University/BALALAB_SHARED/RCN Mentoring/Public/Comparative genomics/2_BLAST_genome/")
#this code is heavily based on the code from the Otter Genomes Project https://github.com/LohmuellerLab/OtterGenomeProject 
# specifically this code: https://github.com/LohmuellerLab/OtterGenomeProject/blob/master/MolecularEvolution/olfactionAnalyses/FunctionalOlfactoryReceptorGene_Analysis/ferret_ORGs/mfur/step_1_c_CleanedByEvalue_ABScript/Step_1_c_CleanByEvaluePrune_Annabel.mfur.R
blast_results<- read.csv("C1KXAEZK01R-Alignment-HitTable.csv")

#based on the length of OPN3 exons in Zontotrichia albicollis I will pick a minimum aa length based on that (minus some for variation in exon length)
blast_results<- blast_results[blast_results$alignmentlength>70,]

#filter out things with low e-value
e<-1e-20
blast_results<- blast_results[blast_results$evalue<e,]

# want to add a unique index to each hit so that I can grab things out later.
blast_results$uniqueHitLabel <- seq(1:dim(blast_results)[1])



#remove hits that are within 100bp of eachother (overlapping), and choosing the hit with the highest e-value. 

####INSTALL GENOMIC RANGES#####
#BiocManager::install("GenomicRanges") #See here if you don't have Bioconductor installed in R (https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
library(GenomicRanges)
#Turn the reading frame information into strandedness for GenomicRanges
#Look at the results and you'll see that the start of the result sequence is AFTER the end.
blast_results$strand<- ifelse(blast_results$subjectframes<0, "-", "+")
#IRanges also requires the start to be smaller than the end, so let's fix it so every start value is smaller than the end.
for(i in 1:nrow(blast_results)){
  sub<- blast_results[i,]
  if(sub$s.end<sub$s.start){
    end=sub$s.start
    start=sub$s.end
    sub$s.start=start
    sub$s.end=end
    blast_results[i,]<- sub
  }else{}
}



grobj <- GRanges(seqnames=blast_results$subjectacc.ver,ranges=IRanges(blast_results$s.start,end=blast_results$s.end),strand=blast_results$strand,bit.score=blast_results$bitscore,evalue=blast_results$evalue,query=blast_results$queryacc.ver,hit.length=blast_results$alignmentlength,uniqueHitLabel=blast_results$uniqueHitLabel)
overlaps <- findOverlaps(grobj,maxgap=100L,minoverlap=0L,drop.self=F,drop.redundant=F)

#this function keeps only those blast hits that DO NOT OVERLAP at all, and have the highest e-value of the options. 
cleanByEValue_prune <- function(grobj,overlaps) {
  queryNums <- unique(queryHits(overlaps))
  keep <- GRanges()
  while(length(queryNums) > 0){
    q <- queryNums[1]
    print(q)
    subs_q <- subjectHits(overlaps[queryHits(overlaps)==q]) # don't need to add the query number in because I have the hits to itself kept into overlaps now     
    grobj_subs <- grobj[subs_q]
    keepThisHit <- which.min(grobj_subs$evalue)
    keep <- c(keep,grobj_subs[keepThisHit])
    queryNums <- queryNums[!(queryNums %in% subs_q)]
    print(paste("dropped:",subs_q))
    print(paste("left to go: ", length(queryNums)))
  }
  return(unique(keep))
}

keep_prune_eval <- cleanByEValue_prune(grobj,overlaps)

# make a .bed file to intersect with your regions

pullOut <- blast_results[blast_results$uniqueHitLabel %in% keep_prune_eval$uniqueHitLabel,] # pull out the original entries 

pullOut$bed_HitStart0 <- pullOut$s.start - 1 # change the start site because bed is zero based (don't change end bc it is non inclusive)

# $species\_$chr\_$start\_$end\_$string\ *** NOTE THAT BLAST IS ONE BASED *** so start end are 1 based
pullOut$bedname <- paste("SYBO",pullOut$subjectacc.ver,pullOut$s.start,pullOut$s.end,pullOut$subjectframes,sep="_")
bed <- pullOut[,c("subjectacc.ver","bed_HitStart0","s.end","bedname","evalue","strand")]
bed<-bed[order(bed$bed_HitStart0,decreasing = TRUE),] #put in DECREASING order by start site because this is on the REVERSE strand. 

write.table(bed,"SYBO_OPN3_Hits.bed",row.names = F,col.names = F,quote=F,sep="\t")
