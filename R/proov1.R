## (1) set the input fasta file name. 

library(VirFinder)
inFaFile <- system.file("data", "contigs.fa", package="VirFinder")
inFaFile <- system.file("data", "assembly_2.fa")
inFaFile <- system.file("data", "assembly_3.fa")
inFaFile <- system.file("data", "assembly_4.fa")
inFaFile <- system.file("data", "crAssphage.fasta", package="VirFinder")
VF.pred(inFaFile)



## (2) prediction
predResult <- VF.pred("fasta/assembly_3.fa")
predResult <- VF.pred("fasta/assembly_5.fa")
predResult <- VF.pred(inFaFile)
predResult

## (A) # Koos Eukaryotide viirustega ennustamine:
## specify the directory of the new model file "VF.modEPV_k8.rda", and load the new model to the work space
modFile <- "VF.modEPV_k8.rda"
load(modFile)
## specify the fasta file containing contigs for prediction
#inFaFile <- "<path_to_the_input_fasta_file>/input.fasta"
#inFaFile <- system.file("data", "crAssphage.fasta", package="VirFinder")

## predict the contigs using the new model
#predResultUser <- VF.pred.user(inFaFile, modEPV)
predResultUser <- VF.pred.user("fasta/assembly_7.fa", modEPV)


#### (2.1) sort sequences by p-value in ascending order
library(tidyverse)
#predResult[order(predResult$pvalue),]
arrange(predResult, pvalue)
arrange(predResultUser, pvalue)

#### (2.2) estimate q-values (false discovery rates) based on p-values
#predResult$qvalue <- VF.qvalue(predResult$pvalue)
predResult <- mutate(predResult, qvalue = VF.qvalue(pvalue))
predResultUser <- mutate(predResultUser, qvalue = VF.qvalue(pvalue))
predResult
predResultUser

#### (2.3) sort sequences by q-value in ascending order
predResult[order(predResult$qvalue),]
predResultUser[order(predResultUser$qvalue),]
arrange(predResult, qvalue) %>% as_data_frame()
arrange(predResultUser, qvalue) %>% as_data_frame()

#10 esimest nime
#nende j2rjetused
#need uude fasta faili

# 10 first ids of interest
ids <- as.character(predResultUser$name[1:10])

# load subset fasta function
source("R/subset-fasta.R")

subset_fasta(ids, "fasta/assembly_7.fa", "output/assembly_7_top10.fa")

