#FROM BIOINFORMATICS R COOKBOOK PAGES 69-72 
#NUCLEOTIDE SEQUENCE ALIGNMENT 
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('Biostrings')#install Biostrings, a package in R

#source("http://bioconductor.org/biocLite.R") #source is a function that brings
#this url into your workspace makes you access
#biocLite("Biostrings") #biocLite, a function helps you install Biostrings
#old version
library(Biostrings) #library, a function is required to call a function
#align two sequence nucleotide sequences of different sizes first define
#to find similarities, identities, mismatches, match, gaps, one sequence is long
#the other is short; to make align i need to insert gaps
sequence1 <- "GAATTCGGCTA"
class(sequence1)
sequence2 <- "GATTACCTA"
#calculate a score will be made up of Number of matches, number mismatches,
#gaps; each a value,
#gaps, two kind: open a gap, extend the gap that i opened, 4nt is more 
#dangerous than 3nt, if coding sequence extending in 3s might one aa, 
#frameshift, nonsense, opening gaps cost more than extending, negative
#SCORE=sum(matches, mismatches, gapopen, gapclose) gap penalty
#there are multiple ways to align two sequence that they are homologous
#a score that i will calculate, match, mismatch, gap
# scoring matrix needs to created by the function below, provide scores
#score means a value made up of (total number matches)*match score
#-(total number of mismatches)*mismatch score-(total number of gapopen)
#-(total number gapextend)
myScoringMat <- nucleotideSubstitutionMatrix(match = 1, 
                                             mismatch = -1, 
                                             baseOnly = FALSE) 
myScoringMat
?nucleotideSubstitutionMatrix
gapOpen <- 4 #gap penalty in my score function it will be substracted
gapExtend <- 1 
myAlignment <- pairwiseAlignment(sequence1, sequence2, 
                                 substitutionMatrix = myScoringMat, 
                                 gapOpening = gapOpen, 
                                 gapExtension = gapExtend, 
                                 type="global", 
                                 scoreOnly = FALSE) 
myAlignment 

myAlignment <- pairwiseAlignment(sequence1, sequence2, 
                                 substitutionMatrix = myScoringMat, 
                                 gapOpening = gapOpen, 
                                 gapExtension = gapExtend, 
                                 type="local", 
                                 scoreOnly = FALSE) 
myAlignment 

#PROTEIN SEQUENCE ALIGNMENT 
data(package="Biostrings") 
data(BLOSUM62) 
subMat <- "BLOSUM62" 
subMat 

sequence1 <- "PAWHEAE" 
sequence2 <- "HEAGAWGHE" 
myAlignProt <- pairwiseAlignment(sequence1, sequence2, 
                                 substitutionMatrix = subMat, 
                                 gapOpening = gapOpen, 
                                 gapExtension = gapExtend, 
                                 type="global", 
                                 scoreOnly = FALSE) 
myAlignProt 

myAlignProt <- pairwiseAlignment(sequence1, sequence2, 
                                 substitutionMatrix = subMat, 
                                 gapOpening = gapOpen, 
                                 gapExtension = gapExtend, 
                                 type="local", scoreOnly = FALSE) 
myAlignProt 
