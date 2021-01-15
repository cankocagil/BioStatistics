################################
## Author     : Can Kocagil  ###
################################

################################
####  --- Dependencies ---  ####
#### 1) seginr              ####
#### 2) msa                 ####
#### 3) Biostrings          ####
#### 3) ape                 ####
################################

#--------------------------------#
##################################
##################################
####    PLEASE READ !!        ####
# IF THE FUNCTÝON DOES NOT WORK  # 
# DUE TO CONNECTÝNG PROBLEM LÝKE #
# "Unable to get any answer from # 
# socket after  10  trials"      #
# RUN THE CODE THAT STARTS LÝNE  #
# 240 ONE BY ONE, ÝF SOCKET      #
# PROBLEM, CHOOSE DATABASE AGAÝN #
##################################
##################################
#--------------------------------#


##################################################################################################################
##################################################################################################################
#########################                              To Do's :                         #########################


# 0) Uploading necessary R packages to proceed
# 1) Choosing a database that holds specific information about the gene and sequence with assigned access number 
# 2) Determining the Substitution Matrix (i.e., BLOSUM50,BLOSUM45,BLOSUM62,BLOSUM100)
# 3) Query to retrieve the necessary information
# 4) Getting sequences from the determined query
# 5) Repeat 3 and 4 until we have enough species and their sequences
# 6) Converting sequences to .fasta format
# 7) Multiple Sequence Alignment
# 8) Calculating the matrix of pairwise distances from the aligned proteins 
# 9) Applying the Neighbor Joining Algorithm
# 10) Plotting the phylogenetic tree

##################################################################################################################
##################################################################################################################

# ------- Step 0 ------- #
##########################

# Uploading necessary packages :
#BiocManager::install("seqinr")
#BiocManager::install("msa")
#BiocManager::install("Biostrings")
#BiocManager::install("ape")

library("seqinr")
library("msa")
library("Biostrings")
library("ape")

## ---- Step 0 Finish ---- ##


phyloGeneticTree <- function(prot_database ,species_list,n_species,
                             access_numbers,fasta_file,substMat,gapOpening,gapExtension)
{
  
  myDataBase = prot_database 
  # Connecting database :
  choosebank(myDataBase)
  
  for (i in 1:n_species) {
    
    if (i ==1){
      # Query to retrieve information for zebrafish 
      myQuery1 = query(species_list[i],access_numbers[i])
      zebrafishSeq = getSequence(myQuery1$req[[1]],as.string=T)
      
    }else if (i==2){
      # Query to retrieve information for human
      myQuery2 = query(species_list[i],access_numbers[i])
      humanSeq = getSequence(myQuery2$req[[1]],as.string=T)
      
    }else if (i==3){
      # Query to retrieve information for mouseearcress
      myQuery3 = query(species_list[i],access_numbers[i])
      mouseearcressSeq = getSequence(myQuery3$req[[1]],as.string=T)
      
    }else if (i==4){
      
      # Query to retrieve information
      myQuery3 = query(species_list[i],access_numbers[i])
      mouseearcressSeq = getSequence(myQuery3$req[[1]],as.string=T)
      
    }else if (i==5){
      # Query to retrieve information
      myQuery4 = query(species_list[i],access_numbers[i])
      bovineSeq = getSequence(myQuery4$req[[1]],as.string=T)
      
    }else if (i==6){
      
      # Query to retrieve information
      myQuery5 = query(species_list[i],access_numbers[i])
      mouseSeq = getSequence(myQuery5$req[[1]], as.string=T)
      
    }else if (i==7){
      
      myQuery6 = query(species_list[i],access_numbers[i])
      frogSeq = getSequence(myQuery6$req[[1]],as.string=T)
      
    }else if (i==8){
      choosebank(myDataBase)
      # Query to retrieve information
      myQuery7 = query(species_list[i],access_numbers[i])
      ratSeq = getSequence(myQuery7$req[[1]],as.string=T)
      
    }else if (i==9){
      choosebank(myDataBase)
      # Query to retrieve information
      myQuery8 = query(species_list[i],access_numbers[i])
      stSeq = getSequence(myQuery8$req[[1]], as.string=T)
      
    }else if (i==10){
      
      # Query to retrieve information
      myQuery9 = query(species_list[i],access_numbers[i])
      eleSeq = getSequence(myQuery9$req[[1]],as.string=T)
      
    }else if (i==11){
      # Query to retrieve information
      myQuery10 = query(species_list[i],access_numbers[i])
      ddSeq = getSequence(myQuery10$req[[1]],as.string=T)
      
    }else if (i==12){
      
      # Query to retrieve information
      myQuery11 = query(species_list[i],access_numbers[i])
      paSeq = getSequence(myQuery11$req[[1]],as.string=T)
      
    }else if (i==13){
      
      # Query to retrieve information
      myQuery12 = query(species_list[i],access_numbers[i])
      chicSeq = getSequence(myQuery12$req[[1]],as.string=T)
      
    }
    
    
  }
  # First converting list :
  list_seq = as.list(c(zebrafishSeq,humanSeq,mouseearcressSeq,
                       bovineSeq,mouseSeq,frogSeq,ratSeq,
                       stSeq,eleSeq,ddSeq,paSeq,chicSeq,dmSeq))
  
  # Converting our sequences in .fasta format :
  write.fasta(sequences = list_seq,names = vect_names,
                      as.string = F,file.out = fasta_file)
  
  # Retrieving sequences in .fasta format :
  seq = readAAStringSet(filepath = fasta_file)
  
  
  # Multiple Sequence Alignment
  myProtAlign <- msa(seq, method = "ClustalW")
  
  print(myProtAlign, show="complete")
  
  myProtAlign2 <- msaConvert(myProtAlign,type="seqinr::alignment")
  
  dist <- dist.alignment(myProtAlign2, "similarity")
  as.matrix(dist)
  
  # Applying the Neighbor Joining Algorithm
  phyloGenTree <- nj(dist)
  
  # Plotting the phylogenetic tree
  layout(matrix(data = 1, nrow = 1, ncol = 1), width=c(1.4))
  plot(phyloGenTree,main = 'Phylogenetic Tree',sub = 'Orthologs of Zebrafish')


}


prot_database = "swissprot" 


# Names of leaves of the tree :
species_list = c('Zebra fish','Human','Mouse arc cress',
                 'bovine','Mus musculus','Xenopus tropicalis',
                 'Rattus norvegicus','Saccharomyces cerevisiae',
                 'Caenorhabditis elegans','Dictyostelium discoideum',
                 'Pongo abelii','Gallus gallus','Drosophila melanogaste')



n_species = length(species_list)

access_numbers = c("AC=Q7ZV68","AC=Q9UBQ0","AC=Q9STT2","AC=Q3T0M0",
                   "AC=Q9QZ88","AC=Q6DEU3","AC=B2RZ78","AC=P38759",
                   "AC=Q9XVX5","AC=Q54IF7","AC=Q5R9Z1","AC=Q5ZIL2","AC=Q9VPX5") 

fasta_file = "myMidtermSeq"

phyloGeneticTree(prot_database,species_list,n_species,access_numbers,fasta_file)





#############################################################################
#############################################################################
#############################################################################
##############     RUN BELOW IF SOCKET PROBLEM ONE BY ONE      ##############
#############################################################################
#############################################################################
#############################################################################

# ------- Step 1 ------- #
##########################

# Choosing our database :

# - UniProt/Swiss-Prot
list_databases = seqinr::choosebank();list_databases

# Choosing UniProt/Swiss-Prot database :
myDataBase = "swissprot" ;myDataBase

# Connecting database :
choosebank(myDataBase)

## ---- Step 1 Finish ---- ##

# ------- Step 2 ------- #
##########################

# Determining Substitution matrix :
# I am going to use default substitution matrix in the following parts

## ---- Step 2 Finish ---- ##

# ------- Step 3,4 and 5 ------- #
##################################

# Query to retrieve the necessary information
# Getting sequences from the determined query

## ----------------- Zebrafish ----------------- ##
# ---- https://www.uniprot.org/uniprot/Q7ZV68 ----#
###################################################
###################################################

# Accesing database number :
access_number1 = "AC=Q7ZV68"

# Organism :
organism1 = "Danio rerio"

# Query to retrieve information
myQuery1 = seqinr::query(listname = organism1,
                         query = access_number1)

zebrafishSeq = seqinr::getSequence(myQuery1$req[[1]],
                                   as.string=T)
zebrafishSeq


###################################################
###################################################

## ----------------- Ortholog1 ----------------- ##
# ---- https://www.uniprot.org/uniprot/Q9UBQ0 ----#
###################################################
###################################################

# Accesing database number :
access_number2 = "AC=Q9UBQ0"

# Organism :
organism2 = "Homo sapiens"

# Query to retrieve information
myQuery2 = seqinr::query(listname = organism2,
                         query = access_number2)

humanSeq = seqinr::getSequence(myQuery2$req[[1]],
                                   as.string=T)
humanSeq

###################################################
###################################################

## ----------------- Ortholog2 ----------------- ##
# ---- https://www.uniprot.org/uniprot/Q9STT2 ----#
###################################################
###################################################

# Accesing database number :
access_number3 = "AC=Q9STT2"

# Organism :
organism3 = "Arabidopsis thaliana"

# Query to retrieve information
myQuery3 = seqinr::query(listname = organism3,
                         query = access_number3)

mouseearcressSeq = seqinr::getSequence(myQuery3$req[[1]],
                                   as.string=T)
mouseearcressSeq

###################################################
###################################################

## ----------------- Ortholog3 ----------------- ##
# ---- https://www.uniprot.org/uniprot/Q3T0M0 ----#
###################################################
###################################################

# Accesing database number :
access_number4 = "AC=Q3T0M0"

# Organism :
organism4 = "Bos taurus"

# Query to retrieve information
myQuery4 = seqinr::query(listname = organism4,
                         query = access_number4)

bovineSeq = seqinr::getSequence(myQuery4$req[[1]],
                                   as.string=T)
bovineSeq

###################################################
###################################################
## ----------------- Ortholog4 ----------------- ##
# ---- https://www.uniprot.org/uniprot/Q9UTI5 ----#
###################################################
###################################################

# Accesing database number :
access_number5 = "AC=Q9QZ88"

# Organism :
organism5 = "Mus musculus"

# Query to retrieve information
myQuery5 = seqinr::query(listname = organism5,
                         query = access_number5)


mouseSeq = seqinr::getSequence(myQuery5$req[[1]],
                                   as.string=TRUE)
mouseSeq

###################################################
###################################################

## ----------------- Ortholog5 ----------------- ##
# ---- https://www.uniprot.org/uniprot/Q6DEU3 ----#
###################################################
###################################################

# Accesing database number :
access_number6 = "AC=Q6DEU3"

# Organism :
organism6 = "Xenopus tropicalis"

# Query to retrieve information
myQuery6 = seqinr::query(listname = organism6,
                         query = access_number6)


frogSeq = seqinr::getSequence(myQuery6$req[[1]],
                                 as.string=TRUE)
frogSeq

###################################################
###################################################

###################################################
###################################################

## ----------------- Ortholog6 ----------------- ##
# ---- https://www.uniprot.org/uniprot/B2RZ78 ----#
###################################################
###################################################
choosebank(myDataBase)
# Accesing database number :
access_number7 = "AC=B2RZ78"

# Organism :
organism7 = "Rattus norvegicus"

# Query to retrieve information
myQuery7 = seqinr::query(listname = organism7,
                         query = access_number7)


ratSeq = seqinr::getSequence(myQuery7$req[[1]],
                             as.string=TRUE)
ratSeq

###################################################
###################################################

## ----------------- Ortholog7 ----------------- ##
# ---- https://www.uniprot.org/uniprot/P38759 ----#
###################################################
###################################################

# Accesing database number :
access_number8 = "AC=P38759"

# Organism :
organism8 = "Saccharomyces cerevisiae"

# Query to retrieve information
myQuery8 = seqinr::query(listname = organism8,
                         query = access_number8)


stSeq = seqinr::getSequence(myQuery8$req[[1]],
                             as.string=TRUE)
stSeq

###################################################
###################################################

###################################################
###################################################

## ----------------- Ortholog8 ----------------- ##
# ---- https://www.uniprot.org/uniprot/P38759 ----#
###################################################
###################################################

# Accesing database number :
access_number9 = "AC=Q9XVX5"

# Organism :
organism9 = "Caenorhabditis elegans"

# Query to retrieve information
myQuery9 = seqinr::query(listname = organism9,
                         query = access_number9)


eleSeq = seqinr::getSequence(myQuery9$req[[1]],
                            as.string=TRUE)
eleSeq

###################################################
###################################################


## ----------------- Ortholog9 ----------------- ##
# ---- https://www.uniprot.org/uniprot/P38759 ----#
###################################################
###################################################

# Accesing database number :
access_number10 = "AC=Q54IF7"

# Organism :
organism10 = "Dictyostelium discoideum"

# Query to retrieve information
myQuery10 = seqinr::query(listname = organism10,
                         query = access_number10)


ddSeq = seqinr::getSequence(myQuery10$req[[1]],
                             as.string=TRUE)
ddSeq

###################################################
###################################################


## ----------------- Ortholog10 ----------------- ##
# ---- https://www.uniprot.org/uniprot/P38759 ----#
###################################################
###################################################

# Accesing database number :
access_number11 = "AC=Q5R9Z1"

# Organism :
organism11 = "Pongo abelii"

# Query to retrieve information
myQuery11 = seqinr::query(listname = organism11,
                          query = access_number11)


paSeq = seqinr::getSequence(myQuery11$req[[1]],
                            as.string=TRUE)
paSeq

###################################################
###################################################

## ----------------- Ortholog11 ----------------- ##
# ---- https://www.uniprot.org/uniprot/P38759 ----#
###################################################
###################################################

# Accesing database number :
access_number12 = "AC=Q5ZIL2"

# Organism :
organism12 = "Gallus gallus"

# Query to retrieve information
myQuery12 = seqinr::query(listname = organism12,
                          query = access_number12)


chicSeq = seqinr::getSequence(myQuery12$req[[1]],
                            as.string=TRUE)
chicSeq

###################################################
###################################################

## ----------------- Ortholog12 ----------------- ##
# ---- https://www.uniprot.org/uniprot/P38759 ----#
###################################################
###################################################

# Accesing database number :
access_number13 = "AC=Q9VPX5"
# Organism :
organism13 = "Drosophila melanogaster"

# Query to retrieve information
myQuery13 = seqinr::query(listname = organism13,
                          query = access_number13)


dmSeq = seqinr::getSequence(myQuery13$req[[1]],
                              as.string=TRUE)
dmSeq

###################################################
###################################################


## ---- Step 3,4,5 Finish ---- ##

# ------- Step 6 ------- #
##########################

# Converting sequences to .fasta format :

# Names of leaves of the tree :
vect_names = c('Zebra fish','Human','Mouse arc cress',
               'bovine','Mus musculus','Xenopus tropicalis',
               'Rattus norvegicus','Saccharomyces cerevisiae',
               'Caenorhabditis elegans','Dictyostelium discoideum',
               'Pongo abelii','Gallus gallus','Drosophila melanogaste')

# First converting list :
list_seq = as.list(c(zebrafishSeq,humanSeq,mouseearcressSeq,
                     bovineSeq,mouseSeq,frogSeq,ratSeq,
                     stSeq,eleSeq,ddSeq,paSeq,chicSeq,dmSeq))

# Converting our sequences in .fasta format :
seqinr::write.fasta(sequences = list_seq,
                    names = vect_names,
                    as.string = F,
                    file.out = "myMidtermSeq")

# Retrieving sequences in .fasta format :
seq = Biostrings::readAAStringSet(filepath = "myMidtermSeq")

## ---- Step 6 Finish ---- ##

# ------- Step 7 ------- #
##########################

# Multiple Sequence Alignment
myProtAlign <- msa::msa(inputSeqs = seq)

print(myProtAlign, show="complete")

myProtAlign2 <- msaConvert(myProtAlign,
                               type="seqinr::alignment")


## ---- Step 7 Finish ---- ##

# ------- Step 8 ------- #
##########################

# Calculating the matrix of pairwise distances from the aligned proteins 

dist <- dist.alignment(myProtAlign2, "similarity")
as.matrix(dist)

?dist.alignment
## ---- Step 8 Finish ---- ## 

# ------- Step 9 ------- #
##########################

# Applying the Neighbor Joining Algorithm
phyloGenTree <- ape::nj(dist)

## ---- Step 9 Finish ---- ##

# ------- Step 10 ------- #
##########################

# Plotting the phylogenetic tree
layout(matrix(data = 1, nrow = 1, ncol = 1), width=c(1.4))
plot(phyloGenTree,main = 'Phylogenetic Tree')


## ---- Step 10  Finish ---- ##