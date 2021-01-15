################################
## Author     : Can Kocagil  ###
################################


X <- 'TTCATA'
Y <- 'TGCTCGTA'

seq.x <- unlist(strsplit(X, ''))
seq.y <- unlist(strsplit(Y, ''))

seq.x <- c(0,seq.x)
seq.y <- c(0,seq.y)
seq.x
seq.y

match <- 5
mismatch <- -2
indel <- -6

## initial the score matrix
score <- matrix(NA, length(seq.x), length(seq.y))
score
score[,1] <- sapply(1:length(seq.x)-1, function(x) x * indel)
score[1,] <- sapply(1:length(seq.y)-1, function(x) x * indel)
score
?sapply
showMethods("sapply")
## The dynamic programming, global alignment recursion
for (i in 2:length(seq.x)) {
  for (j in 2:length(seq.y)){
    # seq.x[i] , seq.y[j] are aligned
    if ( seq.x[i] == seq.y[j]) {
      score[i,j] <- score[i-1, j-1] + match
    } else {
      score[i,j] <- score[i-1, j-1] + mismatch
    }
    # seq.x[i] aligned to -
    sc <- score[i-1,j] + indel
    if (sc > score[i,j])
      score[i,j] = sc
    # seq.y[j] aligned to -
    sc <- score[i,j-1] + indel
    if (sc > score[i,j])
      score[i,j] = sc
  }
}
score
## Traceback
i <- length(seq.x)
j <- length(seq.y)
ax <- character()
ax
ay <- character()
ay
while (i > 1 && j >1){
  ## case 1: best was seq.x[i] aligned to seq.y[j]
  sc <- score[i-1,j-1]
  if (seq.x[i] == seq.y[j]) {
    sc <- sc + match
  } else {
    sc <- sc + mismatch
  }
  if (sc == score[i,j]) {
    ax <- c(seq.x[i], ax)
    ay <- c(seq.y[j], ay)
    i <- i -1
    j <- j-1
    next
  }
  ## case 2: best was seq.x[i] aligned to -
  if ((score[i-1,j] + indel) == score[i,j]) {
    ax <- c(seq.x[i], ax)
    ay <- c("-", ay)
    i <- i-1
    next
  }
  ## case 3: best was seq.y[j] aligned to -
  if ((score[i,j-1] + indel) == score[i,j]) {
    ax <- c("-", ax)
    ay <- c(seq.y[j], ay)
    j <- j-1
    next
  }
}

cat ("Sequence X: ", X,"\n")
cat ("Sequence Y: ", Y,"\n")
cat ("Scoring system: ", match, " for match; ", mismatch, " for mismatch; ", indel, " for gap", "\n\n")

cat ("Dynamic programming matrix:\n")
print (score)
?cat
cat ("\nAlignment:\n")
cat (paste(ax, collapse=''), "\n")
cat (paste(ay, collapse=''),"\n\n")
cat ("Optimum alignment score: ", score[length(score)],"\n")




######################################
#########  QUESTÝON-I  ###############
######################################

## Q1: Do you think the code is efficient enough? Why or Why not? 
## If not what do you think can be done explain your idea.

## To compare the results of the above code with Biostrings package, here is the process:

## Required R package
## It has substitution matrices,
## pairwise alignment tools and so on :
library(Biostrings)

##  Necessary Calculation for Scoring 
##  (total number matches)*match score
## -(total number of mismatches)*mismatch score-(total number of gapopen)
## -(total number gapextend)

## The Weight of Match :
match <- 5 

## The Weight of Mismatch :
mismatch <- -2

## Gap Opening Penalty :
gapOpen <- 6 

## Gap Extension Penalty :
gapExtend <- 6

## Creating Substitution Matrix :
myScoringMat <- nucleotideSubstitutionMatrix(match = match, 
                                             mismatch = mismatch, 
                                             baseOnly = FALSE) 

## Main algorithm for pairwise alignment :

myAlignment <- pairwiseAlignment(X, Y, 
                                 substitutionMatrix = myScoringMat, 
                                 gapOpening = gapOpen, 
                                 gapExtension = gapExtend, 
                                 type="global", 
                                 scoreOnly = FALSE)


myAlignment


## So, two results are basically same, it implies that the dynamic programming code
## is works well, but I try to explain the Space and Run time complexity because they determine
## the efficiency of the code

## We see that there is nested loop that (generally) implies that O(n^2) runtime
## complexity, it is may be okey for this solutions. 

######################################
#########  QUESTÝON-II  ##############
######################################

## Q2: If you were to do "local alignment" what step would you
## change/add please explain and give it a try.
## If I were to do 'Local Alignment', I change the one step with
## respect to Global Alignment here is the explanation :

## The main algorithm for global alignment is described here:

##############################
## --- Global Alignment --- ##
##############################
##  Initialization:         ##
##	Top right = 0           ##	
##  Update Rule:            ##
##  A(i,j) = max{           ##
##	A(i-1,j) -gap           ##
##	A(i,j-1) -gap           ##
##	A(i-1,j-1)              ##
##  }                       ##
##  Termination:            ##
##	 Bottom Right           ##
##############################

## The main algorithm for Local alignment is described here:

##############################
## --- Local Alignment ---  ##
##############################
##  Initialization:         ##
##	Top right = 0           ##	
##  Update Rule:            ##
##  A(i,j) = max{           ##
##	A(i-1,j) -gap           ##
##	A(i,j-1)-gap            ##
##	A(i-1,j-1)              ##
##      0                   ##
##  }                       ##
##  Termination:            ##
##	 Anywhere               ##
##############################


## Therefore, according to the update rule, I should add one conditions that
## 'Any number in score matrix cannot be negative' . The rest is same.
