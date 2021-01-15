################################
## Author     : Can Kocagil  ###
################################

################################
####  --- Dependencies ---  ####
#### 1) ggplot2             ####
#### 2) dplyr               ####
#### 3) multest             ####
#### 4) cluster             ####
#### 5) factoextra          ####
#### 6) pheatmap            ####
#### 7) RColorBrewer        ####
#### 7) tidyverse           ####
#### 7) cowplot             ####
################################


#install.packages("pheatmap")
#install.packages("RColorBrewer")
#BiocManager::install("cluster")
#BiocManager::install("factoextra")

# Load necessary packages :
library("ggplot2")
library("multtest")
library("cluster")    
library("factoextra") 
library("pheatmap") 
library("RColorBrewer")
library("tidyverse")
library("dplyr")
library("cowplot")

# Q: What is most differentially expressed genes means?

# A: A gene is declared differentially expressed if an observed difference
# or change in read counts or expression levels between two experimental
# conditions is statistically significant.
# To identify differentially 
# expressed genes between two conditions, it is important to find statistical
# distributional property of the data to approximate the nature of differential genes.

# Then, first step is to to find a raw p-value, we can follows lot's of different ways.
# I am going to start by applying t-test# to all golub data then, calculate manually their p-values.
# Then, I select the most
# 100 differential genes among findings. After that, I am going to select randomly
# selected samples and apply hierarchical clustering.

# Apply statistical analysis technique test statistics including
# These functions provide a convenient way to compute test statistics,
# e.g., two-sample Welch t-statistics, Wilcoxon statistics,
# F-statistics, paired t-statistics, block F-statistics,
# for each row of a data frame
test_statistics <- mt.teststat(golub,classes)

# Looking test statistics
str(test_statistics)

# Note that lower the p value it is less likely the difference is due to chance 
# combine test statistics in teststat with raw p values
# calculate raw p-values
raw_p.value <- 2 * (1 - pnorm(abs(test_statistics)))

# Let's transfer the p-values more compact form:
df_p.value <- as.data.frame(raw_p.value)

# Look at the head of the p-values :
df_p.value %>% head()

# Note that the actual p-values representing confidence for differential
# expression are raw values. If they were to be corrected for multiple
# hypothesis testing (since there are lots of t-tests), they'd be much higher.

# Then, we need to find most differentially expressed 100 gene's p-values:
most_diff_genes_p.values <- df_p.value %>% arrange(raw_p.value)

# Looking the genes raw p-values:
most_diff_genes_p.values %>% head()

# Then, let's add this column as a new column for golub data:
new_golub <- cbind(golub,df_p.value$raw_p.value)

# Finally, we can find the most 100 differential gene by:
most_diff_gene_indices <- apply(t(new_golub[,39]), 1, order)[1:100]

# let's look the size(100):
dim_diff_exp <- length(most_diff_gene_indices);dim_diff_exp


# Taking random indexes of size 10
random_index_size_1 = 10
random_indexes_new_1 = sample(1:dim_diff_exp,random_index_size_1)

# Finally, indexing to get random gene's from the most differential genes:
feed_index_1 = most_diff_gene_indices[random_indexes_new_1]


#hierarchical clustering and plottings:
dist_mat = dist(t(golub[feed_index_1,]))

# Plottings:

p8 <- heatmap(t(golub[random_indexes_new_1,]), scale = "row")
p9 <- pheatmap(t(golub[random_indexes_new_1,]), color=brewer.pal(9,"Blues"))

# Taking random indexes of size 50
random_index_size_2 = 50
random_indexes_new_2 = sample(1:dim_diff_exp,random_index_size_2)

# Finally, indexing to get random gene's from the most differential genes:
feed_index_2 = most_diff_gene_indices[random_indexes_new_2]

#hierarchical clustering and plottings:
dist_mat_2 = dist(t(golub[feed_index_2,]))


p11 <- heatmap(t(golub[random_indexes_new_1,]), scale = "row")
p12 <- pheatmap(t(golub[random_indexes_new_1,]), color=brewer.pal(9,"Blues"))

par(mfrow=c(2,1))
plot(hclust(dist_mat,method="complete"),
     label=golub.cl,
     main = 'Cluster dendograms with 10 and 50 samples that are chosen randomly in most differentiable 100 genes',
     xlab = 'With 10 samples')
plot(hclust(dist_mat_2,method="complete"),label=golub.cl,xlab = 'With 50 samples')



# Comments:

# We see that part (ii), i.e., randomly selected 50 genes in most differentiable ones
# performed better in clustering because when we look at the cluster dendograms, we see that
# all 0's is in same cluster and all 1's is in same cluster that means it performed better
# than the first one since this is not the case in first part.


