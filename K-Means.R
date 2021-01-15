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

# Retrieving data :
data(golub)

# Take a look at the data :
head(golub)


# Take a look at class vectors:
classes = golub.cl;unique(classes)

# Leukemia (ALL and AML) is the most common form of childhood cancer. 
# It affects the tissues of the body which make the blood cells.


set.seed(218)

###################################################################################################
###################################################################################################
#########################                  Select 50 random integers             ##################
###################################################################################################


# Taking random numbers between 1 and the row number of the Golub data
first_size = 40
random_indexes_1 = sample(1:dim_golub[1],first_size)

# silhouette plot to find best number of cluster k
p1 <- fviz_nbclust(t(golub[random_indexes_1,]), kmeans, method = "silhouette")


###################################################################################################
###################################################################################################
#########################                  Select 400 random integers             #################
###################################################################################################


# Taking random numbers between 1 and the row number of the Golub data
second_size = 400
random_indexes_2 = sample(1:dim_golub[1],second_size)


# silhouette plot to find best number of cluster k
p2 <- fviz_nbclust(t(golub[random_indexes_2,]), kmeans, method = "silhouette")

###################################################################################################
###################################################################################################
#########################                  Select 1000 random integers             ################
###################################################################################################

# Taking random numbers between 1 and the row number of the Golub data
third_size = 1000
random_indexes_3 = sample(1:dim_golub[1],third_size)

# silhouette plot to find best number of cluster k
p3 <- fviz_nbclust(t(golub[random_indexes_3,]), kmeans, method = "silhouette")

# Plot in one figure:
plot_grid(p1,p2,p3)


# Therefore, when we  change the size of the input to the silhouette plot, the
# number of centers are changed. 

# To explain this, let's talk about the silhouette plot,
# Silhoette plot takes input and clustering algorithm, try to find best possiblle k value.
# is used for clustering, in our
# case, kmeans that is unsupervised learning algorith. k is 
# hyperparameter of the algorithm such that we are give possible k values and test it.
# if we happy with the result then we can stick with this k. But if we are not
# okey with this k, we try different k's to optimize the model.

# In our case,
# since we have 2 classes ALL and AML, we think that we should have 2 clusters,
# but it is not always the case due the possible variance of the distribution
# of the features that are characterizing classes. In other means, there may be such
# spread in the data that result in, in some cases, more than 2 cluster. When, we
# look to our results, we see that, in some cases we have 2 cluster, and in another
# cases we have not. It is not surprising thing.


# In our case,
# We see that when the number of examples increases, the number of clusters increases. The reason behind that 
# is increment the spread or variance of the samples so that their distances between itself
# is increased. But note that it may not be the case for all situation, i.e., these correlation
# may or may not be hold. But in our case, it holds.

# Repeate again:
# Taking random numbers between 1 and the row number of the Golub data
set.seed(0)

first_size = 40
random_indexes_1 = sample(1:dim_golub[1],first_size)
second_size = 400
random_indexes_2 = sample(1:dim_golub[1],second_size)
third_size = 1000
random_indexes_3 = sample(1:dim_golub[1],third_size)

# silhouette plot to find best number of cluster k
p4 <- fviz_nbclust(t(golub[random_indexes_1,]), kmeans, method = "silhouette")

# silhouette plot to find best number of cluster k
p5 <- fviz_nbclust(t(golub[random_indexes_2,]), kmeans, method = "silhouette")

# silhouette plot to find best number of cluster k
p6 <- fviz_nbclust(t(golub[random_indexes_3,]), kmeans, method = "silhouette")

plot_grid(p4,p5,p6)



# Finally,
# When we look at the differences between,
# first and second figure, we see that there are such differences since random selection.
