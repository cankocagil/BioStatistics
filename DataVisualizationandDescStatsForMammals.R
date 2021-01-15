################################
## Author     : Can Kocagil  ###
################################

################################
####  --- Dependencies ---  ####
#### 1) ggplot2             ####
#### 1) cowplot             ####
################################





# 0) Uploading necessary R packages to proceed
# 1) Explanatory data analysis to familiarize the dataset(EDA)
# 2) Extra plotting for EDA
# 3) Statistical features of the dataset (Descriptive statistics)
# 4) Question for Hypothesis testing
# 5) Trying to answer question by analyzing data

##################################################################################################################
##################################################################################################################

# ------- Step 0 ------- #
##########################

# Uploading necessary packages :
# install.packages('ggplot2')
library(ggplot2)

data(msleep)

## ---- Step 0 Finish ---- ##

# ------- Step 1 ------- #
##########################

# Preliminary Explanatory data analysis to familiarize the dataset(EDA)

# Getting the data
df <- msleep

# Looking first rows from beginning :
head(df)

# Loaking last rows from end :
tail(df)

ls(df)
# Getting column names :
df_col_names = colnames(df)
df_col_names

# Reading the written features of the data
?msleep

# here are the explanations of the features :
#############################################

# name
# ---- common name

# genus
# vore
# --- carnivore, omnivore or herbivore?
  
# order
# conservation
# --- the conservation status of the animal

# sleep_total
# --- total amount of sleep, in hours

# sleep_rem
# --- rem sleep, in hours

# sleep_cycle
# --- length of sleep cycle, in hours

# awake
# --- amount of time spent awake, in hours

# brainwt
# --- brain weight in kilograms

# bodywt
# --- body weight in kilograms

#############################################

# Dimension/size of the dataset
size_df = dim(df) ;size_df

# So we have 83 rows with 11 features

# Checking the missing values of the dataset
df_na = is.na.data.frame(df) ;df_na


##########################
## ---- Step 1 Finish ---- ##

# ------- Step 2 ------- #
##########################

# Extra plotting for Explanatory Data Analysis :
# In this section, I am going to visualize the data to gain insight of it

# the columns 6-9 are selected since they are not categorical(numerical) :
pairs(df[ ,c(6:9)]) 

# Plot the distribution of vores in the msleep dataset
pie(table(df$vore))

# Histogram total sleep grouped by vores
ggplot(data = df,aes(df$sleep_total)) + geom_histogram() +facet_grid(~vore)

# Histogram of the number of hours per day mammals are asleep (using 20 bins along the x axis)
ggplot(data=df, aes(x=sleep_total)) +       
       geom_histogram(bins=20) +                     
       labs(x='Total sleep time (h)', 
       y='Frequency') + theme_bw()

# Explanatory data analysis for body weights and sleeping time grouped by vore
ggplot(df)+ geom_point(aes(log10(df$bodywt),sleep_total))+facet_wrap(~vore)

# Explanatory data analysis for brain weights and sleeping time grouped by vore
ggplot(df)+ geom_point(aes(log10(df$brainwt),sleep_total))+facet_wrap(~vore)


# Labels for following plotting :
labels = c('herbi'='Herbivore','carni'='Carnivore','omni'='Omnivore','insecti'='Insectivore')

# Visuzalization of the vores w.r.t. number of them
ggplot(data=df,aes(x=vore,fill=vore)) +                        
       geom_bar() +                                 
       scale_x_discrete(name = 'Feeding Type',labels=labels) + 
       scale_fill_discrete(name = 'Feeding Type',    
       labels=labels) + 
       labs(x='Feeding Type',y='# of Species') 
       
       

# Plot REM sleep against total sleep in units of hour per day :
ggplot(data=df,aes(x=sleep_total,y=sleep_rem)) +
       geom_point(na.rm=TRUE) +                  
       labs(x='Total sleep',y='REM sleep')+
       geom_smooth(method='lm',na.rm=T)
       


# ------- Step 3 ------- #
##########################

# Statistical features of the dataset :


# Basic central tendency analysis :
# Measures of central tendency for sleep_total
# The argument na.rm=TRUE removes any missing values 

# Expected value for sleep time
mean(df$sleep_total, na.rm=T)    

# Medium value for sleep time
median(df$sleep_total, na.rm=T) 

# Now looking the spread of the data :
sd(df$sleep_total,na.rm = T)

# Mean absolute deviation :
mad(df$sleep_total, na.rm=T) 

# Inter quartile range :
IQR(df$sleep_total, na.rm=T)  

# Getting statistical information about the dataset
summary(df)

# Maximum sleep time :
max(df$sleep_total, na.rm=T) 

# Minimum sleep time :
min(df$sleep_total, na.rm=T)

# Range of sleep time : (Maximum - Minimum)
range(df$sleep_total, na.rm=T)  



# Correlation part :

# I am going to use Pearson correlation coefficient :
# Correlations close to -1 are very strong negative relationships
# Correlations close to 0 are weak relationships
# Correlations close to 1 are very strong positive relationships

# Finding the correlation between body weight and total sleep time :
cor((df$bodywt),msleep$sleep_total,use = "complete.obs",method = 'pearson')

# Finding the correlation between body weight and total sleep time in log10 bases :
cor(log10(df$bodywt),msleep$sleep_total,use = "complete.obs",method = 'pearson')

# Finding the correlation between brain weight and total sleep time :
cor((df$brainwt),msleep$sleep_total,use = "complete.obs",method = 'pearson')

# Finding the correlation between brain weight and total sleep time in log10 bases :
cor(log10(df$brainwt),msleep$sleep_total,use = "complete.obs",method = 'pearson')

# Qualitative analysis :
table(df$vore)
table(df$sleep_total,df$vore)

# I am going to further Qualitative analysis part in the visualization part

##########################


# ------- Step 4 ------- #
##########################
# Question for Hypothesis testing :
# I am going to find an answer to the question :
# Does body weight or brain weight, which one have more
# impact on with total rem sleep time?

# Expected value for sleep time
expected = mean(df$sleep_rem, na.rm=T)  ;expected

# Medium value for sleep time
medium = median(df$sleep_rem, na.rm=T) ;medium 

# Now looking the spread of the data :
std = sd(df$sleep_rem,na.rm = T);std

# Maximum sleep time :
max_val =max(df$sleep_rem, na.rm=T) ;max_val

# Minimum sleep time :
min_val = min(df$sleep_rem, na.rm=T);min_val


# Finding the correlation between body weight and total sleep time :
corr = cor((df$bodywt),df$sleep_rem,use = "complete.obs",method = 'pearson')

# Finding the correlation between body weight and total sleep time in log10 bases :
corrlog10 = cor(log10(df$bodywt),df$sleep_rem,use = "complete.obs",method = 'pearson');corrlog10

info = cat('Mean :', round(expected,2), '\nMedium :', round(medium,2), '\nStandart Deviation :' ,round(std));info
##########################

# ------- Step 5 ------- #
##########################

# Trying to answer question by analyzing data :
library(cowplot)


theme <- theme(plot.title = element_text(face = "bold",size = 14,color = 'deepskyblue1',hjust = 0.5),
      plot.subtitle = element_text(face = 'bold',color = 'greenyellow',hjust = 0.5),
      plot.caption = element_text(color = "grey52", face = "italic",hjust = 0.5),
      legend.position = 'right',
      legend.title = element_text(colour="black", size=10, 
                                  face="bold"),
      legend.background = element_rect(fill="lightblue",
                                       size=0.5, linetype="solid", 
                                       colour ="darkblue"),
      axis.title.x = element_text(color="darkblue", size=14, face="bold"),
      axis.title.y = element_text(color="darkblue", size=14, face="bold"))

labels = c('herbi'='Herbivore','carni'='Carnivore','omni'='Omnivore','insecti'='Insectivore')
caption1 = 'Descriptive Statistics of Total Rem Sleep Time \n
Mean : 1.87 Medium : 1.5 Standart Deviation : 1.29\n Max : 6.6 Min : 0.1
       Pearson Correlation : -0.32
       Data source : Msleep '
caption2 = 'Descriptive Statistics of Total Rem Sleep Time \nMean : 10.43 Medium : 10.1 
Standart Deviation : 4.45\n Max : 19.9 Min : 1.9
       Pearson Correlation : -0.56
       Data source : Msleep' 

p1 <- ggplot(data=msleep,aes(x=log10(df$bodywt),y=sleep_rem,fill = vore,colour = log10(bodywt))) +
  geom_point(na.rm=TRUE) +
  geom_smooth(method='loess',na.rm=T)+
  facet_wrap(~vore)+
  labs(x='Total Rem Sleep(hour per day)',y='Body Weight(kg)',caption = caption1,
       colour = 'Body Weight',fill = 'Vores')+
  ggtitle('The Effect of Body Weight on Total Rem Sleep Time\n',subtitle = 'Mammals Grouped by Vores')+theme
  
p2 <- ggplot(data=msleep,aes(x=log10(df$brainwt),y=sleep_rem,fill = vore,colour = log10(brainwt))) +
  geom_point(na.rm=TRUE) +  geom_smooth(method='loess',na.rm=T)+
  facet_wrap(~vore)+  labs(x='Total Rem Sleep(hour per day)',y='Brain Weight(kg)',caption = caption2,
   colour = 'Brain Weight',fill = 'Vores')+
  ggtitle('The Effect of Brain Weight on Total Rem Sleep Time\n',subtitle = 'Mammals Grouped by Vores')+  theme

p3 <-ggplot(data=df,aes(x=vore,fill=vore)) +geom_bar() +scale_x_discrete(name = 'Feeding Type',labels=labels) + 
  scale_fill_discrete(name = 'Feeding Type',labels=labels)+labs(x='Feeding Type',y='# of Species')


p4 <-ggplot(data=df,aes(x=sleep_total,y=sleep_rem))+geom_point(na.rm=T) +labs(x='Total sleep',y='REM sleep')+
     geom_smooth(method='loess',na.rm=T)

plot_grid(p1,p2,p3,p4)

