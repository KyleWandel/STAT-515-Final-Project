# Code for Data modeling and explortation of breast cancer data from 
# University of Wisconsin Hospitals, Madison from Dr. William H. Wolberg.

# Install Packages if needed # 
packages <- c( "randomForest","tree","dplyr", "tidyverse", "plotly", "psych", "ISLR", "leaps", "reshape2","readr")
# Check if packages are not installed, then install
for (package in packages) {
  if (!(package %in% installed.packages())) {
    install.packages(package, dependencies = TRUE)
  }
}
# Load libraries
library(dplyr)
library(tidyverse)
library(plotly)
library(psych)
library(ISLR)
library(leaps)
library(reshape2)
library(readr)
library(tree)
library(randomForest)
# Create theme for plots
hw <- theme_gray()+ theme(
  plot.title=element_text(hjust=0.5),
  plot.subtitle=element_text(hjust=0.5),
  plot.caption=element_text(hjust=-.5),
  
  strip.text.y = element_blank(),
  strip.background=element_rect(fill=rgb(.9,.95,1),
                                colour=gray(.5), linewidth =.2),
  
  panel.border=element_rect(fill=FALSE,colour=gray(.70)),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.spacing.x = unit(0.10,"cm"),
  panel.spacing.y = unit(0.05,"cm"),
  
  # axis.ticks.y= element_blank()
  axis.ticks=element_blank(),
  axis.text=element_text(colour="black"),
  axis.text.y=element_text(margin=margin(0,3,0,3)),
  axis.text.x=element_text(margin=margin(-1,0,3,0))
)
# Pull data from Github Respository
url <- "https://raw.githubusercontent.com/KyleWandel/STAT-515-Final-Project/main/breast-cancer-wisconsin.csv"
df <- read.table(url, header = TRUE, sep = ",")

head(df)
# remove ID column 
df <- subset(df, select = -id)
str(df)
df$barenuclei <- as.integer(df$barenuclei)
colSums(is.na(df))
# Our data has 16 rows that are missing data for the barenuclei column. For now we will omit these from the dataset
df_clean <- na.omit(df)
str(df_clean)

# Transform benormal column to malignant (4) = 1 and benign (2) = 0
df_clean$benormal <- ifelse(df_clean$benormal == 4, 1, ifelse(df_clean$benormal == 2, 0, df_clean$benormal))

# Evaluating the predictor variables
# Histrograms
create_histograms <- function(dataframe) {
  # Initialize an empty list to store plots
  plots_list <- list()
  
  # Loop through each column
  for (col in names(dataframe)) {
    # Skip columns that are not numeric
    if (!is.numeric(dataframe[[col]])) {
      next
    }
    
    # Create histogram for the column
    hist_plot <- ggplot(data = dataframe, aes(x = .data[[col]])) +
      geom_histogram(fill = "skyblue", color = "black", bins = 10) +
      labs(title = paste("Histogram of", col), x = col, y = "Frequency") +
      hw
    
    # Store the plot in the list
    plots_list[[col]] <- hist_plot
  }
  
  return(plots_list)
}

# Create side-by-side histograms
histograms <- create_histograms(df_clean)

# Print histograms
for (i in seq_along(histograms)) {
  print(histograms[[i]])
}

# Correlations
# VIF number comparison testing for mutlicollinearity
# Chisq test on predictor variables 
corr.test(df_clean)
# There are some variables that show a high correlation to each other, but not of them were significant
pairs.panels(df_clean)


# Scatters of relationship to response variable use the logistic regression notes
# loop through 


# For our dataset we want to predict if benormal = 1, therefore we will initially be using a logisitc regression model
model_1 <- glm(benormal ~ ., data = df_clean, family = binomial)
# Summary of our initial Model
summary(model_1)
model_1
# Explain the model (Hypothesis test 1?)
# VIF number comparison testing for mutlicollinearity
# Can we make this better through variable transformation and selection?
# A few of the predictor variables are exhibiting a right skew. We will tranform these variables using the log() function.
variables_to_log <- c("mitoses", "normalnucleoli", "blandchromatin", "epithelial", "margadhesion", "uniformcellshape", "uniformcellsize","clumpthickness")
df_log <- df_clean
df_log[variables_to_log] <- lapply(df_log[variables_to_log], log)
# Run the model Again
model_2 <- glm(benormal ~ ., data = df_log, family = binomial)
# Summary of our initial Model
summary(model_2)
model_2
# Compare the models
AIC(model_1, model_2)
# Hypothesis test 2
anova(model_1, model_2, test = "Chisq")
# No significance difference between the models in fact, non-log was better.

# We have tried to log some of the columns, lets now try a new modeling approach

# Using Decision Trees and Random Forest
set.seed(1)
train = sample(1:nrow(df_clean), nrow(df_clean)/2)
tree.df_clean=tree(benormal~.,df_clean,subset=train)
summary(tree.df_clean)

# Plot the model 
plot(tree.df_clean)
text(tree.df_clean,pretty=0)

# Checking if Pruning will improve the model
cv.df_clean=cv.tree(tree.df_clean)
plot(cv.df_clean$size,cv.df_clean$dev,type='b')
# Tree is maxed at 4, 3 is vert close to 4 
prune.df_clean=prune.tree(tree.df_clean,best=3)
plot(prune.df_clean)
text(prune.df_clean,pretty=0)
# Testing models
yhat=predict(tree.df_clean,newdata=df_clean[-train,])
df_clean.test=df_clean[-train,"benormal"]#testing y values

mean((yhat-df_clean.test)^2) #MSE testing

# Creating a random forest model
set.seed(2)
rf.df_clean=randomForest(benormal~.,data=df_clean,subset=train,mtry=6,importance=TRUE)
yhat.rf = predict(rf.df_clean,newdata=df_clean[-train,])
mean((yhat.rf-df_clean.test)^2) #test set MSE
# plotting the forest
plot(rf.df_clean)
# What are the most important factors?
varImpPlot(rf.df_clean)

# Do a PCA on the models reduced dataset for question 3
states=row.names(df_clean)
apply(df_clean, 2, mean)
apply(df_clean, 2, var)
# There isn't much difference in the mean and variance, so we do not need to scale/standardize the variables.

# Calculate the principal components
pr.out=prcomp(df_clean, scale=TRUE)
names(pr.out)
pr.out$center
head(pr.out$x,n = 10)
pairs.panels(pr.out$x)

pr.out$rotation=-pr.out$rotation
pr.out$x=-pr.out$x
biplot(pr.out, scale=0)

(pr.var=pr.out$sdev^2) #variance of each PC
pve=pr.var/sum(pr.var)
pve

# Scree plot
par(mfrow =c(1,2))
plot(pve, xlab="Principal Component", 
     ylab="Proportion of Variance Explained", ylim=c(0,1),type='b')

plot(cumsum(pve), xlab="Principal Component",
     ylab="Cumulative Proportion of Variance Explained", 
     ylim=c(0,1),type='b')
