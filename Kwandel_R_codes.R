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
hw <- theme_gray() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0),
    
    strip.text.y = element_blank(),
    strip.background = element_rect(fill = rgb(0.9, 0.95, 1), colour = gray(0.5), size = 0.2),
    
    panel.border = element_rect(fill = FALSE, colour = gray(0.70)),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing.x = unit(0.10, "cm"),
    panel.spacing.y = unit(0.05, "cm"),
    
    axis.ticks = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.text.y = element_text(margin = margin(0, 3, 0, 3)),
    axis.text.x = element_text(margin = margin(-1, 0, 3, 0))
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

summary(df_clean)
table(df_clean$benormal)
# Pairs Panel to show correlations, histograms and scatter plots
pairs.panels(df_clean)
# There are some variables that show a high correlation to each other, but not of them were significant
# Some of the variables show a skewness to the left.
# So we should log the predictor variables to make them more normally distrubted.
variables_to_log <- c("mitoses", "normalnucleoli", "blandchromatin", "epithelial", "margadhesion", "uniformcellshape", "uniformcellsize","clumpthickness")
df_log <- df_clean
df_log[variables_to_log] <- lapply(df_log[variables_to_log], log)


# First model, logistic regression with no variable changing
# For our dataset we want to predict if benormal = 1, therefore we will initially be using a logisitc regression model
model_1 <- glm(benormal ~ ., data = df_clean, family = binomial)
# Summary of our initial Model
summary(model_1)
model_1
# First model, logistic regression with logging the predictor variables
model_2 <- glm(benormal ~ ., data = df_log, family = binomial)
# Summary of our initial Model
summary(model_2)
model_2
# Compare the models
AIC(model_1, model_2)
# Hypothesis test 2
anova(model_1, model_2, test = "Chisq")
# No significance difference between the models in fact, non-log was better.
# Now lets make the model final by removing some of the insignificant varaibles. The less variables needed to explain the data/results the better. 
model_3 <- glm(benormal ~ clumpthickness + margadhesion + barenuclei + blandchromatin, data = df_clean, family = binomial)
summary(model_3)
model_3
anova(model_1, model_3, test = "Chisq")
# P-value is <.05 so we can conclude the more complex model is significantly better than the simpler model.
# For logisitic regression modeling the best model is:
coef(model_1)

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
# Tree is maxed at 4, 3 is very close to 4 
prune.df_clean=prune.tree(tree.df_clean,best=3)
plot(prune.df_clean)
text(prune.df_clean,pretty=0)
# Testing models
yhat=predict(tree.df_clean,newdata=df_clean[-train,])
df_clean.test=df_clean[-train,"benormal"]#testing y values

mean((yhat-df_clean.test)^2) #MSE testing

# Creating a random forest model
set.seed(1)
train = sample(1:nrow(df_clean), nrow(df_clean)/2)
cancer_rf=randomForest(benormal~.,data=df_clean,subset=train,mtry=3,importance=TRUE)
yhat.rf = predict(cancer_rf,newdata=df_clean[-train,])
df_clean.test=df_clean[-train,"benormal"]
mean((yhat.rf-df_clean.test)^2) #test set MSE
# plotting the forest
plot(cancer_rf)
# What are the most important factors?
varImpPlot(cancer_rf)
cancer_rf


# Do a PCA on the models reduced dataset for question 3
df_clean$benormal <- as.integer(df_clean$benormal)
states=row.names(df_clean)
apply(df_clean, 2, mean)
apply(df_clean, 2, var)
# There isn't much difference in the mean and variance, so we do not need to scale/standardize the variables.
# Calculate the principal components
pr.out=prcomp(df_clean, scale=TRUE)
names(pr.out)
pr.out$center
head(pr.out$x,n = 10)

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
