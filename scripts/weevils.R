# if git is ahead by X commits do this: git reset --soft HEAD~1 (8=# of commits)

## ---- analysis ----
library(dplyr)
library(ggplot2)
library(tidyverse)
library(flextable)
library(officer)
library(kableExtra)
library(knitr)
library(officedown)
library(tidybayes)
library(brms)
library(broom.mixed)
library(parameters)
library(posterior)





#### read in data file ####
behaviour <- read.csv(file="data/raw/female_earwig_personality_data.csv", header=TRUE, sep=",", dec=".") %>%
  as.data.frame()

## read in family (id_mere) data and join to behaviour file
family <- read.csv(file="data/raw/earwig_data.csv", header=TRUE, sep=",", dec=".") %>%
  as.data.frame()














## ---- end