
#install.packages("readr")
#install.packages("janitor")
#install.packages("ggplot2")
#install.packages("gdplyr")
#install.packages("forcats")
library(readr)
library(janitor)
library(ggplot2)
library(dplyr)
library(forcats)


# Question 5.1 

#load data from the study
dsdna_virus_data <- read_csv("question-5-data/Cui_etal2014.csv")

#shows the number of columns and rows within table
View(dsdna_virus_data)



# Question 5.2 

#Need to visualize the data on virion volume and genome length
names(dsdna_virus_data)

#Column names for these are "Virion volume (nm×nm×nm)" and "Genome length (kb)" 
    #These need to cleaned so can be analysed 
clean_virus_data <- clean_names(dsdna_virus_data)
names(clean_virus_data)

#visualize data for virion volume 
ggplot(clean_virus_data, aes(x = virion_volume_nm_nm_nm)) +
  geom_histogram()

#visualize data for genome length
ggplot(clean_virus_data, aes(x = genome_length_kb)) +
  geom_histogram()

#data for both virion volume and genome length appear exponential so can be log-transformed
logged_virus_data <- clean_virus_data %>%
  mutate(log_virion_volume_nm_nm_nm = log(virion_volume_nm_nm_nm),
         log_genome_length_kb = log(genome_length_kb))



# Question 5.3

#produce linear model from the logged data 
linear_model <- lm(log_virion_volume_nm_nm_nm ~ log_genome_length_kb, data = logged_virus_data)

#visualise model output: intercept estimate = log(B), log_genome_length_kb estimate = a, p-values also shown
summary(linear_model)

#tranform log(B) to B
exp(7.0748)



# Question 5.4 

ggplot(logged_virus_data, aes(x = log_genome_length_kb, y = log_virion_volume_nm_nm_nm)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  labs(x = expression(bold("log [Virion volume (nm3)")), y = expression(bold("log [Genome length (kb)]")))


# Question 5.5 

#function to calculate volume
calculate_vol<- function(B, L, a) {
  V <- B * L^a
  return(V)
}

#calculate V when L = 300
vol <-  calculate_vol(1181.807, 300, 1.5152)
vol



