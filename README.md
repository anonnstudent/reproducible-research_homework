# Reproducible research: version control and R

### Question 1-3: Link to "logistic_growth" repo

[<https://github.com/annonpers/logistic_growth>](https://github.com/anonnstudent/logistic_growth)

### Question 4

**4.1** - *A script for simulating a random_walk is provided in the `question-4-code` folder of this repo. Execute the code to produce the paths of two random walks. What do you observe? (10 points)*

The two graphs have different axes and the plots span almost the entire graph. Plot1 (left) axis are x = -7.5 to 0.5 and y = -6 to 1, and Plot2 (right) axis are x = -3.5 to 3.5 and y = -1.5 to 3.5. Due to this, it can be concluded that these two random walks span different areas relative to each other. Both walks start at the same point (0,0) but Plot1 finishes further from this point (-5.25,-2) in comparison to Plot2 (-1.95,0.1). From this, it can also be seen that both these plots finish on the left from their original starting point. However, Plot1 mainly stays within the (-ve,-ve) areas of the plot apart from a small number of time points within the first 100t. Plot2 is centred around (0,0) more in comparison but spends most of the walk around (+ve/-ve, +ve). The walk in Plot1 appears to be more spread out as there are fewer overlapping routes that circle over itself in comparison to Plot2. However, both plots have multiple time points where the walk stays within the area such as those centred around (-6.2,-2.5) for Plot1 and (-1.5, 2.25) for Plot2.

Below are the plots discussed. 

<img width="988" alt="Screenshot 2023-12-07 at 09 38 27" src="https://github.com/anonnstudent/reproducible-research_homework/assets/150166627/4dd0eaf0-6f17-431a-8238-f48b97020ae7">



**4.2** - *Investigate the term **random seeds**. What is a random seed and how does it work? (5 points)*

A random seed is used to make generating random numbers reproducible. In R, it is an integer vector that stores the state of the random number generator (RNG). RNGs cannot be truly random as they use deterministic algorithms to produce numbers that appear statistically random. Setting a random seed allows us to run an RNG that uses the same algorithm each time the code is run, meaning the same random numbers are always produced. If a random seed is not used, then each time the code for the RNG is run the output will be different, making it impossible for anyone to replicate your results. In R, the function set.seed( ) is used to specify the random seed used for RNGs within the code it is run for. This allows the same numbers to always be produced by this code making it reproducible.

**4.3** - *Edit the script to make a reproducible simulation of Brownian motion. Commit the file and push it to your forked `reproducible-research_homework` repo. (10 points)*

See folder `question-4-code` for script.

**4.4** - *Go to your commit history and click on the latest commit. Show the edit you made to the code in the comparison view.*

<img width="1392" alt="Screenshot 2023-12-07 at 14 21 32" src="https://github.com/anonnstudent/reproducible-research_homework/assets/150166627/ecdf84a3-286d-4117-931a-13a8b149e7fd">


### Question 5

Code below is also in a script in folder `question-5-answers`

**5.1** - *Import the data for double-stranded DNA (dsDNA) viruses taken from the Supplementary Materials of the original paper into Posit Cloud (the csv file is in the `question-5-data` folder). How many rows and columns does the table have? (3 points)*

There are 33 rows and 13 columns in this table.

```{r}
#install.packages("readr")
library(readr)

#load data from the study
dsdna_virus_data <- read_csv("question-5-data/Cui_etal2014.csv")

#shows the number of columns and rows within table
View(dsdna_virus_data)
```

**5.2** - *What transformation can you use to fit a linear model to the data? Apply the transformation. (3 points)*

As the equation V = $\beta$ L\^($\alpha$) is allometric (nonlinear) this data can be log-transformed to make it linear. To ensure this transformation would be appropriate for the data I also visualised the data for both virion volume and genome length.

The log transformation can be seen in the code, where columns showing logged values are added to the data. Now a linear model can be fitted to the data: ln(V) = ln($\beta$) + $\alpha$\*ln(L).

```{r}
#install.packages("janitor")
#install.packages("ggplot2")
#install.packages("gdplyr")
library(janitor)
library(ggplot2)
library(dplyr)

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

```

**5.3** - *Find the exponent (*$\alpha$*) and scaling factor (*$\beta$*) of the allometric law for dsDNA viruses and write the p-values from the model you obtained, are they statistically significant? Compare the values you found to those shown in **Table 2** of the paper, did you find the same values? (10 points)*

To find the exponent and scaling factor a linear regression can be run on the log-transformed data, producing a linear model. The output of this model shows ln($\beta$) as the intercept estimate and $\alpha$ as the slope estimate, ln($\beta$) must be transformed to $\beta$. Their values and corresponding p-values are as follows: exponent($\alpha$) = 1.5152 (p = 6.44e-10) scaling factor ($\beta$) = 1181.807 (p = 2.28e-10) As both these p-values are below p = 0.05 (and p = 0.001), the values obtained from the model are statistically significant. In the paper, they found exponent = 1.52 and scaling factor = 1182 for dsDNA viruses. These are the same values that were found in this analysis, rounded to 3sf and 4sf respectively.

```{r}
#produce linear model from the logged data 
linear_model <- lm(log_virion_volume_nm_nm_nm ~ log_genome_length_kb, data = logged_virus_data)

#visualise model output: intercept estimate = log(B), log_genome_length_kb estimate = a, p-values also shown
summary(linear_model)

#tranform log(B) to B
exp(7.0748)
```

**5.4** - *Write the code to reproduce the figure shown below. (10 points)*

```{r}
ggplot(logged_virus_data, aes(x = log_genome_length_kb, y = log_virion_volume_nm_nm_nm)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  labs(x = expression(bold("log [Virion volume (nm3)")), y = expression(bold("log [Genome length (kb)]")))

```
<img width="797" alt="Screenshot 2023-12-07 at 14 39 13" src="https://github.com/anonnstudent/reproducible-research_homework/assets/150166627/86b49860-d96b-4971-a81f-7609a95c1678">


**5.5** - *What is the estimated volume of a 300 kb dsDNA virus? (4 points)*

The estimated virion volume for a 300 kb dsDNA virus = 6697006 nm\^3

```{r}
#function to calculate volume
calculate_vol<- function(B, L, a) {
  V <- B * L^a
  return(V)
}

#calculate V when L = 300
vol <-  calculate_vol(1181.807, 300, 1.5152)
vol
```

### Bonus Question

**(10 points)** - *Explain the difference between reproducibility and replicability in scientific research. How can git and GitHub be used to enhance the reproducibility and replicability of your work? what limitations do they have?*

Reproducibility is the ability of another scientific researcher to obtain the same results using the methods or procedures used within the original study. If a study is reproducible then it is ensured that the results of the original study were valid and not artificial or an error.
Replicability is the ability for similar results to be obtained within a new study using different methods. If similar results are replicated by multiple studies, then this indicates they are reliable. 
Git and GitHub can enhance both within someone’s work for multiple reasons. Most importantly, they allow you to track all the changes made to your work (by yourself or others) and build a history of the project so it can be more easily reproduced in the future (with access to this information). A detailed description of the project can also be available in README files which help people understand it in more detail, again aiding its reproducibility. GitHub allows individuals to fork other people’s repositories to create their own copy that they can modify as they would like without affecting the original. This can benefit both reproducibility (as the project itself can be run) but also replicability, as different methods could be applied to see if similar results are obtained. 
One limitation of git and GitHub is that they requires a good understanding by all collaborators within a project. If some people within a research group are not as efficient at using them or if work and plans are not communicated effectively then this could lead to errors occurring, which may be difficult to correct. 
Furthermore, they are limited to coding only. There are other ways of carrying out scientific investigation, other than coding, that are not, or not easily, accessible through git and GitHub. Practical aspects of studies cannot be shared on this platform but can be shared on others such as protocols.io. As a result of this, the usefullness of git and GitHub is limited to coding, modelling, data analysis and similar computational science.




## Instructions

The homework for this Computer skills practical is divided into 5 questions for a total of 100 points (plus an optional bonus question worth 10 extra points). First, fork this repo and make sure your fork is made **Public** for marking. Answers should be added to the \# INSERT ANSWERS HERE \# section above in the **README.md** file of your forked repository.

Questions 1, 2 and 3 should be answered in the **README.md** file of the `logistic_growth` repo that you forked during the practical. To answer those questions here, simply include a link to your logistic_growth repo.

**Submission**: Please submit a single **PDF** file with your candidate number (and no other identifying information), and a link to your fork of the `reproducible-research_homework` repo with the completed answers. All answers should be on the `main` branch.

## Assignment questions

1)  (**10 points**) Annotate the **README.md** file in your `logistic_growth` repo with more detailed information about the analysis. Add a section on the results and include the estimates for $N_0$, $r$ and $K$ (mention which \*.csv file you used).

2)  (**10 points**) Use your estimates of $N_0$ and $r$ to calculate the population size at $t$ = 4980 min, assuming that the population grows exponentially. How does it compare to the population size predicted under logistic growth?

3)  (**20 points**) Add an R script to your repository that makes a graph comparing the exponential and logistic growth curves (using the same parameter estimates you found). Upload this graph to your repo and include it in the **README.md** file so it can be viewed in the repo homepage.

4)  (**30 points**) Sometimes we are interested in modelling a process that involves randomness. A good example is Brownian motion. We will explore how to simulate a random process in a way that it is reproducible:

    -   A script for simulating a random_walk is provided in the `question-4-code` folder of this repo. Execute the code to produce the paths of two random walks. What do you observe? (10 points)
    -   Investigate the term **random seeds**. What is a random seed and how does it work? (5 points)
    -   Edit the script to make a reproducible simulation of Brownian motion. Commit the file and push it to your forked `reproducible-research_homework` repo. (10 points)
    -   Go to your commit history and click on the latest commit. Show the edit you made to the code in the comparison view (add this image to the **README.md** of the fork). (5 points)

5)  (**30 points**) In 2014, Cui, Schlub and Holmes published an article in the *Journal of Virology* (doi: <https://doi.org/10.1128/jvi.00362-14>) showing that the size of viral particles, more specifically their volume, could be predicted from their genome size (length). They found that this relationship can be modelled using an allometric equation of the form $`V = \beta L^{\alpha}`$, where $`V`$ is the virion volume in nm<sup>3</sup> and $`L`$ is the genome length in nucleotides.

    -   Import the data for double-stranded DNA (dsDNA) viruses taken from the Supplementary Materials of the original paper into Posit Cloud (the csv file is in the `question-5-data` folder). How many rows and columns does the table have? (3 points)
    -   What transformation can you use to fit a linear model to the data? Apply the transformation. (3 points)
    -   Find the exponent ($\alpha$) and scaling factor ($\beta$) of the allometric law for dsDNA viruses and write the p-values from the model you obtained, are they statistically significant? Compare the values you found to those shown in **Table 2** of the paper, did you find the same values? (10 points)
    -   Write the code to reproduce the figure shown below. (10 points)

<p align="center">

<img src="https://github.com/josegabrielnb/reproducible-research_homework/blob/main/question-5-data/allometric_scaling.png" width="600" height="500"/>

</p>

-   What is the estimated volume of a 300 kb dsDNA virus? (4 points)

**Bonus** (**10 points**) Explain the difference between reproducibility and replicability in scientific research. How can git and GitHub be used to enhance the reproducibility and replicability of your work? what limitations do they have? (e.g. check the platform [protocols.io](https://www.protocols.io/)).
