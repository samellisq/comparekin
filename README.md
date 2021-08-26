# comparekin

# To install:

```r
install.packages("devtools")
devtools::install_github("samellisq/comparekin")
library(comparekin)
```
# Introduction

comparekin is a package to calculate pairwise relatedness from pedigree data
while (1) controlling pedigree depth and (2) idenitfying unknown pairwise 
relatedness's. 

Briefly: 

(1) Pedigrees in different systems can be of different depths. Within populations
 this can also be true with individuals first observed at the start of the study 
 will have shallower pedigrees than individuals seen later in the study. This
 package limits the pedigree depth for each indivdiual when calcualting relatedness
 controlling for this vairablity both within and between populations. 

(2) the explicit or implicit assumption of many functions for calculating 
pairwise relatedness from pedigree data is that individuals who do not share a
relative have a relatedness of 0. In biological pedigree data, which often has
many missing maternties and/or paternities this assumption is often unlikely to 
be valid. This package offers a method to identify pairs of indivdiual's whose 
relatedness cannot be established given the pedigree.

The method works on the basis of treating the pedigree as a directed network.

These confounding factors and the methodology used here to avoid them will be described 
in more detail in a future publication.

Below are two examples outlining the functionality and uses of this package. 


# To get a pairwise relatedness matrix

Simple example. This snippet of code just generates some code and then calculates a pairwise relatedness matrix

```r
#Example data taken from kinship2::kinship()
test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
                    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
                    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
                    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))
#some renaming
names(test1)[2] = "mother"
names(test1)[3] = "father"
test1$sex = ifelse(test1$sex ==1, "F", "M")

test1$id = as.character(test1$id)
test1$mother = ifelse(test1$mother!=0, as.character(test1$mother), "UNK")
test1$father = ifelse(test1$father!=0, as.character(test1$father), "UNK")
test1

r.mat =
  get__pairwise_relatednesses_matrix(
    ids = test1$id ,
    population.data = test1,
    parents.to.use = "all.avaliable",
    r.degree = 2,
    kinship.network = NULL,
    r.type = "coef.of.relationship"
  )
r.mat

```
# To get all relatives to the nth degree

This is just a simple example of how to use this package to get all of an individuals related
in the nth degree to a given individual (the "ego.id"" or "from.id")

```r
#Example data taken from kinship2::kinship()
test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
                    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
                    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
                    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))
#some renaming
names(test1)[2] = "mother"
names(test1)[3] = "father"
test1$sex = ifelse(test1$sex ==1, "F", "M")

test1$id = as.character(test1$id)
test1$mother = ifelse(test1$mother!=0, as.character(test1$mother), "UNK")
test1$father = ifelse(test1$father!=0, as.character(test1$father), "UNK")
test1


## USE PACKAGE
r.df = make_relation_df(test1)
net = make_kinship_networks(r.df)

get_nth_degree_relatives(network = population.network, from.id = test1$id[5], r.degree = 2)

ego.network = get_nth_degree_ego_network(network = population.network,
                                          ego.id = test1$id[5],
                                          r.degree = 2,
                                          seperate.lineages = TRUE
                                          )

plot_ego_network(ego.network, census.data = test1)

```
