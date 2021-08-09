# comparekin

```r
install.packages("devtools")
devtools::install_github("samellisq/comparekin")
library(socmixmods)
```

# To get a pairwise relatedness matrix

Simple example. This snippet of code just generates some code and then calculates a pairwise relatedness matrix

```r
rm(list = ls())

#This just the test data from the kinship2 package
test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
                    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
                    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
                    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))
#some renaming
names(test1)[2] = "mother"
names(test1)[3] = "father"
test1$sex = ifelse(test1$sex ==1, "F", "M")

#change the names from numbers: my funtions need characters and also no zeros
library(babynames)
Mnames = unique(filter(babynames, sex == "M", prop >0.0005)$name)
Fnames = unique(filter(babynames, sex == "F", prop >0.0005)$name)
Mnames = sample(Mnames, 200)
Fnames = sample(Fnames, 200)
allnames = c(rbind(Mnames,Fnames))
allnames

test1$id = allnames[test1$id]
test1$mother = ifelse(test1$mother!=0, allnames[test1$mother], "UNK")
test1$father = ifelse(test1$father!=0, allnames[test1$father], "UNK")
test1


## OK we're ready. Calculate the pairwise relatedness matrix
r.mat = 
  get__pairwise_relatednesses_matrix(
                                   ids = test1$id , # IDs between which to calcualte relatedness, best just to do what I've done here and copy the id column in. But you could also only put in ids present in a certain year, for example. 
                                   population.data = test1, # this is the usual input data (slightly renamed) a dataframe with columns: id, mother, father, sex
                                   parents.to.use = "mother.only", # other options: "father.only", "all.avaliable" [default]
                                   r.degree = 2 #, # pedigree depth
                                   # kinship.network = NULL, # commented out becaasue I recommend leave this as null but might change answer if it is really slow to run
                                   # r.type = "coef.of.relationship" # again leave this as it is. OTher option is: kinship.coef which is just coef.of.r /2
)


r.mat
```
# To get all relatives to the nth degree

This is just a simple example of how to use this package to get all of an individuals related
in the nth degree to a given individual (the "ego.id"" or "from.id")

```r
##USE PACKAGE
rm(list = ls())
test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
                    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
                    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
                    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))
names(test1)[2] = "mother"
names(test1)[3] = "father"
test1$sex = ifelse(test1$sex ==1, "F", "M")
library(babynames)
Mnames = unique(filter(babynames, sex == "M", prop >0.0005)$name)
Fnames = unique(filter(babynames, sex == "F", prop >0.0005)$name)
Mnames = sample(Mnames, 200)
Fnames = sample(Fnames, 200)
allnames = c(rbind(Mnames,Fnames))
allnames

test1$id = allnames[test1$id]
test1$mother = ifelse(test1$mother!=0, allnames[test1$mother], "UNK")
test1$father = ifelse(test1$father!=0, allnames[test1$father], "UNK")
test1


## USE PACKAGE
population.network = make_kinship_network(make_relation_df(test1, parents.to.use = "mother.only"))
#note the parents.to.use term. Can be "mothers.only" =  maternal relatendess, "fathers.only" = paternal relateness or "all.avalaible" = both 
test1

get_nth_degree_relatives(network = population.network, from.id = test1$id[5], r.degree = 2)

ego.network = get_nth_degree_ego_network(network = population.network,
                                         ego.id = test1$id[5],
                                         r.degree = 2)
plot_ego_network(ego.network, census.data = test1)


```
