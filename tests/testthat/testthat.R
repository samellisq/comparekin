library(testthat)
library(comparekin)

# test_check("comparekin")

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

r.mat =
  get__pairwise_relatednesses_matrix(
    ids = test1$id ,
    population.data = test1,
    parents.to.use = "all.avaliable",
    r.degree = 2,
    kinship.network = NULL,
    r.type = "coef.of.relationship"
  )

test_that("Output of get__pairwise_relatednesses_matrix is correct",{

  expect_equal(length(test1$id), dim(r.mat)[1])
  expect_equal(length(test1$id), dim(r.mat)[2])


})

r.df = make_relation_df(test1)
net = make_kinship_network(r.df)


test_that("Output of make_kinship_networks is a network", {

  expect_equal(igraph::vcount(net),length(test1$id) )

})
