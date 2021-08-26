#' Simplify data frame to pedigree information
#'
#'Simple function to prepare a data frame of population.data data frame to be
#'used in the \code{\link{make_kinship_network}} function.
#'
#'The function checks column names and trims a data frame to only include the
#'necessary columns. If population.data already only has the necessary columns
#'named correctly then this function is unnecessary, but will not do any harm.
#'
#'@param population.data Data frame containing columns named id and one or both
#'of 'mother' and 'father'.
#'@param parents.to.use Either \code{"all.avaliable"}, \code{"mother.only"} or
#'\code{father.only}. This defines the type of  network that will be created and
#'ultimately the type of relatedness calculated made either biparental, maternal
#'only or paternal only. Option \code{"all.avaliable"} will use both mother and
#'father if both columns are avalaible, otherwise will only use one or the other
#'as are present in the population.data.
#'
#'@return Returns a dataframe with one column named id, and columns mother,
#'father or both as appropriate.
#'
#'@examples
#'#Example data taken from kinship2::kinship()
#'test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#'                    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
#'                    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
#'                    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))
#'#some renaming
#'names(test1)[2] = "mother"
#'names(test1)[3] = "father"
#'test1$sex = ifelse(test1$sex ==1, "F", "M")

#'test1$id = as.character(test1$id)
#'test1$mother = ifelse(test1$mother!=0, as.character(test1$mother), "UNK")
#'test1$father = ifelse(test1$father!=0, as.character(test1$father), "UNK")
#'test1
#'
#'make_relation_df(test1, parents.to.use = "all.avaliable")
#'make_relation_df(test1, parents.to.use = "mother.only")
#'make_relation_df(test1, parents.to.use = "father.only")
#'
#'
#'@export
make_relation_df = function(population.data, parents.to.use = "all.avaliable"){
  if(parents.to.use == "all.avaliable"){
    if("mother" %in% names(population.data) & "father" %in% names(population.data)){
      relation.df = dplyr::select(population.data, id, mother, father)
    }
    if("mother" %in% names(population.data) & !("father" %in% names(population.data))){
      relation.df = dplyr::select(population.data, id, mother)
    }
    if("father" %in% names(population.data) & !("mother" %in% names(population.data))){
      relation.df = dplyr::select(population.data, id, father)
    }
    if(!("mother" %in% names(population.data)) & !("father" %in% names(population.data))){
      stop("Population.data must contain columns named one or both of 'mother' and 'father'")
    }
  }

  if(parents.to.use == "mother.only"){
    if("mother" %in%  names(population.data)){
      relation.df = dplyr::select(population.data, id, mother)
    } else {
      stop("population.df must contain a column called 'mother'")
    }

  }

  if(parents.to.use == "father.only"){
    if("father" %in%  names(population.data)){
      relation.df = dplyr::select(population.data, id, father)
    } else {
      stop("population.df must contain a column called 'father'")
    }

  }

  return(relation.df)
}

#' Convert pedigree to network
#'
#'\code{make_kinship_network} Creates a kinship network from a pedigree data
#'frame created by the function \code{\link{make_relation_df}}.
#'
#'The object returned is the population pedigree represented as a network.
#'This network underlies much of the subsequent relatedness analysis, and is
#'used to both limit the depth of the pedigree and find pairs where relatedness
#'is unknown.
#'
#'@param relation.data.frame a data frame produced by the
#'\code{\link{make_relation_df}} function.
#'
#'@return An igraph network object with nodes as all ids in the population
#'(including parents) and directed edges representing parent --> offspring
#'relationships.
#'
#'@examples
#'# Example data taken from kinship2::kinship()
#'test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#'                    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
#'                    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
#'                    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))
#'#some renaming
#'names(test1)[2] = "mother"
#'names(test1)[3] = "father"
#'test1$sex = ifelse(test1$sex ==1, "F", "M")

#'test1$id = as.character(test1$id)
#'test1$mother = ifelse(test1$mother!=0, as.character(test1$mother), "UNK")
#'test1$father = ifelse(test1$father!=0, as.character(test1$father), "UNK")
#'test1
#'
#'r.df = make_relation_df(test1)
#'net = make_kinship_network(r.df)
#'plot(net)
#'
#'@export
make_kinship_network = function(relation.data.frame){
  if(ncol(relation.data.frame) == 2){
    #little fiddle to make sure everything is in the right order
    idcol = which(names(relation.data.frame) == "id")
    parentcol = which(names(relation.data.frame) != "id")
    relation.data.frame = data.frame(relation.data.frame[, parentcol],
                                     relation.data.frame[,idcol]
    )
    names(relation.data.frame)= c("parent", "id")
    edge.type = "A"
    edgelist = as.matrix(relation.data.frame)
  }

  if(ncol(relation.data.frame) == 3){
    relation.data.frame = data.frame(from = c(relation.data.frame$mother, relation.data.frame$father),
                                     to = c(relation.data.frame$id, relation.data.frame$id)
    )
    names(relation.data.frame)= c("parent", "id")
    edge.type = "B"
    edgelist = as.matrix(relation.data.frame)
  }


  all.ids = unique(c(edgelist[,1],edgelist[,2]))


  adj.mat = matrix(0,length(all.ids), length(all.ids))
  rownames(adj.mat) = all.ids
  colnames(adj.mat) = all.ids



  adj.mat[edgelist] = 1
  if(sum(colnames(adj.mat) == "UNK") >0){
    adj.mat= adj.mat[-which(colnames(adj.mat)=="UNK"),-which(rownames(adj.mat)=="UNK")] # remove UNK from the matrix
  }


  network = igraph::graph_from_adjacency_matrix(adj.mat)

  return(network)

}



#' Create ego-pedigree
#'
#'\code{get_nth_degree_ego_network} trims the full population kinship-network
#'down to only the ego network for a given id up to the nth relatedness degree.
#'
#' @param network A kinship-network derived from the
#' \code{\link{make_kinship_network}} function. The type of network input -
#' biparental, maternal only, paternal only - defines which side(s) of relatives
#' are found.
#' @param ego.id Character. id of the ego for whom the ego network will be
#' drawn.
#' @param r.degree degree of relationship to which to calculate relatedness (
#' \url{https://en.wikipedia.org/wiki/Coefficient_of_relationship}. Degree is
#' calculated separately for each lineage i.e. half-sisters are considered
#' equivalent to full-sisters in count (but not in relatedness calculation).
#' @param seperate.lineages If TRUE each node in the ego network is given a
#' property \code{lineage} labeling each node as part of the lineage A, B or
#' both. Where A and B are paternal and maternal lineages (order determined by
#' the data).
#'
#' @return igraph network object showing showing all relatives to the
#' \code{r.degree} of the \code{ego.id}.
#'
#' @examples
#'#Example data taken from kinship2::kinship()
#'test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#'                    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
#'                    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
#'                    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))
#'#some renaming
#'names(test1)[2] = "mother"
#'names(test1)[3] = "father"
#'test1$sex = ifelse(test1$sex ==1, "F", "M")

#'test1$id = as.character(test1$id)
#'test1$mother = ifelse(test1$mother!=0, as.character(test1$mother), "UNK")
#'test1$father = ifelse(test1$father!=0, as.character(test1$father), "UNK")
#'test1
#'
#'r.df = make_relation_df(test1)
#'net = make_kinship_network(r.df)
#'
#'ego.network = get_nth_degree_ego_network(network = population.network,
#'ego.id = test1$id[5],
#'r.degree = 2,
#'seperate.lineages = TRUE)
#'
#'plot_ego_network(ego.network, census.data = test1)
#'
#'@export
#'
get_nth_degree_ego_network = function(network,
                                      ego.id,
                                      r.degree = 2,
                                      seperate.lineages = FALSE){

  ego.in.edges = igraph::incident(network, ego.id, mode = "in")

  if(length(ego.in.edges) >2){
    stop(paste("inidividual ", ego.id, " has more the 2 parents", sep =""))
  }

  if(seperate.lineages == TRUE & length(ego.in.edges) ==1){
    warning("Only one lineage avaliable. Lineages will not be seperated")
    seperate.lineages = FALSE
  }

  relatives = list()

  if(seperate.lineages == FALSE){
    relatives[[1]] = get_nth_degree_relatives(network, ego.id, r.degree)
  }

  if(seperate.lineages == TRUE){
    for(i in 1:2){
      Xlineal.net = igraph::delete_edges(network, ego.in.edges[i])
      relatives[[i]] = get_nth_degree_relatives(Xlineal.net, ego.id, r.degree)
    }
  }

  all.relatives = unique(unlist(relatives))
  ids = igraph::V(network)$name
  not.ego.ids = ids[!(ids %in% c(all.relatives, ego.id))]
  ego.network = igraph::delete_vertices(network, not.ego.ids)


  if(seperate.lineages == TRUE){
    igraph::V(ego.network)$lineage = ifelse(igraph::V(ego.network)$name %in% relatives[[1]], "A", "B")
    igraph::V(ego.network)$lineage = ifelse(igraph::V(ego.network)$name %in% relatives[[1]] & igraph::V(ego.network)$name %in% relatives[[2]], "both", igraph::V(ego.network)$lineage)
  } else {
    igraph::V(ego.network)$lineage = rep.int("A", igraph::vcount(ego.network))
  }
  igraph::V(ego.network)$lineage = ifelse(igraph::V(ego.network)$name == ego.id, "focal", igraph::V(ego.network)$lineage)


  return(ego.network)
}



#' Find all nth degree relatives from a pedigree-network
#'
#' \code{get_nth_degree_relatives} derives the ids of all relatives within the
#' nth degree for a given individual from a population kinship-network.
#'
#' @param network A kinship-network derived from the
#' \code{\link{make_kinship_network}} function. The type of network input -
#' biparental, maternal only, paternal only - defines which side(s) of relatives
#' are found.
#' @param from.id Character. The id of the individual for whom nth degree
#' relatives will be found
#' @param r.degree degree of relationship to which to calculate relatedness (
#' \url{https://en.wikipedia.org/wiki/Coefficient_of_relationship}. Degree is
#' calculated separately for each lineage i.e. half-sisters are considered
#' equivalent to full-sisters in count (but not in relatedness calculation).
#'
#'@return Returns a character vector of the ids of all relatives up to the
#'\code{r.degree} degree of the \code{from.id}. For example, where
#'\code{r.degree = 2} this will return a vector of all 1st and 2nd degree
#'maternal and paternal relatives (if they are both present in the network).
#'Vector also includes \code{from.id}.
#'
#' @examples
#' #Example data taken from kinship2::kinship()
#'test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#'                    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
#'                    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
#'                    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))
#'#some renaming
#'names(test1)[2] = "mother"
#'names(test1)[3] = "father"
#'test1$sex = ifelse(test1$sex ==1, "F", "M")

#'test1$id = as.character(test1$id)
#'test1$mother = ifelse(test1$mother!=0, as.character(test1$mother), "UNK")
#'test1$father = ifelse(test1$father!=0, as.character(test1$father), "UNK")
#'test1
#'
#'population.network = make_kinship_network(make_relation_df(test1))
#'plot(population.network)
#'get_nth_degree_relatives(network = population.network, from.id = test1$id[5],
#'r.degree = 2)

#'
#' @export

# Defintions.
# Indirect relatives: relatives who are not within the degree path length but are rather related to individuals to whom you are directly related. This defitionis is strictly orthodox but serves a purpose (it does not line up with lienal and collatoral realtives)
# e.g. Nibblings are 3 away, but sibs are within 2.

#### NOTE to make this only materernal or paternal sides will just need to add a wrapper function to remove the focals mother/father edge and then run the function.
get_nth_degree_relatives = function(network, from.id, r.degree){
  ids = igraph::V(network)$name

  # all.distances = distances(network, v = from.id, to = ids, mode = "all")
  # ids = ids[which(ids == colnames(all.distances))]
  # all.nth.distances.rs = ids[all.distances <= r.degree& all.distances >0]

  all.in.distances = igraph::distances(network, v = from.id, to = ids, mode = "in")
  ids = ids[which(ids == colnames(all.in.distances))]
  all.nth.in.distances.rs = ids[all.in.distances <= r.degree& all.in.distances >0]

  all.out.distances = igraph::distances(network, v = from.id, to = ids, mode = "out")
  ids = ids[which(ids == colnames(all.out.distances))]
  all.nth.out.distances.rs = ids[all.out.distances <= r.degree& all.out.distances >0]

  # within.range.relatives = all.nth.distances.rs #gets those within 'degree' steps (will get direct and those at close angles e.g. sibs)
  direct.relatives = c(all.nth.in.distances.rs, all.nth.out.distances.rs)
  #print(sort(within.range.relatives))

  indirect.relatives = indirect_relative_iterator(network, from.id, r.degree) # this gets those at an angle (e.g. aunts, nephews)
  #print(sort(indirect.relatives))


  all.relatives.within.degree = unique(c(direct.relatives, indirect.relatives)) # this is the list of relatives

  # not.ego.ids = ids[!(ids %in% c(all.relatives.within.degree, from.id))]
  # ego.network = delete_vertices(network, not.ego.ids)

  return(all.relatives.within.degree)
}

#' Find indirect relatives from a network
#'
indirect_relative_iterator = function(network, from.id, r.degree){
  allids = igraph::V(network)$name
  in.distances = igraph::distances(network, v = from.id, to = allids, mode = "in")

  n.ins.iterator = r.degree
  n.out.iterator = 1
  output = list()
  while(n.ins.iterator >0){
    output[[n.out.iterator]] = list()
    for(j in n.out.iterator:1){
      reordered.ids = allids[which(allids==colnames(in.distances))]
      ids.at.in = reordered.ids[in.distances == n.ins.iterator]
      if(length(ids.at.in) > 0){
        indirect.relatives = list()
        for(i in 1:length(ids.at.in)){
          in.id.distances = igraph::distances(network, v = ids.at.in[i], to = reordered.ids, mode = "out")
          reordered.ids = allids[which(allids==colnames(in.id.distances))]
          indirect.relatives[[i]] = reordered.ids[in.id.distances == j]
        }
        indirect.relatives = unlist(indirect.relatives)
      } else {
        indirect.relatives = NULL
      }

      output[[n.out.iterator]][[j]] = indirect.relatives
    }


    n.ins.iterator = n.ins.iterator -1
    n.out.iterator = n.out.iterator+1
  }
  return(unlist(output))
}

#' Plot an ego-pedigree
#'
#' \code{plot_ego_network} plots the ego network produced by the
#' \code{get_nth_degree_ego_network} function
#'
#' @param ego_network. igraph network object produced by the
#' \code{get_nth_degree_ego_network} function.
#' @param census.data Option to include a population.data dataframe which if
#' included will add more features to the plotted graph.
#'
#' @return A network plot. Colour representing lineage: orange = lineage A, blue
#' = lineage B, purple = both lineages, black = the ego. Shape represents sex,
#' circle = female, square = male.
#'
#'@examples
#'#Example data taken from kinship2::kinship()
#'test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#'                    mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
#'                    dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
#'                    sex =c(0, 1, 0, 1, 0, 1, 0, 1, 0,  0,  0,  1,  1,  1))
#'#some renaming
#'names(test1)[2] = "mother"
#'names(test1)[3] = "father"
#'test1$sex = ifelse(test1$sex ==1, "F", "M")

#'test1$id = as.character(test1$id)
#'test1$mother = ifelse(test1$mother!=0, as.character(test1$mother), "UNK")
#'test1$father = ifelse(test1$father!=0, as.character(test1$father), "UNK")
#'test1
#'
#'r.df = make_relation_df(test1)
#'net = make_kinship_network(r.df)
#'
#'ego.network = get_nth_degree_ego_network(network = population.network,
#'ego.id = test1$id[5],
#'r.degree = 2,
#'seperate.lineages = TRUE)
#'
#'plot_ego_network(ego.network, census.data = test1)
#'
#' @export

plot_ego_network = function(ego_network, census.data = NULL){
  igraph::V(ego_network)$colour = ifelse(igraph::V(ego_network)$lineage == "A", "darkorange1", "deepskyblue1")
  igraph::V(ego_network)$colour = ifelse(igraph::V(ego_network)$lineage == "focal", "black", igraph::V(ego_network)$colour)
  igraph::V(ego_network)$colour = ifelse(igraph::V(ego_network)$lineage == "both", "darkorchid1", igraph::V(ego_network)$colour)
  igraph::V(ego_network)$label.colour = "black"
  igraph::V(ego_network)$label.colour = ifelse(igraph::V(ego_network)$lineage == "focal", "white", igraph::V(ego_network)$label.colour)

  if(!is.null(census.data)){
    sex.df = data.frame(id = igraph::V(ego_network)$name)
    sex.df =suppressMessages(
      suppressWarnings(
        dplyr::left_join(sex.df, dplyr::select(census.data, id, sex))
      )
    ) # census data must contain columns "sex" and "id
    if("mother" %in% names(census.data)){
      sex.df$sex = ifelse(sex.df$id %in% census.data$mother, "F", sex.df$sex)}
    if("father" %in% names(census.data)){
      sex.df$sex = ifelse(sex.df$id %in% census.data$father, "M", sex.df$sex)}
    igraph::V(ego_network)$sex = sex.df$sex
  }

  if(is.null(census.data)){
    egoplot = igraph::plot.igraph(ego_network,
                   vertex.color = igraph::V(ego_network)$colour,
                   vertex.size = 5,
                   vertex.label.color = igraph::V(ego_network)$label.colour
    )
  } else {
    igraph::V(ego_network)$shape = rep.int("crectangle", igraph::vcount(ego_network))
    igraph::V(ego_network)$shape = ifelse(igraph::V(ego_network)$sex == "M" & !is.na(igraph::V(ego_network)$sex), "square", igraph::V(ego_network)$shape)
    igraph::V(ego_network)$shape = ifelse(igraph::V(ego_network)$sex == "F" & !is.na(igraph::V(ego_network)$sex), "circle", igraph::V(ego_network)$shape)
    egoplot = igraph::plot.igraph(ego_network,
                   vertex.color = igraph::V(ego_network)$colour,
                   vertex.size = 5,
                   vertex.label.color = igraph::V(ego_network)$label.colour
    )
  }

  return(egoplot)
}

