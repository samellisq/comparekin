#' Simplify data frame to pedigree information
#'
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
get_nth_degree_ego_network = function(network, ego.id, r.degree = 2, seperate.lineages = FALSE){

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



#' Find relatives from a pedigree-network

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
    egoplot = igraph::plot(ego_network,
                   vertex.color = igraph::V(ego_network)$colour,
                   vertex.size = 5,
                   vertex.label.color = igraph::V(ego_network)$label.colour
    )
  } else {
    igraph::V(ego_network)$shape = rep.int("crectangle", igraph::vcount(ego_network))
    igraph::V(ego_network)$shape = ifelse(igraph::V(ego_network)$sex == "M" & !is.na(igraph::V(ego_network)$sex), "square", igraph::V(ego_network)$shape)
    igraph::V(ego_network)$shape = ifelse(igraph::V(ego_network)$sex == "F" & !is.na(igraph::V(ego_network)$sex), "circle", igraph::V(ego_network)$shape)
    egoplot = igraph::plot(ego_network,
                   vertex.color = igraph::V(ego_network)$colour,
                   vertex.size = 5,
                   vertex.label.color = igraph::V(ego_network)$label.colour
    )
  }

  return(egoplot)
}

