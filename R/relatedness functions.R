#' Get pairwise relatedness from an ego-pedigree

get_relatedness = function(from.id, to.id, fromid.ego.network, population.data, r.type = "coef.of.relationship"){

  population.data$id = as.character(population.data$id)
  population.data$sex = as.character(population.data$sex)
  population.data$mother = as.character(population.data$mother)
  population.data$father = as.character(population.data$father)

  if(!(r.type %in% c("coef.of.relationship", "kinship.coef" ))){
    stop("r.type must be either 'coef.of.relationship' or 'kinship.coef'")
  }

  ids.in.range = igraph::V(fromid.ego.network)$name
  ids.population.data = dplyr::filter(population.data, id %in% ids.in.range)

  # Add mothers and fathers where missing
  if("mother" %in% names(population.data) & !("father" %in% names(population.data))){
    ids.population.data$father = paste("Dadid" ,ids.population.data$id, sep =".")
  }

  if("father" %in% names(population.data) & !("mother" %in% names(population.data))){
    ids.population.data$mother = paste("Mumid" ,ids.population.data$id, sep =".")
  }

  ids.population.data$mother =
    ifelse(ids.population.data$mother == "UNK", paste("Mumid", ids.population.data$id, sep = "."), ids.population.data$mother)
  ids.population.data$father =
    ifelse(ids.population.data$father == "UNK", paste("Dadid", ids.population.data$id, sep = "."), ids.population.data$father)

  # add dummy rows for the fake mothers and fathers and those at the ends of the degree distribution
  unk.mothers = unique(ids.population.data$mother[!(ids.population.data$mother %in% ids.population.data$id)])
  unk.fathers = unique(ids.population.data$father[!(ids.population.data$father %in% ids.population.data$id)])

  unk.parents.df = data.frame(
    id = c(unk.mothers, unk.fathers),
    sex = c(rep.int("F", length(unk.mothers)), rep.int("M", length(unk.fathers))),
    mother = NA,
    father = NA
  )

  ids.population.data = rbind(dplyr::select(ids.population.data, id, sex, mother, father), unk.parents.df)
  ids.population.data$sex = ifelse(ids.population.data$sex == "M", "male", ids.population.data$sex)
  ids.population.data$sex = ifelse(ids.population.data$sex == "F", "female", ids.population.data$sex)
  ids.population.data$sex = ifelse(ids.population.data$sex == "UNK", "unknown", ids.population.data$sex)


  # calculate all kinships
  ped = kinship2::pedigree(id = ids.population.data$id, dadid = ids.population.data$father, momid = ids.population.data$mother, sex = ids.population.data$sex)
  if(r.type == "coef.of.relationship"){
    full.kin.mat = kinship2::kinship(ped)*2 # kinship automatically calcualted kinship coefficent
  } else {
    full.kin.mat = kinship2::kinship(ped)
  }

  return(full.kin.mat[to.id,from.id])

}



#' Tests if an ego network is complete

is_complete = function(ego.network, focal, r.degree, n.lineages = 2){
  if(!(n.lineages %in% c(1,2))){
    stop("Too many or too few lineages. Must be 1 or 2.")
  }

  i = 1
  complete = TRUE
  current.level.ids = focal

  if(n.lineages == 1){
    while(i <= r.degree & complete == TRUE){
      in.neighbours = igraph::adjacent_vertices(ego.network, current.level.ids, mode = "in")
      in.neighbours = igraph::V(ego.network)$name[unlist(in.neighbours)]
      if(length(in.neighbours) > n.lineages){stop(paste("One of ",current.level.ids ," has too many parents", sep = ""))}
      if(length(in.neighbours)== n.lineages){ # no exponetial growth in number of relatives if only one lienage, hence seperation
        i = i+1
        current.level.ids = in.neighbours
      } else {
        complete = FALSE
      }
    } # close while
  } # close if

  if(n.lineages == 2){
    while(i <= r.degree & complete == TRUE){
      in.neighbours = igraph::adjacent_vertices(ego.network, current.level.ids, mode = "in")
      in.neighbours = igraph::V(ego.network)$name[unlist(in.neighbours)]
      if(length(in.neighbours) > n.lineages^i){stop(paste("One of ",current.level.ids ," has too many parents", sep = ""))}
      if(length(in.neighbours)== n.lineages^i){
        i = i+1
        current.level.ids = in.neighbours
      } else {
        complete = FALSE
      }
    } # close while
  } # close if


  return(complete)
}



#' Calculate pairwise relatedness

calculate_pairwise_relatendess = function(from.id,
                                          to.id,
                                          population.data,
                                          r.degree = 2,
                                          population.network = NULL,
                                          fromid.ego.network = NULL,
                                          toid.ego.network = NULL,
                                          r.type = "coef.of.relationship",
                                          plot.ego.networks = FALSE
                                          ){


  if(!(from.id %in% population.data$id)){
    stop("from.id is not present in this population data")
  }
  if(!(to.id %in% population.data$id)){
    stop("to.id is not present in this population data")
  }


  if(is.null(population.network)){
    population.network = make_kinship_network(make_relation_df(population.data))
  }

  if(!(from.id %in% igraph::V(population.network)$name)){
    stop("from.id is not present in this population network")
  }
  if(!(to.id %in% igraph::V(population.network)$name)){
    stop("to.id is not present in this population network")
  }

  if(is.null(fromid.ego.network)){
    fromid.ego.network = get_nth_degree_ego_network(network = population.network,
                                                    ego.id = from.id,
                                                    r.degree = r.degree)
  }




  #Is to.id in from.id's ego network? ie are they related in the Xth degree according to the data
  if(to.id %in% igraph::V(fromid.ego.network)$name) {
    # If they are related
    relatedness = get_relatedness(from.id, to.id, fromid.ego.network, population.data, r.type = r.type)

    if(plot.ego.networks == TRUE){
      plot_ego_network(fromid.ego.network, population.data)
    }


  } else { # if they are not related


    if(is.null(toid.ego.network)){
      toid.ego.network = get_nth_degree_ego_network(network = population.network,
                                                    ego.id = to.id,
                                                    r.degree = r.degree)
    }


    if("mother" %in% names(population.data) & "father" %in% names(population.data)){
      n.lineages = 2
    } else{
      n.lineages = 1
    }


    fromid.complete = is_complete(fromid.ego.network, focal = from.id, r.degree = r.degree, n.lineages = n.lineages)
    toid.complete = is_complete(toid.ego.network, focal = to.id, r.degree = r.degree, n.lineages = n.lineages)

    if(plot.ego.networks == TRUE){
      par(mfrow = c(2,1))
      plot_ego_network(fromid.ego.network, population.data)
      plot_ego_network(toid.ego.network, population.data)
      par(mfrow = c(1,1))
    }


    if(fromid.complete == TRUE & toid.complete == TRUE){ # if both networks are complete we know they are not related
      relatedness = 0
    } else {
      # if we wanted a bit that said "<0.5" here is where it would go.
      relatedness = NA
    }


  }

  return(relatedness)


}


#' Get pairwise relatedness between multiple individuals

get__pairwise_relatednesses_matrix = function(ids = "ALL",# the ids of the individuals to calcualte pairwise relatedness between
                                          population.data,
                                          parents.to.use = "all.avaliable",
                                          r.degree = 2,
                                          kinship.network = NULL, # recommend leave this as null but might speed up really big populations
                                          r.type = "coef.of.relationship"){

  if(ids[1] == "ALL"){
    ids = population.data$id
  }

  if((is.null(kinship.network))){
    #get the over all network
    relation.df = make_relation_df(population.data, parents.to.use = parents.to.use)
    kinship.network = make_kinship_network(relation.df)
  }

  message(paste("Creating ", length(ids), " pedigree ego networks", sep = ""))
  all.ego.networks = lapply(as.list(ids), get_nth_degree_ego_network, network = kinship.network, r.degree = r.degree)


  r.matrix = matrix (100, length(ids), length(ids))
  rownames(r.matrix) = ids
  colnames(r.matrix) = ids

  j = 0
  message(paste("Calculating all ", (ceiling(length(ids)^2)/2) ," pairwise relatednesses. This can take some time if the number of individuals is large", sep= ""))
  for(i in 1:length(ids)){
    j = i
    while(j <= length(ids)){
      r.matrix[i,j] = calculate_pairwise_relatendess(from.id = rownames(r.matrix)[i],
                                                     to.id = colnames(r.matrix)[j],
                                                     population.data = population.data,
                                                     r.degree = r.degree,
                                                     population.network = kinship.network,
                                                     r.type = r.type,
                                                     plot.ego.networks = FALSE,
                                                     fromid.ego.network = all.ego.networks[[i]],
                                                     toid.ego.network = all.ego.networks[[j]]
                                                     )
      j= j+1
    }
  }

  r.matrix = Matrix::forceSymmetric(r.matrix, uplo = "U")

  if( sum(r.matrix == 100, na.rm = TRUE) >0){
    warning("not all pairwsie relatednesses were calcualated. Uncalculated rs returned as 100. Calculate pairs manually?")
  }

  #r.matrix = as.numeric(r.matrix)
  return(r.matrix)
}



#' Get local relatedness of an individual from a kinship matrix
#'
#' Might not include in the final model
get_local_relatedness = function(id, kinship.matrix, group.members = NULL, local.r.to = "all"){
  id.row = kinship.matrix[,id]
  id.row = id.row[names(id.row) != id]

  if(is.null(group.members)){
    group.members = names(id.row)
  }


  id.row = id.row[names(id.row) %in% group.members]

  local.r.df =   data.frame(local.r = mean(id.row, na.rm = TRUE),
                            n.r.known = sum(!is.na(id.row)),
                            n.r.NA = sum(is.na(id.row))
  )

  local.r.df$local.r = ifelse(is.nan(local.r.df$local.r), NA, local.r.df$local.r)

  return(local.r.df)

}
