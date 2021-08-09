###############################################
# Getting the ego networks
#############################################

require(igraph)


make_relation_df = function(population.data, parents.to.use = "all.avaliable"){
  if(parents.to.use == "all.avaliable"){
    if("mother" %in% names(population.data) & "father" %in% names(population.data)){
      relation.df = select(population.data, id, mother, father)
    } 
    if("mother" %in% names(population.data) & !("father" %in% names(population.data))){
      relation.df = select(population.data, id, mother)
    }
    if("father" %in% names(population.data) & !("mother" %in% names(population.data))){
      relation.df = select(population.data, id, father)
    }
    if(!("mother" %in% names(population.data)) & !("father" %in% names(population.data))){
      stop("Population.data must contain columns named one or both of 'mother' and 'father'")
    }
  }
  
  if(parents.to.use == "mother.only"){
    if("mother" %in%  names(population.data)){
      relation.df = select(population.data, id, mother)
    } else {
      stop("population.df must contain a column called 'mother'")
    }

  }
  
  if(parents.to.use == "father.only"){
    if("father" %in%  names(population.data)){
      relation.df = select(population.data, id, father)
    } else {
      stop("population.df must contain a column called 'father'")
    }
    
  }
  
  return(relation.df)
}

#Make kinship network
#relation.data.frame = data frame with 2 or 3 columns. column names must be 'id' + either/or both of 'mother', 'father'
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

  
  network = graph_from_adjacency_matrix(adj.mat)
  
  return(network)
  
}




get_nth_degree_ego_network = function(network, ego.id, r.degree = 2, seperate.lineages = FALSE){
  
  ego.in.edges = incident(network, ego.id, mode = "in")
  
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
      Xlineal.net = delete_edges(network, ego.in.edges[i])
      relatives[[i]] = get_nth_degree_relatives(Xlineal.net, ego.id, r.degree)
    }
  }

  all.relatives = unique(unlist(relatives))
  ids = V(network)$name
  not.ego.ids = ids[!(ids %in% c(all.relatives, ego.id))]
  ego.network = delete_vertices(network, not.ego.ids)
  
  
  if(seperate.lineages == TRUE){
    V(ego.network)$lineage = ifelse(V(ego.network)$name %in% relatives[[1]], "A", "B")
    V(ego.network)$lineage = ifelse(V(ego.network)$name %in% relatives[[1]] & V(ego.network)$name %in% relatives[[2]], "both", V(ego.network)$lineage)
  } else {
    V(ego.network)$lineage = rep.int("A", vcount(ego.network))
  }
  V(ego.network)$lineage = ifelse(V(ego.network)$name == ego.id, "focal", V(ego.network)$lineage)

  
  return(ego.network)
}



# Defintions. 
# Indirect relatives: relatives who are not within the degree path length but are rather related to individuals to whom you are directly related. This defitionis is strictly orthodox but serves a purpose (it does not line up with lienal and collatoral realtives)
# e.g. Nibblings are 3 away, but sibs are within 2. 

#### NOTE to make this only materernal or paternal sides will just need to add a wrapper function to remove the focals mother/father edge and then run the function. 
get_nth_degree_relatives = function(network, from.id, r.degree){
  ids = V(network)$name
  
  # all.distances = distances(network, v = from.id, to = ids, mode = "all")
  # ids = ids[which(ids == colnames(all.distances))]
  # all.nth.distances.rs = ids[all.distances <= r.degree& all.distances >0]
  
  all.in.distances = distances(network, v = from.id, to = ids, mode = "in")
  ids = ids[which(ids == colnames(all.in.distances))]
  all.nth.in.distances.rs = ids[all.in.distances <= r.degree& all.in.distances >0]
  
  all.out.distances = distances(network, v = from.id, to = ids, mode = "out")
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


indirect_relative_iterator = function(network, from.id, r.degree){
  allids = V(network)$name
  in.distances = distances(network, v = from.id, to = allids, mode = "in")
  
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
          in.id.distances = distances(network, v = ids.at.in[i], to = reordered.ids, mode = "out")
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

plot_ego_network = function(ego_network, census.data = NULL){ 
  V(ego_network)$colour = ifelse(V(ego_network)$lineage == "A", "darkorange1", "deepskyblue1")
  V(ego_network)$colour = ifelse(V(ego_network)$lineage == "focal", "black", V(ego_network)$colour)
  V(ego_network)$colour = ifelse(V(ego_network)$lineage == "both", "darkorchid1", V(ego_network)$colour)
  V(ego_network)$label.colour = "black"
  V(ego_network)$label.colour = ifelse(V(ego_network)$lineage == "focal", "white", V(ego_network)$label.colour)
  
  if(!is.null(census.data)){
    sex.df = data.frame(id = V(ego_network)$name)
    sex.df =suppressMessages(
      suppressWarnings(
        left_join(sex.df, select(census.data, id, sex))
        )
      ) # census data must contain columns "sex" and "id
    if("mother" %in% names(census.data)){
    sex.df$sex = ifelse(sex.df$id %in% census.data$mother, "F", sex.df$sex)}
    if("father" %in% names(census.data)){
    sex.df$sex = ifelse(sex.df$id %in% census.data$father, "M", sex.df$sex)}
    V(ego_network)$sex = sex.df$sex
  }
  
  if(is.null(census.data)){
    egoplot = plot(ego_network,
                   vertex.color = V(ego_network)$colour,
                   vertex.size = 5,
                   vertex.label.color = V(ego_network)$label.colour
    )
  } else {
    V(ego_network)$shape = rep.int("crectangle", vcount(ego_network))
    V(ego_network)$shape = ifelse(V(ego_network)$sex == "M" & !is.na(V(ego_network)$sex), "square", V(ego_network)$shape)
    V(ego_network)$shape = ifelse(V(ego_network)$sex == "F" & !is.na(V(ego_network)$sex), "circle", V(ego_network)$shape)
    egoplot = plot(ego_network,
                   vertex.color = V(ego_network)$colour,
                   vertex.size = 5,
                   vertex.label.color = V(ego_network)$label.colour
    )
  }

  return(egoplot)
}


###########################################
# Calculating kinship
#############################################
require(kinship2)


get_relatedness = function(from.id, to.id, fromid.ego.network, population.data, r.type = "coef.of.relationship"){
  if(!(r.type %in% c("coef.of.relationship", "kinship.coef" ))){
    stop("r.type must be either 'coef.of.relationship' or 'kinship.coef'")
  }

  ids.in.range = V(fromid.ego.network)$name
  ids.population.data = filter(population.data, id %in% ids.in.range)
  
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
  
  ids.population.data = rbind(select(ids.population.data, id, sex, mother, father), unk.parents.df)
  ids.population.data$sex = ifelse(ids.population.data$sex == "M", "male", ids.population.data$sex)
  ids.population.data$sex = ifelse(ids.population.data$sex == "F", "female", ids.population.data$sex)
  ids.population.data$sex = ifelse(ids.population.data$sex == "UNK", "unknown", ids.population.data$sex)
  
  
  # calculate all kinships
  ped = pedigree(id = ids.population.data$id, dadid = ids.population.data$father, momid = ids.population.data$mother, sex = ids.population.data$sex)
  if(r.type == "coef.of.relationship"){
    full.kin.mat = kinship(ped)*2 # kinship automatically calcualted kinship coefficent
  } else {
    full.kin.mat = kinship(ped)
  }

  return(full.kin.mat[to.id,from.id])

}
 


######################################################
# Dealing with missing nodes
#############################################################

is_complete = function(ego.network, focal, r.degree, n.lineages = 2){
  if(!(n.lineages %in% c(1,2))){
    stop("Too many or too few lineages. Must be 1 or 2.")
  }
  
  i = 1
  complete = TRUE
  current.level.ids = focal
  
  if(n.lineages == 1){
    while(i <= r.degree & complete == TRUE){
      in.neighbours = adjacent_vertices(ego.network, current.level.ids, mode = "in")
      in.neighbours = V(ego.network)$name[unlist(in.neighbours)]
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
      in.neighbours = adjacent_vertices(ego.network, current.level.ids, mode = "in")
      in.neighbours = V(ego.network)$name[unlist(in.neighbours)]
      if(length(in.neighbours) > n.lineages*i){stop(paste("One of ",current.level.ids ," has too many parents", sep = ""))}
      if(length(in.neighbours)== n.lineages*i){
        i = i+1
        current.level.ids = in.neighbours
      } else {
        complete = FALSE
      }
    } # close while
  } # close if
  
  
  return(complete)
}

##########################
# Overall
#############################
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
  
  if(!(from.id %in% V(population.network)$name)){
    stop("from.id is not present in this population network")
  }
  if(!(to.id %in% V(population.network)$name)){
    stop("to.id is not present in this population network")
  }
  
  if(is.null(fromid.ego.network)){ 
    fromid.ego.network = get_nth_degree_ego_network(network = population.network,
                                                    ego.id = from.id, 
                                                    r.degree = r.degree)
  }

  
  

  #Is to.id in from.id's ego network? ie are they related in the Xth degree according to the data
  if(to.id %in% V(fromid.ego.network)$name) {
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


########################################
# Next Step out: Running several together
#########################################

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



##################################################################
# Get local relatedness
#################################################################
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
