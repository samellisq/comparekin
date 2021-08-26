#small fix just to remove the 'no visible binding for global variable warning' which is caused by tidyverse
utils::globalVariables(c("id", "sex", "mother", "father"))


#' Check Population Data
#'
check_population_df = function(population.data){


  ##Title names
  if(!("id" %in% names(population.data))){
    stop("population.data must contain a column named id")
  }
  if(!("mother" %in% names(population.data))){
    stop("population.data must contain a column named mother,
         (even if maternal relatedness is not being calculated)")
  }
  if(!("father" %in% names(population.data))){
    stop("population.data must contain a column named father,
         (even if maternal relatedness is not being calculated)")
  }
  if(!("sex" %in% names(population.data))){
    stop("population.data must contain a column named sex")
  }


  ##Check data types

  ##actually don't need this code converts it itself
  # if(typeof(population.data$id) != "character"){
  #   stop("population.data column id must be a character vector")
  # }
  # if(typeof(population.data$mother) != "character"){
  #   stop("population.data column mother must be a character vector")
  # }
  # if(typeof(population.data$father) != "character"){
  #   stop("population.data column father must be a character vector")
  # }
  # if(typeof(population.data$sex) != "character"){
  #   stop("population.data column sex must be a character vector")
  # }

  ## Contents
  if(any(is.na(population.data$id))){
    stop("population.data column id must not contain any NAs")
  }
  if(any(is.na(population.data$mother))){
    stop("population.data column mother must not contain any NAs. Unknown
    mothers should be labelled as UNK)")
  }
  if(any(is.na(population.data$father))){
    stop("population.data column father must not contain any NAs. Unknown
    fathers should be labelled as UNK)")
  }
  if(any(is.na(population.data$sex))){
    stop("population.data column sex must not contain any NAs. Unknown
    sex should be labelled as UNK)")
  }
  acceptable.sex.labels = c("M", "F", "male", "female", "UNK")
  if(!any(population.data$sex %in% acceptable.sex.labels)){
    stop("population.data sex can only be M, F, male, female or UNK")
  }

}
