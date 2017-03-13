#============ stratfold3d function - data stratification and partitioning for crossvalidation ============
#
# Creates stratified folds in three steps:
#   1. Profiles are clustered using k-means clustering according to spatial location
#   2. Each cluster is split to folds, stratified according to profile depth and weighted mean of observed target variable in the profile
#   3. Each final crossvalidation fold is obtained by merging one fold from each cluster
#
# Arguments:
# data        - input data matrix with "ID", "hdepth" and coordinate columns
# target.name - name of target variable
# num.folds   - number of folds
# num.means   - number of centers for k-means clustering
#
# Return value:
# data              - input data plus column of fold indices
# profile.fold.list - list containing a vector of profile IDs that constitute each fold
# obs.fold.list     - list containing a vector of observation IDs that constitute each fold

#target.name = target.name; other.names = main.effect.names; seed = seed; data = training.data; num.folds = num.folds; num.means = num.means 

stratfold3d <- function(target.name, other.names, data, num.folds, num.means, seed) {
  
  # Profiles are clustered using k-means clustering according to spatial location
  target.data <- ddply(data, .(ID), here(summarize), target = weighted.mean(eval(parse(text=target.name)), hdepth))
  other.data <- ddply(data[,c("ID", other.names)], .(ID), colwise(max))
  clustering.data <- cbind(target.data, other.data[,-1])
  clustering.data <- clustering.data[complete.cases(clustering.data),]
  set.seed(seed)
  km <- kmeans(clustering.data[,-1],  centers = num.means)
  clustering.data$cluster <- as.factor(km$cluster)
  
  # Creating empty list to contain stratified folds
  cluster.list <- as.list(rep(NA,length(unique(clustering.data$cluster))))
  names(cluster.list) <- paste("cluster",c(1:length(cluster.list)),sep="")
  
  # Each cluster is split to folds, stratified according to profile depth and weighted mean of observed target variable in the profile
  for(i in 1:length(cluster.list)){
    set.seed(seed)
    cluster.list[[i]] <- createFolds(clustering.data[which(clustering.data$cluster == levels(clustering.data$cluster)[i]),"target"], k = num.folds)
    if(length(cluster.list[[i]]) < num.folds){stop(paste("There is no enough data in cluster", i, sep = " "))}
    for(j in 1:num.folds){
      cluster.list[[i]][[j]] <- clustering.data[which(clustering.data$cluster == levels(clustering.data$cluster)[i]),"ID"][cluster.list[[i]][[j]]]
    }
  }
  
  # List containing a vector of profile IDs that constitute each fold
  profile.fold.list <- as.list(rep(NA,num.folds))
  names(profile.fold.list) <- paste("fold", c(1:num.folds), sep = "")
  for(i in 1:num.folds){
    profile.fold.list[[i]] <- do.call(c,lapply(cluster.list,function(x) x[[i]]))
    names(profile.fold.list[[i]]) <- NULL
  }
  
  # List containing a vector of observation IDs that constitute each fold
  obs.fold.list <- as.list(rep(NA,num.folds))
  names(obs.fold.list) <- paste("fold",c(1:num.folds),sep = "")
  for(i in 1:num.folds){
    obs.fold.list[[i]] <- which(data$ID %in% profile.fold.list[[i]])
  }
  
  # Each row is augmented by fold number
  data.with.foldIDs <- data.frame()
  for(i in 1:length(obs.fold.list)){
    tmp <- data[obs.fold.list[[i]],]
    tmp$fold <- paste("fold",i,sep="")
    data.with.foldIDs <- rbind(tmp,data.with.foldIDs)
  }
  data <- data.with.foldIDs  
  data$fold <- factor(data$fold)
  
  out <- list(data = data, profile.fold.list = profile.fold.list, obs.fold.list = obs.fold.list)
  return(out)
}
