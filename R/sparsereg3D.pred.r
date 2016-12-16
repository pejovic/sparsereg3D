#=========== sparsereg3D.pred - function for prediction on grids =======
#
# Arguments:
# model.info  - output from sparsereg3D.sel function
# depths      - depths at which prediction is performed
# grids       - SpatialPixelsDataframe containing grids of covariates
# chunk.size  - number of pixels in each chunk obtained by splitting grids in order to do parlallel processing
#
# Return value:
# grids with predicted values

#model.info <- sp.reg.som; depths = c(-0.1); grids = gridmaps.sm2D; chunk.size = 20000

sparsereg3D.pred <- function(model.info, depths, grids, chunk.size) {
  
  "%ni%" <- Negate("%in%")
  
  # Extracting data from model.info object
  std.par <- model.info$std.par
  target.name <- model.info$model$target.name
  main.effect.names <- model.info$model$main.effect.names
  depth.int.names <- model.info$model$depth.int.names
  use.interactions <- model.info$model$use.interactions
  use.hier <- model.info$model$use.hier
  poly.deg <- model.info$model$poly.deg
  base.model <- model.info$model$base.model
  
  # Creating 3list of grids, each corresponds to different depth
  grids.3D <- sp3D(grids, stdepths = depths) 
  
  # Number of cores has to be 1, because function mclapply does not work on many cores
  num.cores = 1 
  
  # Converting grids into data frames and renaming altitude column to depth column in each grid, besause function mclapply uses word altitude to denote depth
  covs.data <- mclapply(grids.3D, function(x) as.data.frame(x), mc.cores = num.cores)
  covs.data <- lapply(covs.data, function(x) {plyr::rename(x, replace=c("altitude" = "depth"))}) 
  
  # Adding polynomial depth variables in each grid
  if(poly.deg > 1) { 
    covs.data <- lapply(covs.data, function(x) x <- cbind(x, poly(x$depth, poly.deg, raw=TRUE, simple=TRUE)[,-1]))
    covs.data <- lapply(covs.data, function(x) {names(x) <- c(names(x)[1:(length(names(x))-(poly.deg))], c("depth",paste("depth",c(2:poly.deg),sep="")));return(x)})
  }
  
  # Applying standardization parameters for dummy coding to transform categorical variables in each grid
  covs.data <- mclapply(covs.data, function(x) subset(x, select = all.vars(base.model)[-1], drop=FALSE), mc.cores = num.cores) %>% mclapply(., function(x) predict(std.par$dummy.par, newdata=x),mc.cores = num.cores) %>% mclapply(., function(x) {colnames(x) <- gsub( "\\_|/|\\-|\"|\\s" , "." , colnames(x) ); return(x)}, mc.cores = num.cores) 
  
  # Splitting grids into chunks, computing interactions in parallel on these chunks
  if(use.interactions){
    n <- nrow(covs.data[[1]]) 
    r <- rep(1:ceiling(n/chunk.size),each = chunk.size)[1:n] 
    covs.int.data <- lapply(covs.data, function(x) as.data.frame(x)) %>% lapply(., function(x) split(x,r))
    
    m.cores <- detectCores() 
    registerDoParallel(cores = m.cores)
    
    for(i in 1:length(covs.int.data)){ 
      covs.int.data[[i]] <- foreach(j = covs.int.data[[i]],.combine="rbind") %dopar% {hierNet::compute.interactions.c(as.matrix(j),diagonal = FALSE)}
    }
  }
  
  # Combining main effect grids and grids with computed interactions
  for( i in 1:length(covs.int.data)) {
    covs.data[[i]] <- cbind(covs.data[[i]],covs.int.data[[i]])
  }
  
  # Applying standardization parameters to continual variables 
  covs.data <- mclapply(covs.data, function(x) predict(std.par$cnt.par, newdata=x), mc.cores = num.cores)   # %>% mclapply(.,function(x) subset(x[,which(colnames(x) %in% c(main.effect.names,depth.int.names))]))
  
  # In hierarchical setting, interactions other than with depth must be zero
  if(use.hier){
    for(i in 1:length(covs.data)) {
      covs.int.data[[i]][,colnames(covs.int.data[[i]]) %ni% depth.int.names ] <- 0
      covs.int.data[[i]][,depth.int.names] <- covs.data[[i]][,depth.int.names]
      covs.data[[i]] <- subset(covs.data[[i]][,which(colnames(covs.data[[i]]) %in% main.effect.names)])
    } 
  }else{
    covs.data <- mclapply(covs.data, function(x) subset(x[,which(colnames(x) %in% c(main.effect.names,depth.int.names))]))
  }
  
  
  #Final prediction
  for(i in 1:length(grids.3D)){
    if(!use.hier){
      grids.3D[[i]]$pred <- as.numeric(predict(model.info$model$model, s = model.info$model$model$lambda.min, newx = (covs.data[[i]])))
    }else{
      grids.3D[[i]]$pred <- as.numeric(predict(model.info$model$model, newx = covs.data[[i]], newzz = covs.int.data[[i]]))
    }
    grids.3D[[i]]$pred <- pmax(grids.3D[[i]]$pred, 0) # TODO Ako bolje radi sa 0, staviti 0
    grids.3D[[i]] <- grids.3D[[i]][,"pred"]
  }
  
  return(grids.3D)
}
