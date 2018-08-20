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

#' Spatial prediction at multiple soil depths.
#'
#' @export
#'
#' @param model.info output from \code{sparsereg3D.sel} function
#' @param depths depths at which prediction is performed
#' @param grids SpatialPixelsDataframe containing grids of covariates
#' @param chunk.size number of pixels in each chunk obtained by splitting grids in order to do parlallel processing
#'
#'
#' @return List of prediction grids
#'
#'  @keywords Spatial prediction



sparsereg3D.pred <- function(model.info, depths, grids, chunk.size) {

  "%ni%" <- Negate("%in%")

  if(length(model.info) < 4){
    stop('Final model must be fitted by lasso')
  }

  # Extracting data from model.info object
  std.par <- model.info$std.par
  if(model.info$comp){
    target.name <- model.info$model$target.name
  }else{
    target.name <- all.vars(model.info$model$base.model)[1]
    target.min <- min(profiles[,target.name])
    use.hier <- model.info$model$use.hier
  }
  main.effect.names <- model.info$model$main.effect.names
  depth.int.names <- model.info$model$depth.int.names
  use.interactions <- model.info$model$use.interactions
  poly.deg <- model.info$model$poly.deg
  base.model <- model.info$model$base.model

  grids <- grids[complete.cases(grids@data[,names(grids@data)]), names(grids@data)]

  # creating list of classes of each categorical variables that were used in training
  factor.lvls <- std.par$dummy.par$lvls

  # identifying of classes that exist in cov.maps but not used in training
  for(i in names(std.par$dummy.par$lvls)){
    factor.lvls[[i]] <- (levels(grids@data[,i])[(levels(grids@data[,i]) %ni% std.par$dummy.par$lvls[[i]])])
  }

  # identifying which categorical variables has such classes
  ind <- as.numeric(which(lapply(factor.lvls, function(x) length(x)) != 0))

  # masking of such classes by the second most dominant class
  for(j in names(std.par$dummy.par$lvls[ind])){
    if((sum(as.character(grids@data[,j]) %in% factor.lvls[[j]]) > 0)){
      #grids@data[,j] <- as.character(grids@data[,j])
      grids@data[as.character(grids@data[,j]) %in% factor.lvls[[j]], j] <- attr(sort(summary(grids@data[,j]), decreasing = TRUE)[2], "names")
      grids@data[,j] <- factor(grids@data[,j])
    }
  }

  # Creating 3list of grids, each corresponds to different depth
  grids.3D <- sp3D(grids, stdepths = depths)

  # Number of cores has to be 1, because function mclapply does not work on many cores
  num.cores = 1

  # Converting grids into data frames and renaming altitude column to depth column in each grid, besause function mclapply uses word altitude to denote depth
  covs.data <- mclapply(grids.3D, function(x) as.data.frame(x), mc.cores = num.cores)
  coord.names <- tail(colnames(covs.data[[1]]),3); coord.names.ind <- which(names(covs.data[[1]]) %in% coord.names)
  covs.data <- lapply(covs.data, function(x) {names(x)[coord.names.ind] <- c("x", "y", "depth"); x})

  # Adding polynomial depth variables in each grid
  if(poly.deg > 1) {
    covs.data <- lapply(covs.data, function(x) x <- cbind(x, poly(x$depth, poly.deg, raw=TRUE, simple=TRUE)[,-1]))
    covs.data <- lapply(covs.data, function(x) {names(x) <- c(names(x)[1:(length(names(x))-(poly.deg))], c("depth",paste("depth",c(2:poly.deg),sep="")));return(x)})
    covs.data <- lapply(covs.data, function(x) {x <- x[complete.cases(x),]; return(x)})
    covs.data <- lapply(covs.data, function(x) {catcolwise(as.factor); return(x)})
  }

  if(model.info$comp){
    # Applying standardization parameters for dummy coding to transform categorical variables in each grid
    covs.data <- mclapply(covs.data, function(x) subset(x, select = all.vars(base.model), drop=FALSE), mc.cores = num.cores) %>% mclapply(., function(x) predict(std.par$dummy.par, newdata=x),mc.cores = num.cores) %>% mclapply(., function(x) {colnames(x) <- gsub( "\\_|/|\\-|\"|\\s" , "." , colnames(x) ); return(x)}, mc.cores = num.cores)
  }else{
    covs.data <- mclapply(covs.data, function(x) subset(x, select = all.vars(base.model)[-1], drop=FALSE), mc.cores = num.cores) %>% mclapply(., function(x) predict(std.par$dummy.par, newdata=x),mc.cores = num.cores) %>% mclapply(., function(x) {colnames(x) <- gsub( "\\_|/|\\-|\"|\\s" , "." , colnames(x) ); return(x)}, mc.cores = num.cores)
  }

  # Splitting grids into chunks, computing interactions in parallel on these chunks
  if(use.interactions){
    n <- nrow(covs.data[[1]])
    chunk.size <- round_any(n/10,1000)
    r <- rep(1:ceiling(n/chunk.size),each = chunk.size)[1:n]
    covs.int.data <- lapply(covs.data, function(x) as.data.frame(x)) %>% lapply(., function(x) split(x,r))

    m.cores <- detectCores()
    registerDoParallel(cores = m.cores)


    f <- as.formula(~ .^2)
    #interactions <- hierNet::compute.interactions.c(as.matrix(profiles[,-c(which(colnames(profiles) %in% c("ID",target.name,"hdepth",coord.names)))]),diagonal=FALSE)

      for(i in 1:length(covs.int.data)){
         covs.int.data[[i]] <- foreach(j = 1:length(covs.int.data[[i]]),.combine="rbind") %dopar% {model.matrix(f, (covs.int.data[[i]][[j]]))}
         covs.int.data[[i]] <- covs.int.data[[i]][,-which(colnames(covs.int.data[[i]]) %in% c("(Intercept)", main.effect.names))]

         #covs.int.data[[i]] <- foreach(j = 1:length(covs.int.data[[i]]),.combine="rbind") %dopar% {hierNet::compute.interactions.c(as.matrix(covs.int.data[[i]][[j]]),diagonal = FALSE)}
      }
  }
  stopImplicitCluster()
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

  registerDoParallel(cores = m.cores)
  #Final prediction
  for(i in 1:length(grids.3D)) {
    if(!use.hier){
      if(model.info$comp){
        comp.pred <- as.numeric(predict(model.info$model$model, s = model.info$model$model$lambda.min, newx = (covs.data[[i]]), type = "response")[,,1])
        comp.pred <- data.frame(vars = rep(model.info$model$target.name, each = length(comp.pred)/length(model.info$model$target.name)), prediction = comp.pred)
        comp.pred <- comp.pred %>% group_by(vars) %>% mutate(grouped_id = row_number()) %>% spread(vars, prediction) %>% dplyr::select(-grouped_id)
        grids.3D[[i]]@data <- cbind(grids.3D[[i]]@data, comp.pred)
      }else{
        grids.3D[[i]]$pred <- as.numeric(predict(model.info$model$model, s = model.info$model$model$lambda.min, newx = (covs.data[[i]])))
      }

    }else{
      grids.3D[[i]]$pred <- as.numeric(predict(model.info$model$model, newx = covs.data[[i]], newzz = covs.int.data[[i]]))
    }
    grids.3D[[i]]$pred <- pmax(grids.3D[[i]]$pred, 0) # TODO Ako bolje radi sa 0, staviti 0
    grids.3D[[i]] <- grids.3D[[i]][,"pred"]
  }

  stopImplicitCluster()
  return(grids.3D)
}
