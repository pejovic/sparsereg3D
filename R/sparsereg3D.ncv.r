#=========== sparsereg3D.ncv - function for model assessment via nested crossvalidation =======
#
# Arguments:
# sparse.reg - output from pre.sparsereg3D function
# lambda     - vector of regularization parameter values for lasso regression
# seed       - random number generator seed
#
# Return value:
# RMSE
# R squared


#pre.som <- pre.sparsereg3D(base.model = SOM.fun, use.hier = TRUE, profiles = bor.profs, cov.grids = gridmaps.sm2D)
#sparse.reg = pre.som
#lambda = seq(0,5,0.1)

sparsereg3D.ncv <- function(sparse.reg, lambda, seed = 321){
  
  # Extracting data from sparse.reg object
  profile.fold.list <- sparse.reg$folds$profile.fold.list 
  obs.fold.list <- sparse.reg$folds$obs.fold.list 
  profiles <- sparse.reg$profiles
  num.folds <- sparse.reg$num.folds
  num.means <- sparse.reg$num.means
  target.name <- all.vars(sparse.reg$model$base.model)[1]
  target.min <- min(profiles[,target.name]) # TODO Ako nula daje bolje rezultate, staviti nulu
  use.interactions = sparse.reg$model$use.interactions 
  use.hier = sparse.reg$model$use.hier
  main.effect.names = sparse.reg$model$main.effect.names 
  poly.deg = sparse.reg$model$poly.deg
  depth.int.names = sparse.reg$model$depth.int.names
  all.int.names = sparse.reg$model$all.int.names
  
  # Preparing empty data frames which will contain the results of procedure.
  test.prediction <- data.frame()
  
  # Outer loop of nested crossvalidation
  for(i in 1:length(profile.fold.list)){
    test.data <- profiles[which(profiles$ID %in% profile.fold.list[[i]]),]
    training.obs.ind <- which(profiles$ID %in% do.call(c, profile.fold.list[-i]))
    training.data <- profiles[training.obs.ind,]
    
    # Inner crossvalidation partitioning
    tmp <- stratfold3d(target.name, training.data, num.folds = num.folds, seed = seed, num.means = num.means)
    inner.profile.fold.list <- tmp$profile.fold.list
    inner.obs.fold.list <- tmp$obs.fold.list
    
    if(!use.hier){
      # Keep only columns relevant for model training
      if(use.interactions){
        training.data <- subset(training.data, select = c(target.name, main.effect.names, depth.int.names)) 
        test.data <- subset(test.data, select = c(target.name, main.effect.names, depth.int.names))
      } else {
        training.data <- subset(training.data, select = c(target.name, main.effect.names)) 
        test.data <- subset(test.data, select = c(target.name, main.effect.names))
      }
      
      # Convert list of vectors of indices per fold to vector of fold indices for each observation as glmnet requires
      # list(c(1,4),c(2,5),c(3,6)) -> c(1,2,3,1,2,3)
      inner.fold.indices <- rep(NA,dim(training.data)[1])
      for(i in 1:length(inner.obs.fold.list)){
        inner.fold.indices[inner.obs.fold.list[[i]]] <- i
      }
      
      # Inner crossvalidation loop with model selection
      lasso.cv <- cv.glmnet(as.matrix(training.data[,-1]), training.data[,1], alpha = 1,lambda = lambda, foldid = inner.fold.indices, type.measure = "mse")
      
      # Prediction on test set
      test.pred <- predict(lasso.cv, s = lasso.cv$lambda.min, newx = as.matrix(test.data[,-1]))
      test.pred <- pmax(test.pred, target.min/3)
    }else{
      # Hierarchical setting requires the separation of main effects and interaction effects
      training.main.effects <- as.matrix(training.data[,main.effect.names]) 
      test.main.effects <- as.matrix(test.data[,main.effect.names])
      training.int.effects <- as.matrix(training.data[,all.int.names])
      test.int.effects <- as.matrix(test.data[,all.int.names])
      training.target <- training.data[,target.name]
      test.target <- test.data[,target.name]
      
      # Inner crossvalidation loop with model selection
      hier.path <- hierNet.path(training.main.effects, training.target, zz = training.int.effects, diagonal = FALSE, strong = TRUE, trace = 0, stand.main = FALSE, stand.int = FALSE)
      hier.lasso.cv <- hierNet.cv(hier.path, training.main.effects, training.target, folds = inner.obs.fold.list, trace=0)
      final.hier.model <- hierNet(training.main.effects, training.target, zz = training.int.effects, diagonal=FALSE, strong=TRUE, lam = hier.lasso.cv$lamhat, center = TRUE, stand.main = FALSE, stand.int = FALSE)
      
      # Prediction on test set
      test.pred <- predict(final.hier.model, newx = test.main.effects, newzz = test.int.effects) 
      test.pred <- pmax(test.pred,target.min/3)
    }
    
    # Assembling predictions for final model assessment
    test.obs.pred <- data.frame(obs = test.data[,target.name], pred = as.numeric(test.pred))
    test.prediction <- rbind(test.prediction,test.obs.pred)
  }
  
  # Storing results
  out <- list(RMSE = defaultSummary(test.prediction)[1], Rsquared = defaultSummary(test.prediction)[2])
  return(out)
}
