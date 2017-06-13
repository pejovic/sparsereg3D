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

sparsereg3D.ncv <- function(sparse.reg, lambda, w = NULL){
  
  # Extracting data from sparse.reg object
  seed <- sparse.reg$folds$seed
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
  kmean.vars = sparse.reg$model$kmean.vars
  cum.prop = sparse.reg$model$cum.prop
  mae <- function (actual, predicted) mean(ae(actual, predicted))
  ae <- function (actual, predicted) abs(actual-predicted)
  
  # Preparing empty data frames which will contain the results of procedure.
  test.prediction <- data.frame()
  accuracy <- list()
  data.prediction <- data.frame()
  models.ncv <- as.list(rep(NA,num.folds))
  
  # Outer loop of nested crossvalidation
  for(i in 1:length(profile.fold.list)){
    test.data <- profiles[which(profiles$ID %in% profile.fold.list[[i]]),]
    training.obs.ind <- which(profiles$ID %in% do.call(c, profile.fold.list[-i]))
    training.data <- profiles[training.obs.ind,]
    weight.data <- training.data[,c("mid.depth","hdepth")]
    
    # Inner crossvalidation partitioning
    tmp <- stratfold3d(target.name = target.name, other.names = kmean.vars , data = training.data, num.folds = num.folds, seed = seed, num.means = num.means, cum.prop = cum.prop)
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
      for(j in 1:length(inner.obs.fold.list)){
        inner.fold.indices[inner.obs.fold.list[[j]]] <- j
      }
## MLADEN      
      if(is.null(w)){
        grid <- expand.grid(1:num.folds, lambda)
        train.pred.df <- do.call(rbind,(apply(grid, 1, function(x) data.frame(obs = as.numeric(training.data[which(inner.fold.indices == x[1]), 1]), pred = as.numeric(predict(glmnet(as.matrix(training.data[-which(inner.fold.indices == x[1]), -1]), training.data[-which(inner.fold.indices == x[1]), 1], alpha = 1, lambda = x[2], weights = 1/(1 + 0*weight.data[-which(inner.fold.indices == x[1]),"hdepth"]/100 + 0*abs(weight.data[-which(inner.fold.indices == x[1]),"mid.depth"]))), newx = as.matrix(training.data[which(inner.fold.indices == x[1]), -1]), s = x[2]))))))
        train.pred.df <- cbind(train.pred.df, lambda = (c(rep(lambda, each = length(inner.fold.indices)))))
        train.pred.df$pred <- pmax(train.pred.df$pred, target.min/3)
        lambda.rmse <- ddply(train.pred.df, .(lambda), function(x) cv.error = RMSE(x$pred, x$obs))
        min.cv.error <- min(lambda.rmse[,2])
        lambda.min <- lambda.rmse[which(lambda.rmse[,2] == min.cv.error),"lambda"]
        lasso <- glmnet(as.matrix(training.data[, -1]), training.data[, 1], alpha = 1, lambda = lambda.min, weights = 1/(1 + 0*weight.data[,"hdepth"]/100 + 0*abs(weight.data[,"mid.depth"])))
        weight = NA
      }else{
        train.cv.errors <- matrix(NA, nrow = length(w), ncol = length(lambda))
        for(j in 1:length(w)){
          grid <- expand.grid(1:num.folds, lambda)
          train.pred.df <- do.call(rbind,(apply(grid, 1, function(x) data.frame(obs = as.numeric(training.data[which(inner.fold.indices == x[1]), 1]), pred = as.numeric(predict(glmnet(as.matrix(training.data[-which(inner.fold.indices == x[1]), -1]), training.data[-which(inner.fold.indices == x[1]), 1], alpha = 1, lambda = x[2], weights = 1/(1 + 0*weight.data[-which(inner.fold.indices == x[1]),"hdepth"]/100 + w[j]*abs(weight.data[-which(inner.fold.indices == x[1]),"mid.depth"]))), newx = as.matrix(training.data[which(inner.fold.indices == x[1]), -1]), s = x[2]))))))
          train.pred.df <- cbind(train.pred.df, lambda = (c(rep(lambda, each = length(inner.fold.indices)))))
          train.pred.df$pred <- pmax(train.pred.df$pred, target.min/3)
          lambda.rmse <- ddply(train.pred.df, .(lambda), function(x) cv.error = RMSE(x$pred, x$obs))
          train.cv.errors[j,] <- lambda.rmse[,2]
          }
        min.cv.error <- min(train.cv.errors)
        min.ind <- which(train.cv.errors == min(train.cv.errors), arr.ind=TRUE) # ova funkcija daje kolonu i red minimalne greske
        lambda.min <- sort(lambda, decreasing = FALSE)[min.ind[1,2]]
        lasso <- glmnet(as.matrix(training.data[, -1]), training.data[, 1], alpha = 1, lambda = lambda.min, weights = 1/(1 + 0*weight.data[,"hdepth"]/100 + w[min.ind[1,1]]*abs(weight.data[,"mid.depth"])))
        weight = w[min.ind[1,1]]
      }
## MLADEN      
      
      # Inner crossvalidation loop with model selection
      coef.list <- coef(lasso, s = lambda.min)
      models.ncv[[i]] <- as.list(c(coefficients = coef.list, lambda = lambda.min, weight = weight))
      # Prediction on test set
      test.pred <- predict(lasso, s = lambda.min, newx = as.matrix(test.data[,-1]))
      test.pred <- pmax(test.pred, target.min/3)
      
      training.pred <- predict(lasso, s = lambda.min, newx = as.matrix(training.data[,-1]))
      training.pred <- pmax(training.pred, target.min/3)
      
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
      
      # Extracting the coefficients of each model
      if(poly.deg == 1){ 
        int.coeff <- as.matrix(hier.path$th[,,which(hier.lasso.cv$lamhat==hier.path$lamlist)][,length(main.effect.names)]) 
      } else {
        int.coeff <- as.matrix(hier.path$th[,,which(hier.lasso.cv$lamhat==hier.path$lamlist)][,(length(main.effect.names)-poly.deg+1):length(main.effect.names)]) 
      }
      main.coeff <- hier.path$bp[,which(hier.lasso.cv$lamhat==hier.path$lamlist), drop = F] - hier.path$bn[,which(hier.lasso.cv$lamhat==hier.path$lamlist), drop = F]
      
      coef.list <- data.frame(cov.name=colnames(training.main.effects),main.coeff,int.coeff)
      
      models.ncv[[i]] <- coef.list
      
      # Prediction on test set
      test.pred <- predict(final.hier.model, newx = test.main.effects, newzz = test.int.effects) 
      test.pred <- pmax(test.pred,target.min/3)
      
      training.pred <- predict(final.hier.model, newx = training.main.effects, newzz = training.int.effects)
      training.pred <- pmax(training.pred, target.min/3)
    }
    
    # Assembling predictions for final model assessment
    test.obs.pred <- data.frame(obs = test.data[,target.name], pred = as.numeric(test.pred))
    test.prediction <- rbind(test.prediction,test.obs.pred)
    
    training.obs.pred <- data.frame(obs = training.data[,target.name], pred = as.numeric(training.pred))
    
    if(!use.hier){accuracy[[i]] <- c(test = defaultSummary(test.obs.pred)[1], test = defaultSummary(test.obs.pred)[2], RMSE.train = (min.cv.error))}else{accuracy[[i]] <- defaultSummary(test.obs.pred)}
    

    obs.pred <- data.frame(test.data, pred = as.numeric(test.pred))
    data.prediction <- rbind(data.prediction,obs.pred)
  }
  
  # Storing results
  out <- list(RMSE = defaultSummary(test.prediction)[1], Rsquared = defaultSummary(test.prediction)[2], NCV.Models = accuracy, Data.and.Prediction = data.prediction, models = models.ncv)
  return(out)
}
