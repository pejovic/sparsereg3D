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

#' Model Evaluation with Nested Cross-Validation
#'
#' @export
#'
#' @param sparse.reg output from \code{pre.sparsereg3D} function
#' @param lambda vector of regularization parameter values for lasso regression
#' @param ols logical. If TRUE model will be fitted with OLS insted of using lasso
#' @param step logical. If TRUE stepwise procedure will be used when fitting OLS model
#' @param seed random number generator
#' @param all logical. If TRUE a detailed output will be prepared.
#' @param lambda.1se logical. If TRUE one sigma lambda rule will be used (largest lambda value with cv.err less than or equal to min(cv.err)+ SE).
#' @param w weighted parameter (positive number) that controls the level of using observations from deeper layers as less informative. General weighted model is $w=1/(1+w*depth)$. This option is still under development, so it is not included in the function for model selection.
#'
#'
#' @return List of objects including:
#' \itemize{
#'  \item \code{RMSE}:  Root Mean Squared Error
#'  \item \code{R squared}:   Coefficient of determination
#'  \item Accuracy (if \code{all = TRUE}). Accuracy measures computed for each step in outer cross-valudation loop within nested cross-validation.
#'  \item data.prediction (if \code{all = TRUE}) Data frame containing design matrix extended with predictions obtain through nested cross-validation.
#'  \item models (if \code{all = TRUE}) list of length \code{nfolds} containing models that corresponds to each step in outer cross-valudation loop within nested cross-validation.
#' }
#'
#'  @keywords Model evaluation

sparsereg3D.ncv <- function(sparse.reg, lambda, step = FALSE, ols = FALSE, all = FALSE, lambda.1se = FALSE, w = NULL){

  if(step){
    ols = TRUE
  }

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
      if(!ols){
        if(is.null(w)){
          train.cv <- cv.glmnet(as.matrix(training.data[,-1]), training.data[,1], alpha = 1, lambda = lambda, foldid = inner.fold.indices, type.measure = "mse", weights = 1/(1 + 0*weight.data[,"hdepth"]/100 + 0*abs(weight.data[,"mid.depth"]))) #
          lasso <- train.cv$glmnet.fit
          if(!lambda.1se){
            lambda.min <- train.cv$lambda.min
          }else{
            lambda.min <- train.cv$lambda.1se
          }
          min.cv.error <- min(train.cv$cvm)
          weight = NA
        }else{
          train.cv.errors <- matrix(NA, nrow = length(w), ncol = length(lambda))
          for(j in 1:length(w)){
            train.cv.errors[j,] <- cv.glmnet(as.matrix(training.data[,-1]), training.data[,1], alpha = 1,lambda = lambda, foldid = inner.fold.indices, type.measure = "mse", weights = 1/(1 + 0*weight.data[,"hdepth"]/100 + (w[j])*abs(weight.data[,"mid.depth"])))$cvm #
          }
          min.cv.error <- min(train.cv.errors)
          min.ind <- which(train.cv.errors == min(train.cv.errors), arr.ind=TRUE) # ova funkcija daje kolonu i red minimalne greske
          lambda.min <- sort(lambda, decreasing = TRUE)[min.ind[2]]
          train.cv <- cv.glmnet(as.matrix(training.data[,-1]), training.data[,1], alpha = 1, lambda = lambda, foldid = inner.fold.indices, type.measure = "mse", weights = 1/(1 + 0*weight.data[,"hdepth"]/100 + (w[min.ind[1]])*abs(weight.data[,"mid.depth"])) )
          lasso <- train.cv$glmnet.fit # Ovde sam opet koristio cv.glmnet i ako sam trebao samo glmnet, samo iz ocaja...da bude potpuno isto kao i za slucaj bez tezina.
          weight = w[min.ind[1]]
          if(!lambda.1se){
            lambda.min <- train.cv$lambda.min
          }else{
            lambda.min <- train.cv$lambda.1se
          }
        }

        # Inner crossvalidation loop with model selection
        coef.list <- coef(lasso, s = lambda.min)
        models.ncv[[i]] <- as.list(c(coefficients = coef.list, lambda = lambda.min, weight = weight))
        # Prediction on test set
        test.pred <- predict(lasso, s = lambda.min, newx = as.matrix(test.data[,-1]))
        test.pred <- pmax(test.pred, target.min/3)

        training.pred <- predict(lasso, s = lambda.min, newx = as.matrix(training.data[,-1]))
        training.pred <- pmax(training.pred, target.min/3)

      }else{if(!step){
        ols.model <- lm(formula = as.formula(paste(target.name,"~", paste(names(training.data[, -1]), collapse="+"))), data = training.data)
      }else{
        ols.model <- stepAIC(lm(formula = as.formula(paste(target.name,"~", paste(names(training.data[, -1]), collapse="+"))), data = training.data))
      }

        # Inner crossvalidation loop with model selection
        coef.list <- coefficients(ols.model)
        models.ncv[[i]] <- as.list(coef.list)
        # Prediction on test set
        test.pred <- predict(ols.model, newdata = test.data[,-1])
        test.pred <- pmax(test.pred, target.min/3)

        training.pred <- predict(ols.model, newdata = training.data[,-1])
        training.pred <- pmax(training.pred, target.min/3)
        min.cv.error <- NA
      }
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
      if(!lambda.1se){
        lambda.min <- hier.lasso.cv$lamhat
      }else{
        lambda.min <- hier.lasso.cv$lamhat.1se
      }

      final.hier.model <- hierNet(training.main.effects, training.target, zz = training.int.effects, diagonal = FALSE, strong=TRUE, lam = lambda.min, center = TRUE, stand.main = FALSE, stand.int = FALSE)

      # Extracting the coefficients of each model
      if(poly.deg == 1){
        int.coeff <- as.matrix(hier.path$th[,,which(hier.path$lamlist == lambda.min)][,length(main.effect.names)])
      } else {
        int.coeff <- as.matrix(hier.path$th[,,which(hier.path$lamlist == lambda.min)][,(length(main.effect.names)-poly.deg+1):length(main.effect.names)])
      }
      main.coeff <- hier.path$bp[,which(hier.path$lamlist == lambda.min), drop = F] - hier.path$bn[,which(hier.path$lamlist == lambda.min), drop = F]

      coef.list <- data.frame(cov.name = colnames(training.main.effects), main.coeff, int.coeff)

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

    if(!use.hier){accuracy[[i]] <- c(defaultSummary(test.obs.pred), RMSEcv = sqrt(min.cv.error))}else{accuracy[[i]] <- defaultSummary(test.obs.pred)}


    obs.pred <- data.frame(test.data, pred = as.numeric(test.pred))
    data.prediction <- rbind(data.prediction,obs.pred)
  }

  # Storing results
  if(all){
    out <- list(RMSE = defaultSummary(test.prediction)[1], Rsquared = defaultSummary(test.prediction)[2], accuracy = accuracy, data.prediction = data.prediction, models = models.ncv)
  }else{
    out <- list(RMSE = defaultSummary(test.prediction)[1], Rsquared = defaultSummary(test.prediction)[2])
  }
  return(out)
}
