#=========== sparsereg3D.sel - function for sparse model selection =======
#
# Arguments:
# sparse.reg - output from pre.sparsereg3D function
# lambda     - vector of regularization parameter values for lasso regression
# seed       - random number generator
#
# Return value:
# Model
# Regularization parameter value
# Model coefficients
# Standardization parameters


#' Model Selection based on Cross-Validation
#'
#' @export
#'
#' @param sparse.reg output from \code{pre.sparsereg3D} function
#' @param lambda vector of regularization parameter values for lasso regression
#' @param ols logical. If TRUE model will be fitted with OLS insted of using lasso
#' @param step logical. If TRUE stepwise procedure will be used when fitting OLS model
#' @param seed random number generator
#' @param lambda.1se logical. If TRUE one sigma lambda rule will be used (largest lambda value with cv.err less than or equal to min(cv.err)+ SE).
#' @param grouped logical. If TRUE, it will ensure that the multinomial coefficients for a variable are all in or out together.
#'
#'
#' @return List of objects including:
#' \itemize{
#'  \item \code{model}:  list of objects containing model description
#'  \item \code{lambda}:   Regularization parameter value for lasso models
#'  \item \code{coefficients}:  Model coefficients
#'  \item \code{std.param}:   Standardization parameters
#'  \item \code{comp}: logical. Indicate whether the compositional variable is modeled or not.
#' }
#'
#'  @keywords Model selection



sparsereg3D.sel <- function(sparse.reg, lambda = 0, ols = FALSE, step = FALSE, lambda.1se = FALSE, grouped = TRUE){

  if(step){
    ols = TRUE
  }

  # Extracting data from sparse.reg object
  seed <- sparse.reg$folds$seed
  profile.fold.list <- sparse.reg$folds$profile.fold.list
  obs.fold.list <- sparse.reg$folds$obs.fold.list
  profiles <- sparse.reg$profiles

  if(sparse.reg$comp){
    target.name <- sparse.reg$model$response
  }else{
    target.name <- all.vars(sparse.reg$model$base.model)[1]
    target.min <- min(profiles[,target.name])
  }

  use.interactions <- sparse.reg$model$use.interactions
  use.hier <- sparse.reg$model$use.hier
  main.effect.names <- sparse.reg$model$main.effect.names
  poly.deg <- sparse.reg$model$poly.deg
  depth.int.names <- sparse.reg$model$depth.int.names
  all.int.names <- sparse.reg$model$all.int.names
  base.model <- sparse.reg$model$base.model
  std.par <- sparse.reg$std.par


  # Convert list of vectors of indices per fold to vector of fold indices for each observation as glmnet requires
  # list(c(1,4),c(2,5),c(3,6)) -> c(1,2,3,1,2,3)
  fold.indices <- rep(NA,dim(profiles)[1])
  for(i in 1:length(obs.fold.list)){
    fold.indices[obs.fold.list[[i]]] <- i
  }

  if(!sparse.reg$comp){

  # Lasso training
  if(!ols){
    if(!use.hier){
      if(use.interactions){
        training.data <- subset(profiles, select = c(target.name, main.effect.names, depth.int.names))
      }else{
        training.data <- subset(profiles, select = c(target.name, main.effect.names))
      }

      lasso.cv <- cv.glmnet(as.matrix(training.data[,-1]), training.data[,1], alpha = 1,lambda = lambda, foldid = fold.indices, type.measure = "mse")
      if(!lambda.1se){
        lambda.min <- lasso.cv$lambda.min
      }else{
        lambda.min <- lasso.cv$lambda.1se
      }
      coef.list <- predict(lasso.cv, type="coefficients", s = lambda.min)
      prediction <- predict(lasso.cv, newx = as.matrix(training.data[,-1]), type = "response", s = lambda.min)
    }else{
      # Hierarchical setting requires the separation of main effects and interaction effects
      training.main.effects <- as.matrix(profiles[,main.effect.names])
      training.int.effects <- as.matrix(profiles[,all.int.names])
      training.target <- (profiles[,target.name])

      hier.path = hierNet.path(training.main.effects,training.target, zz = training.int.effects, diagonal = FALSE, strong = TRUE, trace = 0, stand.main = FALSE, stand.int = FALSE)
      hier.lasso.cv = hierNet.cv(hier.path, training.main.effects, training.target, folds = obs.fold.list, trace=0)
      if(!lambda.1se){
        lambda.min <- hier.lasso.cv$lamhat
      }else{
        lambda.min <- hier.lasso.cv$lamhat.1se
      }
      hier.lasso.final <- hierNet(training.main.effects,training.target, zz = training.int.effects, diagonal = FALSE, strong = TRUE, lam = lambda.min, center = TRUE, stand.main = FALSE, stand.int = FALSE)
      prediction <- predict(hier.lasso.final, newx = training.main.effects, zz = training.int.effects)

      # Extracting the coefficients of final model
      if(poly.deg == 1){
        int.coeff <- as.matrix(hier.path$th[,,which(hier.path$lamlist == lambda.min)][,length(main.effect.names)])
      } else {
        int.coeff <- as.matrix(hier.path$th[,,which(hier.path$lamlist == lambda.min)][,(length(main.effect.names)-poly.deg+1):length(main.effect.names)])
      }
      main.coeff <- hier.path$bp[,which(hier.path$lamlist == lambda.min), drop = F] - hier.path$bn[,which(hier.path$lamlist == lambda.min), drop = F]

      coef.list <- data.frame(cov.name = colnames(training.main.effects), main.coeff,int.coeff)
    }
  }else{
    # selecting model with stepwise regression or plain ols
    if(!step){
      ols.model <- lm(formula = as.formula(paste(target.name,"~", paste(names(profiles[, c(main.effect.names, depth.int.names)]), collapse="+"))), data = profiles[,c(target.name, main.effect.names, depth.int.names)])
    }else{
      ols.model <- stepAIC(lm(formula = as.formula(paste(target.name,"~", paste(names(profiles[, c(main.effect.names, depth.int.names)]), collapse="+"))), data = profiles[,c(target.name, main.effect.names, depth.int.names)]))
    }

    # Inner crossvalidation loop with model selection
    coef.list <- coefficients(ols.model)
    prediction <-  predict(ols.model, newdata = profiles[,c(main.effect.names, depth.int.names)])
    model.info <- list(data = data.frame(profiles, prediction = prediction), model = ols.model, coefficients = coef.list)
  }


  # Regression summary list containing the final model, final lambda, model coefficients, and standardization parameters
  if(!ols){
    if(!use.hier){
      model.info <- list(data = data.frame(profiles, prediction = prediction), model = list(model = lasso.cv, lambda = lasso.cv$lambda.min, target.name = target.name, main.effect.names = main.effect.names, depth.int.names = depth.int.names, use.interactions = use.interactions, use.hier = use.hier, poly.deg = poly.deg, base.model = base.model), coefficients = coef.list, std.par = std.par, comp = FALSE)
    }else{
      model.info <- list(data = data.frame(profiles, prediction = prediction), model = list(model = hier.lasso.final, lambda = hier.path$lamlist[which(hier.lasso.cv$lamhat == hier.path$lamlist)], target.name = target.name, main.effect.names = main.effect.names, depth.int.names = depth.int.names, use.interactions = use.interactions, use.hier = use.hier, poly.deg = poly.deg, base.model = base.model), coefficients = coef.list, std.par = std.par, comp = FALSE)
    }
  }
  return(model.info)

  }else{
    if(grouped){grouped = "grouped"}
    # Lasso training
    if(use.interactions){
      training.data <- subset(profiles, select = c(target.name, main.effect.names, depth.int.names))
    }else{
      training.data <- subset(profiles, select = c(target.name, main.effect.names))
    }

    lasso.cv <- cv.glmnet(as.matrix(training.data[,-c(1:3)]), as.matrix(training.data[, c(1:3)]), family="multinomial", alpha = 1,lambda = lambda, foldid = fold.indices, type.multinomial = grouped)
    coef.list <- predict(lasso.cv, type="coefficients", s=lasso.cv$lambda.min)
    prediction <- predict(lasso.cv, newx = as.matrix(training.data[,-c(1:3)]), type = "response", s=lasso.cv$lambda.min)


    # Regression summary list containing the final model, final lambda, model coefficients, and standardization parameters

    model.info <- list(data = data.frame(profiles, prediction = prediction), model = list(model = lasso.cv, lambda = lasso.cv$lambda.min, target.name = target.name, main.effect.names = main.effect.names, depth.int.names = depth.int.names, use.interactions = use.interactions, use.hier = use.hier, poly.deg = poly.deg, base.model = base.model), coefficients = coef.list, std.par = std.par, comp = TRUE)

    return(model.info)

  }
}
