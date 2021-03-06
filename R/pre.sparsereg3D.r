#=========== pre.sparsereg3D - data preprocessing ===========
#
# Arguments:
# base.model       - model description (class "formula" of form "target.variable ~ covariates + depth")
# profiles         - observations of target variable (class "SoilProfileCollection")
# cov.grids        - covariates (class SpatialPixelsDataFrame)
# poly.deg         - degree of polynomial depth function
# num.folds        - number of folds in crossvalidation
# num.means        - number of clusters in k-means clustering
# use.interactions - logical. If TRUE interactions will be included in the model.
# standardize      - logical. If TRUE standardization will be performed.
# use.hier         - logical. If TRUE hierarchy constraints will be enforced.
#
# Return value:
# base.model        - model description extended by polynomial terms if poly.deg is TRUE
# poly.deg          - degree of polynomial depth function
# use.interactions  - should interactions be used in the model (binary)
# use.hier          - should hierarachy constraints be enforced (binary)
# main.effect.names - names of main effect covariates
# depth.int.names   - names of interactions with depth (e.g. slope.depth)
# all.int.names     - names of all interactions (not only including depth)

#' Prepearing data for model selection (\code{sparsereg3D.sel} function) and model evaluation (\code{sparsereg3D.ncv} function)
#'
#' @export
#'
#' @param response Only for compositional data modeling. Default is NA. "response" specifies the names of variables in compositions. For example "response = c("sand","silt","clay")"
#' @param base.model model description (class "formula" of form "target.variable ~ covariates + depth" or in case of compositional data modeling just "~ covariates + depth", without target variable specified in the formula.)
#' @param profiles observations of target variable (class "SoilProfileCollection")
#' @param cov.grids covariates (class SpatialPixelsDataFrame)
#' @param poly.deg degree of polynomial depth function
#' @param num.folds number of folds in crossvalidation
#' @param num.means number of clusters in k-means clustering
#' @param use.interactions logical. If TRUE interactions will be included in the model.
#' @param standardize logical. If TRUE standardization will be performed.
#' @param use.hier logical. If TRUE hierarchy constraints will be enforced.
#' @param s Integer. Only for compositional data modeling. "s" specifies which variable in composition will be used as target variable for stratification. For example, s = 1 in composition "sand, silt, clay"  specifies "sand".
#'
#' @return List of objects including:
#' \itemize{
#'  \item \code{base.model}:  model description extended by polynomial terms if poly.deg is TRUE
#'  \item \code{poly.deg}:   degree of polynomial depth function
#'  \item \code{use.interactions}:  logical. If interactions are included in the model.
#'  \item \code{use.hier}:   logical. If hierarchy constraints  are enforced.
#'  \item \code{main.effect.names}:  names of main effect covariates
#'  \item \code{depth.int.names}:   names of interactions with depth (e.g. slope.depth)
#'  \item \code{all.int.names}:  names of all interactions (not only including depth)
#' }
#'
#'  @keywords preprocessing


# TODO Uraditi transformaciju jedinica u cm ako vec nije

pre.sparsereg3D <- function(response, base.model, profiles, cov.grids, use.hier = FALSE, poly.deg = 1, num.folds = 10, num.means = 3, use.interactions = TRUE, standardize = TRUE, coord.trend = TRUE, seed=321, kmean.vars = NULL, cum.prop = 0.90, s = 1){

  "%ni%" <- Negate("%in%")

  if(is.na(response[1])){

  # check if the depth units are cm
  if(profiles@metadata$depth_unit != "cm"){
    stop("Depth units must be cm")
  }

  if(use.hier == TRUE){use.interactions = TRUE}
  if(use.interactions == FALSE){use.hier = FALSE}

  if(length(all.vars(base.model)[-1]) != 1){
    if(sum(names(cov.grids) %in% all.vars(base.model)[-1])==0){
      stop("None of the covariates in the 'formulaString' matches the names in the 'covariates' object")
    }
  }else{
    if(all.vars(base.model)[-1] != "depth"){
      stop("Depth must be included in the model")
    } else {coord.trend = TRUE}
  }


  # check if all covariates are available:


  p4s = proj4string(profiles)

  # creating data frame with all variables
  profiles <- plyr::join(profiles@horizons,data.frame(data.frame(profiles@site),data.frame(profiles@sp)),type="inner")

  # Names
  target.name <- all.vars(base.model)[1]
  coord.names <- tail(names(profiles),2)

  # Adding depth and horizon depth
  profiles$depth <- - (profiles$Top / 100 + (( profiles$Bottom - profiles$Top ) / 2) / 100)
  profiles$hdepth <- profiles$Bottom - profiles$Top
  profiles$mid.depth <- profiles$depth

  # Spatial overlay
  profiles <- profiles[complete.cases(profiles[,c("ID",coord.names,"hdepth","mid.depth","depth",target.name)]),c("ID",target.name,"hdepth","mid.depth",coord.names,"depth")]
  coordinates(profiles) <- coord.names
  proj4string(profiles) <- p4s
  profiles <- spTransform(profiles, proj4string(cov.grids))
  ov <- sp::over(profiles, cov.grids)

  # Extracting the names of categorical variagles and removing empty classes
  factor.names <- ov %>% subset(., select=which(sapply(., is.factor))) %>% names()
  for(i in factor.names){
    ov[,i] <- factor(ov[,i])
  }

  # Preparing data input matrix with following columns: "ID", target.name, "hdepth", coord.names, "mid.depth", sp.cov.names, "depth"
  sp.cov.names <- names(ov[,which(names(ov) %in% c(all.vars(base.model)))])
  profiles <- cbind(as.data.frame(profiles), ov[,sp.cov.names])
  profiles <- profiles[complete.cases(profiles[,all.vars(base.model)]),c("ID",target.name,"hdepth","mid.depth",coord.names, sp.cov.names, "depth")]

  if(coord.trend){base.model <- as.formula(paste(target.name,"~", paste(c(sp.cov.names,coord.names,"depth"), collapse="+"))); sp.cov.names <- c(sp.cov.names, coord.names)}

  # Adding polynomial depth terms in input data matrix, only if poly.deg > 1
  if(poly.deg > 1){
    profiles <- cbind(profiles,poly(profiles$depth, poly.deg,raw=TRUE,simple=TRUE)[,-1])
    names(profiles) <- c(names(profiles)[1:(length(names(profiles))-(poly.deg-1))],(paste("depth", c(2:poly.deg),sep="")))
    base.model <- as.formula(paste(target.name,"~", paste(c(all.vars(base.model)[-1],paste("depth", c(2:poly.deg),sep="")), collapse="+")))
  }

  # Dummy coding
  dummy.par <- dummyVars(as.formula(paste("~", paste(c(all.vars(base.model))[-1], collapse="+"))), profiles, levelsOnly = FALSE)
  if(coord.trend){profiles <- cbind(profiles[, which(colnames(profiles) %in% c("ID","hdepth", "mid.depth", target.name))], predict(dummy.par, newdata = profiles))
                  }else{profiles <- cbind(profiles[, which(colnames(profiles) %in% c("ID","hdepth","mid.depth", target.name, coord.names))], predict(dummy.par, newdata = profiles))
                  }


  # Names
  colnames(profiles) <- gsub( "\\_|/|\\-|\"|\\s" , "." , colnames(profiles) )
  if(coord.trend){main.effect.names <- colnames(profiles)[-c(1:4)]}else{main.effect.names <- colnames(profiles)[-c(1:6)]}

  #main.effect.names <- gsub( "\\_|/|\\-|\"|\\s" , "." , main.effect.names )


  # Computing interactions in input data matrix
  if (use.interactions == TRUE){
    f <- as.formula(~ .^2)
    #interactions <- hierNet::compute.interactions.c(as.matrix(profiles[,-c(which(colnames(profiles) %in% c("ID",target.name,"hdepth",coord.names)))]),diagonal=FALSE)
    interactions <- model.matrix(f, profiles[,main.effect.names]) %>% subset(., select = -(which(colnames(.) %in% c("(Intercept)", main.effect.names))))
    # Names of interactions including depth
    if(poly.deg > 1){
      depth.int.names <- colnames(interactions[,do.call(c,lapply(strsplit(colnames(interactions),":"), function(x) x[2] %in% c("depth",paste("depth",c(2:poly.deg),sep="")) & x[1] %ni% c("depth",paste("depth",c(2:poly.deg),sep=""))))])
    } else {
      depth.int.names <- (interactions %>% as.data.frame() %>% subset(., select = grep("depth", names(.), value=TRUE)) %>% colnames())
    }

    # In hierarchical setting, interactions other than with depth must be zero
    if(use.hier == TRUE) {
      interactions[,colnames(interactions) %ni% depth.int.names ] <- 0
      profiles <- cbind(profiles,interactions)
    } else {
      profiles <- cbind(profiles,interactions[,colnames(interactions) %in% depth.int.names])
    }

  }

  # Standardization of input data
  if(standardize == TRUE) {
    if(use.interactions == TRUE) {
      cnt.par <- as.data.frame(profiles) %>% subset(., select = c(main.effect.names,depth.int.names)) %>% preProcess(.,method=c("center", "scale"))
    } else {
      cnt.par <- as.data.frame(profiles) %>% subset(., select = c(main.effect.names)) %>% preProcess(.,method=c("center", "scale"))
    }
    profiles <- predict(cnt.par,newdata = profiles)
  }

  #TODO Ne sme da osta ne x i y

  # Data stratification
  profiles <- as.data.frame(profiles)
  if(is.null(kmean.vars)){kmean.vars <- coord.names
    } else {
    kmean.vars.ind <- which(apply(sapply(t(kmean.vars), function(x) grepl(x, main.effect.names)), 1, sum) == 1)
    kmean.vars <- c(coord.names, main.effect.names[kmean.vars.ind])
    }

  tmp <- stratfold3d(target.name = target.name, other.names = kmean.vars, seed = seed, data = profiles, num.folds = num.folds, num.means = num.means, cum.prop = cum.prop)
  profile.fold.list <- tmp$profile.fold.list
  obs.fold.list <- tmp$obs.fold.list
  DataWithFolds <- tmp$data

  # Output object creation
  if(use.interactions == TRUE){
    model <- list(base.model = base.model, poly.deg = poly.deg, use.interactions = use.interactions, use.hier = use.hier, cum.prop = cum.prop, main.effect.names = main.effect.names, kmean.vars = kmean.vars, depth.int.names = depth.int.names, all.int.names = colnames(interactions))
  } else {
    model <- list(base.model = base.model, poly.deg = poly.deg, use.interactions = use.interactions, use.hier = use.hier, cum.prop = cum.prop, main.effect.names = main.effect.names, kmean.vars = kmean.vars, depth.int.names = c(), all.int.names = c())
  }
  folds <- list(profile.fold.list = profile.fold.list, obs.fold.list = obs.fold.list, seed = seed)
  std.par <- list(dummy.par = dummy.par, cnt.par = cnt.par)
  out <- list(comp = TRUE, profiles = profiles, DataWithFolds = DataWithFolds, cov.grids = cov.grids, model = model, num.folds = num.folds, num.means = num.means, std.par = std.par, folds = folds)

  return(out)

  }else{

    if(length(response) < 2) stop("Compisition must have at least two variables.")

    # check if the depth units are cm
    if(profiles@metadata$depth_unit != "cm"){
      stop("Depth units must be cm")
    }


    # check if all covariates are available:
    if(sum(names(cov.grids) %in% all.vars(base.model)[-1]) == 0){
      stop("None of the covariates in the 'formulaString' matches the names in the 'covariates' object")
    }

    p4s = proj4string(profiles)

    # creating data frame with all variables
    profiles <- plyr::join(profiles@horizons,data.frame(data.frame(profiles@site),data.frame(profiles@sp)),type="inner")

    # Names
    target.name <- response  #all.vars(base.model)[1]
    coord.names <- tail(names(profiles),2)

    # Adding depth and horizon depth
    profiles$depth <- - (profiles$Top / 100 + (( profiles$Bottom - profiles$Top ) / 2) / 100)
    profiles$hdepth <- profiles$Bottom - profiles$Top
    profiles$mid.depth <- profiles$depth

    # Spatial overlay
    profiles <- profiles[complete.cases(profiles[,c("ID",coord.names,"hdepth","mid.depth","depth",target.name)]),c("ID",target.name,"hdepth","mid.depth",coord.names,"depth")]
    coordinates(profiles) <- coord.names
    proj4string(profiles) <- p4s
    profiles <- spTransform(profiles, proj4string(cov.grids))
    ov <- sp::over(profiles, cov.grids)

    # Extracting the names of categorical variagles and removing empty classes
    factor.names <- ov %>% subset(., select=which(sapply(., is.factor))) %>% names()
    for(i in factor.names){
      ov[,i] <- factor(ov[,i])
    }

    # Preparing data input matrix with following columns: "ID", target.name, "hdepth", coord.names, "mid.depth", sp.cov.names, "depth"
    sp.cov.names <- names(ov[,which(names(ov) %in% c(all.vars(base.model)))])
    profiles <- cbind(as.data.frame(profiles), ov[,sp.cov.names])
    profiles <- profiles[complete.cases(profiles[,all.vars(base.model)]),c("ID",target.name,"hdepth","mid.depth",coord.names, sp.cov.names, "depth")]

    if(coord.trend){base.model <- as.formula(paste("~", paste(c(sp.cov.names,coord.names,"depth"), collapse="+"))); sp.cov.names <- c(sp.cov.names, coord.names)}

    # Different for compositional data
    # Adding polynomial depth terms in input data matrix, only if poly.deg > 1
    if(poly.deg > 1){
      profiles <- cbind(profiles,poly(profiles$depth, poly.deg,raw=TRUE,simple=TRUE)[,-1])
      names(profiles) <- c(names(profiles)[1:(length(names(profiles))-(poly.deg-1))],(paste("depth", c(2:poly.deg),sep="")))
      base.model <- as.formula(paste("~", paste(c(all.vars(base.model)[-1],paste("depth", c(2:poly.deg),sep="")), collapse="+")))
    }

    # Different for compositional data
    # Dummy coding
    dummy.par <- dummyVars(as.formula(paste("~", paste(c(all.vars(base.model)), collapse="+"))), profiles, levelsOnly = FALSE)
    if(coord.trend){profiles <- cbind(profiles[, which(colnames(profiles) %in% c("ID","hdepth", "mid.depth", target.name))], predict(dummy.par, newdata = profiles))
    }else{profiles <- cbind(profiles[, which(colnames(profiles) %in% c("ID","hdepth","mid.depth", target.name, coord.names))], predict(dummy.par, newdata = profiles))
    }

    # Different for compositional data
    # Names
    colnames(profiles) <- gsub( "\\_|/|\\-|\"|\\s" , "." , colnames(profiles) )
    if(coord.trend){main.effect.names <- colnames(profiles)[-c(1:6)]}else{main.effect.names <- colnames(profiles)[-c(1:8)]}

    #main.effect.names <- gsub( "\\_|/|\\-|\"|\\s" , "." , main.effect.names )

    # Different for compositional data
    # Computing interactions in input data matrix
    if (use.interactions == TRUE){
      f <- as.formula(~ .^2)
      #interactions <- hierNet::compute.interactions.c(as.matrix(profiles[,-c(which(colnames(profiles) %in% c("ID",target.name,"hdepth",coord.names)))]),diagonal=FALSE)
      interactions <- model.matrix(f, profiles[,main.effect.names]) %>% subset(., select = -(which(colnames(.) %in% c("(Intercept)", main.effect.names))))
      # Names of interactions including depth
      if(poly.deg > 1){
        depth.int.names <- colnames(interactions[,do.call(c,lapply(strsplit(colnames(interactions),":"), function(x) x[2] %in% c("depth",paste("depth",c(2:poly.deg),sep="")) & x[1] %ni% c("depth",paste("depth",c(2:poly.deg),sep=""))))])
      } else {
        depth.int.names <- (interactions %>% as.data.frame() %>% subset(., select = grep("depth", names(.), value=TRUE)) %>% colnames())
      }

      profiles <- cbind(profiles,interactions[,colnames(interactions) %in% depth.int.names])
    }

    # Standardization of input data
    if(standardize == TRUE) {
      if(use.interactions == TRUE) {
        cnt.par <- as.data.frame(profiles) %>% subset(., select = c(main.effect.names,depth.int.names)) %>% preProcess(.,method=c("center", "scale"))
      } else {
        cnt.par <- as.data.frame(profiles) %>% subset(., select = c(main.effect.names)) %>% preProcess(.,method=c("center", "scale"))
      }
      profiles <- predict(cnt.par,newdata = profiles)
    }

    #TODO Ne sme da osta ne x i y

    # Different for compositional data
    # Data stratification
    profiles <- as.data.frame(profiles)
    if(is.null(kmean.vars)){kmean.vars <- coord.names
    } else {
      kmean.vars.ind <- which(apply(sapply(t(kmean.vars), function(x) grepl(x, main.effect.names)), 1, sum) == 1)
      kmean.vars <- c(coord.names, main.effect.names[kmean.vars.ind])
    }

    tmp <- stratfold3d(target.name = target.name[s], other.names = kmean.vars, seed = seed, data = profiles, num.folds = num.folds, num.means = num.means, cum.prop = cum.prop)
    profile.fold.list <- tmp$profile.fold.list
    obs.fold.list <- tmp$obs.fold.list
    DataWithFolds <- tmp$data

    # Output object creation
    if(use.interactions == TRUE){
      model <- list(response = target.name, base.model = base.model, poly.deg = poly.deg, use.interactions = use.interactions, cum.prop = cum.prop, main.effect.names = main.effect.names, kmean.vars = kmean.vars, depth.int.names = depth.int.names, all.int.names = colnames(interactions))
    } else {
      model <- list(response = target.name, base.model = base.model, poly.deg = poly.deg, use.interactions = use.interactions, cum.prop = cum.prop, main.effect.names = main.effect.names, kmean.vars = kmean.vars, depth.int.names = c(), all.int.names = c())
    }
    folds <- list(profile.fold.list = profile.fold.list, obs.fold.list = obs.fold.list, seed = seed)
    std.par <- list(dummy.par = dummy.par, cnt.par = cnt.par)
    out <- list(comp = TRUE, profiles = profiles, DataWithFolds = DataWithFolds, cov.grids = cov.grids, model = model, num.folds = num.folds, num.means = num.means, std.par = std.par, folds = folds)

    return(out)


  }

}
