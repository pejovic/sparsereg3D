#=========== pre.sparsereg3D - data preprocessing ===========
#
# Arguments:
# base.model       - model description (class "formula" of form "target.variable ~ covariates + depth")
# profiles         - observations of target variable (class "SoilProfileCollection")
# cov.grids        - covariates (class SpatialPixelsDataFrame)
# poly.deg         - degree of polynomial depth function
# num.folds         - number of folds in crossvalidation 
# num.means        - number of clusters in k-means clustering
# use.interactions - should interactions be included in the model (binary)
# standardize      - should standardization be performed (binary)
# use.hier         - should hierarchy constraints be enforced (binary)
#
# Return value:
# base.model        - model description extended by polynomial terms if poly.deg is TRUE
# poly.deg          - degree of polynomial depth function
# use.interactions  - should interactions be used in the model (binary)
# use.hier          - should hierarachy constraints be enforced (binary)
# main.effect.names - names of main effect covariates
# depth.int.names   - names of interactions with depth (e.g. slope.depth)
# all.int.names     - names of all interactions (not only including depth)


# TODO Uraditi transformaciju jedinica u cm ako vec nije

pre.sparsereg3D <- function(base.model, profiles, cov.grids, use.hier=FALSE, poly.deg = 1, num.folds = 10, num.means = 3, use.interactions = TRUE, standardize=TRUE, seed=321){
  
  "%ni%" <- Negate("%in%")
  
  # check if the depth units are cm
  if(profiles@metadata$depth_unit != "cm"){
    stop("Depth units must be cm")
  }
  
  if(use.hier == TRUE){use.interactions = TRUE}
  if(use.interactions == FALSE){use.hier = FALSE}
  
  
  # check if all covariates are available:
  if(sum(names(cov.grids) %in% all.vars(base.model)[-1])==0){
    stop("None of the covariates in the 'formulaString' matches the names in the 'covariates' object")
  }
  
  p4s = proj4string(profiles)
  
  # creating data frame with all variables
  profiles <- join(profiles@horizons,data.frame(data.frame(profiles@site),data.frame(profiles@sp)),type="inner") 
  
  # Names
  target.name <- all.vars(base.model)[1]
  coord.names <- tail(names(profiles),2)
  
  # Adding depth and horizon depth
  profiles$depth <- - (profiles$Top / 100 + (( profiles$Bottom - profiles$Top ) / 2) / 100) 
  profiles$hdepth <- profiles$Bottom - profiles$Top
  
  # Spatial overlay
  profiles <- profiles[complete.cases(profiles[,c("ID",coord.names,"hdepth","depth",target.name)]),c("ID",target.name,"hdepth",coord.names,"depth")]
  coordinates(profiles) <- ~ x + y
  proj4string(profiles) <- p4s
  profiles <- spTransform(profiles, proj4string(cov.grids))
  ov <- over(profiles, cov.grids)
  
  # Extracting the names of categorical variagles and removing empty classes
  factor.names <- ov %>% subset(., select=which(sapply(., is.factor))) %>% names() 
  for(i in factor.names){
    ov[,i] <- factor(ov[,i])
  }
  
  # Preparing data input matrix with following columns: "ID", target.name, "hdepth", coord.names, sp.cov.names, "depth" 
  sp.cov.names <- names(ov[,which(names(ov) %in% c(all.vars(base.model)))])
  profiles <- cbind(as.data.frame(profiles), ov[,sp.cov.names])
  profiles <- profiles[complete.cases(profiles[,all.vars(base.model)]),c("ID",target.name,"hdepth",coord.names, sp.cov.names, "depth")]
  
  # Adding polynomial depth terms in input data matrix, only if poly.deg > 1  
  if(poly.deg > 1){
    profiles <- cbind(profiles,poly(profiles$depth,poly.deg,raw=TRUE,simple=TRUE)[,-1])
    names(profiles) <- c(names(profiles)[1:(length(names(profiles))-(poly.deg-1))],(paste("depth",c(2:poly.deg),sep="")))
    base.model <- as.formula(paste(target.name,"~", paste(c(all.vars(base.model)[-1],paste("depth",c(2:poly.deg),sep="")), collapse="+")))
  }
  
  # Dummy coding
  dummy.par <- dummyVars(as.formula(paste("~", paste(c(all.vars(base.model))[-1], collapse="+"))),profiles,levelsOnly=FALSE) 
  profiles <- cbind(profiles[,which(colnames(profiles) %in% c("ID","hdepth",target.name,coord.names))], predict(dummy.par, newdata = profiles)) 
  
  # Names
  colnames(profiles) <- gsub( "\\_|/|\\-|\"|\\s" , "." , colnames(profiles) ) 
  main.effect.names <- colnames(profiles)[-which(colnames(profiles) %in% c("ID","hdepth",target.name,coord.names))]
  
  # Computing interactions in input data matrix
  if (use.interactions == TRUE){
    
    interactions <- hierNet::compute.interactions.c(as.matrix(profiles[,-c(which(colnames(profiles) %in% c("ID",target.name,"hdepth",coord.names)))]),diagonal=FALSE)
    
    # Names of interactions including depth
    if(poly.deg > 1){ 
      depth.int.names <- colnames(interactions[,do.call(c,lapply(strsplit(colnames(interactions),":"), function(x) x[2] %in% c("depth",paste("depth",c(2:poly.deg),sep="")) & x[1] %ni% c("depth",paste("depth",c(2:poly.deg),sep=""))))]) 
    } else {
      depth.int.names <- (interactions %>% as.data.frame() %>% subset(., select=grep("depth", names(.), value=TRUE)) %>% colnames())
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
  profiles <- plyr::rename(profiles, replace=c("x" = "longitude", "y" = "latitude"))
  tmp <- stratfold3d(target.name = target.name, seed = seed, data = profiles, num.folds = num.folds, num.means = num.means)
  profile.fold.list <- tmp$profile.fold.list
  obs.fold.list <- tmp$obs.fold.list
  
  # Output object creation
  if(use.interactions == TRUE){
    model <- list(base.model = base.model, poly.deg = poly.deg, use.interactions = use.interactions, use.hier = use.hier, main.effect.names = main.effect.names, depth.int.names = depth.int.names, all.int.names = colnames(interactions))    
  } else {
    model <- list(base.model = base.model, poly.deg = poly.deg, use.interactions = use.interactions, use.hier = use.hier, main.effect.names = main.effect.names, depth.int.names = c(), all.int.names = c())
  }
  folds <- list(profile.fold.list = profile.fold.list, obs.fold.list = obs.fold.list, seed = seed)
  std.par <- list(dummy.par = dummy.par, cnt.par = cnt.par)
  out <- list(profiles = profiles, cov.grids = cov.grids, model = model, num.folds = num.folds, num.means = num.means, std.par = std.par, folds = folds)
  
  return(out)
  
}  
