
library(rgdal)
library(GSIF)
library(gdalUtils)
library(raster)
library(plyr)
library(aqp)
library(psych)
library(mda)
library(classInt)
library(caret)
library(MASS)
library(splines)
library(glmnet)
library(hierNet)
library(magrittr)
library(doParallel)
library(foreach)
library(stargazer)
library(gstat) 

fun.path <- paste(getwd(),"R", sep = "/")
source(paste(fun.path,"stratfold3d.r",sep="/"))
source(paste(fun.path,"pre.sparsereg3D.r",sep="/"))
source(paste(fun.path,"sparsereg3D.ncv.r",sep="/"))
source(paste(fun.path,"sparsereg3D.pred.r",sep="/"))
source(paste(fun.path,"sparsereg3D.sel.r",sep="/"))



data(edgeroi)

## load the 250 m grids:
con <- url("http://gsif.isric.org/lib/exe/fetch.php?media=edgeroi.grids.rda")
load(con)
gridded(edgeroi.grids) <- ~x+y
proj4string(edgeroi.grids) <- CRS("+init=epsg:28355")
PMTGEO5.levels <- levels(edgeroi.grids$PMTGEO5)

stack250m <- stack(edgeroi.grids)

## load the 100 m grids:
con2 <- url("http://gsif.isric.org/lib/exe/fetch.php?media=edgeroi.grids100.rda")
load(con2)
str(edgeroi.grids100)
gridded(edgeroi.grids100) <- ~x+y
proj4string(edgeroi.grids100) <- CRS("+init=epsg:28355")

LNUABS6.levels <- levels(edgeroi.grids100$LNUABS6)

# Categorical grids
cat.100.stack <- stack(edgeroi.grids100[c("LNUABS6")])
# Continual grids
con.100.stack <- stack(edgeroi.grids100[c("MVBSRT6", "TI1LAN6","TI2LAN6", "PCKGAD6", "RUTGAD6" ,"PCTGAD6")])

e <- extent(edgeroi.grids)

cat.100.stack <- crop(cat.100.stack, e)
con.100.stack <- crop(con.100.stack, e)

con.250.stack <- resample(con.100.stack, stack250m, method = "bilinear")
cat.250.stack <- resample(cat.100.stack, stack250m, method = "ngb")

grids250 <- stack(con.250.stack, cat.250.stack, stack250m)
names(grids250)

cov.maps <- as(grids250, "SpatialPixelsDataFrame")

factors <- c("PMTGEO5", "LNUABS6")
f <- colwise(as.factor, .cols = factors)

cov.maps@data[,factors] <- f(cov.maps@data[,factors])

levels(cov.maps$PMTGEO5) <- PMTGEO5.levels

str(cov.maps)

# Data
edgeroi$horizons <- rename(edgeroi$horizons, c("UHDICM"="Top", "LHDICM"="Bottom", "SOURCEID" = "ID"))
edgeroi$horizons$log1pORC <- log1p(edgeroi$horizons$ORCDRC)
edgeroi$sites <- rename(edgeroi$sites, c("SOURCEID" = "ID"))
sites <- edgeroi$sites
coordinates(sites) <- ~ LONGDA94 + LATGDA94
proj4string(sites) <- CRS("+proj=longlat +ellps=GRS80 +datum=WGS84 +no_defs")
sites <- spTransform(sites,proj4string(edgeroi.grids))
edgeroi$sites <- data.frame(sites)
edgeroi$sites <- rename(edgeroi$sites, c("LONGDA94"="x", "LATGDA94"="y"))

edgeroi.spc <- join(edgeroi$horizons, edgeroi$sites, type='inner')
depths(edgeroi.spc) <- ID ~ Top + Bottom
site(edgeroi.spc) <- ~ x + y + TAXGAUC + NOTEOBS
coordinates(edgeroi.spc) <- ~x+y
proj4string(edgeroi.spc) <- CRS(proj4string(edgeroi.grids))

str(edgeroi.spc)

#formula
formulaString <- as.formula(paste(paste("ORCDRC ~ "), paste(c(names(cov.maps),"depth"), collapse="+")))
formulaString

#sparsereg3D
# ORC results
BaseL.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
BaseL.ORC.ncv.time <- system.time(BaseL.ORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.ORC.time <- system.time(BaseL.ORC <- sparsereg3D.sel(sparse.reg = BaseL.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = BaseL.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

BaseP.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
BaseP.ORC.ncv.time <- system.time(BaseP.ORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.ORC.time <- system.time(BaseP.ORC <- sparsereg3D.sel(sparse.reg = BaseP.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.p.pred <- sparsereg3D.pred(model.info = ORC.p.sel, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntL.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntL.ORC.ncv.time <- system.time(IntL.ORC.ncv <- sparsereg3D.ncv(sparse.reg = IntL.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.ORC.time <- system.time(IntL.ORC <- sparsereg3D.sel(sparse.reg = IntL.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntL.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntP.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntP.ORC.ncv.time <- system.time(IntP.ORC.ncv <- sparsereg3D.ncv(sparse.reg = IntP.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.ORC.time <- system.time(IntP.ORC <- sparsereg3D.sel(sparse.reg = IntP.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntP.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntHL.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntHL.ORC.ncv.time <- system.time(IntHL.ORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHL.ORC.time <- system.time(IntHL.ORC <- sparsereg3D.sel(sparse.reg = IntHL.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntHL.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntHP.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntHP.ORC.ncv.time <- system.time(IntHP.ORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.ORC.time <- system.time(IntHP.ORC <- sparsereg3D.sel(sparse.reg = IntHP.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntHP.ORC, chunk.size = 20000, grids = edgeroi.grids, depths = c(-0.1,-0.2,-0.3))

#Results
ll <- length(IntL.ORC$coefficients)
pp <- length(IntHL.ORC$coefficients[,1])+1

cmL.ORC <- data.frame(variable=IntHL.ORC$coefficients[,1], BaseL.ORC.me=BaseL.ORC$coefficients[2:pp], IntL.ORC.me=IntL.ORC$coefficients[2:pp],IntL.ORC.ie=c(IntL.ORC$coefficients[(pp+1):ll],0),IntHL.ORC.me=IntHL.ORC$coefficients[,2],IntHL.ORC.ie=IntHL.ORC$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.ORC$coefficients)
p <- length(IntHP.ORC$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.ORC <- data.frame(variable=IntHP.ORC$coefficients[,1], BaseP.ORC.me=BaseP.ORC$coefficients[2:p], IntP.ORC.me=IntP.ORC$coefficients[2:p],IntP.ORC.ie1=c(IntP.ORC$coefficients[(p+1):l][i1],0,0,0),IntP.ORC.ie2=c(IntP.ORC$coefficients[(p+1):l][i2],0,0,0),IntP.ORC.ie3=c(IntP.ORC$coefficients[(p+1):l][i3],0,0,0),IntHP.ORC.me=IntHP.ORC$coefficients[,2],IntHP.ORC.ie1=IntHP.ORC$coefficients[,3],IntHP.ORC.ie2=IntHP.ORC$coefficients[,4],IntHP.ORC.ie3=IntHP.ORC$coefficients[,5] )
cmORC <- cmP.ORC[,c(1,7:10)]

# Models comparison
ORC.ncv <- data.frame(rbind(BaseL = BaseL.ORC.ncv, BaseP = BaseP.ORC.ncv, IntL = IntL.ORC.ncv, IntP = IntP.ORC.ncv, IntHL = IntHL.ORC.ncv, IntHP = IntHP.ORC.ncv))
ORC.ncv <- data.frame(Model = rownames(ORC.ncv), ORC.ncv)
rownames(ORC.ncv) <- NULL
names(ORC.ncv) <- c("Model","RMSE","R squared")
stargazer(ORC.ncv, summary = FALSE, digits = 2, type = "text")

#Time table

ORC.ncv.time <- rbind(BaseL.ORC.ncv.time, BaseP.ORC.ncv.time, IntL.ORC.ncv.time, IntP.ORC.ncv.time, IntHL.ORC.ncv.time, IntHP.ORC.ncv.time)
ORC.time <- rbind(BaseL.ORC.time, BaseP.ORC.time, IntL.ORC.time, IntP.ORC.time, IntHL.ORC.time, IntHP.ORC.time)

#Number of coefficients

n.coeffs <- function(l.coeffs,p.coeffs){
  BaseL <- data.frame(l.coeffs[,2])
  IntL <- l.coeffs[,3:4]
  IntHL <- l.coeffs[,5:6]
  BaseP <- data.frame(p.coeffs[,2])
  IntP <- p.coeffs[,3:6]
  IntHP <- p.coeffs[,7:10]
  
  n.BaseL <-  data.frame(total = prod(dim(BaseL)), selected = sum(BaseL!=0 ), "main effects" = sum(BaseL!=0 ), "interaction effects" = 0 )
  n.IntL <-   data.frame(total = prod(dim(IntL)),  selected = sum(apply(IntL,2, function(y) sum(y!=0))), "main effects" = apply(IntL,2, function(y) sum(y!=0))[1], "interaction effects" = apply(IntL,2, function(y) sum(y!=0))[2])
  n.IntHL <-   data.frame(total = prod(dim(IntHL)),  selected = sum(apply(IntHL,2, function(y) sum(y!=0))), "main effects" = apply(IntHL,2, function(y) sum(y!=0))[1], "interaction effects" = apply(IntHL,2, function(y) sum(y!=0))[2])
  
  n.BaseP <-  data.frame(total = prod(dim(BaseP)), selected = sum(BaseP!=0 ), "main effects" = sum(BaseP!=0 ), "interaction effects" = 0 )
  n.IntP <-   data.frame(total = prod(dim(IntP)),  selected = sum(apply(IntP,2, function(y) sum(y!=0))), "main effects" = apply(IntP,2, function(y) sum(y!=0))[1], "interaction effects" = sum(apply(IntP,2, function(y) sum(y!=0))[2:4]))
  n.IntHP <-   data.frame(total = prod(dim(IntHP)),  selected = sum(apply(IntHP,2, function(y) sum(y!=0))), "main effects" = apply(IntHP,2, function(y) sum(y!=0))[1], "interaction effects" = sum(apply(IntHP,2, function(y) sum(y!=0))[2:4]))
  
  
  total <- data.frame(rbind(n.BaseL,n.BaseP,n.IntL,n.IntHL,n.IntP,n.IntHP))
  rownames(total) <- NULL
  total <- cbind(model = c("BaseL","BaseP", "IntL","IntHL", "IntP","IntHP"), total)
  return(total)
}

ORC.n.coeffs <- n.coeffs(l.coeffs = cmL.ORC, p.coeffs = cmP.ORC)

save.image(file = "D:/R_projects/edgeroi_ORC.RData")
load(file = "D:/R_projects/edgeroi_ORC.RData")

#formula
formulaString <- as.formula(paste(paste("PHIHO5 ~ "), paste(c(names(cov.maps),"depth"), collapse="+")))
formulaString

# pH results
BaseL.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
BaseL.pH.ncv.time <- system.time(BaseL.pH.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.pH.time <- system.time(BaseL.pH <- sparsereg3D.sel(sparse.reg = BaseL.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = BaseL.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

BaseP.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
BaseP.pH.ncv.time <- system.time(BaseP.pH.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.pH.time <- system.time(BaseP.pH <- sparsereg3D.sel(sparse.reg = BaseP.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.p.pred <- sparsereg3D.pred(model.info = pH.p.sel, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntL.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntL.pH.ncv.time <- system.time(IntL.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntL.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.pH.time <- system.time(IntL.pH <- sparsereg3D.sel(sparse.reg = IntL.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = IntL.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntP.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntP.pH.ncv.time <- system.time(IntP.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntP.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.pH.time <- system.time(IntP.pH <- sparsereg3D.sel(sparse.reg = IntP.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = IntP.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntHL.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntHL.pH.ncv.time <- system.time(IntHL.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHL.pH.time <- system.time(IntHL.pH <- sparsereg3D.sel(sparse.reg = IntHL.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = IntHL.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntHP.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntHP.pH.ncv.time <- system.time(IntHP.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.pH.time <- system.time(IntHP.pH <- sparsereg3D.sel(sparse.reg = IntHP.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = IntHP.pH, chunk.size = 20000, grids = edgeroi.grids, depths = c(-0.1,-0.2,-0.3))

#Results
ll <- length(IntL.pH$coefficients)
pp <- length(IntHL.pH$coefficients[,1])+1

cmL.pH <- data.frame(variable=IntHL.pH$coefficients[,1], BaseL.pH.me=BaseL.pH$coefficients[2:pp], IntL.pH.me=IntL.pH$coefficients[2:pp],IntL.pH.ie=c(IntL.pH$coefficients[(pp+1):ll],0),IntHL.pH.me=IntHL.pH$coefficients[,2],IntHL.pH.ie=IntHL.pH$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.pH$coefficients)
p <- length(IntHP.pH$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.pH <- data.frame(variable=IntHP.pH$coefficients[,1], BaseP.pH.me=BaseP.pH$coefficients[2:p], IntP.pH.me=IntP.pH$coefficients[2:p],IntP.pH.ie1=c(IntP.pH$coefficients[(p+1):l][i1],0,0,0),IntP.pH.ie2=c(IntP.pH$coefficients[(p+1):l][i2],0,0,0),IntP.pH.ie3=c(IntP.pH$coefficients[(p+1):l][i3],0,0,0),IntHP.pH.me=IntHP.pH$coefficients[,2],IntHP.pH.ie1=IntHP.pH$coefficients[,3],IntHP.pH.ie2=IntHP.pH$coefficients[,4],IntHP.pH.ie3=IntHP.pH$coefficients[,5] )
cmpH <- cmP.pH[,c(1,7:10)]

# Models comparison
pH.ncv <- data.frame(rbind(BaseL = BaseL.pH.ncv, BaseP = BaseP.pH.ncv, IntL = IntL.pH.ncv, IntP = IntP.pH.ncv, IntHL = IntHL.pH.ncv, IntHP = IntHP.pH.ncv))
pH.ncv <- data.frame(Model = rownames(pH.ncv), pH.ncv)
rownames(pH.ncv) <- NULL
names(pH.ncv) <- c("Model","RMSE","R squared")
stargazer(pH.ncv, summary = FALSE, digits = 2, type = "text")

#Time table

pH.ncv.time <- rbind(BaseL.pH.ncv.time, BaseP.pH.ncv.time, IntL.pH.ncv.time, IntP.pH.ncv.time, IntHL.pH.ncv.time, IntHP.pH.ncv.time)
pH.time <- rbind(BaseL.pH.time, BaseP.pH.time, IntL.pH.time, IntP.pH.time, IntHL.pH.time, IntHP.pH.time)

#Number of coefficients

n.coeffs <- function(l.coeffs,p.coeffs){
  BaseL <- data.frame(l.coeffs[,2])
  IntL <- l.coeffs[,3:4]
  IntHL <- l.coeffs[,5:6]
  BaseP <- data.frame(p.coeffs[,2])
  IntP <- p.coeffs[,3:6]
  IntHP <- p.coeffs[,7:10]
  
  n.BaseL <-  data.frame(total = prod(dim(BaseL)), selected = sum(BaseL!=0 ), "main effects" = sum(BaseL!=0 ), "interaction effects" = 0 )
  n.IntL <-   data.frame(total = prod(dim(IntL)),  selected = sum(apply(IntL,2, function(y) sum(y!=0))), "main effects" = apply(IntL,2, function(y) sum(y!=0))[1], "interaction effects" = apply(IntL,2, function(y) sum(y!=0))[2])
  n.IntHL <-   data.frame(total = prod(dim(IntHL)),  selected = sum(apply(IntHL,2, function(y) sum(y!=0))), "main effects" = apply(IntHL,2, function(y) sum(y!=0))[1], "interaction effects" = apply(IntHL,2, function(y) sum(y!=0))[2])
  
  n.BaseP <-  data.frame(total = prod(dim(BaseP)), selected = sum(BaseP!=0 ), "main effects" = sum(BaseP!=0 ), "interaction effects" = 0 )
  n.IntP <-   data.frame(total = prod(dim(IntP)),  selected = sum(apply(IntP,2, function(y) sum(y!=0))), "main effects" = apply(IntP,2, function(y) sum(y!=0))[1], "interaction effects" = sum(apply(IntP,2, function(y) sum(y!=0))[2:4]))
  n.IntHP <-   data.frame(total = prod(dim(IntHP)),  selected = sum(apply(IntHP,2, function(y) sum(y!=0))), "main effects" = apply(IntHP,2, function(y) sum(y!=0))[1], "interaction effects" = sum(apply(IntHP,2, function(y) sum(y!=0))[2:4]))
  
  
  total <- data.frame(rbind(n.BaseL,n.BaseP,n.IntL,n.IntHL,n.IntP,n.IntHP))
  rownames(total) <- NULL
  total <- cbind(model = c("BaseL","BaseP", "IntL","IntHL", "IntP","IntHP"), total)
  return(total)
}

pH.n.coeffs <- n.coeffs(l.coeffs = cmL.pH, p.coeffs = cmP.pH)

#save.image(file = "D:/R_projects/edgeroi_pH.RData")

#================= log ORC ==========================================================


#formula
formulaString <- as.formula(paste(paste("log1pORC ~ "), paste(c(names(cov.maps),"depth"), collapse="+")))
formulaString

#sparsereg3D
# logORC results
BaseL.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
BaseL.logORC.ncv.time <- system.time(BaseL.logORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.logORC.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.logORC.time <- system.time(BaseL.logORC <- sparsereg3D.sel(sparse.reg = BaseL.logORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#logORC.l.pred <- sparsereg3D.pred(model.info = BaseL.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

BaseP.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
BaseP.logORC.ncv.time <- system.time(BaseP.logORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.logORC.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.logORC.time <- system.time(BaseP.logORC <- sparsereg3D.sel(sparse.reg = BaseP.logORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#logORC.p.pred <- sparsereg3D.pred(model.info = logORC.p.sel, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntL.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntL.logORC.ncv.time <- system.time(IntL.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntL.logORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.logORC.time <- system.time(IntL.logORC <- sparsereg3D.sel(sparse.reg = IntL.logORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#logORC.l.pred <- sparsereg3D.pred(model.info = IntL.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntP.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntP.logORC.ncv.time <- system.time(IntP.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntP.logORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.logORC.time <- system.time(IntP.logORC <- sparsereg3D.sel(sparse.reg = IntP.logORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#logORC.l.pred <- sparsereg3D.pred(model.info = IntP.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntHL.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntHL.logORC.ncv.time <- system.time(IntHL.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.logORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHL.logORC.time <- system.time(IntHL.logORC <- sparsereg3D.sel(sparse.reg = IntHL.logORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#logORC.l.pred <- sparsereg3D.pred(model.info = IntHL.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntHP.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntHP.logORC.ncv.time <- system.time(IntHP.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.logORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.logORC.time <- system.time(IntHP.logORC <- sparsereg3D.sel(sparse.reg = IntHP.logORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#logORC.l.pred <- sparsereg3D.pred(model.info = IntHP.logORC, chunk.size = 20000, grids = edgeroi.grids, depths = c(-0.1,-0.2,-0.3))

#Results
ll <- length(IntL.logORC$coefficients)
pp <- length(IntHL.logORC$coefficients[,1])+1

cmL.logORC <- data.frame(variable=IntHL.logORC$coefficients[,1], BaseL.logORC.me=BaseL.logORC$coefficients[2:pp], IntL.logORC.me=IntL.logORC$coefficients[2:pp],IntL.logORC.ie=c(IntL.logORC$coefficients[(pp+1):ll],0),IntHL.logORC.me=IntHL.logORC$coefficients[,2],IntHL.logORC.ie=IntHL.logORC$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.logORC$coefficients)
p <- length(IntHP.logORC$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.logORC <- data.frame(variable=IntHP.logORC$coefficients[,1], BaseP.logORC.me=BaseP.logORC$coefficients[2:p], IntP.logORC.me=IntP.logORC$coefficients[2:p],IntP.logORC.ie1=c(IntP.logORC$coefficients[(p+1):l][i1],0,0,0),IntP.logORC.ie2=c(IntP.logORC$coefficients[(p+1):l][i2],0,0,0),IntP.logORC.ie3=c(IntP.logORC$coefficients[(p+1):l][i3],0,0,0),IntHP.logORC.me=IntHP.logORC$coefficients[,2],IntHP.logORC.ie1=IntHP.logORC$coefficients[,3],IntHP.logORC.ie2=IntHP.logORC$coefficients[,4],IntHP.logORC.ie3=IntHP.logORC$coefficients[,5] )
cmlogORC <- cmP.logORC[,c(1,7:10)]

# Models comparison
logORC.ncv <- data.frame(rbind(BaseL = BaseL.logORC.ncv, BaseP = BaseP.logORC.ncv, IntL = IntL.logORC.ncv, IntP = IntP.logORC.ncv, IntHL = IntHL.logORC.ncv, IntHP = IntHP.logORC.ncv))
logORC.ncv <- data.frame(Model = rownames(logORC.ncv), logORC.ncv)
rownames(logORC.ncv) <- NULL
names(logORC.ncv) <- c("Model","RMSE","R squared")
stargazer(logORC.ncv, summary = FALSE, digits = 2, type = "text")

#Time table

logORC.ncv.time <- rbind(BaseL.logORC.ncv.time, BaseP.logORC.ncv.time, IntL.logORC.ncv.time, IntP.logORC.ncv.time, IntHL.logORC.ncv.time, IntHP.logORC.ncv.time)
logORC.time <- rbind(BaseL.logORC.time, BaseP.logORC.time, IntL.logORC.time, IntP.logORC.time, IntHL.logORC.time, IntHP.logORC.time)

#Number of coefficients

n.coeffs <- function(l.coeffs,p.coeffs){
  BaseL <- data.frame(l.coeffs[,2])
  IntL <- l.coeffs[,3:4]
  IntHL <- l.coeffs[,5:6]
  BaseP <- data.frame(p.coeffs[,2])
  IntP <- p.coeffs[,3:6]
  IntHP <- p.coeffs[,7:10]
  
  n.BaseL <-  data.frame(total = prod(dim(BaseL)), selected = sum(BaseL!=0 ), "main effects" = sum(BaseL!=0 ), "interaction effects" = 0 )
  n.IntL <-   data.frame(total = prod(dim(IntL)),  selected = sum(apply(IntL,2, function(y) sum(y!=0))), "main effects" = apply(IntL,2, function(y) sum(y!=0))[1], "interaction effects" = apply(IntL,2, function(y) sum(y!=0))[2])
  n.IntHL <-   data.frame(total = prod(dim(IntHL)),  selected = sum(apply(IntHL,2, function(y) sum(y!=0))), "main effects" = apply(IntHL,2, function(y) sum(y!=0))[1], "interaction effects" = apply(IntHL,2, function(y) sum(y!=0))[2])
  
  n.BaseP <-  data.frame(total = prod(dim(BaseP)), selected = sum(BaseP!=0 ), "main effects" = sum(BaseP!=0 ), "interaction effects" = 0 )
  n.IntP <-   data.frame(total = prod(dim(IntP)),  selected = sum(apply(IntP,2, function(y) sum(y!=0))), "main effects" = apply(IntP,2, function(y) sum(y!=0))[1], "interaction effects" = sum(apply(IntP,2, function(y) sum(y!=0))[2:4]))
  n.IntHP <-   data.frame(total = prod(dim(IntHP)),  selected = sum(apply(IntHP,2, function(y) sum(y!=0))), "main effects" = apply(IntHP,2, function(y) sum(y!=0))[1], "interaction effects" = sum(apply(IntHP,2, function(y) sum(y!=0))[2:4]))
  
  
  total <- data.frame(rbind(n.BaseL,n.BaseP,n.IntL,n.IntHL,n.IntP,n.IntHP))
  rownames(total) <- NULL
  total <- cbind(model = c("BaseL","BaseP", "IntL","IntHL", "IntP","IntHP"), total)
  return(total)
}

logORC.n.coeffs <- n.coeffs(l.coeffs = cmL.logORC, p.coeffs = cmP.logORC)

save.image(file = "D:/R_projects/edgeroi_logORC_pH.RData")


load(file = "D:/R_projects/edgeroi_logORC_pH.RData")
load(file = "D:/R_projects/edgeroi_ORC.RData")
