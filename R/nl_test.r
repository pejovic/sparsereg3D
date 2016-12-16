
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

# load functions
fun.path <- "D:/_R projects/sparsereg3D/R"
source(paste(fun.path,"stratfold3d.r",sep="/"))
source(paste(fun.path,"pre.sparsereg3D.r",sep="/"))
source(paste(fun.path,"sparsereg3D.ncv.r",sep="/"))
source(paste(fun.path,"sparsereg3D.pred.r",sep="/"))
source(paste(fun.path,"sparsereg3D.sel.r",sep="/"))


#list modis files

modis.list<- dir(path=paste(getwd(), "NL250m_covs","MODIS", sep = "/"), pattern=glob2rx("*.tif"), full.names=FALSE)

modis <- readGDAL(paste(getwd(), "NL250m_covs","MODIS", modis.list[1],sep="/"))
names(modis)[1]<-sub(".tif","",modis.list[1])

for(i in modis.list[-1]){
  modis@data[sub(".tif","",i[1])] <- readGDAL(paste(getwd(), "NL250m_covs","MODIS",paste(i),sep="/"))$band1
}

# Other grids

#list files
grid.list<- dir(path=paste(getwd(), "NL250m_covs", sep = "/"), pattern=glob2rx("*.tif"), full.names=FALSE)

# Read grids into R:
gridmaps <- readGDAL(paste(getwd(), "NL250m_covs",grid.list[1],sep="/"))
names(gridmaps)[1]<-sub(".tif","",grid.list[1])

for(i in grid.list[-c(1,8)]){
  gridmaps@data[sub(".tif","",i[1])] <- readGDAL(paste(getwd(), "NL250m_covs",paste(i),sep="/"))$band1
}

proj4string(gridmaps) <- CRS(proj4string(modis))
str(gridmaps)

# DEM 

srtm <- readGDAL(paste(getwd(), "NL250m_covs",grid.list[8],sep="/"))
names(srtm)[1] <- sub(".tif","",grid.list[8])
proj4string(srtm) <- CRS(proj4string(modis))
str(srtm)

# Rasample
stackmodis <- stack(modis)
# Categorical grids
stack.cat.grids <- stack(gridmaps[c("geomorfology", "landcover1970", "landcover1992", "landcover2004", "soilmap", "groundwater")])
# Continual grids
stack.con.grids <- stack(gridmaps[c("relativeElevation")])

stacksrtm <- stack(srtm)

e <- extent(stack.con.grids)

stackmodis <- crop(stackmodis, e)
stacksrtm <- crop(stacksrtm, e)


stack.con.grids <- resample(stack.con.grids, stackmodis, method = "bilinear")
stack.cat.grids <- resample(stack.cat.grids, stackmodis, method = "ngb")

grids250 <- stack(stack.con.grids, stack.cat.grids, srtm, stackmodis)
names(grids250)

cov.maps <- as(grids250, "SpatialPixelsDataFrame")

factors <- c("geomorfology", "landcover1970", "landcover1992", "landcover2004", "soilmap", "groundwater")
f <- colwise(as.factor, .cols = factors)

cov.maps@data[,factors] <- f(cov.maps@data[,factors])

str(cov.maps)
#Data

nl.profiles <- read.csv("nl_data.csv", header = TRUE)
names(nl.profiles) <- c("ID", "SAMPLEID", "x", "y", "Top",   "Bottom",   "ORCDRC" ,  "BLD")
#nl.profiles <- transform(nl.profiles, ID=as.numeric(factor(x)))
nl.profiles$ORCDRC <- pmax(nl.profiles$ORCDRC, 0.1)
nl.profiles$logORCDRC <- log(nl.profiles$ORCDRC)

coordinates(nl.profiles) <- ~ x + y
proj4string(nl.profiles) <- CRS("+init=epsg:4326")
sites <- spTransform(nl.profiles, proj4string(cov.maps))
nl.profiles <- data.frame(sites)


depths(nl.profiles) <- ID ~ Top + Bottom
site(nl.profiles) <- ~ x + y
coordinates(nl.profiles) <- ~ x + y
proj4string(nl.profiles) <- CRS(proj4string(cov.maps))


str(nl.profiles)


formulaString <- as.formula(paste(paste("logORCDRC ~ "), paste(c(names(cov.maps),"depth"), collapse="+")))
formulaString


# ORC results
BaseL.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
BaseL.ORC.ncv.time <- system.time(BaseL.ORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.ORC.time <- system.time(BaseL.ORC <- sparsereg3D.sel(sparse.reg = BaseL.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = BaseL.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

BaseP.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
BaseP.ORC.ncv.time <- system.time(BaseP.ORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.ORC.time <- system.time(BaseP.ORC <- sparsereg3D.sel(sparse.reg = BaseP.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.p.pred <- sparsereg3D.pred(model.info = ORC.p.sel, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntL.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntL.ORC.ncv.time <- system.time(IntL.ORC.ncv <- sparsereg3D.ncv(sparse.reg = IntL.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.ORC.time <- system.time(IntL.ORC <- sparsereg3D.sel(sparse.reg = IntL.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntL.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntP.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntP.ORC.ncv.time <- system.time(IntP.ORC.ncv <- sparsereg3D.ncv(sparse.reg = IntP.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.ORC.time <- system.time(IntP.ORC <- sparsereg3D.sel(sparse.reg = IntP.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntP.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntHL.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntHL.ORC.ncv.time <- system.time(IntHL.ORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.ORC.time <- system.time(IntHL.ORC <- sparsereg3D.sel(sparse.reg = IntHL.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntHL.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntHP.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntHP.ORC.ncv.time <- system.time(IntHP.ORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.ORC.time <- system.time(IntHP.ORC <- sparsereg3D.sel(sparse.reg = IntHP.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntHP.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))





