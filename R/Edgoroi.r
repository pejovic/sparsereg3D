
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

fun.path <- "D:/_R projects/sparsereg3D/R"
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

stack250m <- stack(edgeroi.grids)
str(stack250m)

## load the 100 m grids:
con2 <- url("http://gsif.isric.org/lib/exe/fetch.php?media=edgeroi.grids100.rda")
load(con2)
str(edgeroi.grids100)
gridded(edgeroi.grids100) <- ~x+y
proj4string(edgeroi.grids100) <- CRS("+init=epsg:28355")

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

str(cov.maps)

# Data
edgeroi$horizons <- rename(edgeroi$horizons, c("UHDICM"="Top", "LHDICM"="Bottom", "SOURCEID" = "ID"))
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
BaseL.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 1, num.folds = 10, num.means = 3, cov.grids = cov.maps)    
BaseL.ORC.ncv.time <- system.time(BaseL.ORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.ORC.time <- system.time(BaseL.ORC <- sparsereg3D.sel(sparse.reg = BaseL.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = BaseL.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

BaseP.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 2, num.folds = 10, num.means = 3, cov.grids = cov.maps)    
BaseP.ORC.ncv.time <- system.time(BaseP.ORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.ORC.time <- system.time(BaseP.ORC <- sparsereg3D.sel(sparse.reg = BaseP.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.p.pred <- sparsereg3D.pred(model.info = ORC.p.sel, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntL.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 10, num.means = 3, cov.grids = cov.maps)    
IntL.ORC.ncv.time <- system.time(IntL.ORC.ncv <- sparsereg3D.ncv(sparse.reg = IntL.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.ORC.time <- system.time(IntL.ORC <- sparsereg3D.sel(sparse.reg = IntL.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntL.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntP.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 2, num.folds = 10, num.means = 3, cov.grids = cov.maps)    
IntP.ORC.ncv.time <- system.time(IntP.ORC.ncv <- sparsereg3D.ncv(sparse.reg = IntP.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.ORC.time <- system.time(IntP.ORC <- sparsereg3D.sel(sparse.reg = IntP.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntP.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntHL.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 10, num.means = 3, cov.grids = cov.maps)    
IntHL.ORC.ncv.time <- system.time(IntHL.ORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.ORC.time <- system.time(IntHL.ORC <- sparsereg3D.sel(sparse.reg = IntHL.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntHL.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntHP.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 2, num.folds = 10, num.means = 3, cov.grids = cov.maps)    
IntHP.ORC.ncv.time <- system.time(IntHP.ORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.ORC.time <- system.time(IntHP.ORC <- sparsereg3D.sel(sparse.reg = IntHP.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntHP.ORC, chunk.size = 20000, grids = edgeroi.grids, depths = c(-0.1,-0.2,-0.3))



