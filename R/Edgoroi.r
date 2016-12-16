library(rgdal)
library(aqp)
library(sp)
library(GSIF)
library(plyr)



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

stack100m <- stack(edgeroi.grids100)
str(stack100m)

stack250new <- resample(stack100m, stack250m, method = "ngb")

edgeroi.grids100 <- as(stack250new,"SpatialPixelsDataFrame")
str(edgeroi.grids100)

edgeroi.grids100$LNUABS6 <- factor(edgeroi.grids100$LNUABS6)

edgeroi.grids@data <- data.frame(edgeroi.grids@data, edgeroi.grids100@data) 





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

glm.formulaString <- as.formula(paste(paste("PHIHO5 ~ "), paste(c(names(edgeroi.grids),"depth"), collapse="+")))
glm.formulaString


# SOM results
BaseL.SOM.preproc <- pre.sparsereg3D(base.model = glm.formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = edgeroi.grids)    
BaseL.SOM.ncv.time <- system.time(BaseL.SOM.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.SOM.time <- system.time(BaseL.SOM <- sparsereg3D.sel(sparse.reg = BaseL.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = BaseL.SOM, chunk.size = 20000, grids = edgeroi.grids, depths = c(-0.1,-0.2,-0.3))

BaseP.SOM.preproc <- pre.sparsereg3D(base.model = glm.formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = edgeroi.grids)    
BaseP.SOM.ncv.time <- system.time(BaseP.SOM.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.SOM.time <- system.time(BaseP.SOM <- sparsereg3D.sel(sparse.reg = BaseP.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.p.pred <- sparsereg3D.pred(model.info = SOM.p.sel, chunk.size = 20000, grids = edgeroi.grids, depths = c(-0.1,-0.2,-0.3))


IntL.SOM.preproc <- pre.sparsereg3D(base.model = glm.formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = edgeroi.grids)    
IntL.SOM.ncv.time <- system.time(IntL.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntL.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.SOM.time <- system.time(IntL.SOM <- sparsereg3D.sel(sparse.reg = IntL.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = IntL.SOM, chunk.size = 20000, grids = edgeroi.grids, depths = c(-0.1,-0.2,-0.3))

IntP.SOM.preproc <- pre.sparsereg3D(base.model = glm.formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = edgeroi.grids)    
IntP.SOM.ncv.time <- system.time(IntP.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntP.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.SOM.time <- system.time(IntP.SOM <- sparsereg3D.sel(sparse.reg = IntP.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = IntP.SOM, chunk.size = 20000, grids = edgeroi.grids, depths = c(-0.1,-0.2,-0.3))

IntHL.SOM.preproc <- pre.sparsereg3D(base.model = glm.formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = edgeroi.grids)    
IntHL.SOM.ncv.time <- system.time(IntHL.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.SOM.time <- system.time(IntHL.SOM <- sparsereg3D.sel(sparse.reg = IntHL.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = IntHL.SOM, chunk.size = 20000, grids = edgeroi.grids, depths = c(-0.1,-0.2,-0.3))

IntHP.SOM.preproc <- pre.sparsereg3D(base.model = glm.formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = edgeroi.grids)    
IntHP.SOM.ncv.time <- system.time(IntHP.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.SOM.time <- system.time(IntHP.SOM <- sparsereg3D.sel(sparse.reg = IntHP.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = IntHP.SOM, chunk.size = 20000, grids = edgeroi.grids, depths = c(-0.1,-0.2,-0.3))




download.file("http://gsif.isric.org/zipped/NL250m_covs.zip", "NL250m_covs.zip")
library(R.utils)
unzip("NL250m_covs.zip")


