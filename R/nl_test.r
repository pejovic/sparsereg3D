
path <- "C:/Users/Milutin/Dropbox/Extensions of soil 3D trend models/GSIF_competition_NL"

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
fun.path <- paste(getwd(), "R", sep = "/")
source(paste(fun.path,"stratfold3d.r",sep="/"))
source(paste(fun.path,"pre.sparsereg3D.r",sep="/"))
source(paste(fun.path,"sparsereg3D.ncv.r",sep="/"))
source(paste(fun.path,"sparsereg3D.pred.r",sep="/"))
source(paste(fun.path,"sparsereg3D.sel.r",sep="/"))


#list modis files

modis.list<- dir(path=paste(path,"MODIS", sep = "/"), pattern=glob2rx("*.tif"), full.names=FALSE)

modis <- readGDAL(paste(path, "MODIS", modis.list[1],sep="/"))
names(modis)[1]<-sub(".tif","",modis.list[1])

for(i in modis.list[-1]){
  modis@data[sub(".tif","",i[1])] <- readGDAL(paste(path, "MODIS",paste(i),sep="/"))$band1
}

proj4string(modis) <- "+init=epsg:28992"

# 250m grids 
grid.list2 <- dir(path=paste(path, "grids250m", sep = "/"), pattern=glob2rx("*.tif"), full.names=FALSE)

gridmaps250 <- readGDAL(paste(path, "grids250m" ,grid.list2[1],sep="/"))
names(gridmaps250)[1]<-sub(".tif","",grid.list2[1])

for(i in grid.list2[-c(1)]){
  gridmaps250@data[sub(".tif","",i[1])] <- readGDAL(paste(path, "grids250m", paste(i),sep="/"))$band1
}

proj4string(gridmaps250) <- CRS(proj4string(modis))
str(gridmaps250)


# Other grids

#list files
grid.list <- dir(path=paste(path, sep = "/"), pattern=glob2rx("*.tif"), full.names=FALSE)

# Read grids into R:
gridmaps <- readGDAL(paste(path, grid.list[1],sep="/"))
names(gridmaps)[1]<-sub(".tif","",grid.list[1])

for(i in grid.list[-c(1)]){
  gridmaps@data[sub(".tif","",i[1])] <- readGDAL(paste(path, paste(i),sep="/"))$band1
}

proj4string(gridmaps) <- CRS(proj4string(modis))
str(gridmaps)



# Rasample grids

# Categorical grids
stack.cat.grids <- stack(gridmaps[c("geomorfology", "landcover1970", "landcover1992", "landcover2004", "soilmap", "groundwater")])
# Continual grids
stack.con.grids <- stack(gridmaps[c("relativeElevation")])
stack.con.grids250 <- stack(gridmaps250)
stack.modis <- stack(modis)

e <- extent(stack.con.grids)

stack.modis <- crop(stack.modis, e)
stack.con.grids250 <- crop(stack.con.grids250, e)

stack.con.grids <- resample(stack.con.grids, stack.modis, method = "bilinear")
stack.con.grids250 <- resample(stack.con.grids250, stack.modis, method = "bilinear")
stack.cat.grids <- resample(stack.cat.grids, stack.modis, method = "ngb")

grids250 <- stack(stack.con.grids, stack.cat.grids, stack.con.grids250, stack.modis)
names(grids250)
cov.maps <- as(grids250, "SpatialPixelsDataFrame")


factors <- c("geomorfology", "landcover1970", "landcover1992", "landcover2004", "soilmap", "groundwater")
f <- colwise(as.factor, .cols = factors)
cov.maps@data[,factors] <- f(cov.maps@data[,factors])

str(cov.maps)

#Data
SITE <- read.csv(paste("C:/Users/Milutin/Dropbox/Extensions of soil 3D trend models/GSIF_competition_NL", "Alterra_Soil_xy.csv", sep="/"), header = TRUE)
str(SITE)
## 643 profiles
SITE$SOURCEDB = "Alterra-BODEMDATA"
SITE$SOURCEID <- paste("Alterra", SITE$Profile_id, sep="_")

SITE.s <- SITE[!duplicated(SITE$SOURCEID),]
coordinates(SITE.s) <- ~ X + Y
proj4string(SITE.s) <- "+init=epsg:28992"
SITE.s <- as.data.frame(SITE.s)
#View(SITE.s)

horizons <- read.table(paste("C:/Users/Milutin/Dropbox/Extensions of soil 3D trend models/GSIF_competition_NL","Alterra_Soil_data.csv", sep="/"), header =TRUE, na.strings = c("-99","NA"), sep=",")
str(horizons)

horizons$ORCDRC <- horizons$Organic.Matter....*10 /1.724     ## OC in permilles
summary(horizons$ORCDRC)
horizons$logORCDRC <- (log1p(horizons$ORCDRC))
horizons$SOURCEID <- paste("Alterra", horizons$Profile_id, sep="_")
horizons$SAMPLEID <- make.unique(paste(horizons$Profile_id, horizons$Layer, sep="_"))
horizons <- rename(horizons, c("pH.KCl"="PHIKCL", "Sand...50mu....."="SNDPPT", "Clay...2.mu....."="CLYPPT", "Silt..2.50.mu....."="SLTPPT", "CEC..mmol.kg."="CECSUM", "Bulk.density..g.per.cm3."="BLD", "Start..cm."="Top", "End..cm."="Bottom"))
summary(horizons$Top)
summary(horizons$Bottom)
horizons$DEPTH <- horizons$Top + (horizons$Bottom - horizons$Top)/2


SPROPS.Alterra <- join(horizons[,c("SOURCEID","SAMPLEID","Top","Bottom","DEPTH","SNDPPT","CLYPPT","SLTPPT","PHIKCL","ORCDRC","logORCDRC","CECSUM","BLD")], SITE.s[,c("SOURCEID","SOURCEDB","X","Y")], type="left")
SPROPS.Alterra <- rename(SPROPS.Alterra, c("SOURCEID" = "ID", "X"="x", "Y"="y"))
SPROPS.Alterra <- SPROPS.Alterra[!is.na(SPROPS.Alterra$x) & !is.na(SPROPS.Alterra$y) & !is.na(SPROPS.Alterra$DEPTH),]
str(SPROPS.Alterra)
## 2617

plot(SPROPS.Alterra$x, SPROPS.Alterra$y, pch="+")

depths(SPROPS.Alterra) <- SAMPLEID ~ Top + Bottom
site(SPROPS.Alterra) <- ~ ID + x + y
coordinates(SPROPS.Alterra) <- ~ x + y
proj4string(SPROPS.Alterra) <- "+init=epsg:28992"

nl.profiles <- SPROPS.Alterra

str(nl.profiles)

#Aggregation profiles
nl.profiles@horizons <- rename(nl.profiles@horizons, c("PHIKCL"="pH", "ORCDRC"="SOC"))
agg <- slab(nl.profiles, fm= ~ SOC + pH, slab.structure=seq(0,70,5))

## see ?slab for details on the default aggregate function
head(agg)

Figure2 <- xyplot(top ~ p.q50 | variable, data=agg, ylab='Depth',
                  xlab='median bounded by 25th and 75th percentiles',
                  lower=agg$p.q25, upper=agg$p.q75, ylim=c(70,-2),
                  panel=panel.depth_function,
                  alpha=0.25, sync.colors=TRUE,
                  par.settings=list(superpose.line=list(col='RoyalBlue', lwd=2)),
                  prepanel=prepanel.depth_function,
                  cf=agg$contributing_fraction, cf.col='black', cf.interval=5, 
                  layout=c(2,1), strip=strip.custom(bg=grey(0.8)),
                  scales=list(x=list(tick.number=4, alternating=3, relation='free'))
)

class(Figure2)

pdf("NLAgg.pdf",width=6,height=8)
plot(Figure2) # Make plot
dev.off()
#==========================================================


formulaString <- as.formula(paste(paste("ORCDRC ~ "), paste(c(names(cov.maps),"depth"), collapse="+"))) #[-which(names(cov.maps) %in% c("slope","TWI","LS_factor"))]
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
IntHL.ORC.time <- system.time(IntHL.ORC <- sparsereg3D.sel(sparse.reg = IntHL.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntHL.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntHP.ORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntHP.ORC.ncv.time <- system.time(IntHP.ORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.ORC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.ORC.time <- system.time(IntHP.ORC <- sparsereg3D.sel(sparse.reg = IntHP.ORC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORC.l.pred <- sparsereg3D.pred(model.info = IntHP.ORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))


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
  
  
  total <- data.frame(rbind(n.BaseL,n.BaseP,n.IntL,n.IntP,n.IntHL,n.IntHP))
  rownames(total) <- NULL
  total <- cbind(model = c("BaseL","BaseP", "IntL", "IntP","IntHL","IntHP"), total)
  return(total)
}

ORC.n.coeffs <- n.coeffs(l.coeffs = cmL.ORC, p.coeffs = cmP.ORC)
#save.image(file = "D:/R_projects/nl_ORCDRC.RData")

load(file = "D:/R_projects/nl_ORCDRC.RData")

#rm(list=ls())















path <- "C:/Users/Milutin/Dropbox/Extensions of soil 3D trend models/GSIF_competition_NL"

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
fun.path <- paste(getwd(), "R", sep = "/")
source(paste(fun.path,"stratfold3d.r",sep="/"))
source(paste(fun.path,"pre.sparsereg3D.r",sep="/"))
source(paste(fun.path,"sparsereg3D.ncv.r",sep="/"))
source(paste(fun.path,"sparsereg3D.pred.r",sep="/"))
source(paste(fun.path,"sparsereg3D.sel.r",sep="/"))


#list modis files

modis.list<- dir(path=paste(path,"MODIS", sep = "/"), pattern=glob2rx("*.tif"), full.names=FALSE)

modis <- readGDAL(paste(path, "MODIS", modis.list[1],sep="/"))
names(modis)[1]<-sub(".tif","",modis.list[1])

for(i in modis.list[-1]){
  modis@data[sub(".tif","",i[1])] <- readGDAL(paste(path, "MODIS",paste(i),sep="/"))$band1
}

proj4string(modis) <- "+init=epsg:28992"

# 250m grids 
grid.list2 <- dir(path=paste(path, "grids250m", sep = "/"), pattern=glob2rx("*.tif"), full.names=FALSE)

gridmaps250 <- readGDAL(paste(path, "grids250m" ,grid.list2[1],sep="/"))
names(gridmaps250)[1]<-sub(".tif","",grid.list2[1])

for(i in grid.list2[-c(1)]){
  gridmaps250@data[sub(".tif","",i[1])] <- readGDAL(paste(path, "grids250m", paste(i),sep="/"))$band1
}

proj4string(gridmaps250) <- CRS(proj4string(modis))
str(gridmaps250)


# Other grids

#list files
grid.list <- dir(path=paste(path, sep = "/"), pattern=glob2rx("*.tif"), full.names=FALSE)

# Read grids into R:
gridmaps <- readGDAL(paste(path, grid.list[1],sep="/"))
names(gridmaps)[1]<-sub(".tif","",grid.list[1])

for(i in grid.list[-c(1)]){
  gridmaps@data[sub(".tif","",i[1])] <- readGDAL(paste(path, paste(i),sep="/"))$band1
}

proj4string(gridmaps) <- CRS(proj4string(modis))
str(gridmaps)



# Rasample grids

# Categorical grids
stack.cat.grids <- stack(gridmaps[c("geomorfology", "landcover1970", "landcover1992", "landcover2004", "soilmap", "groundwater")])
# Continual grids
stack.con.grids <- stack(gridmaps[c("relativeElevation")])
stack.con.grids250 <- stack(gridmaps250)
stack.modis <- stack(modis)

e <- extent(stack.con.grids)

stack.modis <- crop(stack.modis, e)
stack.con.grids250 <- crop(stack.con.grids250, e)

stack.con.grids <- resample(stack.con.grids, stack.modis, method = "bilinear")
stack.con.grids250 <- resample(stack.con.grids250, stack.modis, method = "bilinear")
stack.cat.grids <- resample(stack.cat.grids, stack.modis, method = "ngb")

grids250 <- stack(stack.con.grids, stack.cat.grids, stack.con.grids250, stack.modis)
names(grids250)
cov.maps <- as(grids250, "SpatialPixelsDataFrame")


factors <- c("geomorfology", "landcover1970", "landcover1992", "landcover2004", "soilmap", "groundwater")
f <- colwise(as.factor, .cols = factors)
cov.maps@data[,factors] <- f(cov.maps@data[,factors])

str(cov.maps)

#Data
SITE <- read.csv(paste("C:/Users/Milutin/Dropbox/Extensions of soil 3D trend models/GSIF_competition_NL", "Alterra_Soil_xy.csv", sep="/"), header = TRUE)
str(SITE)
## 643 profiles
SITE$SOURCEDB = "Alterra-BODEMDATA"
SITE$SOURCEID <- paste("Alterra", SITE$Profile_id, sep="_")
SITE.s <- SITE[!duplicated(SITE$SOURCEID),]
coordinates(SITE.s) <- ~ X + Y
proj4string(SITE.s) <- "+init=epsg:28992"
SITE.s <- as.data.frame(SITE.s)
#View(SITE.s)
horizons <- read.table(paste("C:/Users/Milutin/Dropbox/Extensions of soil 3D trend models/GSIF_competition_NL","Alterra_Soil_data.csv", sep="/"), header =TRUE, na.strings = c("-99","NA"), sep=",")
str(horizons)

horizons$ORCDRC <- horizons$Organic.Matter....*10 /1.724     ## OC in permilles
summary(horizons$ORCDRC)
horizons$logORCDRC <- (log1p(horizons$ORCDRC))
horizons$SOURCEID <- paste("Alterra", horizons$Profile_id, sep="_")
horizons$SAMPLEID <- make.unique(paste(horizons$Profile_id, horizons$Layer, sep="_"))
horizons <- rename(horizons, c("pH.KCl"="PHIKCL", "Sand...50mu....."="SNDPPT", "Clay...2.mu....."="CLYPPT", "Silt..2.50.mu....."="SLTPPT", "CEC..mmol.kg."="CECSUM", "Bulk.density..g.per.cm3."="BLD", "Start..cm."="Top", "End..cm."="Bottom"))
summary(horizons$Top)
summary(horizons$Bottom)
horizons$DEPTH <- horizons$Top + (horizons$Bottom - horizons$Top)/2


SPROPS.Alterra <- join(horizons[,c("SOURCEID","SAMPLEID","Top","Bottom","DEPTH","SNDPPT","CLYPPT","SLTPPT","PHIKCL","ORCDRC","logORCDRC","CECSUM","BLD")], SITE.s[,c("SOURCEID","SOURCEDB","X","Y")], type="left")
SPROPS.Alterra <- rename(SPROPS.Alterra, c("SOURCEID" = "ID", "X"="x", "Y"="y"))
SPROPS.Alterra <- SPROPS.Alterra[!is.na(SPROPS.Alterra$x) & !is.na(SPROPS.Alterra$y) & !is.na(SPROPS.Alterra$DEPTH),]
str(SPROPS.Alterra)
## 2617

plot(SPROPS.Alterra$x, SPROPS.Alterra$y, pch="+")

depths(SPROPS.Alterra) <- SAMPLEID ~ Top + Bottom
site(SPROPS.Alterra) <- ~ ID + x + y
coordinates(SPROPS.Alterra) <- ~ x + y
proj4string(SPROPS.Alterra) <- "+init=epsg:28992"

nl.profiles <- SPROPS.Alterra

str(nl.profiles)


formulaString <- as.formula(paste(paste("logORCDRC ~ "), paste(c(names(cov.maps),"depth"), collapse="+"))) #[-which(names(cov.maps) %in% c("slope","TWI","LS_factor"))]
formulaString


# logORC results
BaseL.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
BaseL.logORC.ncv.time <- system.time(BaseL.logORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.logORC.preproc, lambda = seq(0,0.2,0.001), seed = 321))
BaseL.logORC.time <- system.time(BaseL.logORC <- sparsereg3D.sel(sparse.reg = BaseL.logORC.preproc ,lambda = seq(0,0.2,0.001), seed = 321))
#logORC.l.pred <- sparsereg3D.pred(model.info = BaseL.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

BaseP.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
BaseP.logORC.ncv.time <- system.time(BaseP.logORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.logORC.preproc, lambda = seq(0,0.2,0.001), seed = 321))
BaseP.logORC.time <- system.time(BaseP.logORC <- sparsereg3D.sel(sparse.reg = BaseP.logORC.preproc ,lambda = seq(0,0.2,0.001), seed = 321))
#logORC.p.pred <- sparsereg3D.pred(model.info = logORC.p.sel, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntL.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntL.logORC.ncv.time <- system.time(IntL.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntL.logORC.preproc, lambda = seq(0,0.2,0.001), seed = 321))
IntL.logORC.time <- system.time(IntL.logORC <- sparsereg3D.sel(sparse.reg = IntL.logORC.preproc ,lambda = seq(0,0.2,0.001), seed = 321))
#logORC.l.pred <- sparsereg3D.pred(model.info = IntL.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntP.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntP.logORC.ncv.time <- system.time(IntP.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntP.logORC.preproc, lambda = seq(0,0.2,0.001), seed = 321))
IntP.logORC.time <- system.time(IntP.logORC <- sparsereg3D.sel(sparse.reg = IntP.logORC.preproc ,lambda = seq(0,0.2,0.001), seed = 321))
#logORC.l.pred <- sparsereg3D.pred(model.info = IntP.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

library(doParallel)

registerDoParallel(cores=4)

pack = (.packages())


result = foreach (j = 1:2, .packages = pack) %dopar% { 
  
  if (j==1){
    IntHL.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
    IntHL.logORC.ncv.time <- system.time(IntHL.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.logORC.preproc, lambda = seq(0,0.2,0.001), seed = 321))
    IntHL.logORC.time <- system.time(IntHL.logORC <- sparsereg3D.sel(sparse.reg = IntHL.logORC.preproc ,lambda = seq(0,0.2,0.001), seed = 321))
    #logORC.l.pred <- sparsereg3D.pred(model.info = IntHL.logORC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))
    list(IntHL.logORC.preproc, IntHL.logORC.ncv.time, IntHL.logORC.ncv, IntHL.logORC.time, IntHL.logORC)
  } else if (j==2) {
    IntHP.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
    IntHP.logORC.ncv.time <- system.time(IntHP.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.logORC.preproc, lambda = seq(0,0.2,0.001), seed = 321))
    IntHP.logORC.time <- system.time(IntHP.logORC <- sparsereg3D.sel(sparse.reg = IntHP.logORC.preproc ,lambda = seq(0,0.2,0.001), seed = 321))
    #logORC.l.pred <- sparsereg3D.pred(model.info = IntHP.logORC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))
    list(IntHP.logORC.preproc, IntHP.logORC.ncv.time, IntHP.logORC.ncv, IntHP.logORC.time, IntHP.logORC)
  }
  
}

stopImplicitCluster()

IntHL.logORC.preproc <- result[[1]][[1]]
IntHL.logORC.ncv.time <- result[[1]][[2]] 
IntHL.logORC.ncv <- result[[1]][[3]]
IntHL.logORC.time <- result[[1]][[4]]
IntHL.logORC <- result[[1]][[5]]

IntHP.logORC.preproc <- result[[2]][[1]]
IntHP.logORC.ncv.time <- result[[2]][[2]] 
IntHP.logORC.ncv <- result[[2]][[3]]
IntHP.logORC.time <- result[[2]][[4]]
IntHP.logORC <- result[[2]][[5]]


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
  
  
  total <- data.frame(rbind(n.BaseL,n.BaseP,n.IntL,n.IntP,n.IntHL,n.IntHP))
  rownames(total) <- NULL
  total <- cbind(model = c("BaseL","BaseP", "IntL", "IntP","IntHL","IntHP"), total)
  return(total)
}

logORC.n.coeffs <- n.coeffs(l.coeffs = cmL.logORC, p.coeffs = cmP.logORC)
save.image(file = "D:/R_projects/nl_logORCDRC.RData")

load(file = "D:/R_projects/nl_logORCDRC.RData")

load(file = "D:/R_projects/result.rda")

rm(list=ls())














path <- "C:/Users/User/Dropbox/Extensions of soil 3D trend models/GSIF_competition_NL"

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
fun.path <- paste(getwd(), "R", sep = "/")
source(paste(fun.path,"stratfold3d.r",sep="/"))
source(paste(fun.path,"pre.sparsereg3D.r",sep="/"))
source(paste(fun.path,"sparsereg3D.ncv.r",sep="/"))
source(paste(fun.path,"sparsereg3D.pred.r",sep="/"))
source(paste(fun.path,"sparsereg3D.sel.r",sep="/"))


#list modis files

modis.list<- dir(path=paste(path,"MODIS", sep = "/"), pattern=glob2rx("*.tif"), full.names=FALSE)

modis <- readGDAL(paste(path, "MODIS", modis.list[1],sep="/"))
names(modis)[1]<-sub(".tif","",modis.list[1])

for(i in modis.list[-1]){
  modis@data[sub(".tif","",i[1])] <- readGDAL(paste(path, "MODIS",paste(i),sep="/"))$band1
}

proj4string(modis) <- "+init=epsg:28992"

# 250m grids 
grid.list2 <- dir(path=paste(path, "grids250m", sep = "/"), pattern=glob2rx("*.tif"), full.names=FALSE)

gridmaps250 <- readGDAL(paste(path, "grids250m" ,grid.list2[1],sep="/"))
names(gridmaps250)[1]<-sub(".tif","",grid.list2[1])

for(i in grid.list2[-c(1)]){
  gridmaps250@data[sub(".tif","",i[1])] <- readGDAL(paste(path, "grids250m", paste(i),sep="/"))$band1
}

proj4string(gridmaps250) <- CRS(proj4string(modis))
str(gridmaps250)


# Other grids

#list files
grid.list <- dir(path=paste(path, sep = "/"), pattern=glob2rx("*.tif"), full.names=FALSE)

# Read grids into R:
gridmaps <- readGDAL(paste(path, grid.list[1],sep="/"))
names(gridmaps)[1]<-sub(".tif","",grid.list[1])

for(i in grid.list[-c(1)]){
  gridmaps@data[sub(".tif","",i[1])] <- readGDAL(paste(path, paste(i),sep="/"))$band1
}

proj4string(gridmaps) <- CRS(proj4string(modis))
str(gridmaps)



# Rasample grids

# Categorical grids
stack.cat.grids <- stack(gridmaps[c("geomorfology", "landcover1970", "landcover1992", "landcover2004", "soilmap", "groundwater")])
# Continual grids
stack.con.grids <- stack(gridmaps[c("relativeElevation")])
stack.con.grids250 <- stack(gridmaps250)
stack.modis <- stack(modis)

e <- extent(stack.con.grids)

stack.modis <- crop(stack.modis, e)
stack.con.grids250 <- crop(stack.con.grids250, e)

stack.con.grids <- resample(stack.con.grids, stack.modis, method = "bilinear")
stack.con.grids250 <- resample(stack.con.grids250, stack.modis, method = "bilinear")
stack.cat.grids <- resample(stack.cat.grids, stack.modis, method = "ngb")

grids250 <- stack(stack.con.grids, stack.cat.grids, stack.con.grids250, stack.modis)
names(grids250)
cov.maps <- as(grids250, "SpatialPixelsDataFrame")


factors <- c("geomorfology", "landcover1970", "landcover1992", "landcover2004", "soilmap", "groundwater")
f <- colwise(as.factor, .cols = factors)
cov.maps@data[,factors] <- f(cov.maps@data[,factors])

str(cov.maps)

#Data
SITE <- read.csv(paste("C:/Users/User/Dropbox/Extensions of soil 3D trend models/GSIF_competition_NL", "Alterra_Soil_xy.csv", sep="/"), header = TRUE)
str(SITE)
## 643 profiles
SITE$SOURCEDB = "Alterra-BODEMDATA"
SITE$SOURCEID <- paste("Alterra", SITE$Profile_id, sep="_")
SITE.s <- SITE[!duplicated(SITE$SOURCEID),]
coordinates(SITE.s) <- ~ X + Y
proj4string(SITE.s) <- "+init=epsg:28992"
SITE.s <- as.data.frame(SITE.s)
#View(SITE.s)
horizons <- read.table(paste("C:/Users/User/Dropbox/Extensions of soil 3D trend models/GSIF_competition_NL","Alterra_Soil_data.csv", sep="/"), header =TRUE, na.strings = c("-99","NA"), sep=",")
str(horizons)

horizons$ORCDRC <- horizons$Organic.Matter....*10 /1.724     ## OC in permilles
summary(horizons$ORCDRC)
horizons$logORCDRC <- (log1p(horizons$ORCDRC))
horizons$SOURCEID <- paste("Alterra", horizons$Profile_id, sep="_")
horizons$SAMPLEID <- make.unique(paste(horizons$Profile_id, horizons$Layer, sep="_"))
horizons <- rename(horizons, c("pH.KCl"="PHIKCL", "Sand...50mu....."="SNDPPT", "Clay...2.mu....."="CLYPPT", "Silt..2.50.mu....."="SLTPPT", "CEC..mmol.kg."="CECSUM", "Bulk.density..g.per.cm3."="BLD", "Start..cm."="Top", "End..cm."="Bottom"))
summary(horizons$Top)
summary(horizons$Bottom)
horizons$DEPTH <- horizons$Top + (horizons$Bottom - horizons$Top)/2


SPROPS.Alterra <- join(horizons[,c("SOURCEID","SAMPLEID","Top","Bottom","DEPTH","SNDPPT","CLYPPT","SLTPPT","PHIKCL","ORCDRC","logORCDRC","CECSUM","BLD")], SITE.s[,c("SOURCEID","SOURCEDB","X","Y")], type="left")
SPROPS.Alterra <- rename(SPROPS.Alterra, c("SOURCEID" = "ID", "X"="x", "Y"="y"))
SPROPS.Alterra <- SPROPS.Alterra[!is.na(SPROPS.Alterra$x) & !is.na(SPROPS.Alterra$y) & !is.na(SPROPS.Alterra$DEPTH),]
str(SPROPS.Alterra)
## 2617

plot(SPROPS.Alterra$x, SPROPS.Alterra$y, pch="+")

depths(SPROPS.Alterra) <- SAMPLEID ~ Top + Bottom
site(SPROPS.Alterra) <- ~ ID + x + y
coordinates(SPROPS.Alterra) <- ~ x + y
proj4string(SPROPS.Alterra) <- "+init=epsg:28992"

nl.profiles <- SPROPS.Alterra

str(nl.profiles)



# pH 
formulaString <- as.formula(paste(paste("PHIKCL ~ "), paste(c(names(cov.maps),"depth"), collapse="+")))
formulaString


# pH results
BaseL.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
BaseL.pH.ncv.time <- system.time(BaseL.pH.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.pH.time <- system.time(BaseL.pH <- sparsereg3D.sel(sparse.reg = BaseL.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = BaseL.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

BaseP.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
BaseP.pH.ncv.time <- system.time(BaseP.pH.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.pH.time <- system.time(BaseP.pH <- sparsereg3D.sel(sparse.reg = BaseP.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.p.pred <- sparsereg3D.pred(model.info = pH.p.sel, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntL.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntL.pH.ncv.time <- system.time(IntL.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntL.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.pH.time <- system.time(IntL.pH <- sparsereg3D.sel(sparse.reg = IntL.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = IntL.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntP.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntP.pH.ncv.time <- system.time(IntP.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntP.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.pH.time <- system.time(IntP.pH <- sparsereg3D.sel(sparse.reg = IntP.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = IntP.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntHL.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntHL.pH.ncv.time <- system.time(IntHL.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHL.pH.time <- system.time(IntHL.pH <- sparsereg3D.sel(sparse.reg = IntHL.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = IntHL.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntHP.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
IntHP.pH.ncv.time <- system.time(IntHP.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.pH.time <- system.time(IntHP.pH <- sparsereg3D.sel(sparse.reg = IntHP.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = IntHP.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))


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
  
  
  total <- data.frame(rbind(n.BaseL,n.BaseP,n.IntL,n.IntP,n.IntHL,n.IntHP))
  rownames(total) <- NULL
  total <- cbind(model = c("BaseL","BaseP", "IntL", "IntP","IntHL","IntHP"), total)
  return(total)
}

pH.n.coeffs <- n.coeffs(l.coeffs = cmL.pH, p.coeffs = cmP.pH)




save.image(file = "D:/R_projects/nl_pHKCL.RData")

load(file = "D:/R_projects/result.rda")

IntHL.pH.preproc <- result[[1]][[1]]
IntHL.pH.ncv.time <- result[[1]][[2]] 
IntHL.pH.ncv <- result[[1]][[3]]
IntHL.pH.time <- result[[1]][[4]]
IntHL.pH <- result[[1]][[5]]

IntHP.pH.preproc <- result[[2]][[1]]
IntHP.pH.ncv.time <- result[[2]][[2]] 
IntHP.pH.ncv <- result[[2]][[3]]
IntHP.pH.time <- result[[2]][[4]]
IntHP.pH <- result[[2]][[5]]

IntHL.logORC.preproc <- result[[3]][[1]]
IntHL.logORC.ncv.time <- result[[3]][[2]] 
IntHL.logORC.ncv <- result[[3]][[3]]
IntHL.logORC.time <- result[[3]][[4]]
IntHL.logORC <- result[[3]][[5]]

IntHP.logORC.preproc <- result[[4]][[1]]
IntHP.logORC.ncv.time <- result[[4]][[2]] 
IntHP.logORC.ncv <- result[[4]][[3]]
IntHP.logORC.time <- result[[4]][[4]]
IntHP.logORC <- result[[4]][[5]]


