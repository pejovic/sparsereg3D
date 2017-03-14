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


#load("nl.profs.rda")
#load("nl.covs.rda")


formulaString <- as.formula(paste(paste("logORCDRC ~ "), paste(c(names(cov.maps),"depth"), collapse="+"))) #[-which(names(cov.maps) %in% c("slope","TWI","LS_factor"))]
formulaString

# logORC results
BaseL.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)    
BaseL.logORC.ncv.time <- system.time(BaseL.logORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.logORC.preproc, lambda = seq(0,0.2,0.001)))
BaseL.logORC.time <- system.time(BaseL.logORC <- sparsereg3D.sel(sparse.reg = BaseL.logORC.preproc ,lambda = seq(0,0.2,0.001)))
#logORC.l.pred <- sparsereg3D.pred(model.info = BaseL.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

BaseP.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)    
BaseP.logORC.ncv.time <- system.time(BaseP.logORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.logORC.preproc, lambda = seq(0,0.2,0.001)))
BaseP.logORC.time <- system.time(BaseP.logORC <- sparsereg3D.sel(sparse.reg = BaseP.logORC.preproc ,lambda = seq(0,0.2,0.001)))
#logORC.p.pred <- sparsereg3D.pred(model.info = logORC.p.sel, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntL.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)    
IntL.logORC.ncv.time <- system.time(IntL.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntL.logORC.preproc, lambda = seq(0,0.2,0.001)))
IntL.logORC.time <- system.time(IntL.logORC <- sparsereg3D.sel(sparse.reg = IntL.logORC.preproc ,lambda = seq(0,0.2,0.001)))
#logORC.l.pred <- sparsereg3D.pred(model.info = IntL.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

#model.info = IntL.logORC; grids = cov.maps; depths = c(-0.1)

IntP.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)    
IntP.logORC.ncv.time <- system.time(IntP.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntP.logORC.preproc, lambda = seq(0,0.2,0.001), seed = 555))
IntP.logORC.time <- system.time(IntP.logORC <- sparsereg3D.sel(sparse.reg = IntP.logORC.preproc ,lambda = seq(0,0.2,0.001), seed = 555))
#logORC.l.pred <- sparsereg3D.pred(model.info = IntP.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

library(doParallel)

registerDoParallel(cores=4)

pack = (.packages())


SOCresult = foreach (j = 1:2, .packages = pack) %dopar% { 
  
  if (j==1){
    IntHL.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
    IntHL.logORC.ncv.time <- system.time(IntHL.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.logORC.preproc, lambda = seq(0,0.2,0.001)))
    IntHL.logORC.time <- system.time(IntHL.logORC <- sparsereg3D.sel(sparse.reg = IntHL.logORC.preproc ,lambda = seq(0,0.2,0.001)))
    #logORC.l.pred <- sparsereg3D.pred(model.info = IntHL.logORC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))
    list(IntHL.logORC.preproc, IntHL.logORC.ncv.time, IntHL.logORC.ncv, IntHL.logORC.time, IntHL.logORC)
  } else if (j==2) {
    IntHP.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps)    
    IntHP.logORC.ncv.time <- system.time(IntHP.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.logORC.preproc, lambda = seq(0,0.2,0.001)))
    IntHP.logORC.time <- system.time(IntHP.logORC <- sparsereg3D.sel(sparse.reg = IntHP.logORC.preproc ,lambda = seq(0,0.2,0.001)))
    #logORC.l.pred <- sparsereg3D.pred(model.info = IntHP.logORC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))
    list(IntHP.logORC.preproc, IntHP.logORC.ncv.time, IntHP.logORC.ncv, IntHP.logORC.time, IntHP.logORC)
  }
  
}

stopImplicitCluster()

IntHL.logORC.preproc <- SOCresult[[1]][[1]]
IntHL.logORC.ncv.time <- SOCresult[[1]][[2]] 
IntHL.logORC.ncv <- SOCresult[[1]][[3]]
IntHL.logORC.time <- SOCresult[[1]][[4]]
IntHL.logORC <- SOCresult[[1]][[5]]

IntHP.logORC.preproc <- SOCresult[[2]][[1]]
IntHP.logORC.ncv.time <- SOCresult[[2]][[2]]
IntHP.logORC.ncv <- SOCresult[[2]][[3]]
IntHP.logORC.time <- SOCresult[[2]][[4]]
IntHP.logORC <- SOCresult[[2]][[5]]


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

stargazer(cmP.logORC, summary = FALSE, digits = 2, type = "text")

# Models comparison
logORC.ncv <- data.frame(rbind(BaseL = BaseL.logORC.ncv[1:2], BaseP = BaseP.logORC.ncv[1:2], IntL = IntL.logORC.ncv[1:2], IntP = IntP.logORC.ncv[1:2], IntHL = IntHL.logORC.ncv[1:2], IntHP = IntHP.logORC.ncv[1:2]))
logORC.ncv <- data.frame(Model = rownames(logORC.ncv), logORC.ncv)
rownames(logORC.ncv) <- NULL
names(logORC.ncv) <- c("Model","RMSE","R squared")
stargazer(logORC.ncv, summary = FALSE, digits = 2, type = "text")

#Time table

logORC.ncv.time <- rbind(BaseL.logORC.ncv.time, BaseP.logORC.ncv.time, IntL.logORC.ncv.time, IntP.logORC.ncv.time, IntHL.logORC.ncv.time, IntHP.logORC.ncv.time)
logORC.time <- rbind(BaseL.logORC.time, BaseP.logORC.time, IntL.logORC.time, IntP.logORC.time, IntHL.logORC.time, IntHP.logORC.time)

#Number of coefficients


logORC.n.coeffs <- n.coeffs(l.coeffs = cmL.logORC, p.coeffs = cmP.logORC)



# pH 
formulaString <- as.formula(paste(paste("PHIKCL ~ "), paste(c(names(cov.maps),"depth"), collapse="+")))
formulaString
seed = 321

# pH results
BaseL.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)    
BaseL.pH.ncv.time <- system.time(BaseL.pH.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.pH.preproc, lambda = seq(0,0.2,0.001)))
BaseL.pH.time <- system.time(BaseL.pH <- sparsereg3D.sel(sparse.reg = BaseL.pH.preproc ,lambda = seq(0,0.2,0.001)))
#pH.l.pred <- sparsereg3D.pred(model.info = BaseL.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

BaseP.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)    
BaseP.pH.ncv.time <- system.time(BaseP.pH.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.pH.preproc, lambda = seq(0,0.2,0.001)))
BaseP.pH.time <- system.time(BaseP.pH <- sparsereg3D.sel(sparse.reg = BaseP.pH.preproc ,lambda = seq(0,0.2,0.001)))
#pH.p.pred <- sparsereg3D.pred(model.info = pH.p.sel, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntL.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)    
IntL.pH.ncv.time <- system.time(IntL.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntL.pH.preproc, lambda = seq(0,0.2,0.001)))
IntL.pH.time <- system.time(IntL.pH <- sparsereg3D.sel(sparse.reg = IntL.pH.preproc ,lambda = seq(0,0.2,0.001)))
#pH.l.pred <- sparsereg3D.pred(model.info = IntL.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntP.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)    
IntP.pH.ncv.time <- system.time(IntP.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntP.pH.preproc, lambda = seq(0,0.2,0.001)))
IntP.pH.time <- system.time(IntP.pH <- sparsereg3D.sel(sparse.reg = IntP.pH.preproc ,lambda = seq(0,0.2,0.001)))
#pH.l.pred <- sparsereg3D.pred(model.info = IntP.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

library(doParallel)

registerDoParallel(cores=4)

pack = (.packages())


pHresult = foreach (j = 1:2, .packages = pack) %dopar% { 
  
  if (j==1){
    IntHL.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)    
    IntHL.pH.ncv.time <- system.time(IntHL.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.pH.preproc, lambda = seq(0,0.2,0.001)))
    IntHL.pH.time <- system.time(IntHL.pH <- sparsereg3D.sel(sparse.reg = IntHL.pH.preproc ,lambda = seq(0,0.2,0.001)))
    #pH.l.pred <- sparsereg3D.pred(model.info = IntHL.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))
    list(IntHL.pH.preproc, IntHL.pH.ncv.time, IntHL.pH.ncv, IntHL.pH.time, IntHL.pH)
  } else if (j==2) {
    IntHP.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = nl.profiles, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)    
    IntHP.pH.ncv.time <- system.time(IntHP.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.pH.preproc, lambda = seq(0,0.2,0.001)))
    IntHP.pH.time <- system.time(IntHP.pH <- sparsereg3D.sel(sparse.reg = IntHP.pH.preproc ,lambda = seq(0,0.2,0.001)))
    #pH.l.pred <- sparsereg3D.pred(model.info = IntHP.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))
    list(IntHP.pH.preproc, IntHP.pH.ncv.time, IntHP.pH.ncv, IntHP.pH.time, IntHP.pH)
  }
  
}

stopImplicitCluster()

IntHL.pH.preproc <- pHresult[[1]][[1]]
IntHL.pH.ncv.time <- pHresult[[1]][[2]] 
IntHL.pH.ncv <- pHresult[[1]][[3]]
IntHL.pH.time <- pHresult[[1]][[4]]
IntHL.pH <- pHresult[[1]][[5]]

IntHP.pH.preproc <- pHresult[[2]][[1]]
IntHP.pH.ncv.time <- pHresult[[2]][[2]] 
IntHP.pH.ncv <- pHresult[[2]][[3]]
IntHP.pH.time <- pHresult[[2]][[4]]
IntHP.pH <- pHresult[[2]][[5]]



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

stargazer(cmP.pH, summary = FALSE, digits = 2, type = "latex")

# Models comparison
pH.ncv <- data.frame(rbind(BaseL = BaseL.pH.ncv[1:2], BaseP = BaseP.pH.ncv[1:2], IntL = IntL.pH.ncv[1:2], IntP = IntP.pH.ncv[1:2], IntHL = IntHL.pH.ncv[1:2], IntHP = IntHP.pH.ncv[1:2]))
pH.ncv <- data.frame(Model = rownames(pH.ncv), pH.ncv)
rownames(pH.ncv) <- NULL
names(pH.ncv) <- c("Model","RMSE","R squared")
stargazer(pH.ncv, summary = FALSE, digits = 2, type = "text")

#Time table

pH.ncv.time <- rbind(BaseL.pH.ncv.time, BaseP.pH.ncv.time, IntL.pH.ncv.time, IntP.pH.ncv.time, IntHL.pH.ncv.time, IntHP.pH.ncv.time)
pH.time <- rbind(BaseL.pH.time, BaseP.pH.time, IntL.pH.time, IntP.pH.time, IntHL.pH.time, IntHP.pH.time)

#Number of coefficients


pH.n.coeffs <- n.coeffs(l.coeffs = cmL.pH, p.coeffs = cmP.pH)


