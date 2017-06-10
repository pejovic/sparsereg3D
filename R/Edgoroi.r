
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

#edg.kml <- aqp::slice(edgeroi.spc, 2.5 ~ ORCDRC)
#proj4string(edg.kml) <- CRS(proj4string(edgeroi.grids))
#plotKML(edg.kml["ORCDRC"])

#Aggregation profiles
edgeroi.spc@horizons <- rename(edgeroi.spc@horizons, c("PHIHO5"="pH", "log1pORC"="logSOC"))
agg <- slab(edgeroi.spc, fm= ~ logSOC, slab.structure=20)

## see ?slab for details on the default aggregate function
head(agg)

Figure2 <- xyplot(top ~ p.q50 | variable, data=agg, ylab='Depth',
                xlab='',
                lower=agg$p.q25, upper=agg$p.q75, ylim=c(800,-2),
                panel=panel.depth_function,
                alpha=0.25, sync.colors=TRUE,
                par.settings=list(superpose.line=list(col='RoyalBlue', lwd=2)),
                prepanel=prepanel.depth_function,
                #cf=agg$contributing_fraction, cf.col='black', cf.interval=20, 
                layout=c(1,1), strip=strip.custom(bg=grey(0.8)),
                scales=list(x=list(tick.number=4, alternating=3, relation='free'))
)

Figure2

pdf("EdgeroilogSOCAgg.pdf",width=4,height=6)
plot(Figure2) # Make plot
dev.off()
#=================================================================== 



#formula
formulaString <- as.formula(paste(paste("PHIHO5 ~ "), paste(c(names(cov.maps),"depth"), collapse="+")))
seed = 1

# pH results
BaseL.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)#, kmean.vars = all.vars(formulaString), cum.prop = 0.9)    # , kmean.vars = all.vars(formulaString), cum.prop = 0.9
BaseL.pH.ncv.time <- system.time(BaseL.pH.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.pH.preproc, lambda = seq(0,0.2,0.001), w = NULL)) #seq(0.1, 1, 0.1)
BaseL.pH.time <- system.time(BaseL.pH <- sparsereg3D.sel(sparse.reg = BaseL.pH.preproc ,lambda = seq(0,0.1,0.001)))
#pH.l.pred <- sparsereg3D.pred(model.info = BaseL.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

BaseP.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)#, kmean.vars = all.vars(formulaString), cum.prop = 0.9)    
BaseP.pH.ncv.time <- system.time(BaseP.pH.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.pH.preproc, lambda = seq(0,0.2,0.001), w = NULL))
BaseP.pH.time <- system.time(BaseP.pH <- sparsereg3D.sel(sparse.reg = BaseP.pH.preproc ,lambda = seq(0,0.2,0.001)))
#pH.p.pred <- sparsereg3D.pred(model.info = pH.p.sel, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntL.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)#, kmean.vars = all.vars(formulaString), cum.prop = 0.9)    
IntL.pH.ncv.time <- system.time(IntL.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntL.pH.preproc, lambda = seq(0,0.2,0.001), w = NULL))
IntL.pH.time <- system.time(IntL.pH <- sparsereg3D.sel(sparse.reg = IntL.pH.preproc ,lambda = seq(0,0.2,0.001)))
#pH.l.pred <- sparsereg3D.pred(model.info = IntL.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntP.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)#, kmean.vars = all.vars(formulaString), cum.prop = 0.9)    
IntP.pH.ncv.time <- system.time(IntP.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntP.pH.preproc, lambda = seq(0,0.2,0.001), w = NULL))
IntP.pH.time <- system.time(IntP.pH <- sparsereg3D.sel(sparse.reg = IntP.pH.preproc ,lambda = seq(0,0.2,0.001)))
pH.l.pred <- sparsereg3D.pred(model.info = IntP.pH, chunk.size = 20000, grids = cov.maps, depths = c(-0.05, -0.15, -0.30))

##### Writing prediction to file ##########################################
writeGDAL(pH.l.pred[[1]], paste("D:/R_projects","Edg.pH.pred1.tiff",sep="/"),drivername = "GTiff")
writeGDAL(pH.l.pred[[2]], paste("D:/R_projects","Edg.pH.pred2.tiff",sep="/"),drivername = "GTiff")
writeGDAL(pH.l.pred[[3]], paste("D:/R_projects","Edg.pH.pred3.tiff",sep="/"),drivername = "GTiff")


rbind(BaseL.pH.ncv[1:2], BaseP.pH.ncv[1:2], IntL.pH.ncv[1:2], IntP.pH.ncv[1:2])

IntP.pH.ncv[1:2]

library(doParallel)

registerDoParallel(cores=4)

pack = (.packages())


result = foreach (j = 1:2, .packages = pack) %dopar% { 
  
  if (j==1){
    IntHL.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)#, kmean.vars = all.vars(formulaString), cum.prop = 0.90)    
    IntHL.pH.ncv.time <- system.time(IntHL.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.pH.preproc, lambda = seq(0,0.2,0.001)))
    IntHL.pH.time <- system.time(IntHL.pH <- sparsereg3D.sel(sparse.reg = IntHL.pH.preproc ,lambda = seq(0,0.2,0.001)))
    #pH.l.pred <- sparsereg3D.pred(model.info = IntHL.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))
    list(IntHL.pH.preproc, IntHL.pH.ncv.time, IntHL.pH.ncv, IntHL.pH.time, IntHL.pH)
  } else if (j==2) {
    IntHP.pH.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)#, kmean.vars = all.vars(formulaString), cum.prop = 0.90)    
    IntHP.pH.ncv.time <- system.time(IntHP.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.pH.preproc, lambda = seq(0,0.2,0.001)))
    IntHP.pH.time <- system.time(IntHP.pH <- sparsereg3D.sel(sparse.reg = IntHP.pH.preproc ,lambda = seq(0,0.2,0.001)))
    #pH.l.pred <- sparsereg3D.pred(model.info = IntHP.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))
    list(IntHP.pH.preproc, IntHP.pH.ncv.time, IntHP.pH.ncv, IntHP.pH.time, IntHP.pH)
  }
  
}

stopImplicitCluster()

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

#Results
ll <- length(IntL.pH$coefficients)
pp <- length(IntHL.pH$coefficients[,1])+1

cmL.pH <- data.frame(variable=IntHL.pH$coefficients[,1], BaseL.pH.me=BaseL.pH$coefficients[2:pp], IntL.pH.me=IntL.pH$coefficients[2:pp],IntL.pH.ie=c(IntL.pH$coefficients[(pp+1):ll],0),IntHL.pH.me=IntHL.pH$coefficients[,2],IntHL.pH.ie=IntHL.pH$coefficients[,3] )
stargazer(cmL.pH, summary = FALSE, digits = 2, type = "latex")

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

#do.call(cbind, lapply(BaseL.pH.ncv$models, function(x) x[[1]]))

pH.n.coeffs <- n.coeffs(l.coeffs = cmL.pH, p.coeffs = cmP.pH)

#save.image(file = "D:/R_projects/edgeroi_pH.RData")

#================= log ORC ==========================================================


#formula
formulaString <- as.formula(paste(paste("log1pORC ~ "), paste(c(names(cov.maps),"depth"), collapse="+")))
formulaString

seed = 1
source(paste(fun.path,"sparsereg3D.ncv.r",sep="/"))
source(paste(fun.path,"stratfold3d.r",sep="/"))
source(paste(fun.path,"pre.sparsereg3D.r",sep="/"))
#sparsereg3D
# logORC results
BaseL.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)#, kmean.vars = all.vars(formulaString), cum.prop = 0.90)     # , kmean.vars = all.vars(formulaString), cum.prop = 0.90
BaseL.logORC.ncv.time <- system.time(BaseL.logORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.logORC.preproc, lambda = seq(0,0.2,0.001))) #w = seq(0, 1, 0.1)
BaseL.logORC.time <- system.time(BaseL.logORC <- sparsereg3D.sel(sparse.reg = BaseL.logORC.preproc ,lambda = seq(0,0.2,0.001)))
#logORC.l.pred <- sparsereg3D.pred(model.info = BaseL.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

BaseP.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)#, kmean.vars = all.vars(formulaString), cum.prop = 0.90)    
BaseP.logORC.ncv.time <- system.time(BaseP.logORC.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.logORC.preproc, lambda = seq(0,0.2,0.001)))
BaseP.logORC.time <- system.time(BaseP.logORC <- sparsereg3D.sel(sparse.reg = BaseP.logORC.preproc ,lambda = seq(0,2,0.1)))
#logORC.p.pred <- sparsereg3D.pred(model.info = logORC.p.sel, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntL.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)#, kmean.vars = all.vars(formulaString), cum.prop = 0.90)    
IntL.logORC.ncv.time <- system.time(IntL.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntL.logORC.preproc, lambda = seq(0,0.2,0.001), w = seq(0, 1, 0.1)))
IntL.logORC.time <- system.time(IntL.logORC <- sparsereg3D.sel(sparse.reg = IntL.logORC.preproc ,lambda = seq(0,0.2,0.001))) 
#logORC.l.pred <- sparsereg3D.pred(model.info = IntL.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1,-0.2,-0.3))

IntP.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)#, kmean.vars = all.vars(formulaString), cum.prop = 0.90)    
IntP.logORC.ncv.time <- system.time(IntP.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntP.logORC.preproc, lambda = seq(0,0.2,0.001), w = seq(0, 5, 0.2)))
IntP.logORC.time <- system.time(IntP.logORC <- sparsereg3D.sel(sparse.reg = IntP.logORC.preproc ,lambda = seq(0,0.2,0.001)))
#logORC.l.pred <- sparsereg3D.pred(model.info = IntP.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.05, -0.15, -0.3))

logORC.l.pred[[1]]$pred <- expm1(logORC.l.pred[[1]]$pred)
logORC.l.pred[[2]]$pred <- expm1(logORC.l.pred[[2]]$pred)
logORC.l.pred[[3]]$pred <- expm1(logORC.l.pred[[3]]$pred)

##### Writing prediction to file ##########################################
writeGDAL(logORC.l.pred[[1]], paste("D:/R_projects","Edg.logSOC.pred1.tiff",sep="/"),drivername = "GTiff")
writeGDAL(logORC.l.pred[[2]], paste("D:/R_projects","Edg.logSOC.pred2.tiff",sep="/"),drivername = "GTiff")
writeGDAL(logORC.l.pred[[3]], paste("D:/R_projects","Edg.logSOC.pred3.tiff",sep="/"),drivername = "GTiff")

#rbind(BaseL.logORC.ncv[1:2], BaseP.logORC.ncv[1:2], IntL.logORC.ncv[1:2], IntP.logORC.ncv[1:2])

library(doParallel)

registerDoParallel(cores=4)

pack = (.packages())


SOCresult = foreach (j = 1:2, .packages = pack) %dopar% { 
  
  if (j==1){
    IntHL.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)#, kmean.vars = all.vars(formulaString), cum.prop = 0.90)    
    IntHL.logORC.ncv.time <- system.time(IntHL.logORC.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.logORC.preproc, lambda = seq(0,0.2,0.001)))
    IntHL.logORC.time <- system.time(IntHL.logORC <- sparsereg3D.sel(sparse.reg = IntHL.logORC.preproc ,lambda = seq(0,0.2,0.001)))
    #logORC.l.pred <- sparsereg3D.pred(model.info = IntHL.logORC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))
    list(IntHL.logORC.preproc, IntHL.logORC.ncv.time, IntHL.logORC.ncv, IntHL.logORC.time, IntHL.logORC)
  } else if (j==2) {
    IntHP.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = TRUE, profiles = edgeroi.spc, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = seed)#, kmean.vars = all.vars(formulaString), cum.prop = 0.90)    
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

stargazer(cmL.logORC, summary = FALSE, digits = 2, type = "latex")

#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.logORC$coefficients)
p <- length(IntHP.logORC$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.logORC <- data.frame(variable=IntHP.logORC$coefficients[,1], BaseP.logORC.me=BaseP.logORC$coefficients[2:p], IntP.logORC.me=IntP.logORC$coefficients[2:p],IntP.logORC.ie1=c(IntP.logORC$coefficients[(p+1):l][i1],0,0,0),IntP.logORC.ie2=c(IntP.logORC$coefficients[(p+1):l][i2],0,0,0),IntP.logORC.ie3=c(IntP.logORC$coefficients[(p+1):l][i3],0,0,0),IntHP.logORC.me=IntHP.logORC$coefficients[,2],IntHP.logORC.ie1=IntHP.logORC$coefficients[,3],IntHP.logORC.ie2=IntHP.logORC$coefficients[,4],IntHP.logORC.ie3=IntHP.logORC$coefficients[,5] )
cmlogORC <- cmP.logORC[,c(1,7:10)]

stargazer(cmP.logORC, summary = FALSE, digits = 2, type = "latex")

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

n.coeffs <- function(l.coeffs,p.coeffs){
  BaseL <- data.frame(l.coeffs[,2])
  IntL <- l.coeffs[,3:4]
  IntHL <- l.coeffs[,5:6]
  BaseP <- data.frame(p.coeffs[,2])
  IntP <- p.coeffs[,3:6]
  IntHP <- p.coeffs[,7:10]
  
  n.BaseL <-  data.frame(total = prod(dim(BaseL)), selected = sum(BaseL!=0 ), "Number of Variables" = sum(apply(BaseL,1, function(y) sum(abs(y)))!=0), "main effects" = sum(BaseL!=0 ), "interaction effects" = 0 )
  n.IntL <-   data.frame(total = prod(dim(IntL)),  selected = sum(apply(IntL,2, function(y) sum(abs(y)!=0))), "Number of Variables" = sum(apply(IntL,1, function(y) sum(abs(y)))!=0) ,"main effects" = apply(IntL,2, function(y) sum(abs(y)!=0))[1], "interaction effects" = apply(IntL,2, function(y) sum(abs(y)!=0))[2])
  n.IntHL <-   data.frame(total = prod(dim(IntHL)),  selected = sum(apply(IntHL,2, function(y) sum(abs(y)!=0))), "Number of Variables" = sum(apply(IntHL,1, function(y) sum(abs(y)))!=0), "main effects" = apply(IntHL,2, function(y) sum(abs(y)!=0))[1], "interaction effects" = apply(IntHL,2, function(y) sum(abs(y)!=0))[2])
  
  n.BaseP <-  data.frame(total = prod(dim(BaseP)), selected = sum(BaseP!=0 ), "Number of Variables" = sum(apply(BaseP,1, function(y) sum(abs(y)))!=0), "main effects" = sum(BaseP!=0 ), "interaction effects" = 0 )
  n.IntP <-   data.frame(total = prod(dim(IntP)),  selected = sum(apply(IntP,2, function(y) sum(abs(y)!=0))), "Number of Variables" = sum(apply(IntP[-c(dim(IntP)[1]:(dim(IntP)[1]-1)),], 1, function(y) sum(y))!=0), "main effects" = apply(IntP,2, function(y) sum(abs(y)!=0))[1], "interaction effects" = sum(apply(IntP,2, function(y) sum(abs(y)!=0))[2:4]))
  n.IntHP <-   data.frame(total = prod(dim(IntHP)),  selected = sum(apply(IntHP,2, function(y) sum(abs(y)!=0))), "Number of Variables" = sum(apply(IntHP[-c(dim(IntHP)[1]:(dim(IntHP)[1]-1)),], 1, function(y) sum(abs(y)))!=0), "main effects" = apply(IntHP,2, function(y) sum(abs(y)!=0))[1], "interaction effects" = sum(apply(IntHP,2, function(y) sum(abs(y)!=0))[2:4]))
  
  
  total <- data.frame(rbind(n.BaseL,n.BaseP,n.IntL,n.IntP,n.IntHL,n.IntHP))
  rownames(total) <- NULL
  total <- cbind(model = c("BaseL","BaseP", "IntL", "IntP","IntHL","IntHP"), total)
  return(total)
}

logORC.n.coeffs <- n.coeffs(l.coeffs = cmL.logORC, p.coeffs = cmP.logORC)
pH.n.coeffs <- n.coeffs(l.coeffs = cmL.pH, p.coeffs = cmP.pH)

#do.call(cbind, lapply(BaseL.logORC.ncv$models, function(x) x[[1]]))


stargazer(cbind(logORC.n.coeffs, logORC.ncv[,-1], pH.n.coeffs[,-1], pH.ncv[,-1]), summary = FALSE, digits = 2, type = "latex")


#load(file = "D:/R_projects/edgeroi_logORC_pH.RData")



logSOC.formula <- as.formula(paste(paste("log1pORC ~ "), paste(c(names(cov.maps),"depth"), collapse="+")))
pH.formula <- as.formula(paste(paste("PHIHO5 ~ "), paste(c(names(cov.maps),"depth"), collapse="+")))


intbasedif <- function(formula, depth.inc, profiles, cov.grids, poly = 3, num.folds = 5, num.means = 3, poly.deg = poly.deg, seed = 1, lambda.seq = seq(0,0.2,0.001), weights = seq(0,5,0.2)){
  variable <- all.vars(formula)[1]
  seed = seed
  max.depth = max(profiles$Bottom)
  depth.seq = c(0, sort(seq(depth.inc, (round(max.depth/depth.inc)-1)*depth.inc, depth.inc), decreasing = FALSE))
  pH.seq = data.frame()
  IntP.seq = data.frame()
  BaseP.seq = data.frame()
  profiles.i = profiles

  
  for(i in 1:length(depth.seq)){
    profiles.i@horizons <- profiles@horizons[profiles@horizons$Bottom < round(max.depth - depth.seq[i]),]
    
    IntP.preproc <- pre.sparsereg3D(base.model = formula, use.hier = FALSE, profiles = profiles.i, use.interactions = TRUE, poly.deg = poly.deg, num.folds = num.folds, num.means = num.means, cov.grids = cov.grids, seed = seed)
    IntP.seq <- rbind(IntP.seq, sparsereg3D.ncv(sparse.reg = IntP.preproc, lambda = lambda.seq, w = weights)[1:2])
    
    BaseP.preproc <- pre.sparsereg3D(base.model = formula, use.hier = FALSE, profiles = profiles.i, use.interactions = FALSE, poly.deg = poly.deg, num.folds = num.folds, num.means = num.means, cov.grids = cov.maps, seed = seed) 
    BaseP.seq <- rbind(BaseP.seq, sparsereg3D.ncv(sparse.reg = BaseP.preproc, lambda = lambda.seq, w = weights)[1:2])
  }
    IntP.data <- data.frame(Depth = c(max.depth, sort(depth.seq[-1], decreasing = TRUE)), RMSE = IntP.seq$RMSE, variable = rep("IntP",length(depth.seq)))
    BaseP.data <- data.frame(Depth = c(max.depth, sort(depth.seq[-1], decreasing = TRUE)), RMSE = BaseP.seq$RMSE, variable = rep("BaseP",length(depth.seq)))
    
    
    BIP.data <- rbind(BaseP.data,IntP.data)
    BIP.data$variable <- factor(BIP.data$variable)
    plot <- qplot(Depth, RMSE, data = BIP.data, geom = c("point", "line"), color = variable)+theme_bw()+theme(legend.text = element_text(size = 12)) + theme(axis.text = element_text(size=12), axis.title = element_text(size = 14, face = "bold")) + labs(x = "Depth [cm]",y = "RMSE")
    
  return(plot)
}


logSOC.dif <- intbasedif(formula = logSOC.formula, depth.inc = 100, profiles = edgeroi.spc, cov.grids = cov.maps, num.folds = 5, num.means = 3, poly.deg = 3, seed = 1, lambda.seq = seq(0,0.2,0.001), weights = seq(0,5,0.2))
pH.dif <- intbasedif(formula = pH.formula, depth.inc = 100, profiles = edgeroi.spc, cov.grids = cov.maps, num.folds = 5, num.means = 3, poly.deg = 3, seed = 1, lambda.seq = seq(0,0.2,0.001), weights = NULL)



pdf("logSOC.wdif.pdf",width=10,height=6)
logSOC.dif
dev.off()


pdf("pH.dif.pdf",width=10,height=6)
pH.dif
dev.off()

save.image(file = "D:/R_projects/edgeroi_logORC_pH.RData")

load("D:/R_projects/edgeroi_logORC_pH.RData")


