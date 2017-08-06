# Napisati nesto ovde
#
#
path1 <- "C:/Users/Milutin/Dropbox/Extensions of soil 3D trend models/Data and Scripts"



library(rgdal)
library(GSIF)
library(gdalUtils)
library(raster)
library(plyr)
library(aqp)
#library(psych)
#library(mda)
library(classInt)
library(caret)
library(MASS)
library(splines)
library(glmnet)
library(hierNet)
library(magrittr)
library(doParallel)
library(foreach)
#library(stargazer)
library(gstat)

load(paste(path1,"BorData.rda",sep = "/"))
load(paste(path1,"covmaps.rda",sep = "/"))

fun.path <- "D:/_R projects/sparsereg3D/R"
source(paste(fun.path,"stratfold3d.r",sep="/"))
source(paste(fun.path,"pre.sparsereg3D.r",sep="/"))
source(paste(fun.path,"sparsereg3D.ncv.r",sep="/"))
source(paste(fun.path,"sparsereg3D.pred.r",sep="/"))
source(paste(fun.path,"sparsereg3D.sel.r",sep="/"))

#================== Spatial references ===============================================================
gk_7 <- "+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 +y_0=0 +ellps=bessel +towgs84=574.027,170.175,401.545,4.88786,-0.66524,-13.24673,0.99999311067 +units=m"
utm <- "+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m"

#================== Names and abbrevations of covariates =================================================

CovNames <- c("DEM","Aspect","Slope","TWI","ConvInd","CrSectCurv","LongCurv","ChNetBLevel","VDistChNet","NegOp", "PosOp","WEeast","WEnw","clc","SoilType")

#=================== DATA ============================================================
bor.profs <- bor[,c("ID","x","y","Soil.Type","Top","Bottom","As","pH","Co","SOM")]
bor.profs$logSOM <- log(bor.profs$SOM)
bor.profs$logAs <- log(bor.profs$As)
bor.profs$ORCDRC <- bor.profs$SOM*10/1.724
bor.profs$logORCDRC <- log(bor.profs$ORCDRC)

bor.profs$mid <- (bor.profs$Top+bor.profs$Bottom)/2
depths(bor.profs) <- ID ~ Top + Bottom
site(bor.profs) <- ~ Soil.Type + x + y
coordinates(bor.profs) <- ~ x+y
proj4string(bor.profs) <- CRS(utm)
#=====================================================================================

summary(bor.profs@horizons[,c("logORCDRC","ORCDRC","pH","Top","Bottom")])

bor.profs@site$depth <- profileApply(bor.profs, estimateSoilDepth, name='As', top='Top', bottom='Bottom')
bor.profs@site$n <- (ddply(bor.profs@horizons,.(ID), summarize, br = length(As)))[,2]

sum.n.d <- summary(bor.profs@site[,c("depth","n")])
sum.As.SOM.pH <- summary(bor.profs@horizons[,c("As","SOM","pH")])[,1]
rbind(sum.As.SOM.pH,sum.n.d)



# Formulas

ORCDRC.fun <- as.formula(paste("ORCDRC ~", paste(c(CovNames,"depth"), collapse="+")))
log.ORCDRC.fun <- as.formula(paste("logORCDRC ~", paste(c(CovNames,"depth"), collapse="+")))
pH.fun <- as.formula(paste("pH ~", paste(c(CovNames,"depth"), collapse="+")))

# ORCDRC results
BaseL.ORCDRC.preproc <- pre.sparsereg3D(base.model = ORCDRC.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
BaseL.ORCDRC.ncv.time <- system.time(BaseL.ORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.ORCDRC.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.ORCDRC.time <- system.time(BaseL.ORCDRC <- sparsereg3D.sel(sparse.reg = BaseL.ORCDRC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORCDRC.l.pred <- sparsereg3D.pred(model.info = BaseL.ORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

BaseP.ORCDRC.preproc <- pre.sparsereg3D(base.model = ORCDRC.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
BaseP.ORCDRC.ncv.time <- system.time(BaseP.ORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.ORCDRC.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.ORCDRC.time <- system.time(BaseP.ORCDRC <- sparsereg3D.sel(sparse.reg = BaseP.ORCDRC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORCDRC.p.pred <- sparsereg3D.pred(model.info = ORCDRC.p.sel, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


IntL.ORCDRC.preproc <- pre.sparsereg3D(base.model = ORCDRC.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
IntL.ORCDRC.ncv.time <- system.time(IntL.ORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = IntL.ORCDRC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.ORCDRC.time <- system.time(IntL.ORCDRC <- sparsereg3D.sel(sparse.reg = IntL.ORCDRC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORCDRC.l.pred <- sparsereg3D.pred(model.info = IntL.ORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntP.ORCDRC.preproc <- pre.sparsereg3D(base.model = ORCDRC.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
IntP.ORCDRC.ncv.time <- system.time(IntP.ORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = IntP.ORCDRC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.ORCDRC.time <- system.time(IntP.ORCDRC <- sparsereg3D.sel(sparse.reg = IntP.ORCDRC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORCDRC.l.pred <- sparsereg3D.pred(model.info = IntP.ORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHL.ORCDRC.preproc <- pre.sparsereg3D(base.model = ORCDRC.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
IntHL.ORCDRC.ncv.time <- system.time(IntHL.ORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.ORCDRC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHL.ORCDRC.time <- system.time(IntHL.ORCDRC <- sparsereg3D.sel(sparse.reg = IntHL.ORCDRC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORCDRC.l.pred <- sparsereg3D.pred(model.info = IntHL.ORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHP.ORCDRC.preproc <- pre.sparsereg3D(base.model = ORCDRC.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
IntHP.ORCDRC.ncv.time <- system.time(IntHP.ORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.ORCDRC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.ORCDRC.time <- system.time(IntHP.ORCDRC <- sparsereg3D.sel(sparse.reg = IntHP.ORCDRC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#ORCDRC.l.pred <- sparsereg3D.pred(model.info = IntHP.ORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

#============================= Coefficients tables for ORCDRC models ==============================================================================

ll <- length(IntL.ORCDRC$coefficients)
pp <- length(IntHL.ORCDRC$coefficients[,1])+1

cmL.ORCDRC <- data.frame(variable=IntHL.ORCDRC$coefficients[,1], BaseL.ORCDRC.me=BaseL.ORCDRC$coefficients[2:pp], IntL.ORCDRC.me=IntL.ORCDRC$coefficients[2:pp],IntL.ORCDRC.ie=c(IntL.ORCDRC$coefficients[(pp+1):ll],0),IntHL.ORCDRC.me=IntHL.ORCDRC$coefficients[,2],IntHL.ORCDRC.ie=IntHL.ORCDRC$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.ORCDRC$coefficients)
p <- length(IntHP.ORCDRC$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.ORCDRC <- data.frame(variable=IntHP.ORCDRC$coefficients[,1], BaseP.ORCDRC.me=BaseP.ORCDRC$coefficients[2:p], IntP.ORCDRC.me=IntP.ORCDRC$coefficients[2:p],IntP.ORCDRC.ie1=c(IntP.ORCDRC$coefficients[(p+1):l][i1],0,0,0),IntP.ORCDRC.ie2=c(IntP.ORCDRC$coefficients[(p+1):l][i2],0,0,0),IntP.ORCDRC.ie3=c(IntP.ORCDRC$coefficients[(p+1):l][i3],0,0,0),IntHP.ORCDRC.me=IntHP.ORCDRC$coefficients[,2],IntHP.ORCDRC.ie1=IntHP.ORCDRC$coefficients[,3],IntHP.ORCDRC.ie2=IntHP.ORCDRC$coefficients[,4],IntHP.ORCDRC.ie3=IntHP.ORCDRC$coefficients[,5] )
cmORCDRC <- cmP.ORCDRC[,c(1,7:10)]

# Models comparison
ORCDRC.ncv <- data.frame(rbind(BaseL = BaseL.ORCDRC.ncv, BaseP = BaseP.ORCDRC.ncv, IntL = IntL.ORCDRC.ncv, IntP = IntP.ORCDRC.ncv, IntHL = IntHL.ORCDRC.ncv, IntHP = IntHP.ORCDRC.ncv))
ORCDRC.ncv <- data.frame(Model = rownames(ORCDRC.ncv), ORCDRC.ncv)
rownames(ORCDRC.ncv) <- NULL
names(ORCDRC.ncv) <- c("Model","RMSE","R squared")
stargazer(ORCDRC.ncv, summary = FALSE, digits = 2, type = "text")

#Time table

ORCDRC.ncv.time <- rbind(BaseL.ORCDRC.ncv.time, BaseP.ORCDRC.ncv.time, IntL.ORCDRC.ncv.time, IntP.ORCDRC.ncv.time, IntHL.ORCDRC.ncv.time, IntHP.ORCDRC.ncv.time)
ORCDRC.time <- rbind(BaseL.ORCDRC.time, BaseP.ORCDRC.time, IntL.ORCDRC.time, IntP.ORCDRC.time, IntHL.ORCDRC.time, IntHP.ORCDRC.time)

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

ORCDRC.n.coeffs <- n.coeffs(l.coeffs = cmL.ORCDRC, p.coeffs = cmP.ORCDRC)

#ORCDRC.results <- list(ORCDRC.ncv = ORCDRC.ncv, l.coeffs = cmL.ORCDRC, p.coeffs = cmP.ORCDRC, IntP.ORCDRC = IntP.ORCDRC, IntHP.ORCDRC = IntHP.ORCDRC, BaseP.ORCDRC = BaseP.ORCDRC)
#save(ORCDRC.results, file = "ORCDRC.results.rda" )
#==============================================================================================================================================


# logORCDRC results
BaseL.logORCDRC.preproc <- pre.sparsereg3D(base.model = log.ORCDRC.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
BaseL.logORCDRC.ncv.time <- system.time(BaseL.logORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.logORCDRC.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.logORCDRC.time <- system.time(BaseL.logORCDRC <- sparsereg3D.sel(sparse.reg = BaseL.logORCDRC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#logORCDRC.l.pred <- sparsereg3D.pred(model.info = BaseL.logORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

BaseP.logORCDRC.preproc <- pre.sparsereg3D(base.model = log.ORCDRC.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
BaseP.logORCDRC.ncv.time <- system.time(BaseP.logORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.logORCDRC.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.logORCDRC.time <- system.time(BaseP.logORCDRC <- sparsereg3D.sel(sparse.reg = BaseP.logORCDRC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#logORCDRC.p.pred <- sparsereg3D.pred(model.info = logORCDRC.p.sel, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


IntL.logORCDRC.preproc <- pre.sparsereg3D(base.model = log.ORCDRC.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
IntL.logORCDRC.ncv.time <- system.time(IntL.logORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = IntL.logORCDRC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.logORCDRC.time <- system.time(IntL.logORCDRC <- sparsereg3D.sel(sparse.reg = IntL.logORCDRC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#logORCDRC.l.pred <- sparsereg3D.pred(model.info = IntL.logORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntP.logORCDRC.preproc <- pre.sparsereg3D(base.model = log.ORCDRC.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
IntP.logORCDRC.ncv.time <- system.time(IntP.logORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = IntP.logORCDRC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.logORCDRC.time <- system.time(IntP.logORCDRC <- sparsereg3D.sel(sparse.reg = IntP.logORCDRC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#logORCDRC.l.pred <- sparsereg3D.pred(model.info = IntP.logORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHL.logORCDRC.preproc <- pre.sparsereg3D(base.model = log.ORCDRC.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
IntHL.logORCDRC.ncv.time <- system.time(IntHL.logORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.logORCDRC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHL.logORCDRC.time <- system.time(IntHL.logORCDRC <- sparsereg3D.sel(sparse.reg = IntHL.logORCDRC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#logORCDRC.l.pred <- sparsereg3D.pred(model.info = IntHL.logORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHP.logORCDRC.preproc <- pre.sparsereg3D(base.model = log.ORCDRC.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
IntHP.logORCDRC.ncv.time <- system.time(IntHP.logORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.logORCDRC.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.logORCDRC.time <- system.time(IntHP.logORCDRC <- sparsereg3D.sel(sparse.reg = IntHP.logORCDRC.preproc ,lambda = seq(0,5,0.1), seed = 321))
#logORCDRC.l.pred <- sparsereg3D.pred(model.info = IntHP.logORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

#============================= Coefficients tables for logORCDRC models ==============================================================================

ll <- length(IntL.logORCDRC$coefficients)
pp <- length(IntHL.logORCDRC$coefficients[,1])+1

cmL.logORCDRC <- data.frame(variable=IntHL.logORCDRC$coefficients[,1], BaseL.logORCDRC.me=BaseL.logORCDRC$coefficients[2:pp], IntL.logORCDRC.me=IntL.logORCDRC$coefficients[2:pp],IntL.logORCDRC.ie=c(IntL.logORCDRC$coefficients[(pp+1):ll],0),IntHL.logORCDRC.me=IntHL.logORCDRC$coefficients[,2],IntHL.logORCDRC.ie=IntHL.logORCDRC$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.logORCDRC$coefficients)
p <- length(IntHP.logORCDRC$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.logORCDRC <- data.frame(variable=IntHP.logORCDRC$coefficients[,1], BaseP.logORCDRC.me=BaseP.logORCDRC$coefficients[2:p], IntP.logORCDRC.me=IntP.logORCDRC$coefficients[2:p],IntP.logORCDRC.ie1=c(IntP.logORCDRC$coefficients[(p+1):l][i1],0,0,0),IntP.logORCDRC.ie2=c(IntP.logORCDRC$coefficients[(p+1):l][i2],0,0,0),IntP.logORCDRC.ie3=c(IntP.logORCDRC$coefficients[(p+1):l][i3],0,0,0),IntHP.logORCDRC.me=IntHP.logORCDRC$coefficients[,2],IntHP.logORCDRC.ie1=IntHP.logORCDRC$coefficients[,3],IntHP.logORCDRC.ie2=IntHP.logORCDRC$coefficients[,4],IntHP.logORCDRC.ie3=IntHP.logORCDRC$coefficients[,5] )
cmlogORCDRC <- cmP.logORCDRC[,c(1,7:10)]

# Models comparison
logORCDRC.ncv <- data.frame(rbind(BaseL = BaseL.logORCDRC.ncv, BaseP = BaseP.logORCDRC.ncv, IntL = IntL.logORCDRC.ncv, IntP = IntP.logORCDRC.ncv, IntHL = IntHL.logORCDRC.ncv, IntHP = IntHP.logORCDRC.ncv))
logORCDRC.ncv <- data.frame(Model = rownames(logORCDRC.ncv), logORCDRC.ncv)
rownames(logORCDRC.ncv) <- NULL
names(logORCDRC.ncv) <- c("Model","RMSE","R squared")
stargazer(logORCDRC.ncv, summary = FALSE, digits = 2, type = "text")

#Time table

logORCDRC.ncv.time <- rbind(BaseL.logORCDRC.ncv.time, BaseP.logORCDRC.ncv.time, IntL.logORCDRC.ncv.time, IntP.logORCDRC.ncv.time, IntHL.logORCDRC.ncv.time, IntHP.logORCDRC.ncv.time)
logORCDRC.time <- rbind(BaseL.logORCDRC.time, BaseP.logORCDRC.time, IntL.logORCDRC.time, IntP.logORCDRC.time, IntHL.logORCDRC.time, IntHP.logORCDRC.time)

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

logORCDRC.n.coeffs <- n.coeffs(l.coeffs = cmL.logORCDRC, p.coeffs = cmP.logORCDRC)

#logORCDRC.results <- list(logORCDRC.ncv = logORCDRC.ncv, l.coeffs = cmL.logORCDRC, p.coeffs = cmP.logORCDRC, IntP.logORCDRC = IntP.logORCDRC, IntHP.logORCDRC = IntHP.logORCDRC, BaseP.logORCDRC = BaseP.logORCDRC)
#save(logORCDRC.results, file = "logORCDRC.results.rda" )
#==============================================================================================================================================




# pH results
BaseL.pH.preproc <- pre.sparsereg3D(base.model = pH.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
BaseL.pH.ncv.time <- system.time(BaseL.pH.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.pH.time <- system.time(BaseL.pH <- sparsereg3D.sel(sparse.reg = BaseL.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = BaseL.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

BaseP.pH.preproc <- pre.sparsereg3D(base.model = pH.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
BaseP.pH.ncv.time <- system.time(BaseP.pH.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.pH.time <- system.time(BaseP.pH <- sparsereg3D.sel(sparse.reg = BaseP.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.p.pred <- sparsereg3D.pred(model.info = pH.p.sel, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


IntL.pH.preproc <- pre.sparsereg3D(base.model = pH.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
IntL.pH.ncv.time <- system.time(IntL.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntL.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.pH.time <- system.time(IntL.pH <- sparsereg3D.sel(sparse.reg = IntL.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = IntL.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntP.pH.preproc <- pre.sparsereg3D(base.model = pH.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
IntP.pH.ncv.time <- system.time(IntP.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntP.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.pH.time <- system.time(IntP.pH <- sparsereg3D.sel(sparse.reg = IntP.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = IntP.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHL.pH.preproc <- pre.sparsereg3D(base.model = pH.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
IntHL.pH.ncv.time <- system.time(IntHL.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHL.pH.time <- system.time(IntHL.pH <- sparsereg3D.sel(sparse.reg = IntHL.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = IntHL.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHP.pH.preproc <- pre.sparsereg3D(base.model = pH.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)
IntHP.pH.ncv.time <- system.time(IntHP.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.pH.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.pH.time <- system.time(IntHP.pH <- sparsereg3D.sel(sparse.reg = IntHP.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
#pH.l.pred <- sparsereg3D.pred(model.info = IntHP.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


#============================= Coefficients tables for pH models ==============================================================================

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

#==============================================================================================================================================

# Models comparison
pH.ncv <- data.frame(rbind(BaseL = BaseL.pH.ncv, BaseP = BaseP.pH.ncv, IntL = IntL.pH.ncv, IntP = IntP.pH.ncv, IntHL = IntHL.pH.ncv, IntHP = IntHP.pH.ncv))
pH.ncv <- data.frame(Model = rownames(pH.ncv), pH.ncv)
rownames(pH.ncv) <- NULL
names(pH.ncv) <- c("Model","RMSE","R squared")
stargazer(pH.ncv, summary = FALSE, digits = 2, type = "text")

pH.results <- list(pH.ncv = pH.ncv, l.coeffs = cmL.pH, p.coeffs = cmP.pH, IntP.pH = IntP.pH, IntHP.pH = IntHP.pH, BaseP.pH = BaseP.pH)

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


#save(pH.results, file = "pH.results.rda" )



















#============================= Time Table =======================================================================================================

As.ncv.time <- rbind(BaseL.As.ncv.time, BaseP.As.ncv.time, IntL.As.ncv.time, IntP.As.ncv.time, IntHL.As.ncv.time, IntHP.As.ncv.time)
As.time <- rbind(BaseL.As.time, BaseP.As.time, IntL.As.time, IntP.As.time, IntHL.As.time, IntHP.As.time)

SOM.ncv.time <- rbind(BaseL.SOM.ncv.time, BaseP.SOM.ncv.time, IntL.SOM.ncv.time, IntP.SOM.ncv.time, IntHL.SOM.ncv.time, IntHP.SOM.ncv.time)
SOM.time <- rbind(BaseL.SOM.time, BaseP.SOM.time, IntL.SOM.time, IntP.SOM.time, IntHL.SOM.time, IntHP.SOM.time)

pH.ncv.time <- rbind(BaseL.pH.ncv.time, BaseP.pH.ncv.time, IntL.pH.ncv.time, IntP.pH.ncv.time, IntHL.pH.ncv.time, IntHP.pH.ncv.time)
pH.time <- rbind(BaseL.pH.time, BaseP.pH.time, IntL.pH.time, IntP.pH.time, IntHL.pH.time, IntHP.pH.time)

#================================================================================================================================================

model = data.frame(Model = c("BaseL","BaseP","IntL","IntP","IntHL","IntHP"))

ORCDRC.ncv.time <- data.frame(Assessment.ORCDRC = rbind(BaseL.ORCDRC.ncv.time,BaseP.ORCDRC.ncv.time, IntL.ORCDRC.ncv.time, IntP.ORCDRC.ncv.time, IntHL.ORCDRC.ncv.time,IntHP.ORCDRC.ncv.time))
ORCDRC.time <- data.frame(Training.ORCDRC = rbind(BaseL.ORCDRC.time,BaseP.ORCDRC.time, IntL.ORCDRC.time, IntP.ORCDRC.time, IntHL.ORCDRC.time,IntHP.ORCDRC.time))

logORCDRC.ncv.time <- data.frame(Assessment.SOM = rbind(BaseL.SOM.ncv.time,BaseP.SOM.ncv.time, IntL.SOM.ncv.time, IntP.SOM.ncv.time, IntHL.SOM.ncv.time,IntHP.SOM.ncv.time))
SOM.time <- data.frame(Training.SOM = rbind(BaseL.SOM.time,BaseP.SOM.time, IntL.SOM.time, IntP.SOM.time, IntHL.SOM.time,IntHP.SOM.time))

pH.ncv.time <- data.frame(Assessment.pH = rbind(BaseL.pH.ncv.time,BaseP.pH.ncv.time, IntL.pH.ncv.time, IntP.pH.ncv.time, IntHL.pH.ncv.time,IntHP.pH.ncv.time))
pH.time <- data.frame(Training.pH = rbind(BaseL.pH.time,BaseP.pH.time, IntL.pH.time, IntP.pH.time, IntHL.pH.time,IntHP.pH.time))


timeALL <- cbind(model, ORCDRC.time, ORCDRC.ncv.time, SOM.time, SOM.ncv.time, pH.time, pH.ncv.time)
stargazer(timeALL, summary = FALSE, digits = 1, type = 'text')

