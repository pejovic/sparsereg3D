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

CovNames <- c("DEM","Aspect","Slope","TWI","ConvInd","CrSectCurv","LongCurv","ChNetBLevel","VDistChNet","NegOp", "PosOp","WEeast","WEnw","DD","CD","clc","SoilType")

#=================== DATA ============================================================
bor.profs <- bor[,c("ID","x","y","Soil.Type","Top","Bottom","As","pH","Co","SOM")]
bor.profs$logSOM <- log(bor.profs$SOM)
bor.profs$logAs <- log(bor.profs$As)
bor.profs$mid <- (bor.profs$Top+bor.profs$Bottom)/2
depths(bor.profs) <- ID ~ Top + Bottom
site(bor.profs) <- ~ Soil.Type + x + y
coordinates(bor.profs) <- ~ x+y
proj4string(bor.profs) <- CRS(utm)
#=====================================================================================

summary(bor.profs@horizons[,c("As","SOM","pH","Top","Bottom")])

bor.profs@site$depth <- profileApply(bor.profs, estimateSoilDepth, name='As', top='Top', bottom='Bottom')
bor.profs@site$n <- (ddply(bor.profs@horizons,.(ID), summarize, br = length(As)))[,2]

sum.n.d <- summary(bor.profs@site[,c("depth","n")])
sum.As.SOM.pH <- summary(bor.profs@horizons[,c("As","SOM","pH")])[,1]
rbind(sum.As.SOM.pH,sum.n.d)



# Formulas
As.fun <- as.formula(paste("As ~", paste(c(CovNames,"depth"), collapse="+")))
SOM.fun <- as.formula(paste("SOM ~", paste(c(CovNames,"depth"), collapse="+")))
log.As.fun <- as.formula(paste("logAs ~", paste(c(CovNames,"depth"), collapse="+")))
log.SOM.fun <- as.formula(paste("logSOM ~", paste(c(CovNames,"depth"), collapse="+")))
pH.fun <- as.formula(paste("pH ~", paste(c(CovNames,"depth"), collapse="+")))

# SOM results
BaseL.SOM.preproc <- pre.sparsereg3D(base.model = SOM.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
BaseL.SOM.ncv.time <- system.time(BaseL.SOM.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.SOM.time <- system.time(BaseL.SOM <- sparsereg3D.sel(sparse.reg = BaseL.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = BaseL.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

BaseP.SOM.preproc <- pre.sparsereg3D(base.model = SOM.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
BaseP.SOM.ncv.time <- system.time(BaseP.SOM.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.SOM.time <- system.time(BaseP.SOM <- sparsereg3D.sel(sparse.reg = BaseP.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.p.pred <- sparsereg3D.pred(model.info = SOM.p.sel, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


IntL.SOM.preproc <- pre.sparsereg3D(base.model = SOM.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntL.SOM.ncv.time <- system.time(IntL.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntL.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.SOM.time <- system.time(IntL.SOM <- sparsereg3D.sel(sparse.reg = IntL.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = IntL.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntP.SOM.preproc <- pre.sparsereg3D(base.model = SOM.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntP.SOM.ncv.time <- system.time(IntP.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntP.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.SOM.time <- system.time(IntP.SOM <- sparsereg3D.sel(sparse.reg = IntP.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = IntP.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHL.SOM.preproc <- pre.sparsereg3D(base.model = SOM.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntHL.SOM.ncv.time <- system.time(IntHL.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.SOM.time <- system.time(IntHL.SOM <- sparsereg3D.sel(sparse.reg = IntHL.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = IntHL.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHP.SOM.preproc <- pre.sparsereg3D(base.model = SOM.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntHP.SOM.ncv.time <- system.time(IntHP.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.SOM.time <- system.time(IntHP.SOM <- sparsereg3D.sel(sparse.reg = IntHP.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = IntHP.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

#============================= Coefficients tables for SOM models ==============================================================================

ll <- length(IntL.SOM$coefficients)
pp <- length(IntHL.SOM$coefficients[,1])+1

cmL.SOM <- data.frame(variable=IntHL.SOM$coefficients[,1], BaseL.SOM.me=BaseL.SOM$coefficients[2:pp], IntL.SOM.me=IntL.SOM$coefficients[2:pp],IntL.SOM.ie=c(IntL.SOM$coefficients[(pp+1):ll],0),IntHL.SOM.me=IntHL.SOM$coefficients[,2],IntHL.SOM.ie=IntHL.SOM$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.SOM$coefficients)
p <- length(IntHP.SOM$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.SOM <- data.frame(variable=IntHP.SOM$coefficients[,1], BaseP.SOM.me=BaseP.SOM$coefficients[2:p], IntP.SOM.me=IntP.SOM$coefficients[2:p],IntP.SOM.ie1=c(IntP.SOM$coefficients[(p+1):l][i1],0,0,0),IntP.SOM.ie2=c(IntP.SOM$coefficients[(p+1):l][i2],0,0,0),IntP.SOM.ie3=c(IntP.SOM$coefficients[(p+1):l][i3],0,0,0),IntHP.SOM.me=IntHP.SOM$coefficients[,2],IntHP.SOM.ie1=IntHP.SOM$coefficients[,3],IntHP.SOM.ie2=IntHP.SOM$coefficients[,4],IntHP.SOM.ie3=IntHP.SOM$coefficients[,5] )
cmSOM <- cmP.SOM[,c(1,7:10)]

# Models comparison
SOM.ncv <- data.frame(rbind(BaseL = BaseL.SOM.ncv, BaseP = BaseP.SOM.ncv, IntL = IntL.SOM.ncv, IntP = IntP.SOM.ncv, IntHL = IntHL.SOM.ncv, IntHP = IntHP.SOM.ncv))
SOM.ncv <- data.frame(Model = rownames(SOM.ncv), SOM.ncv)
rownames(SOM.ncv) <- NULL
names(SOM.ncv) <- c("Model","RMSE","R squared")
stargazer(SOM.ncv, summary = FALSE, digits = 2, type = "text")

#SOM.results <- list(SOM.ncv = SOM.ncv, l.coeffs = cmL.SOM, p.coeffs = cmP.SOM, IntP.SOM = IntP.SOM, IntHP.SOM = IntHP.SOM, BaseP.SOM = BaseP.SOM)
#save(SOM.results, file = "SOM.results.rda" )
#==============================================================================================================================================


# As results
BaseL.As.preproc <- pre.sparsereg3D(base.model = As.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
BaseL.As.ncv.time <- system.time(BaseL.As.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.As.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.As.time <- system.time(BaseL.As <- sparsereg3D.sel(sparse.reg = BaseL.As.preproc ,lambda = seq(0,5,0.1), seed = 321))
#As.l.pred <- sparsereg3D.pred(model.info = BaseL.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

BaseP.As.preproc <- pre.sparsereg3D(base.model = As.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
BaseP.As.ncv.time <- system.time(BaseP.As.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.As.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.As.time <- system.time(BaseP.As <- sparsereg3D.sel(sparse.reg = BaseP.As.preproc ,lambda = seq(0,5,0.1), seed = 321))
#As.p.pred <- sparsereg3D.pred(model.info = As.p.sel, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


IntL.As.preproc <- pre.sparsereg3D(base.model = As.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntL.As.ncv.time <- system.time(IntL.As.ncv <- sparsereg3D.ncv(sparse.reg = IntL.As.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.As.time <- system.time(IntL.As <- sparsereg3D.sel(sparse.reg = IntL.As.preproc ,lambda = seq(0,5,0.1), seed = 321))
#As.l.pred <- sparsereg3D.pred(model.info = IntL.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntP.As.preproc <- pre.sparsereg3D(base.model = As.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntP.As.ncv.time <- system.time(IntP.As.ncv <- sparsereg3D.ncv(sparse.reg = IntP.As.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.As.time <- system.time(IntP.As <- sparsereg3D.sel(sparse.reg = IntP.As.preproc ,lambda = seq(0,5,0.1), seed = 321))
#As.l.pred <- sparsereg3D.pred(model.info = IntP.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHL.As.preproc <- pre.sparsereg3D(base.model = As.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntHL.As.ncv.time <- system.time(IntHL.As.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.As.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.As.time <- system.time(IntHL.As <- sparsereg3D.sel(sparse.reg = IntHL.As.preproc ,lambda = seq(0,5,0.1), seed = 321))
#As.l.pred <- sparsereg3D.pred(model.info = IntHL.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHP.As.preproc <- pre.sparsereg3D(base.model = As.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntHP.As.ncv.time <- system.time(IntHP.As.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.As.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.As.time <- system.time(IntHP.As <- sparsereg3D.sel(sparse.reg = IntHP.As.preproc ,lambda = seq(0,5,0.1), seed = 321))
#As.l.pred <- sparsereg3D.pred(model.info = IntHP.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


#============================= Coefficients tables for As models ==============================================================================

ll <- length(IntL.As$coefficients)
pp <- length(IntHL.As$coefficients[,1])+1

cmL.As <- data.frame(variable=IntHL.As$coefficients[,1], BaseL.As.me=BaseL.As$coefficients[2:pp], IntL.As.me=IntL.As$coefficients[2:pp],IntL.As.ie=c(IntL.As$coefficients[(pp+1):ll],0),IntHL.As.me=IntHL.As$coefficients[,2],IntHL.As.ie=IntHL.As$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.As$coefficients)
p <- length(IntHP.As$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.As <- data.frame(variable=IntHP.As$coefficients[,1], BaseP.As.me = BaseP.As$coefficients[2:p], IntP.As.me=IntP.As$coefficients[2:p],IntP.As.ie1=c(IntP.As$coefficients[(p+1):l][i1],0,0,0),IntP.As.ie2=c(IntP.As$coefficients[(p+1):l][i2],0,0,0),IntP.As.ie3=c(IntP.As$coefficients[(p+1):l][i3],0,0,0),IntHP.As.me=IntHP.As$coefficients[,2],IntHP.As.ie1=IntHP.As$coefficients[,3],IntHP.As.ie2=IntHP.As$coefficients[,4],IntHP.As.ie3=IntHP.As$coefficients[,5] )
cmSOM <- cmP.As[,c(1,7:10)]

# Models comparison
As.ncv <- data.frame(rbind(BaseL = BaseL.As.ncv, BaseP = BaseP.As.ncv, IntL = IntL.As.ncv, IntP = IntP.As.ncv, IntHL = IntHL.As.ncv, IntHP = IntHP.As.ncv))
As.ncv <- data.frame(Model = rownames(As.ncv), As.ncv)
rownames(As.ncv) <- NULL
names(As.ncv) <- c("Model","RMSE","R squared")
stargazer(As.ncv, summary = FALSE, digits = 2, type = "text")

As.results <- list(As.ncv = As.ncv, l.coeffs = cmL.As, p.coeffs = cmP.As, IntP.As = IntP.As, IntHP.As = IntHP.As, IntL.As = IntL.As)
save(As.results, file = "As.results.rda" )
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
IntL.pH.time <- system.time(IntHL.pH <- sparsereg3D.sel(sparse.reg = IntHL.pH.preproc ,lambda = seq(0,5,0.1), seed = 321))
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
cmSOM <- cmP.pH[,c(1,7:10)]

#==============================================================================================================================================

# Models comparison
pH.ncv <- data.frame(rbind(BaseL = BaseL.pH.ncv, BaseP = BaseP.pH.ncv, IntL = IntL.pH.ncv, IntP = IntP.pH.ncv, IntHL = IntHL.pH.ncv, IntHP = IntHP.pH.ncv))
pH.ncv <- data.frame(Model = rownames(pH.ncv), pH.ncv)
rownames(pH.ncv) <- NULL
names(pH.ncv) <- c("Model","RMSE","R squared")
stargazer(pH.ncv, summary = FALSE, digits = 2, type = "text")

pH.results <- list(pH.ncv = pH.ncv, l.coeffs = cmL.pH, p.coeffs = cmP.pH, IntP.pH = IntP.pH, IntHP.pH = IntHP.pH, BaseP.pH = BaseP.pH)
save(pH.results, file = "pH.results.rda" )



load("As.results.rda")


#============================= Time Table =======================================================================================================

As.ncv.time <- rbind(BaseL.As.ncv.time, BaseP.As.ncv.time, IntL.As.ncv.time, IntP.As.ncv.time, IntHL.As.ncv.time, IntHP.As.ncv.time)
As.time <- rbind(BaseL.As.time, BaseP.As.time, IntL.As.time, IntP.As.time, IntHL.As.time, IntHP.As.time)

SOM.ncv.time <- rbind(BaseL.SOM.ncv.time, BaseP.SOM.ncv.time, IntL.SOM.ncv.time, IntP.SOM.ncv.time, IntHL.SOM.ncv.time, IntHP.SOM.ncv.time)
SOM.time <- rbind(BaseL.SOM.time, BaseP.SOM.time, IntL.SOM.time, IntP.SOM.time, IntHL.SOM.time, IntHP.SOM.time)

pH.ncv.time <- rbind(BaseL.pH.ncv.time, BaseP.pH.ncv.time, IntL.pH.ncv.time, IntP.pH.ncv.time, IntHL.pH.ncv.time, IntHP.pH.ncv.time)
pH.time <- rbind(BaseL.pH.time, BaseP.pH.time, IntL.pH.time, IntP.pH.time, IntHL.pH.time, IntHP.pH.time)

#================================================================================================================================================
  
IntP.coeffs <- list(cmAs = As.results$p.coeffs[,c(1,3:6)], cmpH = pH.results$p.coeffs[,c(1,3:6)], cmSOM = SOM.results$p.coeffs[,c(1,3:6)])

IntP.coeffs <- cbind(IntP.coeffs$cmSOM,IntP.coeffs$cmpH[,-1])
IntP.coeffs[,1] <- as.character(IntP.coeffs[,1])
IntP.coeffs <- cbind(As.results$p.coeffs[,c(1,3:6)],IntP.coeffs[,-1])
names(IntP.coeffs) <- c("variable","As.me","As.ie1","As.ie2","As.ie3","SOM.me","SOM.ie1","SOM.ie2","SOM.ie3","pH.me","pH.ie1","pH.ie2","pH.ie3")

stargazer(IntP.coeffs,summary=FALSE, digits=2, type="latex")


IntHP.coeffs <- list(cmAs = As.results$p.coeffs[,c(1,7:10)], cmpH = pH.results$p.coeffs[,c(1,7:10)], cmSOM = SOM.results$p.coeffs[,c(1,7:10)])

IntHP.coeffs <- cbind(IntHP.coeffs$cmSOM,IntHP.coeffs$cmpH[,-1])
IntHP.coeffs[,1] <- as.character(IntHP.coeffs[,1])
IntHP.coeffs <- cbind(As.results$p.coeffs[,c(1,7:10)],IntHP.coeffs[,-1])
names(IntHP.coeffs) <- c("variable","As.me","As.ie1","As.ie2","As.ie3","SOM.me","SOM.ie1","SOM.ie2","SOM.ie3","pH.me","pH.ie1","pH.ie2","pH.ie3")

stargazer(IntHP.coeffs,summary=FALSE, digits=2, type="latex")


IntL.coeffs <- cbind(As.results$l.coeffs[,-c(2)],SOM.results$l.coeffs[,-c(1,2)],pH.results$l.coeffs[,-c(1,2)])
stargazer(IntL.coeffs,summary=FALSE, digits=2, type="latex")

BaseL <- list(cmAs.l = As.results$l.coeffs[,c(1,2)], cmpH.l = pH.results$l.coeffs[,c(1,2)], cmSOM.l = SOM.results$l.coeffs[,c(1,2)])
BaseL <- cbind(BaseL$cmAs.l, BaseL$cmSOM.l[,-1],BaseL$cmpH.l[,-1])
BaseL[,1] <- as.character(BaseL[,1])
names(BaseL) <- c("variable","As","SOM","pH")

stargazer(BaseL,summary=FALSE, digits=2, type="latex")

BaseP <- list(cmAs.p = As.results$p.coeffs[,c(1,2)], cmpH.p = pH.results$p.coeffs[,c(1,2)], cmSOM.p = SOM.results$p.coeffs[,c(1,2)])
BaseP <- cbind(BaseP$cmAs.p, BaseP$cmSOM.p[,-1],BaseP$cmpH.p[,-1])
BaseP[,1] <- as.character(BaseP[,1])
names(BaseP) <- c("variable","As","SOM","pH")

stargazer(BaseP, summary=FALSE, digits=2, type="latex")


load("As.results.rda")
load("SOM.results.rda")
load("pH.results.rda")
#============================ As coef path plot =================================================

altitude <- data.frame(depth=seq(0.0,0.40,0.01))
depth <- data.frame(depth = bor.profs$mid/100)
pp <-preProcess(depth, method=c("center", "scale"))

#altitude <- data.frame(altitude=c(-0.1,-0.3,-0.4))
altitude.s <- as.numeric(predict(pp , newdata = altitude)[,1])
variables <- c("CD","DD","DEM","Slope","VDistChNet","WEnw")#,"TWI")
variables <- variables[order(variables)]

#10.4895187  5.7506602 -0.4258973000
#cmP.As[cmP.As$variable=="VDistChNet","IntHP.ie1"] <- 5.7506602

cmP.As <- As.results$l.coeffs

coefs.IntL <- (cmP.As[which(as.character(cmP.As$variable) %in% variables),c(1,3:4)])

coefs.IntL <- coefs.IntL[order(coefs.IntL$variable),]

coefs.IntL <- coefs.IntL[,-1]


effects.IntL <- as.numeric()
for(i in 1:dim(coefs.IntL)[1]){
  effects <- coefs.IntL[i,1] + coefs.IntL[i,2]*altitude.s
  effects <- data.frame(depth = altitude[,1], var = effects)
  effects.IntL <- rbind(effects.IntP, effects)
}  

effects.IntL$Variables <- factor(rep(variables, each = length(altitude.s)))
intPlot <- qplot(var, depth, data = effects.IntL, geom = c("point", "line"), color = Variables)+theme_bw()+theme(legend.text=element_text(size=12))+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+labs(x="Coefficients",y="Depth")

pdf("AsIntHPlot.pdf",width=8,height=10)
intHPlot
dev.off()

pdf("AsIntPlot.pdf",width=8,height=10)
intPlot
dev.off()


#============================ SOM coef path plot =================================================

altitude <- data.frame(depth=seq(0.0,0.40,0.01))
depth <- data.frame(depth = bor.profs$mid/100)
pp <-preProcess(depth, method=c("center", "scale"))

#altitude <- data.frame(altitude=c(-0.1,-0.3,-0.4))
altitude.s <- as.numeric(predict(pp , newdata = altitude)[,1])
variables <- c("CD","DD","DEM","Slope","VDistChNet","WEnw")#,"TWI")
variables <- variables[order(variables)]

#10.4895187  5.7506602 -0.4258973000
#cmP.As[cmP.As$variable=="VDistChNet","IntHP.ie1"] <- 5.7506602

cmP.As <- SOM.results$L.coeffs

coefs.IntHP <- (cmP.As[which(as.character(cmP.As$variable) %in% variables),c(1,7:10)])
coefs.IntP <- (cmP.As[which(as.character(cmP.As$variable) %in% variables),c(1,3:6)])

coefs.IntHP <- coefs.IntHP[order(coefs.IntHP$variable),]
coefs.IntP <- coefs.IntP[order(coefs.IntP$variable),]

coefs.IntHP <- coefs.IntHP[,-1]
coefs.IntP <- coefs.IntP[,-1]


effects.IntHP <- as.numeric()
for(i in 1:dim(coefs.IntHP)[1]){
  effects <- coefs.IntHP[i,1] + coefs.IntHP[i,2]*altitude.s + coefs.IntHP[i,3]*altitude.s^2 + coefs.IntHP[i,4]*altitude.s^3
  effects <- data.frame(depth = altitude[,1], var = effects)
  effects.IntHP <- rbind(effects.IntHP, effects)
}  

effects.IntHP$Variables <- factor(rep(variables, each = length(altitude.s)))
intHPlot <- qplot(var, depth, data = effects.IntHP, geom = c("point", "line"), color = Variables)+theme_bw()+theme(legend.text=element_text(size=12))+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+labs(x="Coefficients",y="Depth")


effects.IntP <- as.numeric()
for(i in 1:dim(coefs.IntP)[1]){
  effects <- coefs.IntP[i,1] + coefs.IntP[i,2]*altitude.s + coefs.IntP[i,3]*altitude.s^2 + coefs.IntP[i,4]*altitude.s^3
  effects <- data.frame(depth = altitude[,1], var = effects)
  effects.IntP <- rbind(effects.IntP, effects)
}  

effects.IntP$Variables <- factor(rep(variables, each = length(altitude.s)))
intPlot <- qplot(var, depth, data = effects.IntP, geom = c("point", "line"), color = Variables)+theme_bw()+theme(legend.text=element_text(size=12))+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+labs(x="Coefficients",y="Depth")

pdf("AsIntHPlot.pdf",width=8,height=10)
intHPlot
dev.off()

pdf("AsIntPlot.pdf",width=8,height=10)
intPlot
dev.off()

#============= Number of selected covariates =================================================

#============= IntP models ===================================================================


p.coeffs <- As.results$p.coeffs

np.coeffs <- function(p.coeffs){
  BaseP <- p.coeffs[,2]
  IntP <- p.coeffs[,3:6]
  IntHP <- p.coeffs[,7:10]
  
  n.BaseP <-  c(sum((BaseP!=0 )),0,0,0, sum((BaseP!=0 )))
  n.IntP <-  c(apply(IntP,2, function(y) sum(y!=0)),total = sum(apply(IntP,2, function(y) sum(y!=0))))
  n.IntHP <- c(apply(IntHP,2, function(y) sum(y!=0)),total = sum(apply(IntHP,2, function(y) sum(y!=0))))
  
  total <- data.frame(rbind(n.BaseP,n.IntP,n.IntHP))
  total <- cbind(model = c("BaseP", "IntP", "IntHP"), total)
  rownames(total) <- NULL
  names(total) <- c("model","me","ie(d)","ie(d2)","ie(d3)","total")
  return(total)
}

As.np.coeffs <- np.coeffs(As.results$p.coeffs)
SOM.np.coeffs <- np.coeffs(SOM.results$p.coeffs)
pH.np.coeffs <- np.coeffs(pH.results$p.coeffs)

nl.coeffs <- function(l.coeffs){
  BaseL <- l.coeffs[,2]
  IntL <- l.coeffs[,3:4]
  IntHL <- l.coeffs[,5:6]
  
  n.BaseL <-  c(sum((BaseL!=0 )),0, sum((BaseL!=0 )))
  n.IntL <-  c(apply(IntL,2, function(y) sum(y!=0)),total = sum(apply(IntL,2, function(y) sum(y!=0))))
  n.IntHL <- c(apply(IntHL,2, function(y) sum(y!=0)),total = sum(apply(IntHL,2, function(y) sum(y!=0))))
  
  total <- data.frame(rbind(n.BaseL,n.IntL,n.IntHL))
  total <- cbind(model = c("BaseL", "IntL", "IntHL"), total)
  rownames(total) <- NULL
  names(total) <- c("model","me","ie(d)","total")
  return(total)
}


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


As.n.coeffs <- n.coeffs(l.coeffs = As.results$l.coeffs, p.coeffs = As.results$p.coeffs)
SOM.n.coeffs <- n.coeffs(l.coeffs = SOM.results$l.coeffs, p.coeffs = SOM.results$p.coeffs)
pH.n.coeffs <-n.coeffs(l.coeffs = pH.results$l.coeffs, p.coeffs = pH.results$p.coeffs)

stargazer(cbind(As.n.coeffs,SOM.n.coeffs[,-1],pH.n.coeffs[,-1]), summary = FALSE, type = 'latex')



As.nl.coeffs <- nl.coeffs(As.results$l.coeffs)
SOM.nl.coeffs <- nl.coeffs(SOM.results$l.coeffs)
pH.nl.coeffs <- nl.coeffs(pH.results$l.coeffs)

l.coeffs <- As.results$l.coeffs
p.coeffs <- As.results$p.coeffs



######################### LOG Results ###########################################################
source(paste(getwd(),"sparsereg3D.sel.r",sep="/"))
source(paste(getwd(),"pre.sparsereg3D.r",sep="/"))
# log.SOM results
BaseL.log.SOM.preproc <- pre.sparsereg3D(base.model = log.SOM.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
BaseL.log.SOM.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.log.SOM.preproc, lambda = seq(0,5,0.1), seed = 321)
BaseL.log.SOM <- sparsereg3D.sel(sparse.reg = BaseL.log.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321)
#log.SOM.l.pred <- sparsereg3D.pred(model.info = BaseL.log.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

BaseP.log.SOM.preproc <- pre.sparsereg3D(base.model = log.SOM.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
BaseP.log.SOM.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.log.SOM.preproc, lambda = seq(0,5,0.1), seed = 321)
BaseP.log.SOM <- sparsereg3D.sel(sparse.reg = BaseP.log.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321)
#log.SOM.p.pred <- sparsereg3D.pred(model.info = log.SOM.p.sel, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


IntL.log.SOM.preproc <- pre.sparsereg3D(base.model = log.SOM.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntL.log.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntL.log.SOM.preproc, lambda = seq(0,5,0.1), seed = 321)
IntL.log.SOM <- sparsereg3D.sel(sparse.reg = IntL.log.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321)
#log.SOM.l.pred <- sparsereg3D.pred(model.info = IntL.log.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntP.log.SOM.preproc <- pre.sparsereg3D(base.model = log.SOM.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntP.log.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntP.log.SOM.preproc, lambda = seq(0,5,0.1), seed = 321)
IntP.log.SOM <- sparsereg3D.sel(sparse.reg = IntP.log.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321)
#log.SOM.l.pred <- sparsereg3D.pred(model.info = IntP.log.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHL.log.SOM.preproc <- pre.sparsereg3D(base.model = log.SOM.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntHL.log.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.log.SOM.preproc, lambda = seq(0,5,0.1), seed = 321)
IntHL.log.SOM <- sparsereg3D.sel(sparse.reg = IntHL.log.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321)
#log.SOM.l.pred <- sparsereg3D.pred(model.info = IntHL.log.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHP.log.SOM.preproc <- pre.sparsereg3D(base.model = log.SOM.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntHP.log.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.log.SOM.preproc, lambda = seq(0,5,0.1), seed = 321)
IntHP.log.SOM <- sparsereg3D.sel(sparse.reg = IntHP.log.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321)
#log.SOM.l.pred <- sparsereg3D.pred(model.info = IntHP.log.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

#============================= Coefficients tables for log.SOM models ==============================================================================

ll <- length(IntL.log.SOM$coefficients)
pp <- length(IntHL.log.SOM$coefficients[,1])+1

cmL.log.SOM <- data.frame(variable=IntHL.log.SOM$coefficients[,1], BaseL.log.SOM.me=BaseL.log.SOM$coefficients[2:pp], IntL.log.SOM.me=IntL.log.SOM$coefficients[2:pp],IntL.log.SOM.ie=c(IntL.log.SOM$coefficients[(pp+1):ll],0),IntHL.log.SOM.me=IntHL.log.SOM$coefficients[,2],IntHL.log.SOM.ie=IntHL.log.SOM$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.log.SOM$coefficients)
p <- length(IntHP.log.SOM$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.log.SOM <- data.frame(variable=IntHP.log.SOM$coefficients[,1], BaseP.log.SOM.me=BaseP.log.SOM$coefficients[2:p], IntP.log.SOM.me=IntP.log.SOM$coefficients[2:p],IntP.log.SOM.ie1=c(IntP.log.SOM$coefficients[(p+1):l][i1],0,0,0),IntP.log.SOM.ie2=c(IntP.log.SOM$coefficients[(p+1):l][i2],0,0,0),IntP.log.SOM.ie3=c(IntP.log.SOM$coefficients[(p+1):l][i3],0,0,0),IntHP.log.SOM.me=IntHP.log.SOM$coefficients[,2],IntHP.log.SOM.ie1=IntHP.log.SOM$coefficients[,3],IntHP.log.SOM.ie2=IntHP.log.SOM$coefficients[,4],IntHP.log.SOM.ie3=IntHP.log.SOM$coefficients[,5] )
cmSOM <- cmP.log.SOM[,c(1,7:10)]

# Models comparison
log.SOM.ncv <- data.frame(rbind(BaseL = BaseL.log.SOM.ncv, BaseP = BaseP.log.SOM.ncv, IntL = IntL.log.SOM.ncv, IntP = IntP.log.SOM.ncv, IntHL = IntHL.log.SOM.ncv, IntHP = IntHP.log.SOM.ncv))
log.SOM.ncv <- data.frame(Model = rownames(log.SOM.ncv), log.SOM.ncv)
rownames(log.SOM.ncv) <- NULL
names(log.SOM.ncv) <- c("Model","RMSE","R squared")
stargazer(log.SOM.ncv, summary = FALSE, digits = 2, type = "text")

#SOM.results <- list(SOM.ncv = SOM.ncv, l.coeffs = cmL.SOM, p.coeffs = cmP.SOM, IntP.SOM = IntP.SOM, IntHP.SOM = IntHP.SOM, BaseP.SOM = BaseP.SOM)
#save(SOM.results, file = "SOM.results.rda" )
#==============================================================================================================================================


# log.As results
BaseL.log.As.preproc <- pre.sparsereg3D(base.model = log.As.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
BaseL.log.As.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.log.As.preproc, lambda = seq(0,5,0.1), seed = 321)
BaseL.log.As <- sparsereg3D.sel(sparse.reg = BaseL.log.As.preproc ,lambda = seq(0,5,0.1), seed = 321)
#log.As.l.pred <- sparsereg3D.pred(model.info = BaseL.log.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

BaseP.log.As.preproc <- pre.sparsereg3D(base.model = log.As.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
BaseP.log.As.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.log.As.preproc, lambda = seq(0,5,0.1), seed = 321)
BaseP.log.As <- sparsereg3D.sel(sparse.reg = BaseP.log.As.preproc ,lambda = seq(0,5,0.1), seed = 321)
#log.As.p.pred <- sparsereg3D.pred(model.info = log.As.p.sel, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


IntL.log.As.preproc <- pre.sparsereg3D(base.model = log.As.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntL.log.As.ncv <- sparsereg3D.ncv(sparse.reg = IntL.log.As.preproc, lambda = seq(0,5,0.1), seed = 321)
IntL.log.As <- sparsereg3D.sel(sparse.reg = IntL.log.As.preproc ,lambda = seq(0,5,0.1), seed = 321)
#log.As.l.pred <- sparsereg3D.pred(model.info = IntL.log.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntP.log.As.preproc <- pre.sparsereg3D(base.model = log.As.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntP.log.As.ncv <- sparsereg3D.ncv(sparse.reg = IntP.log.As.preproc, lambda = seq(0,5,0.1), seed = 321)
IntP.log.As <- sparsereg3D.sel(sparse.reg = IntP.log.As.preproc ,lambda = seq(0,5,0.1), seed = 321)
#log.As.l.pred <- sparsereg3D.pred(model.info = IntP.log.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

As.profs <- bor.profs
As.profs@horizons$As.pred <- IntP.log.As$training.prediction

IntHL.log.As.preproc <- pre.sparsereg3D(base.model = log.As.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntHL.log.As.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.log.As.preproc, lambda = seq(0,5,0.1), seed = 321)
IntHL.log.As <- sparsereg3D.sel(sparse.reg = IntHL.log.As.preproc ,lambda = seq(0,5,0.1), seed = 321)
#log.As.l.pred <- sparsereg3D.pred(model.info = IntHL.log.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHP.log.As.preproc <- pre.sparsereg3D(base.model = log.As.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntHP.log.As.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.log.As.preproc, lambda = seq(0,5,0.1), seed = 321)
IntHP.log.As <- sparsereg3D.sel(sparse.reg = IntHP.log.As.preproc ,lambda = seq(0,5,0.1), seed = 321)
#log.As.l.pred <- sparsereg3D.pred(model.info = IntHP.log.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


#============================= Coefficients tables for log.As models ==============================================================================

ll <- length(IntL.log.As$coefficients)
pp <- length(IntHL.log.As$coefficients[,1])+1

cmL.log.As <- data.frame(variable=IntHL.log.As$coefficients[,1], BaseL.log.As.me=BaseL.log.As$coefficients[2:pp], IntL.log.As.me=IntL.log.As$coefficients[2:pp],IntL.log.As.ie=c(IntL.log.As$coefficients[(pp+1):ll],0),IntHL.log.As.me=IntHL.log.As$coefficients[,2],IntHL.log.As.ie=IntHL.log.As$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.log.As$coefficients)
p <- length(IntHP.log.As$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.log.As <- data.frame(variable=IntHP.log.As$coefficients[,1], BaseP.log.As.me = BaseP.log.As$coefficients[2:p], IntP.log.As.me=IntP.log.As$coefficients[2:p],IntP.log.As.ie1=c(IntP.log.As$coefficients[(p+1):l][i1],0,0,0),IntP.log.As.ie2=c(IntP.log.As$coefficients[(p+1):l][i2],0,0,0),IntP.log.As.ie3=c(IntP.log.As$coefficients[(p+1):l][i3],0,0,0),IntHP.log.As.me=IntHP.log.As$coefficients[,2],IntHP.log.As.ie1=IntHP.log.As$coefficients[,3],IntHP.log.As.ie2=IntHP.log.As$coefficients[,4],IntHP.log.As.ie3=IntHP.log.As$coefficients[,5] )
cmSOM <- cmP.log.As[,c(1,7:10)]

# Models comparison
log.As.ncv <- data.frame(rbind(BaseL = BaseL.log.As.ncv, BaseP = BaseP.log.As.ncv, IntL = IntL.log.As.ncv, IntP = IntP.log.As.ncv, IntHL = IntHL.log.As.ncv, IntHP = IntHP.log.As.ncv))
log.As.ncv <- data.frame(Model = rownames(log.As.ncv), log.As.ncv)
rownames(log.As.ncv) <- NULL
names(log.As.ncv) <- c("Model","RMSE","R squared")
stargazer(log.As.ncv, summary = FALSE, digits = 2, type = "text")


#=================== TIME ==========================================================================================================

# SOM results
BaseL.SOM.preproc <- pre.sparsereg3D(base.model = SOM.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
BaseL.SOM.ncv.time <- system.time(BaseL.SOM.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.SOM.time <- system.time(BaseL.SOM <- sparsereg3D.sel(sparse.reg = BaseL.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = BaseL.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

BaseP.SOM.preproc <- pre.sparsereg3D(base.model = SOM.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
BaseP.SOM.ncv.time <- system.time(BaseP.SOM.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.SOM.time <- system.time(BaseP.SOM <- sparsereg3D.sel(sparse.reg = BaseP.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.p.pred <- sparsereg3D.pred(model.info = SOM.p.sel, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


IntL.SOM.preproc <- pre.sparsereg3D(base.model = SOM.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntL.SOM.ncv.time <- system.time(IntL.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntL.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.SOM.time <- system.time(IntL.SOM <- sparsereg3D.sel(sparse.reg = IntL.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = IntL.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntP.SOM.preproc <- pre.sparsereg3D(base.model = SOM.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntP.SOM.ncv.time <- system.time(IntP.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntP.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.SOM.time <- system.time(IntP.SOM <- sparsereg3D.sel(sparse.reg = IntP.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = IntP.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHL.SOM.preproc <- pre.sparsereg3D(base.model = SOM.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntHL.SOM.ncv.time <- system.time(IntHL.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHL.SOM.time <- system.time(IntHL.SOM <- sparsereg3D.sel(sparse.reg = IntHL.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = IntHL.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHP.SOM.preproc <- pre.sparsereg3D(base.model = SOM.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntHP.SOM.ncv.time <- system.time(IntHP.SOM.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.SOM.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.SOM.time <- system.time(IntHP.SOM <- sparsereg3D.sel(sparse.reg = IntHP.SOM.preproc ,lambda = seq(0,5,0.1), seed = 321))
#SOM.l.pred <- sparsereg3D.pred(model.info = IntHP.SOM, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

# As results
BaseL.As.preproc <- pre.sparsereg3D(base.model = As.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
BaseL.As.ncv.time <- system.time(BaseL.As.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.As.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseL.As.time <- system.time(BaseL.As <- sparsereg3D.sel(sparse.reg = BaseL.As.preproc ,lambda = seq(0,5,0.1), seed = 321))
#As.l.pred <- sparsereg3D.pred(model.info = BaseL.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

BaseP.As.preproc <- pre.sparsereg3D(base.model = As.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
BaseP.As.ncv.time <- system.time(BaseP.As.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.As.preproc, lambda = seq(0,5,0.1), seed = 321))
BaseP.As.time <- system.time(BaseP.As <- sparsereg3D.sel(sparse.reg = BaseP.As.preproc ,lambda = seq(0,5,0.1), seed = 321))
#As.p.pred <- sparsereg3D.pred(model.info = As.p.sel, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


IntL.As.preproc <- pre.sparsereg3D(base.model = As.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntL.As.ncv.time <- system.time(IntL.As.ncv <- sparsereg3D.ncv(sparse.reg = IntL.As.preproc, lambda = seq(0,5,0.1), seed = 321))
IntL.As.time <- system.time(IntL.As <- sparsereg3D.sel(sparse.reg = IntL.As.preproc ,lambda = seq(0,5,0.1), seed = 321))
#As.l.pred <- sparsereg3D.pred(model.info = IntL.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntP.As.preproc <- pre.sparsereg3D(base.model = As.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntP.As.ncv.time <- system.time(IntP.As.ncv <- sparsereg3D.ncv(sparse.reg = IntP.As.preproc, lambda = seq(0,5,0.1), seed = 321))
IntP.As.time <- system.time(IntP.As <- sparsereg3D.sel(sparse.reg = IntP.As.preproc ,lambda = seq(0,5,0.1), seed = 321))
#As.l.pred <- sparsereg3D.pred(model.info = IntP.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHL.As.preproc <- pre.sparsereg3D(base.model = As.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntHL.As.ncv.time <- system.time(IntHL.As.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.As.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHL.As.time <- system.time(IntHL.As <- sparsereg3D.sel(sparse.reg = IntHL.As.preproc ,lambda = seq(0,5,0.1), seed = 321))
#As.l.pred <- sparsereg3D.pred(model.info = IntHL.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntHP.As.preproc <- pre.sparsereg3D(base.model = As.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntHP.As.ncv.time <- system.time(IntHP.As.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.As.preproc, lambda = seq(0,5,0.1), seed = 321))
IntHP.As.time <- system.time(IntHP.As <- sparsereg3D.sel(sparse.reg = IntHP.As.preproc ,lambda = seq(0,5,0.1), seed = 321))
#As.l.pred <- sparsereg3D.pred(model.info = IntHP.As, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

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

model = data.frame(Model = c("BaseL","BaseP","IntL","IntP","IntHL","IntHP"))

As.ncv.time <- data.frame(Assessment.As = c(BaseL.As.ncv.time[3],BaseP.As.ncv.time[3], IntL.As.ncv.time[3], IntP.As.ncv.time[3], IntHL.As.ncv.time[3],IntHP.As.ncv.time[3]))
As.time <- data.frame(Training.As = c(BaseL.As.time[3],BaseP.As.time[3], IntL.As.time[3], IntP.As.time[3], IntHL.As.time[3],IntHP.As.time[3]))

SOM.ncv.time <- data.frame(Assessment.SOM = c(BaseL.SOM.ncv.time[3],BaseP.SOM.ncv.time[3], IntL.SOM.ncv.time[3], IntP.SOM.ncv.time[3], IntHL.SOM.ncv.time[3],IntHP.SOM.ncv.time[3]))
SOM.time <- data.frame(Training.SOM = c(BaseL.SOM.time[3],BaseP.SOM.time[3], IntL.SOM.time[3], IntP.SOM.time[3], IntHL.SOM.time[3],IntHP.SOM.time[3]))

pH.ncv.time <- data.frame(Assessment.pH = c(BaseL.pH.ncv.time[3],BaseP.pH.ncv.time[3], IntL.pH.ncv.time[3], IntP.pH.ncv.time[3], IntHL.pH.ncv.time[3],IntHP.pH.ncv.time[3]))
pH.time <- data.frame(Training.pH = c(BaseL.pH.time[3],BaseP.pH.time[3], IntL.pH.time[3], IntP.pH.time[3], IntHL.pH.time[3],IntHP.pH.time[3]))

time <- cbind(model,As.time, As.ncv.time, SOM.time, SOM.ncv.time, pH.time, pH.ncv.time)


stargazer(time, summary = FALSE, digits = 1, type = 'latex')


model = data.frame(Model = c("BaseL","BaseP","IntL","IntP","IntHL","IntHP"))

As.ncv.time <- data.frame(Assessment.As = rbind(BaseL.As.ncv.time,BaseP.As.ncv.time, IntL.As.ncv.time, IntP.As.ncv.time, IntHL.As.ncv.time,IntHP.As.ncv.time))
As.time <- data.frame(Training.As = rbind(BaseL.As.time,BaseP.As.time, IntL.As.time, IntP.As.time, IntHL.As.time,IntHP.As.time))

SOM.ncv.time <- data.frame(Assessment.SOM = rbind(BaseL.SOM.ncv.time,BaseP.SOM.ncv.time, IntL.SOM.ncv.time, IntP.SOM.ncv.time, IntHL.SOM.ncv.time,IntHP.SOM.ncv.time))
SOM.time <- data.frame(Training.SOM = rbind(BaseL.SOM.time,BaseP.SOM.time, IntL.SOM.time, IntP.SOM.time, IntHL.SOM.time,IntHP.SOM.time))

pH.ncv.time <- data.frame(Assessment.pH = rbind(BaseL.pH.ncv.time,BaseP.pH.ncv.time, IntL.pH.ncv.time, IntP.pH.ncv.time, IntHL.pH.ncv.time,IntHP.pH.ncv.time))
pH.time <- data.frame(Training.pH = rbind(BaseL.pH.time,BaseP.pH.time, IntL.pH.time, IntP.pH.time, IntHL.pH.time,IntHP.pH.time))


timeALL <- cbind(model,As.time, As.ncv.time, SOM.time, SOM.ncv.time, pH.time, pH.ncv.time)
stargazer(timeALL, summary = FALSE, digits = 1, type = 'text')

