# Napisati nesto ovde
#
#
path1 <- "C:/Users/User/Dropbox/Extensions of soil 3D trend models/Data and Scripts"



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

fun.path <- "D:/R_projects/sparsereg3D/R"
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
  total <- cbind(model = c("BaseL","BaseP", "IntL","IntP", "IntHL","IntHP"), total)
  return(total)
}


#================== Spatial references ===================================================================
gk_7 <- "+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 +y_0=0 +ellps=bessel +towgs84=574.027,170.175,401.545,4.88786,-0.66524,-13.24673,0.99999311067 +units=m"
utm <- "+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m"

#================== Names and abbrevations of covariates =================================================

CovNames <- c("DEM","Aspect","Slope","TWI","ConvInd","CrSectCurv","LongCurv","ChNetBLevel","VDistChNet","NegOp", "PosOp","WEeast","WEnw","clc","SoilType")

#=================== DATA ================================================================================
bor.profs <- bor[,c("ID","x","y","Soil.Type","Top","Bottom","As","pH","Co","SOM")]
bor.profs$logSOM <- log(bor.profs$SOM)
bor.profs$logAs <- log(bor.profs$As)
bor.profs$ORCDRC <- bor.profs$SOM*10/1.724
bor.profs$logORCDRC <- log1p(bor.profs$ORCDRC)

bor.profs$mid <- (bor.profs$Top+bor.profs$Bottom)/2
depths(bor.profs) <- ID ~ Top + Bottom
site(bor.profs) <- ~ Soil.Type + x + y
coordinates(bor.profs) <- ~ x+y
proj4string(bor.profs) <- CRS(utm)
#=====================================================================================

#Aggregation profiles
bor.profs@horizons <- rename(bor.profs@horizons, c("ORCDRC"="SOC"))
agg <- slab(bor.profs, fm= ~ SOC + pH, slab.structure=seq(0,70,5))

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

pdf("BorAgg.pdf",width=6,height=8)
plot(Figure2) # Make plot
dev.off()
#==========================================================
edg <- edgeroi.spc@horizons[,c("SOC","pH")]
nl <- nl.profiles@horizons[,c("SOC","pH")]
bors <- bor.profs@horizons[,c("SOC","pH")]

mean.and.sd <- function(values) {
  nas <- sum(is.na(values))
  narm.values <- (values[!is.na(values)])
  m <- mean(values, na.rm=TRUE)
  s <- sd(values, na.rm=TRUE)
  s.mean <- s/sqrt(length(narm.values))
  mini <- min(values, na.rm=TRUE)
  maxi <- max(values, na.rm=TRUE)
  q25 <- quantile(values, probs = 0.25, na.rm=TRUE)
  q75 <- quantile(values, probs = 0.75, na.rm=TRUE)
  cvar <- s/m
  med <- median(values,na.rm=TRUE)
  inqr <- IQR(values,na.rm=TRUE)
  n <- length(narm.values)
  res <- c(min=mini, q1=q25, mean=m , median=med, q3=q75 ,max=maxi, sd=s, obs=n, NAs = nas)
  return(res)
}

edg.res <- apply(edg, 2, mean.and.sd)
nl.res <- apply(nl, 2, mean.and.sd)
bor.res <- apply(bors, 2, mean.and.sd)

des.res <- cbind(bor.res, edg.res, nl.res)

stargazer(des.res, summary = FALSE, digits = 2, type = "latex")


#==========================================================

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

#==============================================================================================================================================
seed = 1

log.ORCDRC.fun <- as.formula(paste("logORCDRC ~", paste(c(CovNames,"depth"), collapse="+")))


# logORCDRC results
BaseL.logORCDRC.preproc <- pre.sparsereg3D(base.model = log.ORCDRC.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, seed = seed, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    #, kmean.vars = all.vars(log.ORCDRC.fun), cum.prop = 0.90
BaseL.logORCDRC.ncv.time <- system.time(BaseL.logORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.logORCDRC.preproc, lambda = seq(0,0.2,0.001), w = NULL)) #seq(0.1, 1, 0.1)
BaseL.logORCDRC.time <- system.time(BaseL.logORCDRC <- sparsereg3D.sel(sparse.reg = BaseL.logORCDRC.preproc ,lambda = seq(0,0.2,0.001)))
#logORCDRC.l.pred <- sparsereg3D.pred(model.info = BaseL.logORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

BaseP.logORCDRC.preproc <- pre.sparsereg3D(base.model = log.ORCDRC.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, seed = seed, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)  
BaseP.logORCDRC.ncv.time <- system.time(BaseP.logORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.logORCDRC.preproc, lambda = seq(0,0.2,0.001), w = NULL))
BaseP.logORCDRC.time <- system.time(BaseP.logORCDRC <- sparsereg3D.sel(sparse.reg = BaseP.logORCDRC.preproc ,lambda = seq(0,0.2,0.001)))
#logORCDRC.p.pred <- sparsereg3D.pred(model.info = logORCDRC.p.sel, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


IntL.logORCDRC.preproc <- pre.sparsereg3D(base.model = log.ORCDRC.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, seed = seed, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntL.logORCDRC.ncv.time <- system.time(IntL.logORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = IntL.logORCDRC.preproc, lambda = seq(0,0.2,0.001), w = NULL))
IntL.logORCDRC.time <- system.time(IntL.logORCDRC <- sparsereg3D.sel(sparse.reg = IntL.logORCDRC.preproc ,lambda = seq(0,0.2,0.001)))
#logORCDRC.l.pred <- sparsereg3D.pred(model.info = IntL.logORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntP.logORCDRC.preproc <- pre.sparsereg3D(base.model = log.ORCDRC.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, seed = seed, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
IntP.logORCDRC.ncv.time <- system.time(IntP.logORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = IntP.logORCDRC.preproc, lambda = seq(0,0.2,0.001), w = NULL))
IntP.logORCDRC.time <- system.time(IntP.logORCDRC <- sparsereg3D.sel(sparse.reg = IntP.logORCDRC.preproc ,lambda = seq(0,0.2,0.001)))
#logORCDRC.l.pred <- sparsereg3D.pred(model.info = IntP.logORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

rbind(BaseL.logORCDRC.ncv[1:2], BaseP.logORCDRC.ncv[1:2], IntL.logORCDRC.ncv[1:2], IntP.logORCDRC.ncv[1:2])

library(doParallel)

registerDoParallel(cores=4)

pack = (.packages())


result = foreach (j = 1:2, .packages = pack) %dopar% { 
  
  if (j==1){
    IntHL.logORCDRC.preproc <- pre.sparsereg3D(base.model = log.ORCDRC.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, seed = seed, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
    IntHL.logORCDRC.ncv.time <- system.time(IntHL.logORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.logORCDRC.preproc, lambda = seq(0,5,0.1)))
    IntHL.logORCDRC.time <- system.time(IntHL.logORCDRC <- sparsereg3D.sel(sparse.reg = IntHL.logORCDRC.preproc ,lambda = seq(0,5,0.1)))
    #logORCDRC.l.pred <- sparsereg3D.pred(model.info = IntHL.logORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))
    list(IntHL.logORCDRC.preproc, IntHL.logORCDRC.ncv.time, IntHL.logORCDRC.ncv, IntHL.logORCDRC.time, IntHL.logORCDRC)
  } else if (j==2) {
    IntHP.logORCDRC.preproc <- pre.sparsereg3D(base.model = log.ORCDRC.fun, use.hier = TRUE, profiles = bor.profs, use.interactions = TRUE, seed = seed, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
    IntHP.logORCDRC.ncv.time <- system.time(IntHP.logORCDRC.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.logORCDRC.preproc, lambda = seq(0,5,0.1)))
    IntHP.logORCDRC.time <- system.time(IntHP.logORCDRC <- sparsereg3D.sel(sparse.reg = IntHP.logORCDRC.preproc ,lambda = seq(0,5,0.1)))
    #logORCDRC.l.pred <- sparsereg3D.pred(model.info = IntHP.logORCDRC, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))
    list(IntHP.logORCDRC.preproc, IntHP.logORCDRC.ncv.time, IntHP.logORCDRC.ncv, IntHP.logORCDRC.time, IntHP.logORCDRC)
  }
  
}

stopImplicitCluster()

IntHL.logORCDRC.preproc <- result[[1]][[1]]
IntHL.logORCDRC.ncv.time <- result[[1]][[2]] 
IntHL.logORCDRC.ncv <- result[[1]][[3]]
IntHL.logORCDRC.time <- result[[1]][[4]]
IntHL.logORCDRC <- result[[1]][[5]]

IntHP.logORCDRC.preproc <- result[[2]][[1]]
IntHP.logORCDRC.ncv.time <- result[[2]][[2]] 
IntHP.logORCDRC.ncv <- result[[2]][[3]]
IntHP.logORCDRC.time <- result[[2]][[4]]
IntHP.logORCDRC <- result[[2]][[5]]

#============================= Observed vs Predicted ================================================================================================
#base.model = log.ORCDRC.fun; use.hier = FALSE; profiles = bor.profs; use.interactions = FALSE; poly.deg = 1; num.folds = 20; num.means = 1; cov.grids = gridmaps.sm2D

#df <- BaseP.logORCDRC.ncv$data.prediction

#previous.theme = theme_set(theme_bw()) #set black and white ggplot theme

#Define data to be plotted
#dfplot = ggplot(df, aes(x = logORCDRC, y = pred))

#Assemble a x-axis title. The 'atop' function allows you to have two lines of 
#text in a pasted-together expression
#my.xlab = expression(atop(paste("Observed log(1+SOC)")))
#Assemble a y-axis title using the same method. Two spaces are needed after 
#"Temperature" to make proper spacing on the y-axis title for some reason
#my.ylab = expression(atop(paste("Predicted log(1+SOC)")))

#Open a new png device to print the figure out to (or use tiff, pdf, etc).
#png(filename = "BaseL.Obs.Pred.png", width = 600, height = 600, units = 'px')
#print(dfplot + 
        #Define point shape and set alpha transparency
        #geom_point(alpha = 1/5, color = alpha("red"), size = 3) +
        #draw a dashed line f unity with geom_abline
        #geom_abline(intercept = 0, slope = 1, linetype = 2) +
        #xlab(my.xlab) +  #insert the x-axis title
        #ylab(my.ylab) +  #insert the y-axis title
        #xlim(1,5) +  #set the x-axis limits explicitly
        #ylim(1,5)   #set the y-axis limits explicitly
        #Adjust the plot margins slightly

#) #end of print statement
#dev.off() #close the png device to save the figure. 

#============================= Coefficients tables for logORCDRC models ==============================================================================

ll <- length(IntL.logORCDRC$coefficients)
pp <- length(IntHL.logORCDRC$coefficients[,1])+1

cmL.logORCDRC <- data.frame(variable=IntHL.logORCDRC$coefficients[,1], BaseL.logORCDRC.me=BaseL.logORCDRC$coefficients[2:pp], IntL.logORCDRC.me=IntL.logORCDRC$coefficients[2:pp],IntL.logORCDRC.ie=c(IntL.logORCDRC$coefficients[(pp+1):ll],0),IntHL.logORCDRC.me=IntHL.logORCDRC$coefficients[,2],IntHL.logORCDRC.ie=IntHL.logORCDRC$coefficients[,3] )
stargazer(cmL.logORCDRC, summary = FALSE, digits = 2, type = "latex")

#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.logORCDRC$coefficients)
p <- length(IntHP.logORCDRC$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.logORCDRC <- data.frame(variable=IntHP.logORCDRC$coefficients[,1], BaseP.logORCDRC.me=BaseP.logORCDRC$coefficients[2:p], IntP.logORCDRC.me=IntP.logORCDRC$coefficients[2:p],IntP.logORCDRC.ie1=c(IntP.logORCDRC$coefficients[(p+1):l][i1],0,0,0),IntP.logORCDRC.ie2=c(IntP.logORCDRC$coefficients[(p+1):l][i2],0,0,0),IntP.logORCDRC.ie3=c(IntP.logORCDRC$coefficients[(p+1):l][i3],0,0,0),IntHP.logORCDRC.me=IntHP.logORCDRC$coefficients[,2],IntHP.logORCDRC.ie1=IntHP.logORCDRC$coefficients[,3],IntHP.logORCDRC.ie2=IntHP.logORCDRC$coefficients[,4],IntHP.logORCDRC.ie3=IntHP.logORCDRC$coefficients[,5] )
cmlogORCDRC <- cmP.logORCDRC[,c(1,7:10)]

stargazer(cmP.logORCDRC, summary = FALSE, digits = 2, type = "latex")

# Models comparison
logORCDRC.ncv <- data.frame(rbind(BaseL = BaseL.logORCDRC.ncv[1:2], BaseP = BaseP.logORCDRC.ncv[1:2], IntL = IntL.logORCDRC.ncv[1:2], IntP = IntP.logORCDRC.ncv[1:2], IntHL = IntHL.logORCDRC.ncv[1:2], IntHP = IntHP.logORCDRC.ncv[1:2]))
logORCDRC.ncv <- data.frame(Model = rownames(logORCDRC.ncv), logORCDRC.ncv)
rownames(logORCDRC.ncv) <- NULL
names(logORCDRC.ncv) <- c("Model","RMSE","R squared")
stargazer(logORCDRC.ncv, summary = FALSE, digits = 2, type = "text")

#Time table
logORCDRC.ncv.time <- rbind(BaseL.logORCDRC.ncv.time, BaseP.logORCDRC.ncv.time, IntL.logORCDRC.ncv.time, IntP.logORCDRC.ncv.time, IntHL.logORCDRC.ncv.time, IntHP.logORCDRC.ncv.time)
logORCDRC.time <- rbind(BaseL.logORCDRC.time, BaseP.logORCDRC.time, IntL.logORCDRC.time, IntP.logORCDRC.time, IntHL.logORCDRC.time, IntHP.logORCDRC.time)

#Number of coefficients

logORCDRC.n.coeffs <- n.coeffs(l.coeffs = cmL.logORCDRC, p.coeffs = cmP.logORCDRC)

stargazer(cbind(logORCDRC.n.coeffs, logORCDRC.ncv[,-1], pH.n.coeffs[,-1], pH.ncv[,-1]), summary = FALSE, digits = 2, type = "latex")

#logORCDRC.results <- list(logORCDRC.ncv = logORCDRC.ncv, l.coeffs = cmL.logORCDRC, p.coeffs = cmP.logORCDRC, IntP.logORCDRC = IntP.logORCDRC, IntHP.logORCDRC = IntHP.logORCDRC, BaseP.logORCDRC = BaseP.logORCDRC)
#save(logORCDRC.results, file = "logORCDRC.results.rda" )
#==============================================================================================================================================


seed = 1

pH.fun <- as.formula(paste("pH ~", paste(c(CovNames,"depth"), collapse="+")))

# pH results
BaseL.pH.preproc <- pre.sparsereg3D(base.model = pH.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, seed = seed, poly.deg = 1, num.folds = 5, num.means = 5, cov.grids = gridmaps.sm2D) 
BaseL.pH.ncv.time <- system.time(BaseL.pH.ncv <- sparsereg3D.ncv(sparse.reg = BaseL.pH.preproc, lambda = seq(0,0.2,0.001)))
BaseL.pH.time <- system.time(BaseL.pH <- sparsereg3D.sel(sparse.reg = BaseL.pH.preproc ,lambda = seq(0,0.2,0.001)))
#pH.l.pred <- sparsereg3D.pred(model.info = BaseL.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

BaseP.pH.preproc <- pre.sparsereg3D(base.model = pH.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = FALSE, seed = seed, poly.deg = 3, num.folds = 5, num.means = 5, cov.grids = gridmaps.sm2D)   
BaseP.pH.ncv.time <- system.time(BaseP.pH.ncv <- sparsereg3D.ncv(sparse.reg = BaseP.pH.preproc, lambda = seq(0,0.2,0.001)))
BaseP.pH.time <- system.time(BaseP.pH <- sparsereg3D.sel(sparse.reg = BaseP.pH.preproc ,lambda = seq(0,0.2,0.001)))
#pH.p.pred <- sparsereg3D.pred(model.info = pH.p.sel, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))


IntL.pH.preproc <- pre.sparsereg3D(base.model = pH.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, seed = seed, poly.deg = 1, num.folds = 5, num.means = 5, cov.grids = gridmaps.sm2D)   
IntL.pH.ncv.time <- system.time(IntL.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntL.pH.preproc, lambda = seq(0,0.2,0.001)))
IntL.pH.time <- system.time(IntL.pH <- sparsereg3D.sel(sparse.reg = IntL.pH.preproc ,lambda = seq(0,0.2,0.001)))
#pH.l.pred <- sparsereg3D.pred(model.info = IntL.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

IntP.pH.preproc <- pre.sparsereg3D(base.model = pH.fun, use.hier = FALSE, profiles = bor.profs, use.interactions = TRUE, seed = seed, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)#, kmean.vars = all.vars(log.ORCDRC.fun), cum.prop = 0.90)    
IntP.pH.ncv.time <- system.time(IntP.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntP.pH.preproc, lambda = seq(0,0.2,0.001)))
IntP.pH.time <- system.time(IntP.pH <- sparsereg3D.sel(sparse.reg = IntP.pH.preproc ,lambda = seq(0,0.2,0.001)))
#pH.l.pred <- sparsereg3D.pred(model.info = IntP.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))

rbind(BaseL.pH.ncv[1:2], BaseP.pH.ncv[1:2], IntL.pH.ncv[1:2], IntP.pH.ncv[1:2])

library(doParallel)

registerDoParallel(cores=4)

pack = (.packages())


pH.result = foreach (j = 1:2, .packages = pack) %dopar% { 
  
  if (j==1){
    IntHL.pH.preproc <- pre.sparsereg3D(base.model = pH.fun, use.hier = TRUE, seed = seed, profiles = bor.profs, use.interactions = TRUE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
    IntHL.pH.ncv.time <- system.time(IntHL.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntHL.pH.preproc, lambda = seq(0,5,0.1)))
    IntHL.pH.time <- system.time(IntHL.pH <- sparsereg3D.sel(sparse.reg = IntHL.pH.preproc ,lambda = seq(0,5,0.1)))
    #pH.l.pred <- sparsereg3D.pred(model.info = IntHL.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))
    list(IntHL.pH.preproc, IntHL.pH.ncv.time, IntHL.pH.ncv, IntHL.pH.time, IntHL.pH)
  } else if (j==2) {
    IntHP.pH.preproc <- pre.sparsereg3D(base.model = pH.fun, use.hier = TRUE, seed = seed, profiles = bor.profs, use.interactions = TRUE, poly.deg = 3, num.folds = 5, num.means = 3, cov.grids = gridmaps.sm2D)    
    IntHP.pH.ncv.time <- system.time(IntHP.pH.ncv <- sparsereg3D.ncv(sparse.reg = IntHP.pH.preproc, lambda = seq(0,5,0.1)))
    IntHP.pH.time <- system.time(IntHP.pH <- sparsereg3D.sel(sparse.reg = IntHP.pH.preproc ,lambda = seq(0,5,0.1)))
    #pH.l.pred <- sparsereg3D.pred(model.info = IntHP.pH, chunk.size = 20000, grids = gridmaps.sm2D, depths = c(-0.1,-0.2,-0.3))
    list(IntHP.pH.preproc, IntHP.pH.ncv.time, IntHP.pH.ncv, IntHP.pH.time, IntHP.pH)
  }
  
}

stopImplicitCluster()

IntHL.pH.preproc <- pH.result[[1]][[1]]
IntHL.pH.ncv.time <- pH.result[[1]][[2]] 
IntHL.pH.ncv <- pH.result[[1]][[3]]
IntHL.pH.time <- pH.result[[1]][[4]]
IntHL.pH <- pH.result[[1]][[5]]

IntHP.pH.preproc <- pH.result[[2]][[1]]
IntHP.pH.ncv.time <- pH.result[[2]][[2]] 
IntHP.pH.ncv <- pH.result[[2]][[3]]
IntHP.pH.time <- pH.result[[2]][[4]]
IntHP.pH <- pH.result[[2]][[5]]


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

stargazer(cmP.pH, summary = FALSE, digits = 2, type = "latex")
#==============================================================================================================================================

# Models comparison
pH.ncv <- data.frame(rbind(BaseL = BaseL.pH.ncv[1:2], BaseP = BaseP.pH.ncv[1:2], IntL = IntL.pH.ncv[1:2], IntP = IntP.pH.ncv[1:2], IntHL = IntHL.pH.ncv[1:2], IntHP = IntHP.pH.ncv[1:2]))
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
  
  
  total <- data.frame(rbind(n.BaseL,n.BaseP,n.IntL,n.IntP,n.IntHL,n.IntHP))
  rownames(total) <- NULL
  total <- cbind(model = c("BaseL","BaseP", "IntL","IntP", "IntHL","IntHP"), total)
  return(total)
}

pH.n.coeffs <- n.coeffs(l.coeffs = cmL.pH, p.coeffs = cmP.pH)


save.image(file = "D:/R_projects/Bor_results_final.RData" )

#load("D:/R_projects/Bor_results_final.RData")

#============================== Depth Accuracy ======================================================

ORCDRC.fun <- as.formula(paste("ORCDRC ~", paste(c(CovNames,"depth"), collapse="+")))
log.ORCDRC.fun <- as.formula(paste("logORCDRC ~", paste(c(CovNames,"depth"), collapse="+")))
pH.fun <- as.formula(paste("pH ~", paste(c(CovNames,"depth"), collapse="+")))

formula = log.ORCDRC.fun; depth.inc = 20; profiles = bor.profs; cov.grids = gridmaps.sm2D; num.folds = 5; num.means = 3; poly.deg = 3; seed = 1; lambda.seq = seq(0,0.2,0.001); weights = NULL

intbasedif <- function(formula, depth.inc, profiles, cov.grids, poly = 3, num.folds = 5, num.means = 3, poly.deg = poly.deg, seed = 1, lambda.seq = seq(0,0.2,0.001), weights = NULL){
  variable <- all.vars(formula)[1]
  seed = seed
  max.depth = max(profiles$Bottom)
  depth.seq = c(0, sort(seq(depth.inc, (round(max.depth/depth.inc)-1)*depth.inc, depth.inc), decreasing = FALSE))
  pH.seq = data.frame()
  IntP.seq = data.frame()
  BaseP.seq = data.frame()
  profiles.i = profiles
  
  
  for(i in 1:length(depth.seq)){
    profiles.i@horizons <- profiles@horizons[profiles@horizons$Bottom <= round(max.depth - depth.seq[i]),]
    
    IntP.preproc <- pre.sparsereg3D(base.model = formula, use.hier = FALSE, profiles = profiles.i, use.interactions = TRUE, poly.deg = poly.deg, num.folds = num.folds, num.means = num.means, cov.grids = cov.grids, seed = seed)
    IntP.seq <- rbind(IntP.seq, sparsereg3D.ncv(sparse.reg = IntP.preproc, lambda = lambda.seq, w = weights)[1:2])
    
    BaseP.preproc <- pre.sparsereg3D(base.model = formula, use.hier = FALSE, profiles = profiles.i, use.interactions = FALSE, poly.deg = poly.deg, num.folds = num.folds, num.means = num.means, cov.grids = cov.grids, seed = seed) 
    BaseP.seq <- rbind(BaseP.seq, sparsereg3D.ncv(sparse.reg = BaseP.preproc, lambda = lambda.seq, w = weights)[1:2])
  }
  IntP.data <- data.frame(Depth = c(max.depth, sort(depth.seq[-1], decreasing = TRUE)), RMSE = IntP.seq$RMSE, variable = rep("IntP",length(depth.seq)))
  BaseP.data <- data.frame(Depth = c(max.depth, sort(depth.seq[-1], decreasing = TRUE)), RMSE = BaseP.seq$RMSE, variable = rep("BaseP",length(depth.seq)))
  
  
  BIP.data <- rbind(BaseP.data,IntP.data)
  BIP.data$variable <- factor(BIP.data$variable)
  plot <- qplot(Depth, RMSE, data = BIP.data, geom = c("point", "line"), color = variable) + theme_bw() + theme(legend.text = element_text(size = 12)) + theme(axis.text = element_text(size=12), axis.title = element_text(size = 14, face = "bold")) + labs(x = "Depth [cm]",y = "RMSE")
  
  return(plot)
}

logSOC.dif <- intbasedif(formula = log.ORCDRC.fun, depth.inc = 20, profiles = bor.profs, cov.grids = gridmaps.sm2D, num.folds = 5, num.means = 3, poly.deg = 3, seed = 1, lambda.seq = seq(0,0.2,0.001), weights = NULL)
pH.dif <- intbasedif(formula = pH.fun, depth.inc = 20, profiles = bor.profs, cov.grids = gridmaps.sm2D, num.folds = 5, num.means = 3, poly.deg = 3, seed = 1, lambda.seq = seq(0,0.2,0.001), weights = NULL)



pdf("BOR_logSOC.dif.pdf",width=10,height=6)
logSOC.dif
dev.off()


pdf("BOR_pH.dif.pdf",width=10,height=6)
pH.dif
dev.off()




#============================== SOM coef path plot ==================================================
depth <- data.frame(depth=seq(-0.1,-0.30,-0.1))
depth.mean <- as.data.frame(t(IntHP.pH$std.par$cnt.par$mean))[,"depth"]
depth.sd <- as.data.frame(t(IntHP.pH$std.par$cnt.par$std))[,"depth"]

depth.s <- as.numeric(((depth-depth.mean)/depth.sd)[,1])#              as.numeric(predict(IntHP.ORCDRC$std.par , newdata = depth)[,1])
#variables <- IntHP.ORCDRC$std.par$dummy.par$vars[-which(IntHP.ORCDRC$std.par$dummy.par$vars %in% c(IntHP.ORCDRC$std.par$dummy.par$facVars,"depth", "depth2","depth3"))] #
variables <- IntHP.pH$model$main.effect.names[-which(IntHP.pH$model$main.effect.names %in% c("depth", "depth2","depth3"))] #
variables <- variables[order(variables)]

coefs.IntHP <- (cmP.pH[which(as.character(cmP.pH$variable) %in% variables),c(1,7:10)])
coefs.IntP <- (cmP.pH[which(as.character(cmP.pH$variable) %in% variables),c(1,3:6)])

coefs.IntHP <- coefs.IntHP[order(coefs.IntHP$variable),]
coefs.IntP <- coefs.IntP[order(coefs.IntP$variable),]


effects.table <- data.frame()

for(i in 1:dim(coefs.IntHP)[1]){
  effects <- (data.frame(coefs.IntHP[i,1], t(coefs.IntHP[i,2] + coefs.IntHP[i,3]*depth.s + coefs.IntHP[i,4]*depth.s^2 + coefs.IntHP[i,5]*depth.s^3)))
  names(effects) <- c("variable", as.character(as.numeric(depth[,1])))
  effects.table <- rbind(effects.table, effects)
}

effects.table <- effects.table[match(IntHP.ORCDRC$model$main.effect.names[-which(IntHP.ORCDRC$model$main.effect.names %in% c("depth", "depth2","depth3"))], effects.table$variable),]
pH.effects.table <- effects.table
soc.effects.table <- effects.table

effects.table <- data.frame(soc.effects.table, pH.effects.table[,-1])

stargazer(effects.table, summary = FALSE, digits = 2, type = "latex")


#plot
coefs.plot <- coefs.IntP[,-1]
zero.effects <- as.numeric(which(apply(coefs.plot,1,sum) == 0))
plot.variables <- as.character(coefs.IntHP[-zero.effects,"variable"])
coefs.plot <- coefs.plot[-zero.effects,]
effects.plot <- as.numeric()

for(i in 1:dim(coefs.plot)[1]){
  effects <- coefs.plot[i,1] + coefs.plot[i,2]*depth.s + coefs.plot[i,3]*depth.s^2 + coefs.plot[i,4]*depth.s^3
  effects <- data.frame(depth = depth[,1], var = effects)
  names(effects) <- c("depth","var")
  effects.plot <- rbind(effects.plot, effects)
}

effects.plot$Variables <- factor(rep(plot.variables, each = length(depth.s)))
Plot.effects <- qplot(var, depth, data = effects.plot, geom = c("point", "line"), color = Variables)+theme_bw()+theme(legend.text=element_text(size=12))+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+labs(x="Coefficients",y="Depth")

pdf("Plot.effects.pdf",width=8,height=10)
intHPlot
dev.off()

pdf("SOMIntPlot.pdf",width=8,height=10)
intPlot
dev.off()


















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

