Sparsereg3D
================

Introduction
------------

This package implements the methodology proposed in the manuscript "Sparse Regression Interaction Models for Spatial Predictions of Soil Properties in 3D" by Pejovic M., Nikolic M., Heuvelink G., Hengl T., Kilibarda M. and Bajat B., that was subbmited to the journal "Computers and Geosciences" on August 2017. The proposed methodology utilizes lasso regression (standard and hierarchical) for making 3D interaction model of soil variables. Once the 3D model is created, spatial prediction can be made at different soil depths. The term "Interaction" in the title refers to the interaction effects between spatial covariates and soil depth, which are particularly important for modeling varying influences of external factors with soil depth. Broadly, the methodology is centered around two key points: 1) using lasso to select the best predicting linear model in cases where many varaibles (that arise from including interactions) are available; 2) propor model evaluation by using Nested Cross-Validation; and 3) spatial prediction over grids that refer to different soil depths.

How to install the package
--------------------------

The package on this repository is in its "raw" form, often referred to as 'source'. If you have Rtools and LaTeX installed, you will have no problem building the package from source by hand. This has been made easy using package devtools using the following commands:

``` r
install_github("hadley/devtools")
library(devtools)
install_github("pejovic/sparsereg3D")
library(sparsereg3D)
```

How to use 'sparsereg3D' package
--------------------------------

The presented example uses the Edgeroi data set \[<http://plotkml.r-forge.r-project.org/edgeroi.html>\]. This is one of the standard soil data sets used to test soil mapping methods in Australia. It contains 359 soil profiles with soil observations in sites and horizons tables;

``` r
library(GSIF)
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

# Resampling 100m grigs to 250m.

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

# Coercing categorical variables as factor
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

# Defining model
formulaString <- as.formula(paste(paste("log1pORC ~ "), paste(c(names(cov.maps),"depth"), collapse="+")))


# Spatial prediction with sparsereg3D

# Preprocessing the data
MLR.logORC.preproc <- pre.sparsereg3D(base.model = formulaString, use.hier = FALSE, profiles = edgeroi.spc, use.interactions = FALSE, poly.deg = 1, num.folds = 5, num.means = 3, cov.grids = cov.maps, seed = 1)

# Nested cross-validation
MLR.logORC.ncv <- sparsereg3D.ncv(sparse.reg = MLR.logORC.preproc, lambda = seq(0,0.2,0.001))

# Model selection
MLR.logORC <- sparsereg3D.sel(sparse.reg = MLR.logORC.preproc , lambda = seq(0,0,0), step = TRUE)

# Spatial prediction at 0.1m depth
MLR.pred <- sparsereg3D.pred(model.info = MLR.logORC, chunk.size = 20000, grids = cov.maps, depths = c(-0.1))
```
