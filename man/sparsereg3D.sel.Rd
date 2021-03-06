% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparsereg3D.sel.r
\name{sparsereg3D.sel}
\alias{sparsereg3D.sel}
\title{Model Selection based on Cross-Validation}
\usage{
sparsereg3D.sel(sparse.reg, lambda = 0, ols = FALSE, step = FALSE,
  lambda.1se = FALSE)
}
\arguments{
\item{sparse.reg}{output from \code{pre.sparsereg3D} function}

\item{lambda}{vector of regularization parameter values for lasso regression}

\item{ols}{logical. If TRUE model will be fitted with OLS insted of using lasso}

\item{step}{logical. If TRUE stepwise procedure will be used when fitting OLS model}

\item{lambda.1se}{logical. If TRUE one sigma lambda rule will be used (largest lambda value with cv.err less than or equal to min(cv.err)+ SE).}

\item{seed}{random number generator}
}
\value{
List of objects including:
\itemize{
 \item \code{model}:  list of objects containing model description
 \item \code{lambda}:   Regularization parameter value for lasso models
 \item \code{coefficients}:  Model coefficients
 \item \code{std.param}:   Standardization parameters
 \item \code{comp}: logical. Indicate whether the compositional variable is modeled or not.
}

 @keywords Model selection
}
\description{
Model Selection based on Cross-Validation
}
