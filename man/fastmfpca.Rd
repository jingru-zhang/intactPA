\name{fastmfpca}
\alias{fastmfpca}
\title{Fast multilevel functional principal component analysis}

\description{
\code{fastmfpca} performs fast multilevel functional principal component analysis 
for longitudinal or repeated-measures functional data. It decomposes subject-level 
and visit-level variability and provides functional principal component (FPC) 
approximations of the observed curves, even in the presence of missing values.
}

\usage{
fastmfpca(
  Y,
  id,
  visit = NULL,
  twoway = TRUE,
  weight = "obs",
  argvals = NULL,
  pve = 0.99,
  npc = NULL,
  pve2 = NULL,
  npc2 = NULL,
  p = 3,
  m = 2,
  knots = 35,
  lambda.Gt = NULL,
  lambda.Gw = NULL,
  silent = TRUE
)
}

\arguments{
  \item{Y}{A numeric matrix of multilevel functional data observed on a regular grid. 
  Each row corresponds to one visit for one subject. Missing values are allowed and 
  must be labeled as \code{NA}. Required.}

  \item{id}{A vector of subject identifiers corresponding to the rows of \code{Y}. 
  Required.}

  \item{visit}{An optional vector of visit identifiers. If not provided, visits 
  are assumed to be indexed sequentially as 1, 2, \dots{} for each subject.}

  \item{twoway}{Logical; if \code{TRUE}, perform two-way ANOVA decomposition and 
  compute visit-specific means. Defaults to \code{TRUE}.}

  \item{weight}{Specifies how to weight sample covariances. 
  If \code{"obs"}, weights are proportional to the number of observations; 
  if \code{"subj"}, each subject is weighted equally. Defaults to \code{"obs"}.}

  \item{argvals}{Optional numeric vector of observed locations in the functional domain.}

  \item{pve}{Proportion of variance explained used to select the number of principal 
  components for both levels. Defaults to \code{0.99}.}

  \item{npc}{Optional integer giving a pre-specified number of principal components 
  at level 1 (subject level). Overrides \code{pve} if supplied.}

  \item{pve2}{Optional proportion of variance explained used to select the number of 
  components at level 2 (visit level). Overrides \code{pve} if supplied.}

  \item{npc2}{Optional integer giving a pre-specified number of principal components 
  at level 2. Overrides \code{npc} if supplied.}

  \item{p}{Degree of B-spline basis functions used. Defaults to 3.}

  \item{m}{Order of the difference penalty applied to the spline basis. Defaults to 2.}

  \item{knots}{Number of knots or a vector of knot locations used in the B-spline basis. 
  Defaults to 35.}

  \item{lambda.Gt}{Optional tuning parameter for FACE estimation of the total covariance.}

  \item{lambda.Gw}{Optional tuning parameter for FACE estimation of the within-subject covariance.}

  \item{silent}{Logical; if \code{TRUE}, suppress progress messages. Defaults to \code{TRUE}.}
}

\value{
A list containing:

\describe{
  \item{Yhat}{Functional principal component (FPC) approximation of \code{Y}; 
  fitted curves for all subjects and visits.}

  \item{Y.df}{The observed data matrix with missing values imputed using FACE.}

  \item{Thetamat}{The covariance matrix used for spectral decomposition in FACE.}

  \item{funbasis}{The transformed cubic B-spline basis used in the decomposition.}

  \item{lambda_face}{A length-2 list of tuning parameters used by FACE for 
  estimating total and within-subject covariances, respectively.}
}
}

\references{
Zhang, J., Cui, E., Li, H., & Shou, H. (2025). 
\emph{INTACT: A method for integration of longitudinal physical activity data from multiple sources}. 
Biometrics.

Cui, E., Li, R., Crainiceanu, C. M., & Xiao, L. (2023). 
Fast multilevel functional principal component analysis. 
\emph{Journal of Computational and Graphical Statistics}, 32(2), 366--377.
}

\seealso{
\code{\link{INTACT}}
}



