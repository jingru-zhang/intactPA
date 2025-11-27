\name{INTACTforXW}
\alias{INTACTforXW}
\title{Integration of residuals from longitudinal physical activity data across multiple sources}

\description{
\code{INTACTforXW} provides harmonized residuals at the subject level (X part) and 
visit level (W part) from longitudinal physical activity data collected across batches/sources. 
It is a subroutine of the \code{\link{INTACT}} method.
}

\usage{
INTACTforXW(
  mydf,
  featurenames,
  daynum,
  formula,
  loga0 = seq(-1.5, 1.5, 0.5),
  cw = 0.01,
  ro = c(10, 0.9),
  itermax = 100,
  nfold = 1,
  dd = NULL,
  endpara = 1e-3,
  th.value = 0.01,
  wpve = 0.5
)
}

\arguments{
  \item{mydf}{A \code{data.frame} containing the subject ID, covariates for the fixed-effects model, 
  batch label, and scaled original data. Each row corresponds to one observation for one subject 
  on one day; columns are variables.}

  \item{featurenames}{A character vector giving the names of the physical activity features.}

  \item{daynum}{A list where each element is the number of repeated measures for a given subject 
  from each source.}

  \item{formula}{A character string specifying the fixed-effects part of the linear mixed-effects 
  model, written in \pkg{lme4} formula notation. It should include covariates and interactions 
  (e.g., \code{"age + sex"}), but must \emph{not} include the batch variable or random effects.}

  \item{loga0}{A numeric vector on the log scale; its exponential defines the candidate range for 
  the penalty parameter \code{a0}, from which the optimal value is selected via cross-validation.}

  \item{cw}{Regularization parameter for the smoothing term.}

  \item{ro}{Length-2 numeric vector giving the regularization parameters for the proximal term.}

  \item{itermax}{Maximum number of iterations allowed.}

  \item{nfold}{Number of folds for cross-validation.}

  \item{dd}{Optional length-2 numeric vector specifying pre-determined numbers of components 
  for the X and W parts, respectively. If provided, this overrides \code{wpve}.}

  \item{endpara}{Convergence tolerance for terminating the iteration.}

  \item{th.value}{Threshold applied to constrain the selection of the number of components.}

  \item{wpve}{Weight (proportion of variance explained) used to select the number of components 
  when \code{dd} is not specified.}
}

\value{
A list with the following components:

\describe{
  \item{Xint}{Harmonized subject-level (X part) residuals (before scaling back) contained in \code{Eta}.}

  \item{Wint}{Harmonized visit-level (W part) residuals (before scaling back) contained in \code{Eta}.}

  \item{Sx}{A list of subject-level (X part) covariance matrices; each element is a Theta matrix 
  for one source, derived from the \code{fastmfpca} function.}

  \item{Sw}{A list of visit-level (W part) covariance matrices; each element is a Theta matrix 
  for one source, derived from the \code{fastmfpca} function.}

  \item{eps}{A list of noise residuals, with one element per source.}

  \item{predicted}{A list of predicted values (before scaling back) based on covariates and 
  fixed-effects estimation, with one element per source.}

  \item{funbasis}{The transformed cubic B-spline basis obtained from the \code{fastmfpca} function.}

  \item{Eta}{Estimated residual term in the raw data.}

  \item{elapsed.int}{Computation time for \emph{INTACT} (and \emph{INTACT0}). This includes 
  the time for fixed-effects estimation (\code{elapsed.fixedeff}), cross-validation 
  (\code{elapsed.cv}), and common eigenspace estimation (\code{elapsed.intact}).}

  \item{{d}}{Length-2 numeric vector giving the numbers of components for the X- and W-parts, respectively.}

  \item{{Xi.x}}{Estimated common eigenfunctions for the subject-level (X) part.}

  \item{{Xi.w}}{Estimated common eigenfunctions for the visit-level (W) part.}

  \item{{fixeffs}}{Estimated common fixed effects in the linear mixed-effects model.}

  \item{{batch_effects_adjusted}}{Estimated source-specific mean-shift curves in the linear mixed-effects model.}
}
}

\references{
Zhang, J., Cui, E., Li, H., & Shou, H. (2025). 
\emph{INTACT: A method for integration of longitudinal physical activity data from multiple sources}. Biometrics.
}

\seealso{
\code{\link{INTACT}}
}




