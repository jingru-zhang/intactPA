\name{INTACT}
\alias{INTACT}
\title{INTACT: Integration of longitudinal physical activity data from multiple sources}

\description{
Implements the INTACT method to harmonize longitudinal physical activity
features collected across sources while accounting for fixed effects
and within-subject correlation.
}

\usage{
INTACT(
  mydf0,
  Yraw,
  formula,
  needlog = 0,
  comparisons = NULL,
  act = 1,
  wpve = 0.5,
  dd = NULL,
  loga0 = seq(-1.5, 1.5, 0.5),
  cw = 0.01,
  ro = c(10, 0.9),
  itermax = 100,
  nfold = 1,
  endpara = 1e-3,
  th.value = 0.01
)
}

\arguments{
  \item{mydf0}{A \code{data.frame} containing subject ID, covariates for the
  fixed-effects model, and a batch label. Each row is one observation for one
  subject on one day; columns are variables.}

  \item{Yraw}{A numeric matrix of physical activity features to be harmonized. 
  Each row corresponds to one observation for a given subject on a given day, 
  and each column represents a distinct physical activity feature.}

  \item{formula}{A character string specifying the fixed-effects part of the
  linear mixed-effects model in \pkg{lme4} formula syntax (right-hand side only),
  including covariates and interactions, e.g., \code{"age + sex"}.
  Do \emph{not} include the batch variable or random effects.}

  \item{needlog}{A 0/1 indicator: \code{1} to log-transform features via
  \code{log(Yraw + 1)}, \code{0} otherwise.}

  \item{comparisons}{A 0/1 indicator: \code{1} to additionally return results
  from competing methods (raw data, simple standardization, quantile
  normalization, and longitudinal ComBat); \code{0} to return only INTACT and
  INTACT0. Default \code{NULL} behaves as \code{0}.}

  \item{act}{A 0/1 indicator: \code{1} to enforce nonnegativity of the
  harmonized features; \code{0} otherwise.}

  \item{wpve}{Weight (proportion of variance explained) used to select the
  number of components when \code{dd} is not provided.}

  \item{dd}{Optional length-2 numeric vector giving the pre-specified numbers of
  components for the X- and W-parts, respectively. If supplied, overrides
  \code{wpve}.}

  \item{loga0}{A numeric vector of candidate values on the log scale; their
  exponentials define the grid for the penalty parameter \code{a0}.
  The best value is selected via cross-validation.}

  \item{cw}{Regularization parameter for the smoothing term.}

  \item{ro}{Length-2 numeric vector of regularization parameters for the
  proximal term.}

  \item{itermax}{Maximum number of iterations allowed.}

  \item{nfold}{Number of folds for cross-validation.}

  \item{endpara}{Convergence tolerance for terminating the iteration.}

  \item{th.value}{Threshold used to constrain/guard the selection of the number
  of components.}
}

\details{
\code{INTACT} harmonizes longitudinal physical activity features across multiple
sources while preserving the biological and behavioral variation explained by
fixed effects and subject-level structure. The function returns results for both
\emph{INTACT} and \emph{INTACT0}, and, if requested, also provides outputs from 
several competing harmonization methods for comparison.
}

\value{
A named list with method-specific results and estimation results:

\preformatted{
list(
  reintact = ..., reintact0 = ...,
  ## present only if \code{comparisons == 1}:
  reraw = ..., resim = ..., requan = ..., relcombat = ...,
  ## always estimation results:
  d = ..., Xi.x = ..., Xi.w = ..., fixeffs = ..., batch_effects_adjusted = ...
)
}

\strong{Method-specific results}
(\code{reintact}, \code{reintact0}, and, if \code{comparisons == 1}, also
\code{reraw}, \code{resim}, \code{requan}, \code{relcombat}) — each is a list with:
\describe{
  \item{{Y}}{Numeric matrix of harmonized features. Rows are subject–day observations;
  columns are physical-activity features.}
  \item{{Eta}}{Residual term(s) after harmonization.}
  \item{{elapsed}}{Computation time. For \emph{INTACT}/\emph{INTACT0}, a list with
  fixed-effects estimation (\code{elapsed.fixedeff}), cross-validation (\code{elapsed.cv}),
  common-eigenspace estimation (\code{elapsed.intact}), and total time (\code{elapsed.total}).
  For competing methods, a single numeric giving total runtime.}
  \item{{multipar}}{Overall scaling factor \code{c0}; present for \emph{INTACT}/\emph{INTACT0} only.}
}

\strong{Estimation results}:
\describe{
  \item{{d}}{Length-2 numeric vector giving the estimated numbers of components for the X- and W-parts, respectively.}
  \item{{Xi.x}}{Estimated common eigenfunctions for the subject-level (X) part.}
  \item{{Xi.w}}{Estimated common eigenfunctions for the visit-level (W) part.}
  \item{{fixeffs}}{Estimated common fixed effects in the linear mixed-effects model.}
  \item{{batch_effects_adjusted}}{Estimated source-specific mean-shift curves in the linear mixed-effects model.}
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

\examples{
## The "example" dataset provides mydf0, Yraw, and formula.
data(example)

## Yraw does not need log-transform:
needlog <- 0

## Also compute competing methods for comparison:
comparisons <- 1

reall <- INTACT(mydf0, Yraw, formula, needlog = needlog, comparisons = comparisons)
}

