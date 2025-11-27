\name{getMVDV}
\alias{getMVDV}
\title{Quantify estimation accuracy of covariance structures from harmonized data}

\description{
\code{getMVDV} evaluates the accuracy of estimated covariance structures 
derived from harmonized data. It computes differences between the true 
common covariance and its estimates, as well as differences between 
covariance structures of harmonized data across multiple sources.
}

\usage{
getMVDV(Phixfun, Lamx, Phiwfun, Lamw, Y, multipar)
}

\arguments{
  \item{Phixfun}{A \eqn{p \times k_1} matrix of common subject-level 
  eigenfunctions (X part). Each column is one eigenfunction evaluated at 
  \eqn{p} time points.}

  \item{Lamx}{A \eqn{k_1 \times k_1} diagonal matrix of subject-level 
  eigenvalues (X part).}

  \item{Phiwfun}{A \eqn{p \times k_2} matrix of common visit-level 
  eigenfunctions (W part). Each column is one eigenfunction evaluated at 
  \eqn{p} time points.}

  \item{Lamw}{A \eqn{k_2 \times k_2} diagonal matrix of visit-level 
  eigenvalues (W part).}

  \item{Y}{A list of length \eqn{M}, where each element contains the harmonized 
  residuals from one source.}

  \item{multipar}{Overall scaling factor \eqn{c_0} for \emph{INTACT} (and 
  \emph{INTACT0}).}
}

\value{
A list containing:

\describe{
  \item{MV.X}{Squared Frobenius norm of the difference between the true 
  common covariance and the subject-level (X part) sample covariance from 
  harmonized data.}

  \item{MV.W}{Squared Frobenius norm of the difference between the true 
  common covariance and the visit-level (W part) sample covariance from 
  harmonized data.}

  \item{DV.X}{Differences between subject-level (X part) sample covariances 
  of harmonized data across multiple sources.}

  \item{DV.W}{Differences between visit-level (W part) sample covariances 
  of harmonized data across multiple sources.}

  \item{var.X}{A list of subject-level (X part) sample covariance matrices of 
  harmonized data, with one element per source.}

  \item{var.W}{A list of visit-level (W part) sample covariance matrices of 
  harmonized data, with one element per source.}
}
}

\references{
Zhang, J., Cui, E., Li, H., & Shou, H. (2025). 
\emph{INTACT: A method for integration of longitudinal physical activity data from multiple sources}. 
Biometrics.
}

\seealso{
\code{\link{INTACT}}
}


