% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.R
\name{cv_Mean}
\alias{cv_Mean}
\title{Cross-validation to select bandwidth for Mean of FDA}
\usage{
cv_Mean(data, h_seq, k = 5, d = 1)
}
\arguments{
\item{data}{a list containing matrices, each of which represents a sample of the
data in the FDA. Each matrix should have a time column as the first column,
followed by p columns response values.}

\item{h_seq}{a sequence of candidate bandwidths to be evaluated in the
cross-validation.}

\item{k}{the number of folds in the cross-validation.}

\item{d}{the degree of the polynomial trend to be fitted in the mean estimation.}
}
\value{
The value of \code{h} that minimizes the cross-validation mean squared error.
}
\description{
This function uses cross-validation to select the bandwidth parameter \code{h}
for the mean estimation method.
}
\examples{
h_seq <- seq(0.1,1,length=10)
cv_Mean(FDAdata,h_seq,k=5,d=1)

}
