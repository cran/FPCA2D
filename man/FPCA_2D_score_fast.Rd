\name{FPCA_2D_score_fast}
\alias{FPCA_2D_score_fast}
\title{
   Two Dimensional Functional Principal Component Analysis
}
\description{
   Calcualte the two dimensional functional principal component scores by using Fourier Basis
}
\usage{
FPCA_2D_score_fast(X)
}
\arguments{
  \item{X}{
   X is the input three dimensional array.  The first two dimensions are the dimension of each input image.  All the inputs images
   are organized as the third dimension of the input data array.  All the image should be scaled to the rage from 0 to 1 before running
   this function.
}
}
\details{
  Calcualte the two dimensional functional principal component scores by using Fourier Basis%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{eigen_value}{The eigen value can be used to calcualte the proportion of variance that each FPC score can explain.}
  \item{FPC_score}{The output FPC scores.}
  \item{Eigen_vector}{The eigen_vector represents the directions of the lienar transformation in the functional domain.}
}
\references{
Lin N, Jiang J, Guo S, Xiong M (2015) Functional Principal Component Analysis and Randomized Sparse Clustering Algorithm for Medical Image Analysis. PLoS ONE 10(7): e0132945. doi:10.1371/journal.pone.0132945
}

\examples{
   A = array(sample(seq(0,1,0.001),300),dim=c(10,10,3))
   rlt = FPCA_2D_score_fast(A)
}

