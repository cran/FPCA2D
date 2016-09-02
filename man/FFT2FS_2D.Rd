\name{FFT2FS_2D}
\alias{FFT2FS_2D}
\title{
   Two Dimenstional Fourier Expansion
}
\description{
   Conduct the two dimensional Fourier expansion by using fast Fourier transformation.
}
\usage{
FFT2FS_2D(A)
}
\arguments{
  \item{A}{
   A is a two dimensional numerical matrix.   
}
}
\details{
  Conduct the two dimensional Fourier expansion by using fast Fourier transformation.  To make the results comparable, be sure to 
  transform all the intensity values of each image pixel should be between 0 and 1.
}
\value{
 \item{C}{The two dimensional Fourier expansion of the input matrix}
}
\references{
  Lin N, Jiang J, Guo S, Xiong M (2015) Functional Principal Component Analysis and Randomized Sparse Clustering Algorithm for Medical Image Analysis. PLoS ONE 10(7): e0132945. doi:10.1371/journal.pone.0132945
}

\examples{
   \dontrun{
   A = matrix(sample(seq(0,1,0.01),100),nrow=10)
   rlt = FFT2FS_2D(A)
           } 
}
