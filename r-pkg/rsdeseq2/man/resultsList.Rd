\name{resultsList}
\alias{resultsList}
\title{Evaluate several result requests from one fitted object}
\usage{
resultsList(object, requests, ..., batchSize = 32L)
}
\arguments{
\item{object}{A dataset fitted by \code{rsdeseq2::DESeq}, a Rust fit, a primitive Wald/LRT fit, or supported list-like primitive inputs.}
\item{requests}{A list of argument lists for \code{results}. An empty argument list requests default results. Outer names are preserved. Argument names within each request must be unique and exact.}
\item{...}{Common \code{results} arguments. Each request overrides the common values, including an explicit \code{NULL}.}
\item{batchSize}{An integer from 1 to 64 specifying the maximum native inference requests per chunk.}
}
\description{
Validates the fitted dataset once and reuses aligned primitive arrays across
bounded native contrast batches. Every result keeps its own inference options,
Cook's masking, alpha and independent-filtering/BH decisions.
}
\details{
Requests support the selectors and options of \code{results}, including named,
numeric, list and factor contrasts. Primitive inputs retain their contextual
requirements; explicit Cook's filtering requires retained fitted-dataset context.
LRT requests retain their stored full-versus-reduced statistic and p-value.

The return values have the same shape, class, metadata and gene order as the
corresponding individual calls. Native buffers are borrowed and validated per
bounded chunk; numerical scratch is consumed one request at a time. Storage for
the requested returned tables grows with their total size.
}
\value{A list of result tables in request order, preserving the names of \code{requests}. Dataset inputs return the same Bioconductor result classes as \code{results}.}
\examples{
fit <- waldFitRust(c(10, 20),
  cbind(Intercept = c(1, 1), condition = c(2, -1)),
  betaSE = matrix(0.5, 2, 2))
resultsList(fit, list(
  forward = list(name = "condition"),
  reverse = list(contrast = c(0, -1), alpha = 0.05)
))
}
