\name{fit}
\alias{fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fit function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fit(object, ncomp, hyp = "H0", EM.starting.point = "kmeans", gaussian.or.T = "gaussian", model.association = NULL, model.mean = NULL, model.var = NULL, start.mean = NULL, start.var = NULL, start.alpha = NULL, control = list(tol = 1e-05, max.iter = 3000, min.freq = 4), force.replace = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{ncomp}{
%%     ~~Describe \code{ncomp} here~~
}
  \item{hyp}{
%%     ~~Describe \code{hyp} here~~
}
  \item{EM.starting.point}{
%%     ~~Describe \code{EM.starting.point} here~~
}
  \item{gaussian.or.T}{
%%     ~~Describe \code{gaussian.or.T} here~~
}
  \item{model.association}{
%%     ~~Describe \code{model.association} here~~
}
  \item{model.mean}{
%%     ~~Describe \code{model.mean} here~~
}
  \item{model.var}{
%%     ~~Describe \code{model.var} here~~
}
  \item{start.mean}{
%%     ~~Describe \code{start.mean} here~~
}
  \item{start.var}{
%%     ~~Describe \code{start.var} here~~
}
  \item{start.alpha}{
%%     ~~Describe \code{start.alpha} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
  \item{force.replace}{
%%     ~~Describe \code{force.replace} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (object, ncomp, hyp = "H0", EM.starting.point = "kmeans", 
    gaussian.or.T = "gaussian", model.association = NULL, model.mean = NULL, 
    model.var = NULL, start.mean = NULL, start.var = NULL, start.alpha = NULL, 
    control = list(tol = 1e-05, max.iter = 3000, min.freq = 4), 
    force.replace = FALSE) 
standardGeneric("fit"), generic = structure("fit", package = "CNVtools"), package = "CNVtools", group = list(), valueClass = character(0), signature = c("object", 
"ncomp", "hyp", "EM.starting.point", "gaussian.or.T", "model.association", 
"model.mean", "model.var", "start.mean", "start.var", "start.alpha", 
"control", "force.replace"), default = `NULL`, skeleton = function (object, 
    ncomp, hyp = "H0", EM.starting.point = "kmeans", gaussian.or.T = "gaussian", 
    model.association = NULL, model.mean = NULL, model.var = NULL, 
    start.mean = NULL, start.var = NULL, start.alpha = NULL, 
    control = list(tol = 1e-05, max.iter = 3000, min.freq = 4), 
    force.replace = FALSE) 
stop("invalid call in method dispatch to \"fit\" (no default method)", 
    domain = NA)(object, ncomp, hyp, EM.starting.point, gaussian.or.T, 
    model.association, model.mean, model.var, start.mean, start.var, 
    start.alpha, control, force.replace), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
