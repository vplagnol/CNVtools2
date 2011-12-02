\name{CNVtools.multivariate.signal}
\alias{CNVtools.multivariate.signal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a CNVtools object from multi-dimensional probe level data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CNVtools.multivariate.signal(multivariate.signal, trait, batch, model.mean = formula("~ strata(cn)"), model.var = formula("~  strata(cn)"), model.association = formula("~  cn"), IID = NULL, FID = NULL, f.member = NULL, covariates = data.frame())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{multivariate.signal}{
%%     ~~Describe \code{multivariate.signal} here~~
}
  \item{trait}{
%%     ~~Describe \code{trait} here~~
}
  \item{batch}{
%%     ~~Describe \code{batch} here~~
}
  \item{model.mean}{
%%     ~~Describe \code{model.mean} here~~
}
  \item{model.var}{
%%     ~~Describe \code{model.var} here~~
}
  \item{model.association}{
%%     ~~Describe \code{model.association} here~~
}
  \item{IID}{
%%     ~~Describe \code{IID} here~~
}
  \item{FID}{
%%     ~~Describe \code{FID} here~~
}
  \item{f.member}{
%%     ~~Describe \code{f.member} here~~
}
  \item{covariates}{
%%     ~~Describe \code{covariates} here~~
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
function (multivariate.signal, trait, batch, model.mean = formula("~ strata(cn)"), 
    model.var = formula("~  strata(cn)"), model.association = formula("~  cn"), 
    IID = NULL, FID = NULL, f.member = NULL, covariates = data.frame()) 
{
    if (class(multivariate.signal) != "matrix") 
        stop("Multivariate input signal must be a matrix")
    mean.signal <- apply(multivariate.signal, MAR = 1, FUN = mean, 
        na.rm = TRUE)
    res <- new("CNVtools", trait = trait, signal = mean.signal, 
        batch = batch, IID = IID, FID = FID, f.member = f.member, 
        model.mean = model.mean, model.var = model.var, model.association = model.association)
    res@multiprobe.signal <- multivariate.signal
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
