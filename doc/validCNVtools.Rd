\name{validCNVtools}
\alias{validCNVtools}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check whether a CNVtools object is valid
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
validCNVtools(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
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
function (object) 
{
    if (length(object@signal) != object@nsamples) {
        message("Lengths of signal and nsamples do not match")
        stop()
    }
    if (length(object@trait) != object@nsamples) {
        message("Lengths of trait and nsamples do not match")
        stop()
    }
    if (length(object@batch) != object@nsamples) {
        message("Lengths of batch and nsamples do not match")
        stop()
    }
    if ((nrow(object@covariates) > 0) && (nrow(object@covariates) != 
        object@nsamples)) {
        message("The covariates matrix does not have the right size")
        stop()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
