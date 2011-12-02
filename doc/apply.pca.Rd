\name{apply.pca}
\alias{apply.pca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Apply PCA transformation to probe level CNV data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
apply.pca(.Object, multiprobe.signal = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{.Object}{
%%     ~~Describe \code{.Object} here~~
}
  \item{multiprobe.signal}{
%%     ~~Describe \code{multiprobe.signal} here~~
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
structure(function (.Object, multiprobe.signal = NULL) 
{
    standardGeneric("apply.pca")
  }, generic = structure("apply.pca", package = "CNVtools"), package = "CNVtools", group = list(), valueClass = character(0), signature = c(".Object", 
"multiprobe.signal"), default = `NULL`, skeleton = function (.Object, 
    multiprobe.signal = NULL) 
stop("invalid call in method dispatch to \"apply.pca\" (no default method)", 
    domain = NA)(.Object, multiprobe.signal), class = structure("nonstandardGenericFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
