\name{compact.data.frame}
\alias{compact.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Makes a large data frame used by the CNVtools C code more compact and easier to read
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compact.data.frame(full.frame)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{full.frame}{
%%     ~~Describe \code{full.frame} here~~
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
function (full.frame) 
{
    full.frame <- full.frame[order(full.frame$cn), ]
    full.frame.mod <- do.call(rbind.data.frame, split(x = full.frame$posterior, 
        f = full.frame$IID))
    names(full.frame.mod) <- paste("P", c(1:(dim(full.frame.mod)[2])), 
        sep = "")
    full.frame.mod$cn <- apply(full.frame.mod, FUN = which.max, 
        MAR = 1)
    full.frame.mod$IID <- row.names(full.frame.mod)
    full.frame <- subset(full.frame[, c("IID", "batch", "signal", 
        "trait", "proba.not.outlier")], full.frame$cn == 1)
    full.frame <- merge(full.frame, full.frame.mod)
    return(full.frame)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
