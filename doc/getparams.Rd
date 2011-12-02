\name{getparams}
\alias{getparams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the parameters from a CNVtools fit
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getparams(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
%%     ~~Describe \code{d} here~~
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
function (d) 
{
    p <- list()
    p[["ns"]] <- length(levels(d$batch))
    p[["nc"]] <- range(d$cn)[2]
    p[["nind"]] <- dim(d)[1]/p[["nc"]]
    maxLike <- tapply(d$pr, d$IID, max)
    p[["lnL"]] <- sum(maxLike + log(tapply(exp(d$pr - maxLike[d$IID]), 
        FUN = sum, d$IID)))
    p[["alpha"]] <- matrix(0, nrow = p$nc, ncol = p$ns)
    p[["mean"]] <- matrix(0, nrow = p$nc, ncol = p$ns)
    p[["var"]] <- matrix(0, nrow = p$nc, ncol = p$ns)
    p[["nu"]] <- matrix(0, nrow = p$nc, ncol = p$ns)
    p[["pdc"]] <- matrix(0, nrow = p$nc, ncol = p$ns)
    lev <- levels(d$batch)
    for (j in 1:p$ns) {
        for (i in 1:p$nc) {
            p$mean[i, j] <- .Call("get_first_match", nrow(d), 
                as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), 
                as.numeric(j), d$mean)
            p$alpha[i, j] <- .Call("get_first_match", nrow(d), 
                as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), 
                as.numeric(j), d$alpha)
            p$var[i, j] <- .Call("get_first_match", nrow(d), 
                as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), 
                as.numeric(j), d$var)
            p$nu[i, j] <- .Call("get_first_match", nrow(d), as.numeric(d$cn), 
                as.numeric(d$batch), as.numeric(i), as.numeric(j), 
                d$nu)
            p$pdc[i, j] <- .Call("get_first_match", nrow(d), 
                as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), 
                as.numeric(j), d$pdc)
        }
    }
    return(p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
