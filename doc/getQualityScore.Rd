\name{getQualityScore}
\alias{getQualityScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute a quality score for the CNV clustering.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getQualityScore(x, hyp = "H0")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{hyp}{
%%     ~~Describe \code{hyp} here~~
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
function (x, hyp = "H0") 
{
    if (hyp == "H0") 
        posterior <- x@best.fit.H0
    if (hyp == "H1") 
        posterior <- x@best.fit.H1
    Qvec <- c()
    for (coh in levels(posterior$batch)) {
        posterior2 <- subset(posterior, posterior$batch == coh)
        sds <- tapply(posterior2$signal, FUN = sd, INDEX = posterior2$cn)
        means <- tapply(posterior2$signal, FUN = mean, INDEX = posterior2$cn)
        freq <- table(posterior2$cn)/length(posterior2$cn)
        sds <- subset(sds, table(posterior2$cn) > 4)
        means <- subset(means, table(posterior2$cn) > 4)
        freq <- subset(freq, table(posterior2$cn) > 4)
        l <- length(means)
        if (l == 1) {
            return(NA)
        }
        else {
            dmeans <- abs(means[1:(l - 1)] - means[2:l])
            av.sds <- (freq[1:(l - 1)] * sds[1:(l - 1)] + freq[2:l] * 
                sds[2:l])/(freq[1:(l - 1)] + freq[2:l])
            weights <- freq[1:(l - 1)] * freq[2:l]
            Q <- sum(weights * dmeans/av.sds)/sum(weights)
            Qvec <- append(Qvec, values = Q)
        }
    }
    return(Qvec)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
