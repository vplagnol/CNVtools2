\name{CNVtools-class}
\Rdversion{1.1}
\docType{class}
\alias{CNVtools-class}
\alias{apply.pca,CNVtools-method}
\alias{expand,CNVtools-method}
\alias{fit,CNVtools-method}
\alias{initialize,CNVtools-method}
\alias{plot.cnv,CNVtools-method}
\alias{show,CNVtools-method}
\alias{TDT,CNVtools-method}

\title{Class \code{"CNVtools"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CNVtools", signal, trait, batch, model.association, model.mean, model.var, IID, FID, f.member, covariates)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{signal}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{trait}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{batch}:}{Object of class \code{"factor"} ~~ }
    \item{\code{ncomp}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{IID}:}{Object of class \code{"character"} ~~ }
    \item{\code{FID}:}{Object of class \code{"factor"} ~~ }
    \item{\code{f.member}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{covariates}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{expanded.covariates}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{frame.fitting}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{multiprobe.signal}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{nprobes}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nsamples}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cn.calls}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{model.association}:}{Object of class \code{"formula"} ~~ }
    \item{\code{model.mean}:}{Object of class \code{"formula"} ~~ }
    \item{\code{model.var}:}{Object of class \code{"formula"} ~~ }
    \item{\code{design.matrix.association}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{design.matrix.mean}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{design.matrix.var}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{strata.mean}:}{Object of class \code{"integer"} ~~ }
    \item{\code{strata.var}:}{Object of class \code{"integer"} ~~ }
    \item{\code{family}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{parent.offspring}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{binary.trait}:}{Object of class \code{"logical"} ~~ }
    \item{\code{best.fit.H0}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{mle.H0}:}{Object of class \code{"list"} ~~ }
    \item{\code{best.fit.H1}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{mle.H1}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{apply.pca}{\code{signature(.Object = "CNVtools")}: ... }
    \item{expand}{\code{signature(object = "CNVtools")}: ... }
    \item{fit}{\code{signature(object = "CNVtools")}: ... }
    \item{initialize}{\code{signature(.Object = "CNVtools")}: ... }
    \item{plot.cnv}{\code{signature(.Object = "CNVtools")}: ... }
    \item{show}{\code{signature(object = "CNVtools")}: ... }
    \item{TDT}{\code{signature(.Object = "CNVtools")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("CNVtools")
}
\keyword{classes}
