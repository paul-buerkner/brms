\name{cor.brms}
\alias{cor.brms}
\alias{cor.brms-class}

\title{
Correlation structure classes
}

\description{
Classes of correlation structures available in the \code{brms} package.
}

\value{
Avaiblable correlation structures
\describe{
  \item{cor.ar}{autoregressive process of arbitrary order.}
  \item{cor.arma}{autoregressive-moving average process, with arbitrary orders for the 
                  autoregressive and moving average components.}
  \item{cor.ma}{moving average process of arbitrary order.}                
}
}

\seealso{
\code{\link[brms:cor.ar]{cor.ar}, \link[brms:cor.arma]{cor.arma}, \link[brms:cor.ma]{cor.ma}}
}