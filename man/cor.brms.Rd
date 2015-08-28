\name{cor_brms}
\alias{cor_brms}
\alias{cor_brms-class}

\title{
Correlation structure classes
}

\description{
Classes of correlation structures available in the \code{brms} package.
}

\value{
Avaiblable correlation structures
\describe{
  \item{cor_ar}{autoregressive process of arbitrary order.}
  \item{cor_arma}{autoregressive-moving average process, with arbitrary orders for the 
                  autoregressive and moving average components.}
  \item{cor_ma}{moving average process of arbitrary order.}                
}
}

\seealso{
\code{\link[brms:cor_ar]{cor_ar}, \link[brms:cor_arma]{cor_arma}, \link[brms:cor_ma]{cor_ma}}
}