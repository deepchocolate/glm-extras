% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2wrap.R
\name{R2Wrap}
\alias{R2Wrap}
\title{Create a model container from a fitted glm model}
\usage{
R2Wrap(mod, comparison = ~0, cases = NA_integer_, controls = NA_integer_)
}
\arguments{
\item{mod}{The fitted model object}

\item{comparison}{A formula to compare the model \code{mod} with, e.g. ~-exposure to remove \code{exposure} from \code{mod}}

\item{cases}{N of \code{y=1}}

\item{controls}{N of \code{y=0}}
}
\description{
if passing a formula to argument \code{comparison}, R2.scaleLiability will provide
the incremental R2, i.e. the difference between model \code{mod} and \code{mod} with
the updated formula.
}
\details{
The parameters \code{cases} and \code{controls} is relevant when analyzing a binary
phenotype with a linear model (\code{glm(..., family="gaussian")}), in which
\code{p*(1-p)} is used as a denominator for \code{R2.scaleObserved} instead of
\code{var(y)}. For binomial models these arguments can be ignored.
}
