\name{mra}
\alias{mra}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Conducting Mendelian Randomization Analysis
}
\description{
\code{mra} is used to estimate causal effect of a quantitative exposure on a binary outcome in a Mendelian randomization analysis, of which outcome data is collected from a case-control study.
}
\usage{
mra(oformula, odata, eformula, edata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oformula}{
an object of class "\code{\link[stats]{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted on the case-control dataset. More details of model specification are illustrated in 'Details' and 'Examples'.
}
  \item{odata}{
a data frame containing variables specified in \code{oformula}, including outcome (case-control status), instruments, and covariates (if any).
}
  \item{eformula}{
an object of class "\code{\link[stats]{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted on the exposure dataset. More details of model specification are illustrated in 'Details' and 'Examples'.
}
  \item{edata}{
a data frame containing variables specified in \code{eformula}, including exposure, outcome, instruments, and covariates (if any).
}
}
\details{

\bold{\code{oformula}} specifies the model used to fit case-control data (outcome data), including case-control status, instruments, and covariates (if any). \code{mra} relies on a feature supported by the \code{Formula} package, so that the right hand side of \code{oformula} is separated into two parts by \code{|}. In general, format for \code{oformula} is \code{case-control status ~ covariates | instuments}. For example, \code{y ~ x1 + x2 | g1 + g2} specifies an outcome model for binary outcome \code{y}, with two covariates \code{x1} and \code{x2}, and two instruments \code{g1} and \code{g2} fitted in the model. One can use \code{y ~ 1 | g1 + g2} if no covariate to be adjusted. An intercept is always required in \code{oformula}. %Common operators for R formula are supported when specifying covariates. For example, one can use \code{*} to specified interaction. \code{mra} also supports arithmetic expressions (e.g. \code{log}) and function \code{\link[base]{I}}.
\code{mra} will convert character or factor variables into dummy variable, so one does not have to create dummy variables by their own unless they want to specify the baseline.

\bold{\code{eformula}} specifies the model used to fit exposure data, including exposure variable, outcome status (same variable as the case-control status), instruments, and covariates (if any). The right hand side of \code{eformula} is similar to that of \code{oformula}. The left hand side of \code{eformula} also has two parts separated by \code{|}. In general, format for \code{eformula} is \code{exposure | outcome ~ x1 | g1 + g2}. \code{mra} requires to know outcome status for every sample in exposure data.


\bold{Note 1: } Different covariates could be adjusted in \code{oformula} and \code{eformula}, which is quite common in practice as case-control data and exposure data may be collected for different research purpose, therefore, different covariates are measured. The instruments specified in \code{oformula} and \code{eformula} must be the same.

\bold{Note 2: } The case-control data and exposure data may share some subjects. This happens when researcher picks some subjects from the case-control study to measured their exposure based on their criteria, and has another set of exposure data from other sources. As such, a subject can appear in both datasets. Both \code{odata} and \code{edata} should always have a column named as \code{id}, so that \code{mra} can account for the variation due to data overlapping. This column is needed even if your case-control and exposure datasets do not share any subject.

%jkl

}
\value{
\code{mra} returns an object of class "mra".

The function \code{summary} is used to display a summary of the results. Many generic accessor functions are supported in \code{mra} to extract useful information of the value returned by \code{mra}. See 'Note' for more details.

An object of class "mra" is a list containing the following components:

\item{coefficients }{a named vector of coefficients. \code{bet} is the causal effect. \code{alp.} and \code{phi.} are coefficients estimated for the instruments and covariates in the exposure model. \code{a} and \code{gam.} are coefficients estimated for the intercept and covariates in the outcome model. \code{alp0} and \code{c0} are the intercept and estimated variance of random error in the exposure model for subjects without conditions. If there are subjects with conditions in the exposure data, \code{alp1} and \code{c1} are estimated as the intercept and variance of random error in the exposure model for those subjects. Refer to the paper for more details of model parameterization used in \code{mra}. }
\item{residuals }{the residuals for subjects without conditions in exposure data, that is exposure minus fitted values. Covariates (if any) are also adjusted. }
\item{fitted.values }{the fitted mean values for subjects without conditions in exposure data. Covariates (if any) are also adjusted. }
\item{wald }{Wald test. }
\item{lm }{Lagrange multiplier test. Recommanded for confidence interval and hypothesis testing. }
\item{ct }{test for presence of confounders. Available when exposure is also measured for some subjects with conditions. }
\item{tsr }{generalized two-stage regression method. This method is deprecated as it generates a more biased estimate for causal effect with underestimated standard error, a too-narrow confidence interval, and an underpowered test. }
\item{vcov }{variance-covariance matrix of \code{coefficients}. Using the generic function \code{\link[stats]{vcov}} to access this component is recommanded. }
\item{sigma2 }{the estimated variance of the random errors in the exposure model. \code{c0} for subjects without conditions, \code{c1} for subjects with conditions (if any). }
\item{call }{the matched call. }

}
\references{
Zhang, H., Qin, J., Berndt, S.I., Albanes, D., Gail, M.H., Yu, K. (2019) On Mendelian Randomization in Case-Control Studies. Under review.
}
\author{
Han Zhang
}
\note{
The generic functions \code{\link[base]{summary}}, \code{\link[stats]{coef}}, \code{\link[stats]{residuals}}, \code{\link[stats]{fitted}}, \code{\link[stats]{vcov}} are supported.

The generic function \code{\link[graphics]{plot}} is used for model diagnosis to the Assumption 2 in paper: in controls, random error of the exposure model is independent of instruments conditioned on confounders. This is done by examining the "residuals versus fits plot" of residuals (\code{\link[stats]{residuals}}) and fitted values (\code{\link[stats]{fitted}}) on subjects without conditions in the exposure data. A presence of heteroscedastic among those subjects suggests potential violation of this assumption, and caution is needed when interpreting the result.

The generic function \code{\link[stats]{confint}} can be used to compute Wald's confidence intervals for one or more parameters in fitted models (case-control and exposure models). This works fine for all parameters except for the causal effect, as we showed in paper that Wald's confidence interval for causal effect can have coverage probability much lower than its nominal level. We recommand to use the confidence interval derived from the Lagrange multiplier test, which is also more powerful than the Wald test.

}

\examples{

## This example estimates parameters in the
## following underlying models:
## 1. outcome model. A logistic regression model
##    d ~ z + x, of which the coefficient of
##    exposure z is the causal effect of interest;
## 2. exposure model. A quasi-likelihood model
##    z ~ g + x, of which g are used as instruments.
## In Mendelian randomization, those parameters
## could be estimated by fitting two working models
## with special parameterization:
## a. A logistic regression model d ~ g + x
## b. A quasi-likelihood model z ~ d + g + x

data(edata)
data(odata)

fit <- mra(d ~ x1 + x2 | g1 + g2 + g3,
           odata,
           z | d ~ x2 + x3 | g1 + g2 + g3,
           edata)

## summary tables for outcome model and exposure model
## and for testing the presence of confounder (if available)
summary(fit)

## causal effect estimate and its standard error
coef(fit)['bet']
sqrt(vcov(fit)['bet', 'bet'])

## Lagrange multiplier test
fit$lm

## model diagnosis
plot(fit)


}
