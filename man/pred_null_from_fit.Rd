% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_linear_stan.R
\name{pred_null_from_fit}
\alias{pred_null_from_fit}
\title{predict all saliva, neutrophil and lymphocyte data for multiple participants}
\usage{
pred_null_from_fit(
  fit,
  saliva_model,
  gm_model,
  lymphocyte_model,
  N_samples = 100
)
}
\arguments{
\item{fit}{object returned by fit_saliva_gm_stan: list with two elements:
fit: stanfit object
data: data list used to fit stanfit object}

\item{saliva_model}{either a stan model object used to predict saliva curves,
or a string specifying model name to compile this model}

\item{gm_model}{either a stan model object used to predict gm curves,
or a string specifying model name to compile this model}

\item{lymphocyte_model}{either a stan model object used to predict lymphocyte curves,
or a string specifying model name to compile this model}

\item{N_samples}{integer.  number of samples to use in model predictions}
}
\value{
list of three objects:
time: time series for which to predict
N_predictions: list of length N, where each element is for one participants.
Each element is in turn a list of two elements:
Uhat: matrix of model predictions.  Each row is a time and each column is a sample.
Lbhat: matrix of model predictions.  Each row is a time and each column is a sample.
C_predictions: list of length C, where each element is for one cell population.
Each element is a matrix of model predictions.  Each row is a time and each column is a sample.
}
\description{
predict saliva and gm curves from model fit,
when saliva data and gm data for one cell population are fitted
together for a single participant
}
