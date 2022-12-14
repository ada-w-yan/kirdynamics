% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general_functions.R
\name{solve_ODE_error_handling}
\alias{solve_ODE_error_handling}
\title{solve an ODE with error handling}
\usage{
solve_ODE_error_handling(mod, solving_time, tol = 1e-06)
}
\arguments{
\item{mod}{an ode_system object generated by odin}

\item{solving_time}{times at which to solve the ODE}

\item{tol}{numeric vector of length 1: tolerance}
}
\value{
a deSolve object: solution of the ODE
}
\description{
if an error occurs, solve with a smaller time interval
for solution elements except the time vector, set to tol if the solution is below tol
}
