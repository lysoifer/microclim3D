% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/climate_functions.R
\name{altcorrect_point}
\alias{altcorrect_point}
\title{Altcorrect_point
apply altitudinal correction to point climate data
code adapted from microclimf:::.runmodel2Cpp}
\usage{
altcorrect_point(climdata, dtmc_p, elev_p, altcor = 2)
}
\arguments{
\item{climdata}{climate dataframe}

\item{altcor}{numeric; 1 = fixed lapse rate, 2 = humidity-dependent lapse rate}

\item{dtmc}{elevation of input macroclimate data}
}
\value{
altitduinally corrected climate data (adjusts temp and pressure) for input into
micropoint::runmicropoint
}
\description{
Altcorrect_point
apply altitudinal correction to point climate data
code adapted from microclimf:::.runmodel2Cpp
}
