\name{ranefPlot}
\alias{ranefPlot}
\title{Dotplot of random effects from glmmTMB}
\usage{
ranefPlot(
  x,
  component = c("cond", "zi"),
  grpvar = 1,
  term = "(Intercept)",
  cutSe = Inf,
  maxNames = 40,
  col = "Dark2",
  xlab = "x",
  level = 0.95,
  ...
)
}
\arguments{
\item{x}{model output from glmmTMB}

\item{component}{model component}

\item{grpvar}{name of random effect}

\item{cutSe}{only plot levels where conditional SE is above this threshold}

\item{maxNames}{maximum number of labels on the plot}

\item{col}{vector of colours or name of a colorBrewer pallet}

\item{xlab}{x label}
\item{term}{variable to plot}
\item{level}{level of prediction interval}
\item{...}{Additional arguments for plot}
}
\value{
nothing
}
\description{
Creates a nice dotplot with CI's
}
