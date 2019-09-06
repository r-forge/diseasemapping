\name{fishertestGpu}
\alias{fishertestGpu}

\title{Performs Fisher's exact test}

\description{
Provide marginal row total and marginal column total of a contingency table, streams, work group size and a vector, performs Fisher's exact test using Monte Carlo simulation on GPU, returns the two-sided p-value.
}
\usage{fishertestGpu(rowsum, colsum, results, streams, globalSize, localSize)}

\arguments{
\item{rowsum}{A vector of marginal row total.}
\item{colsum}{A vector of marginal column total.}
\item{results}{A vector of test statistics, whose length is the number of replicates used in the Monte Carlo simulation. }
\item{streams}{Used to generate random numbers on GPU by work items.}
\item{globalSize}{Number of work items used in the device domain.}
\item{localSize}{Number of work items in a work group.}
}

\value{The p-value of Fisher's exact test using Monte Carlo simulation.}

\examples{
## Fisher's test for a 3 by 3 contingency table
TeaTasting <-matrix(c(21,17,5,19,38,6,26,12,9),nrow=3)
TeaTasting

## Using the R function 
fisher.test(TeaTasting, alternative = "t",  simulate.p.value = TRUE, B = 1e5)


## Using the GPU function
## Set global size
Nglobal = c(8,64)

## Create streams
stream_num=prod(Nglobal)
streams = clRNG:::cpp_mrg31k3pCreateStreams(stream_num)

NperItem = 30
results = vclVector(length=as.integer(stream_num*NperItem), type="float")
## Number of replicates in simulation
length(results)   

rowsum = as.integer(rowSums(TeaTasting, na.rm = FALSE))
colsum = as.integer(colSums(TeaTasting, na.rm = FALSE))
rowsum<-as.vclVector(rowsum)
colsum<-as.vclVector(colsum)
 
GpuRandom::fishertestGpu(rowsum, colsum, results,  streams, globalSize = Nglobal)

## calculate the p-value
STATISTIC <- -sum(lfactorial(TeaTasting))
almost.1 <- 1 + 64 * .Machine$double.eps
PVAL <- (1 + sum(   as.vector(results) <= STATISTIC/almost.1   )) / (length(results) + 1)
PVAL
}
