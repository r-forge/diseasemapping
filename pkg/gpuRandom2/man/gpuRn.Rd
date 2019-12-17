\name{gpuRn}
\alias{gpuRn}

\title{Generate Uniform or Normal Random Numbers}

\description{Give streams, work group sizes and vector length, generate a vector of Uniform or Normal random numbers on GPU.}


\usage{gpuRn(x, streams, workgroupSize, random_type=c("uniform","normal"))}

\arguments{
\item{x}{A vector of desired length for storing generated random numbers.}
\item{streams}{Used to generate random numbers on GPU by work items.}
\item{workgroupSize}{Global Sizes of device domain.}
\item{random_type}{Either Uniform random numbers or Normal random numbers.}
\item{verbose}{A logical indicating if global and local dimensions should be showed. }
}

\value{A vector containing the desired number of Normal or Uniform random numbers.}

\examples{
## Create some streams
Nglobal = c(4,4)
streams1 = gpuRandom:::cpp_mrg31k3pCreateStreams(prod(Nglobal));
t(streams1)

## Generate eight random uniforms
x = vclVector(length=8L, type=theType)
streams2<-gpuRandom::gpuRn(x, streams=streams1, Nglobal,random_type = "uniform")
## new streams
t(streams2)[1:4,]
as.vector(x)

## Generate a large number of Uniform random numbers in float precison
Nglobal = c(64,4)
streams = gpuRandom:::cpp_mrg31k3pCreateStreams(prod(Nglobal))
if(gpuInfo()$maxAllocatableMem > 1e10) {
  NperItem = 5000
} else {
  NperItem = 500
}

x = vclVector(length=as.integer(prod(Nglobal)*NperItem), type="float")
length(x)

gpuRandom::gpuRn(x,  streams=streams, Nglobal, random_type = "uniform"))

xVec = as.vector(x)
hist(xVec, xlim=c(0,1),breaks=40)


## Generate a large number of Normal random numbers in double precision
Nglobal = c(64,4)
streams = gpuRandom:::cpp_mrg31k3pCreateStreams(prod(Nglobal))

if(gpuInfo()$maxAllocatableMem > 1e10) {NperItem = 5e5}else {NperItem = 500}

z = vclVector(0, length=as.integer(prod(Nglobal)*NperItem), type="double")

## Number of generated random numbers
length(z)

gpuRandom::gpuRn(z, streams=streams, workgroupSize = Nglobal, random_type="normal")

zVec = as.vector(z)
hist(zVec,breaks=40)
}
