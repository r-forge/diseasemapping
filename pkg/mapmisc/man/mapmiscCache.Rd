\name{persistentCache}
\alias{persistentCache}
\alias{mapmiscCache}
\alias{mapmiscCachePath}
\title{Set a persistent cache}
\description{
Sets a cache folder in temporary space
  }
\usage{
	persistentCache(verbose=TRUE)
	}

\arguments{
\item{verbose}{print location of the cache folder}
}
\details{
	The default cache for map images is tempdir()/mapmiscCache, which will be deleted when an R session ends.  Running this function sets a cache in /tmp/mapmiscCache_[username], which will re-use cached data across R sessions.
}
\examples{
# current cache
getOption("mapmiscCachePath")

# set a new cache
myCache = file.path(tempdir(), 'myCache')
dir.create(myCache)
options(mapmiscCachePath = myCache)
getOption("mapmiscCachePath")

# create a persistent cache
\dontrun{
persistentCache(verbose=TRUE)
}
getOption("mapmiscCachePath")
}