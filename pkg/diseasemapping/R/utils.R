# used to have
# importFrom(INLA, inla, inla.models, inla.make.lincombs)

loadInla = function() {
	res = "INLA" %in% rownames(installed.packages())
	if(res) 
		require("INLA", quietly=TRUE,warn.conflicts=FALSE)
	res
}

inla = function(...) {
	if(loadInla()) {
		INLA::inla(...)
	} else {
		return(list(logfile="INLA is not installed. \n install splines, numDeriv, Rgraphviz, graph,\n fields, rgl, mvtnorm, multicore, pixmap,\n splancs, orthopolynom \n then see www.r-inla.org"))
	}
}

inla.models = function(...) {
	if(loadInla()) {
		INLA::inla.models(...)
	} else {
		return(NULL)
	}
}

#inla.model.properties = function(...) {
#	if(loadInla()) {
#		INLA::inla.model.properties(...)
#	} else {
#		return(NULL)
#	}
#}


inla.make.lincombs = function(...) {
	if(loadInla()) {
		INLA::inla.make.lincombs(...)
	} else {
		return(list())
	}
}
inla.make.lincomb = function(...) {
	if(loadInla()) {
		INLA::inla.make.lincomb(...)
	} else {
		return(list(lc=list()))
	}
}


