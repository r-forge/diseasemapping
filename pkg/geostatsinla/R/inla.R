# used to have
# importFrom(INLA, inla, inla.models, inla.make.lincombs)

loadInla = function() {
	res = require("INLA", quietly=TRUE,warn.conflicts=FALSE)
	if(!res) {
		warning("see www.r-inla.org to install INLA. \n First install:\n splines, numDeriv, Rgraphviz, graph, fields,\n rgl, mvtnorm, multicore, pixmap,\n splancs, orthopolynom")
	}
	res
}

inla = function(...) {
	if(loadInla()) {
		INLA::inla(...)
	} else {
		return(NULL)
	}
}

inla.models = function(...) {
	if(loadInla()) {
		INLA::inla.models(...)
	} else {
		return(NULL)
	}
}

inla.make.lincombs = function(...) {
	if(loadInla()) {
		INLA::inla.make.lincombs(...)
	} else {
		return(NULL)
	}
}