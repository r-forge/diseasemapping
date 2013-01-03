getTiles = function(xlim, ...) {
	UseMethod("getTiles")
	
}

getTiles.matrix = function(xlim, ...){
	if(!all(dim(xlim)==2)){
		warning("wrong dimensions, is xlim a bounding box?")
	}
	webmaps::getTiles(xlim[1,], xlim[2,], ...)
}

getTiles.Extent = function(xlim, ...){
	webmaps::getTiles(
			c(xlim@xmin, xlim@xmax), 
			c(xlim@ymin, xlim@ymax), 
			...)
	
}

getTiles.default = function(xlim, ...){
	webmaps::getTiles(xlim, ...)			
}