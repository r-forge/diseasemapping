
cholGpu = function(x,D, 
	control = list(
		workgroupSize = NULL, 
		localWorkgroupSize = NULL, 
		localStorage = NULL, 
		sizeOfDouble = NULL), 
	verbose=FALSE) {

	file <- system.file("CL", "cholGpu.cl", package = "geostatsgpu")
	kernel <- readChar(file, file.info(file)$size)

	if(missing(D)) {
		D = vclVector(0, nrow(x),
			type='double',
			ctx_id = x@.context_index)
	}

	Scontrol = c(
		'workgroupSize', 'localWorkgroupSize', 
		'localStorage',
		'sizeOfDouble')

	missingControl = setdiff(Scontrol, 
		names(control))

	if(length(missingControl) | 
		any(unlist(lapply(control, is.null)))) {
	
		localInfo = gpuNlocal(
			kernel, 
			'cholOffDiag',
			x@.device_index)

		for(Dcontrol in c(
			'localWorkgroupSize', 
			'sizeOfDouble')) {

			if(!length(control[[Dcontrol]])) {
				control[[Dcontrol]] = localInfo[[Dcontrol]]
			}
		}

		if(!length(control$workgroupSize)) {
			control$workgroupSize = localInfo$maxWorkgroupSize
		}

		if(!length(control$localStorage)) {
			control$localStorage = pmax(
				control$localWorkgroupSize,
				floor(
				0.95*localInfo$localMemory / 
				control$sizeofDouble) ) 
		}
	}

	if(verbose){
		message(paste('global work items', 
			control$workgroupSize, 
			'local work items', 
			control$localWorkgroupSize, 
			'local storage', 
			control$localStorage, 
			"N", ncol(x)))
	}


	fromC = cpp_cholGpu(
		x@address, 
		D@address,
		control$workgroupSize, 
		control$localWorkgroupSize, 
		control$localStorage,
		x@.context_index - 1,
		kernel
		)

	res = list(L = x, D=D, logDet = fromC, 
		control = control)

	res
}