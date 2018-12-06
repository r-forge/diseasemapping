
cholGpu = function(x,D, 
	MCtotal, MClocal, 
	localStorage, sizeOfDouble, 
	verbose=FALSE) {

	file <- system.file("CL", "cholGpu.cl", package = "geostatsgpu")
	kernel <- readChar(file, file.info(file)$size)

	if(missing(D)) {
		D = vclVector(0, nrow(x),
			type='double',
			ctx_id = x@.context_index)
	}


	if(missing(MCtotal) | missing(MClocal) | 
		missing(sizeOfDouble) | missing(localStorage)
		) {
		localInfo = gpuNlocal(
			kernel, 
			'sumLog',
			x@.device_index)

		if(missing(MCtotal)) {
			MCtotal = localInfo$maxWorkgroupSize
		}
		if(missing(MClocal)) {
			MClocal = localInfo$localWorkgroupSize
		}
		if(missing(sizeOfDouble)) {
			sizeofDouble = localInfo$sizeOfDouble
		}
		if(missing(localStorage)) {
			localStorage = floor(
				localInfo$localMemory / 
				sizeofDouble)
		}

	}

	if(verbose){
		message(paste('global work items', MCtotal, 
			'local work items', MClocal, 
			'local storage', localStorage))
	}


	fromC = cpp_cholGpu(
		x@address, 
		D@address,
		MCtotal,
		MClocal, 
		localStorage,
		x@.context_index - 1,
		kernel
		)

	res = list(L = x, D=D, logDet = fromC, 
		extra = c(
			MCglobal = MCtotal, 
			MClocal = MClocal, 
			localStorage = localStorage))

	res
}