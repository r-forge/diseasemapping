#' @export
markdownHeader = function(
		title=NULL, 
		author=NULL,
		date=Pmisc::today(),
		biblatex=1,
	 bibliotitle = 'References',
	 bibliostyle = 'authoryear',
		biblatexoptions = c(maxbibnames=20,
				maxcitenames=2,doi='false',
				isbn='false',	giveninits='true',backend=biber),
		startAndEnd='---',
		beamer=FALSE,
		...
) {
	
	result = list(
			startAndEnd,
			paste('title:', title),
			paste('author:', author),
			paste('date:', date)
	)
	
	dots = list(
			biblatex=biblatex,
			'biblio-title'=bibliotitle,
			'biblio-style'=bibliostyle,
			biblatexoptions = biblatexoptions, 
			...)	
	
	dots = dots[
			unlist(lapply(dots, length))>0
			]
	
	names(dots) = gsub(
			"header[[:punct:]]?includes",
			'header-includes', names(dots),
			ignore.case=TRUE)
	
	for(D in names(dots)) {
		if(length(dots[[D]])>1) {
			dots[[D]] = knitr::combine_words(c('',
							paste(
									'  - ', 
									names(dots[[D]]),
									c("",' = ')[1+as.logical(nchar(names(dots[[D]])))],
									unlist(dots[[D]]),
									sep='')
					), 
					sep='\n', and='')
		}
	}
	
	dotsNotHeaders = grep("^header-includes", names(dots), value=TRUE, invert=TRUE)
	
	for(D in dotsNotHeaders) {
		result[[length(result)+1]] = paste(
				D, ": ", dots[[D]], sep=''
		)	
	}
	
	# header includes
	if(beamer | any(names(dots)=='header-includes')) {
		toAdd <- gsub(
				'\n+', '\n', 
				knitr::combine_words(c(
								'header-includes: ', 
								dots[['header-includes']], 
								slideHeaderIncludes[beamer]
						), sep='\n', and='')
	)
	result[[length(result)+1]] = toAdd 
	}
	result[[length(result)+1]] = startAndEnd
	
	result = knitr::combine_words(c('',result,''), 
			sep='\n', and='')
	knitr::asis_output(result)
	
}

#' @export
today = function(...){
	dots = list(...)
	if(!any(names(dots)=='format'))
		dots$format='%A %d %B %Y'
	do.call(format, c(list(x=Sys.time()), dots))
}

slideHeaderIncludes = 
		knitr::asis_output(
				knitr::combine_words(c('',
								paste(
										'  - \\newcommand{\\',
										c(
												'newcol}[1][0.5]{\\column{#1\\textwidth}}',
												'bcol}[1][0.5]{\\begin{columns}\\column{#1\\textwidth}}',
												'ecol}{\\end{columns}}',
												'newcolb}[2][0.5]{\\end{block}\\column{#1\\textwidth}\\begin{block}{#2}}',
												'bcolb}[2][0.5]{\\begin{columns}\\column{#1\\textwidth}\\begin{block}{#2}}',
												'ecolb}{\\end{block}\\end{columns}}'),
										sep='')),
						sep='\n', and='')
		)


