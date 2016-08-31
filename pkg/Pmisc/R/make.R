#' @export
make = function(x, 
		suffix=NULL, beamer=FALSE,
		run=FALSE, ...) {

	if(!is.null(suffix))
		x = paste(
				tools::file_path_sans_ext(basename(x)),
				suffix, sep='.')
	
theString= paste(' -f', 
		system.file('src','knitrMakefile', package='Pmisc'),
		x)

if(beamer)
	theString = paste(theString, 'pandocTo=beamer')

dots = list(...)
for(D in seq(from=1, by=1, len=length(dots)))
	theString = paste(theString, D)

if(run) {
	system(paste('make ', theString))
} else {
	cat(theString)
}

invisible(theString)

}

#' @export
Makefile = function(x, suffix = NULL,
		beamer = FALSE, output='Makefile') {

if(!is.null(suffix))
		x = paste(
				tools::file_path_sans_ext(basename(x)),
				suffix, sep='.')
	
	
before = list(
	pandocTo = c('beamer', 'latex')[2-beamer],
	docxTemplate = system.file('src/template.docx', package='Pmisc'),
 odtTemplate = system.file('src/template.odt', package='Pmisc')
)

before = paste(paste(names(before), before, sep='='), collapse='\n')
before = paste(before, '\n\nall: ', x, '\n\n')

result = paste(
		before,
'.PRECIOUS: %.md %.tex

%.md: %.Rmd
	R -e "knitr::knit(\'$<\', encoding=\'UTF-8\')" $(Rargs)

%.tex: %.md
	pandoc --standalone --smart --biblatex $(pandocArgs) --to=$(pandocTo) --output=$@ $<

%.bcf: %.tex
	xelatex $<

%.bbl: %.bcf
	biber $<	

%.pdf: %.tex %.bbl
	xelatex -interaction=nonstopmode $<;
	xelatex -interaction=nonstopmode $<

%.html: %.md
	pandoc --smart --standalone --mathjax --filter=pandoc-citeproc $(pandocArgs) --to=html5 --output=$@ $<

%.rtf: %.md
	pandoc --smart --standalone --filter=pandoc-citeproc $(pandocArgs) --output=$@ $<

%.odt: %.md
	pandoc --smart --standalone --filter=pandoc-citeproc $(pandocArgs) --reference-docx=$(odtTemplate) --output=$@ $<

%.docx: %.md
	pandoc --smart --standalone --filter=pandoc-citeproc $(pandocArgs) --reference-docx=$(docxTemplate) --output=$@ $<

clean:
	rm *.run.xml* *.blg *.out *.log *.aux *.bcf *.bbl *.nav *.toc *.vrb'
)	
	
	
if(is.function(output)) {
		output(result)
	}
	if(is.character(output)) {
		cat(result, file=output)
	}
	invisible(result)
}