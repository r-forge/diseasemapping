all: slides.pdf document.html document.docx

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz my*.bib *.nav *.xml *.vrb *.snm 

%.md: %.Rmd
	R -e "knitr::knit('$<', '$@', encoding='UTF-8')" 

%.tex: %.md
	pandoc --standalone --smart --biblatex --to=beamer --template=../inst/src/beamer.Pmisc --output=$@ $<

%.bcf: %.tex
	xelatex $<

%.bbl: %.bcf
	biber $<	

%.pdf: %.tex %.bbl
	xelatex -interaction=nonstopmode $<;
	xelatex -interaction=nonstopmode $<

%.html: %.md
	pandoc --standalone --filter=pandoc-citeproc --self-contained --output=$@ $<

%.docx: %.md
	pandoc --filter=pandoc-citeproc --reference-docx=../inst/src/template.docx --output=$@ $<