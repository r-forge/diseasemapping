all: slides.pdf document.html document.docx document.odt
# add anotherGoat.html

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz my*.bib *.nav *.xml *.vrb *.snm 

empty:
	rm *.pdf *.html *.docx *.odt

slides.tex: slides.md
	pandoc --standalone --smart --biblatex --to=beamer --template=../inst/src/beamer.Pmisc --output=$@ $<

document.tex: document.md
	pandoc --standalone --smart --biblatex --output=$@ $<

%.bcf: %.tex
	xelatex $<

%.bbl: %.bcf
	biber $<	

%.pdf: %.tex %.bbl
	xelatex -interaction=nonstopmode $<;
	xelatex -interaction=nonstopmode $<;
	R --slave -e "tools::compactPDF('$@', gs_quality='ebook')"

%.html: %.md
	pandoc --standalone --filter=pandoc-citeproc --self-contained --output=$@ $<

%.docx: %.md
	pandoc --filter=pandoc-citeproc --reference-docx=../inst/src/template.docx --output=$@ $<
	
%.odt: %.md
	pandoc --filter=pandoc-citeproc --reference-odt=../inst/src/template.odt --output=$@ $<
