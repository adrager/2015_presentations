arne.pdf:  arne.tex
	pdflatex arne
	bibtex arne
	latex arne
	pdflatex arne
	pdflatex arne
	rm *.aux *.bbl *.blg *.lof *.log *.lot *.out
