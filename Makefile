cmsdas_2014_CERN.pdf:  cmsdas_2014_CERN.tex
	pdflatex cmsdas_2014_CERN
	bibtex cmsdas_2014_CERN
	latex cmsdas_2014_CERN
	pdflatex cmsdas_2014_CERN
	pdflatex cmsdas_2014_CERN
	rm *.aux *.bbl *.blg *.lof *.log *.lot *.out
