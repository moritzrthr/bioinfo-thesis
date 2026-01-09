all:
	pdflatex main && biber main && pdflatex main && pdflatex main

latex:
	pdflatex main

clean:
	rm *.glo *.aux *.toc *.out *.log *.idx *.blg *.bbl *.pdf *.lof *.lot _region_.tex
