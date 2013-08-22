TEX   = latexmk
FLAGS = -pvc

default: thesis.tex
	$(TEX) $(FLAGS)

clean:
	$(TEX) -C
	$(RM) *.bbl *.loa *.lol *.pyg *.run.xml *.synctex.gz
