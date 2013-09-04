TEX   = latexmk -quiet
FLAGS = -pvc

default: thesis.tex
	TEXINPUTS="../shared:$(TEXINPUTS)" $(TEX) $(FLAGS)

clean:
	$(TEX) -C
	$(RM) *.bbl *.loa *.lol *.pyg *.run.xml *.synctex.gz
