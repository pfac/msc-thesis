ERB   = erb -P
TEX   = latexmk
FLAGS = -pvc

default: main.tex
	$(TEX) $(FLAGS)

%: en/%.tex
	$(TEX) $(FLAGS) tex/$@.tex

clean:
	$(TEX) -C