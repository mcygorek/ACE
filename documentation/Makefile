HTLATEX=$(shell which htlatex)
$(info HTLATEX: $(HTLATEX))

ifeq ($(HTLATEX),)
  $(warning 'htlatex' required to generate html documentation!)
endif


all: documentation.tex
	pdflatex documentation.tex
	htlatex documentation.tex


.PHONY: clean

clean:
	rm documentation.{html,aux,idv,css,lg,log,tmp,4tc,4ct} documentation*x.png 2>/dev/null || true
