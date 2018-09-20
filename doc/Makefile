# Build the manual, then copy it to the current directory.
manual:
	cd src/user_guide ; \
		pdflatex Rayleigh_User_Guide.tex ; \
		pdflatex Rayleigh_User_Guide.tex
	mv src/user_guide/Rayleigh_User_Guide.pdf .

doxygen:
	doxygen options.dox

.PHONY: manual doxygen

