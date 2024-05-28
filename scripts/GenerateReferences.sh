#!/bin/bash
#Generate a RTF file with references provided in CITE
#Use AGU bibliography style and delete labels for easy copy/paste
#Version: 04/28/11

#Genetare a pdf with the list of references provided in $CITE
CITE="Joughin2009,Rignot2002,Frey2001"

#First erase files; 
rm -rf references.rtf

#create Latex file
cat <<END > references.tex
\documentclass{article}
\bibliographystyle{agu}
\begin{document}
\nocite{$CITE}
\bibliography{$JPL_SVN/publications/bibtex/references}
\end{document}
END
#Generate pdf
echo "Compiling document"
latex  -interaction=batchmode -file-line-error references.tex > /dev/null
echo "Running Bibtex"
bibtex references
echo "Removing labels"
cat references.bbl | sed -e "s-ibitem.*\]-DELETE-" | grep -v DELETE > references.bak1 #delete labels that hold in one line
cat references.bak1 | sed -e "/ibitem/,/]/d" > references.bak2                        #delete labels that are on multiple lines
mv references.bak2 references.bbl
echo "Converting to rtf"
$ISSM_DIR/externalpackages/latex2rtf/install/latex2rtf -P $ISSM_DIR/externalpackages/latex2rtf/install/cfg/ references

#Remove all but rtf file
rm -rf references.[!rtf]*
rm -rf references.tex

#open output
#open references.rtf
