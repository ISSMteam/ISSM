#!/bin/bash
#get number of lines of code
cloc $ISSM_DIR/src $ISSM_DIR/m4 --exclude-dir=.svn --exclude-dir=ad  --exclude-ext=exp --exclude-lang=make --out=temp --force-lang="MATLAB",m
cat temp
./cloc2html.py
rm temp

cat $ISSM_DIR/src/dox/issm.dox | sed '/<table/,//d' > input1
cat $ISSM_DIR/src/dox/issm.dox | sed '1,/<\/table>/d' > input2
cat input1 temp.html input2 > $ISSM_DIR/src/dox/issm.dox
#rm input1 input2 temp.html
