#!/bin/bash

# Translates trivial ISSM tests from MATLAB to JavaScript

for INFILE in "$@"; do
    if ! [[ $INFILE == *.m ]]; then
        echo "Error: Invalid input file detected. Only MATLAB files are accepted as input."
        echo "Instead, found $INFILE."
        echo "Exiting..."
        exit
    fi

    OUTFILE=$(dirname $INFILE)/$(basename $INFILE .m).js
    cp $INFILE $OUTFILE

    printf "Translating to $OUTFILE..."
    TRIANGLEMODELPATH=$(grep '^.*triangle.*' $OUTFILE | awk -F "," '{print $(NF-1)}' | cut -d ')' -f 1)

    if [[ -n $TRIANGLEMODELPATH ]]; then
        TRIANGLEMODEL=$(basename $TRIANGLEMODELPATH | sed 's/\..*//' | tr A-Z a-z)
        TRIANGLEVAL=$(grep '^.*triangle.*' $OUTFILE | awk -F "," '{print $NF}' | cut -d ')' -f 1)
        sed -i "s/\s*triangle.*/triangle(md,$TRIANGLEMODEL[0],$TRIANGLEVAL);/" $OUTFILE # triangle statement
    fi


    sed -i 's/.*setmask\(.*\)/setmask\1/' $OUTFILE # setmask

    sed -i 's:%://:g' $OUTFILE # comments
    sed -i 's/{/[/g' $OUTFILE # left curly brace to bracket
    sed -i 's/}/]/g' $OUTFILE # right curly brace to bracket


    sed -i 's/.*parameterize.*/parameterize(md);/' $OUTFILE # parameterize
    sed -i 's/^.*\(extrude\)\(.*\)/md\.\1\2/' $OUTFILE # extrude

    sed -i '/.*solve.*/! s/^md\s*=\s*//' $OUTFILE # remove md= unless it calls solve

    sed -i 's:^md\.cluster=generic.*://&:' $OUTFILE # comment out oshostname

    sed -i 'N;s/\.\.\.//g' $OUTFILE # JavaScript lets you have multiline statments without explicit indicators

    sed -i 's/md\.results\..*Solution/&[0]/' $OUTFILE # actual results are stored in the 0th index of the solution

    sed -i 's/\(.*\)(:)\s*=\s*\([^;]*\);$/for (var i = 0; i < \1.length; ++i) {\n\t\1[i] = \2;\n}/' $OUTFILE

    sed -i 's/\(.*\)(\([0-9]*\):\([0-9]*\))\s*=\s*\([^;]*\)/for (var i = \2; i < \3; ++i) {\n\t\1[i] = \4\n}/' $OUTFILE

    sed -i 's/\(.*\)\s*=\s*find(\(.*\))/\1=[];\nfor (var i = 0; i < \2.length; ++i) {\n\tif (\2[i] !== 0)\n\t\t\1.push(i);\n}/' $OUTFILE # matlab's find creates an array of the indices with nonzero elements

    sed -i 's/[^%]*if\s\(.*\)$/if (\1) {/' $OUTFILE #if statement
    sed -i 's/else/} & {/' $OUTFILE #else

    sed -i 's/end/}/' $OUTFILE #closing block

    if grep 'zeros(.*)' $OUTFILE; then
        sed -i '2s/^/function zeros(...args) {\n\tvar array = [];\n\tfor (var i = 0; i < args[0]; ++i) {\n\t\tarray.push(args.length == 1 ? 0 : zeros(args.slice(1)));\n\t}\n\treturn array;\n}\nvar md = new model();\n/' $OUTFILE
    fi # include zeros function to generate matrices of zeros


    #sed -i "s/\([^']\)NaN\([^']\)/\1null\2/g" $OUTFILE # NaN translates to null

    sed -i 's/function\s\(.*\)\s*=\s*\(.*\)(\(.*\)/this.\2 = function(\3)/' $OUTFILE

    sed -i 's/(self,*\(.*\))/(\1)/' $OUTFILE

    sed -i 's/classdef\s*\(.*\)/function \1() {/' $OUTFILE # Classes not in javascript - implemented as functions

    sed -i 's/properties (\(.*\))/\/\/properties (\1)/' $OUTFILE #properties block
    sed -i 's/^\(\s*\)methods.*/\1\/\/&/' $OUTFILE #methods block declaraion

    sed -i 's/md\s*=\s*checkfield\(.*\)/checkfield\1/' $OUTFILE #checkfield

    sed -i 's/\(\s*\)error(/\1console.error(/' $OUTFILE #error -> console.error

    sed -i 's/switch \(.*\)/switch(\1) {/' $OUTFILE # switch statement
    sed -i 's/case \(.*\)/&:/' $OUTFILE # case x
    sed -i '$!N;/switch/!s/\n\s*case .*$/break\;&/;P;D' $OUTFILE
    sed -i 's/otherwise/default:/' $OUTFILE # default case
    sed -i '$!N;s/\n\s*default:\s*$/break\;&/;P;D' $OUTFILE # add break before default case
    sed -i 'N;s/\(\s*\)\(.*\)break\;/\1\2\n\1break\;/' $OUTFILE # make break its own line

    perl -i -pe 's/(varargin\[)(\d+)/$1.($2-1)/e' $OUTFILE # matlab starts its indices at 1, javascript at 0

    sed -i 's/self/this/g' $OUTFILE # self -> this
    sed -i 's/md\.transient/md.trans/g' $OUTFILE # self -> this

    sed -i 's/function(\([^ ]*\))/function(\1) {/' $OUTFILE

    #sed -i 's/\(\d+\)\^\(\d+\)/Math.pow(\1,\2)/' $OUTFILE
    sed -i 's/\([0-9]*\)\^\([0-9]*\)/Math.pow(\1,\2)/g' $OUTFILE

    #if [[ $OUTFILE == test*.js ]]; then
    sed -i '2s/^/var md = new model();\n/' $OUTFILE # initialize the model
    #fi

    #sed -i 'N;s/properties.*\n\s*\(.*\)\(\s*\)=\s*\(.*\)/this/' $OUTFILE

    #zeros - make it work for all args
    #power function (10+5^-10)
    printf "completed.\n"
done
