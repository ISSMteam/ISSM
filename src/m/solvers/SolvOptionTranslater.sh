PythonFiles=($(find . -name "*.py"))
#PythonFiles=(~/Model/issm_dev/src/m/classes/clusters/generic.py)
#PythonFiles=($(find . -name $1))
operators=(- + \\* \\/ \< \>)
multop=(+= -= \\*= == \>= \<= !=)
compop=(\< \> == \>= \<= !=)
keywords=(axis delimiter fmt shell mode dtype zlib)

#needs:
#   -negative number don't need spaces
#   -equals in options don't need spaces

for PythonName in "${PythonFiles[@]}"; do
    startword=($(awk '{print $1; exit}' "$PythonName" ))

    if [[ "$startword" == "function" ]]; then
	echo treating "$PythonName"
	cp "$PythonName" "$PythonName".bkp

	#first add the imports
	sed -i '1s/^/from collections import OrderedDict\n\n\n/' "$PythonName"
	sed -i '1s/^/from pairoptions import pairoptions\n/' "$PythonName"
	# and the return
	echo "return solverOptions" >> "$PythonName"

	#define the function
	sed -i "s/function solverOptions=\([a-z]\+\)(varargin)/def \1(\*args):/g" "$PythonName"

	#brutal indentation from line 6, but that should work here
	sed -i '6~1s/^/    /' "$PythonName"


	#first deal with equal signs
	#-should e spaced one except when ==, += , -= and if part of options
	sed -i "s/\([][:alnum:]]\+\)=\([[:alnum:]'\"]\+\)/\1 = \2/g"  "$PythonName"
	sed -i  's/ \+= \+/ = /g' "$PythonName"   #reduce number of spaces to one before and after =

	#replace varargin
	sed -i 's/varargin{:}/*args/g' "$PythonName"
	sed -i 's/varargin/*args/g' "$PythonName"

	#change getfieldvalue
	sed -i 's/getfieldvalue(options, /options.getfieldvalue(/g' "$PythonName"

	#replace struct
	sed -i 's/struct()/OrderedDict()/g' "$PythonName"

	#change to dict format
	sed -i "s/\(solverOptions\).\([a-z]\+[_[a-z]\+]*\)/\1['\2']/g" "$PythonName"

	#shift commment command from matlab
	sed -i 's/%/#/g' "$PythonName"

	#Add spaces after coma and limit to one space
	sed -i  's/,/, /g' "$PythonName"
	sed -i  's/, \+/, /g' "$PythonName"

	#two space for inline comments
	sed -i  's/\([[:alnum:]]\) \+#/\1  #/g' "$PythonName"

	#fix equals continuing with brackets
	sed -i 's/=\[/ = \[/g' "$PythonName"
	sed -i 's/ \+= \[/ = \[/g' "$PythonName"
	#same for quote
	sed -i "s/='/ = '/g" "$PythonName"
	sed -i "s/ \+= '/ = '/g" "$PythonName"


	#deal with operators (+, -, /, *)
	for OP in "${operators[@]}"; do
    	    sed -i  's/'"$OP"'/ '"$OP"' /g' "$PythonName"
    	    sed -i  's/ \+'"$OP"' \+/ '"$OP"' /g' "$PythonName"
	done
	sed -i 's/\* args/\*args/g' "$PythonName"

	# get multiple operators back together
	sed -i  's/+ \+=/+=/g' "$PythonName"
	sed -i  's/- \+=/-=/g' "$PythonName"
	sed -i  's/= \+=/==/g' "$PythonName"
	sed -i  's/> \+=/>=/g' "$PythonName"
	sed -i  's/< \+=/<=/g' "$PythonName"
	sed -i  's/! \+=/!=/g' "$PythonName"
	sed -i  's/\* \+\*/\*\*/g' "$PythonName"
	sed -i  's/e - /e-/g' "$PythonName"

	#power operator does not need spaces
	sed -i  's/ \+\*\* \+/\*\*/g' "$PythonName"

	# and fix their spacing
	for MOP in "${multop[@]}"; do
    	    sed -i  's/'"$MOP"'/ '"$MOP"' /g' "$PythonName"
    	    sed -i  's/ \+'"$MOP"' \+/ '"$MOP"' /g' "$PythonName"
	done

	# get comparison operators without space in strings (for checkfield)
	for OP in "${compop[@]}";do
	    sed -i "s/' "$OP" '/'"$OP"'/g" "$PythonName"
	done


	sed -i  's/ \+= \+/ = /g' "$PythonName"   #reduce number of spaces to one before and after =

	# fix path names
	templength=1
	while [ $templength -gt 0 ]; do
    	    sed -i  "s/'\([^ ]*\) \/ /'\1\//g w temp.tmp" "$PythonName"
    	    templength=$(wc -l < temp.tmp)
	done

	#replace tab by 4 spaces
	sed -i 's/\t/    /g' "$PythonName"

	#some extraneous spaces
	sed -i 's/( /(/g' "$PythonName"
	#remove trailing spaces too
	sed -i  's/[ \t]*$//' "$PythonName"

	#remove end of line semicolon and backslash
	sed -i 's/;$//g' "$PythonName"
	sed -i 's/\\$//g' "$PythonName"
    fi
done
