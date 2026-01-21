#!/bin/bash

PATTERN1="print('launching solution sequence on remote cluster')"
PATTERN2="#Execute Queue job"

#only for line with the following pattern (use PATTERN=PATTERN1 to turn off)
PATTERN=$PATTERN1
#PATTERN="printLine"

FILES=$(find ./ -name "*.h" -o -name "*.m" -o -name "*.py" -o -name "*.cpp" | xargs grep "$PATTERN"| sed -e "s/:/ /g" | awk '{print $1}' | sort -u)
for FILE  in $FILES
do 

	echo "modifying $FILE"
	LINES=$(cat $FILE | grep -n "$PATTERN" | sed -e "s/:/ /g" | awk '{print $1}' )

	for LINE in $LINES
	do 
		cat $FILE | sed -e ""$LINE"s/$PATTERN1/$PATTERN2/g" > $FILE$$.bak
		mv $FILE$$.bak $FILE
	done

done
