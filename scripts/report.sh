#!/bin/bash
#generate html report from info.log output file

#style
#{{{
TITLE_STYLE='width="700px" cellpadding="10"'
TITLE_FONT='style="color:#6495ed; font-family:Arial, Verdana, Tahoma; font-weight: bold; font-size:25px;" align="center"'
SECTION_STYLE='width="700px" cellpadding="5"'
SECTION_FONT='style="color:#6495ed; font-family:Arial, Verdana, Tahoma; font-size:20px; font-weight: bold;" align="left"'
TABLE_STYLE='width="680px" rules=none bgcolor="#ffffdd" border=1 bordercolor="#000000" cellpadding="3"'
TABLE_FONT='style="color:#404040; font-family:Arial, Verdana, Tahoma; font-size:14px; font-weight: normal;" align="left"'
CODE_STYLE='width="700px" rules=none'
CODE_FONT='style="color:#404040; font-family:Arial, Verdana, Tahoma; font-size:12px; font-weight: normal;" align="left"'
BODY_STYLE='width="700px"'
BODY_FONT='style="color:#404040; font-family:Arial, Verdana, Tahoma; font-size:14px;"'
BODY_FONTC=$(echo $BODY_FONT | sed -e "s/style=\"/style=\"text-align:center; /g")
BODY_FONTL=$(echo $BODY_FONT | sed -e "s/style=\"/style=\"text-align:left; /g")
FOOTER_STYLE='width="700px"  cellpadding="10"'
FOOTER_FONT='style="color:#404040; font-family:Arial, Verdana, Tahoma; font-size:12px; font-weight: normal;" align="center"'
#}}}

#process info.log
#{{{
if [ ! -f info.log ]; then
	echo "File info.log not found!" >&2   # Error message to stderr.
	exit 1
fi 
NRNAME=$(      cat info.log | grep "name:"           | awk '{$1=""}1')
TODAY=$(       cat info.log | grep "today"           | awk '{printf("%s %s",$2,$3);}')
USER=$(        cat info.log | grep "user"            | awk '{print $2}')
VERSION=$(     cat info.log | grep "version"         | awk '{print $2}')
HOST_NAME=$(   cat info.log | grep "host"            | awk '{print $2}')
OS=$(          cat info.log | grep "OS"              | awk '{print $2}')
RELEASE=$(     cat info.log | grep "release"         | awk '{print $2}')
EL_INSTALL=$(  cat info.log | grep "elapsed_install" | awk '{print $2}')
EL_TOTAL=$(    cat info.log | grep "elapsed_total"   | awk '{print $2}')
IS_MATLAB=$(   cat info.log | grep "is_matlab"       | awk '{print $2}')
IS_PYTHON=$(   cat info.log | grep "is_python"       | awk '{print $2}')
EL_MATLAB=$(   cat info.log | grep "elapsed_matlab"  | awk '{print $2}')
EL_PYTHON=$(   cat info.log | grep "elapsed_python"  | awk '{print $2}')
CRASH_MATLAB=$(cat info.log | grep "matlab_crash:"   | awk '{print $2}')
CRASH_PYTHON=$(cat info.log | grep "python_crash:"   | awk '{print $2}')

#Did installation work?
if [ $(ls -1 $ISSM_DIR/bin | wc -l) -le 1 ]; then
	IS_INSTALL=0
else
	IS_INSTALL=1
fi
#}}}

#1. summary table 
#{{{
rm report.html
cat << END >> report.html
<div align="center">
<table $TITLE_STYLE><tr><td $TITLE_FONT>$NRNAME</td></tr></table>

<table $TABLE_STYLE>
<tr> 
<td $TABLE_FONT>host: $HOST_NAME ($OS)</td>
<td $TABLE_FONT>date: $TODAY</td>
</tr>
<tr>
<td $TABLE_FONT>user: $USER</td>
<td $TABLE_FONT>release: $RELEASE</td>
</tr>
<tr>
<td $TABLE_FONT>total elapsed time: $EL_TOTAL</td>
<td $TABLE_FONT>installation elapsed time: $EL_INSTALL</td>
</tr>
<tr>
<td $TABLE_FONT>svn version: $VERSION</td>
<td $TABLE_FONT></td>
</tr>
</table>
<br><hr width="700px">
END
# }}}

#stop if did not install
#{{{
if [ $IS_INSTALL -eq 0 ]; then
	cat << END >> report.html
	<table $(echo $BODY_STYLE) style="border-collapse:collapse;">
	<tr><td $BODY_FONTC>Status: <span style="color:#ff0000">Installation failed</span></td></tr>
	</table>
	<table $FOOTER_STYLE><tr><td $FOOTER_FONT><a href="http://issm.jpl.nasa.gov" title="ISSM website" target="_blank">ISSM</a> nightly run report</td></tr></table>
	</div>
END
exit 0
fi
#}}}

#2. matlab report
if [ $IS_MATLAB -eq 1 ]; then
#Process matlab_log.log {{{
cat matlab_log.log        | egrep 'ERROR|SUCCESS|FAILURE' | grep -v "PETSC" | sed -e "s/>/\&gt;/g" | sed -e "s/</\&lt;/g" > matlab.log
cat matlab.log        | grep -v "SUCCESS" > matlab_short.log
cat matlab_log.log        | grep "PETSC" | sed -e "s/>/\&gt;/g" | sed -e "s/</\&lt;/g" > petscerror.log
NUM_MATLAB_TOT=$(wc -l matlab.log | awk '{print $1}')
NUM_MATLAB_ERR=$(cat matlab.log | grep 'ERROR'   | grep -v "PETSC" | wc -l)
NUM_MATLAB_SUC=$(cat matlab.log | grep 'SUCCESS' | wc -l)
NUM_MATLAB_FAI=$(cat matlab.log | grep 'FAILURE' | wc -l)
#}}}
#write report {{{
cat << END >> report.html
<table $SECTION_STYLE><tr><td $SECTION_FONT>Matlab tests</td></tr></table>
<table $(echo $BODY_STYLE) style="border-collapse:collapse;">
$(if [ $CRASH_MATLAB -eq 0 ]; then
echo "<tr><td $BODY_FONTL>Status: <span style=\"color:#008000\">all tests have been run</span></td></tr>"
else
	echo "<tr><td $BODY_FONTL>Status: <span style=\"color:#ff0000\">Matlab crashed</span></td></tr>"
fi)
<tr><td $BODY_FONTL>Total execution time: $EL_MATLAB</td></tr>
<tr><td $BODY_FONTL>Number of successes: $NUM_MATLAB_SUC/$NUM_MATLAB_TOT</td></tr>
<tr><td $BODY_FONTL>Number of errors: $NUM_MATLAB_ERR/$NUM_MATLAB_TOT</td></tr>
<tr><td $BODY_FONTL>Number of failures: $NUM_MATLAB_FAI/$NUM_MATLAB_TOT</td></tr>
</table>
END
#}}}
fi

#2. python report
if [ $IS_PYTHON -eq 1 ]; then
#Process python_log.log {{{
cat python_log.log        | egrep 'ERROR|SUCCESS|FAILURE' | grep -v "PETSC" | sed -e "s/>/\&gt;/g" | sed -e "s/</\&lt;/g" > python.log
cat python.log        | grep -v "SUCCESS" > python_short.log
cat python_log.log        | grep "PETSC" | sed -e "s/>/\&gt;/g" | sed -e "s/</\&lt;/g" > petscerror.log
NUM_PYTHON_TOT=$(wc -l python.log | awk '{print $1}')
NUM_PYTHON_ERR=$(cat python.log | grep 'ERROR'   | grep -v "PETSC" | wc -l)
NUM_PYTHON_SUC=$(cat python.log | grep 'SUCCESS' | wc -l)
NUM_PYTHON_FAI=$(cat python.log | grep 'FAILURE' | wc -l)
#}}}
#write report {{{
cat << END >> report.html
<table $SECTION_STYLE><tr><td $SECTION_FONT>Python tests</td></tr></table>
<table $(echo $BODY_STYLE) style="border-collapse:collapse;">
$(if [ $CRASH_PYTHON -eq 0 ]; then
	echo "<tr><td $BODY_FONTL>Status: <span style=\"color:#008000\">all tests have been run</span></td></tr>"
else
	echo "<tr><td $BODY_FONTL>Status: <span style=\"color:#ff0000\">Python crashed</span></td></tr>"
fi)
<tr><td $BODY_FONTL>Total execution time: $EL_PYTHON</td></tr>
<tr><td $BODY_FONTL>Number of successes: $NUM_PYTHON_SUC/$NUM_PYTHON_TOT</td></tr>
<tr><td $BODY_FONTL>Number of errors: $NUM_PYTHON_ERR/$NUM_PYTHON_TOT</td></tr>
<tr><td $BODY_FONTL>Number of failures: $NUM_PYTHON_FAI/$NUM_PYTHON_TOT</td></tr>
</table>
END
#}}}
fi

#3. Matlab and python tables
if [ $IS_MATLAB -eq 1 ]; then
#Matlab{{{
#display table ONLY if installation worked and there has been at leat one FAILURE or ERROR
if [ $IS_INSTALL -eq 1 ] && [ $NUM_MATLAB_TOT -gt 1 ] && [ $NUM_MATLAB_SUC -ne $NUM_MATLAB_TOT ]
then
	cat << END >> report.html
<table $SECTION_STYLE><tr><td $(echo $SECTION_FONT)>List of Matlab tests</td></tr></table>
<table $BODY_STYLE style="border-collapse:collapse;">
<tr> 
<th $BODY_FONT>Result</th> 
<th $BODY_FONT>Tolerance</th> 
<th $BODY_FONT>Test id</th>  
<th $BODY_FONT>Test name</th> 
<th $BODY_FONT>Field checked</th>
</tr>
$(cat matlab_short.log | while read line
  do
	  echo "<tr>"
	  STATUS=`echo $line | awk '{print $1}'`
	  #FAILURE
	  if [ "$STATUS" = "FAILURE" ]
	  then
		  FONTC=$(echo "$BODY_FONTC bgcolor=#ffff00");
		  FONTL=$(echo "$BODY_FONTL bgcolor=#ffff00");
		  echo $line | awk -v FONTC="$FONTC" -v FONTL="$FONTL" '
		  { printf("<td %s id=FAILURE>%s</td>\n<td %s>%s</td>\n<td %s>%s</td>\n<td %s>%s</td>\n<td %s>%s</td>\n\n",FONTL,$1,FONTC,$3,FONTC,$6,FONTL,$9,FONTL,$11);}
		  '; 
	  else
		  #SUCCESS
		  if [ "$STATUS" = "SUCCESS" ]
		  then
			  FONTC=$(echo "$BODY_FONTC bgcolor=#ddffdd")
			  FONTL=$(echo "$BODY_FONTL bgcolor=#ddffdd")
			  #do not write anything
		  #ERROR
		  else
			  FONTC=$(echo "$BODY_FONTC bgcolor=#ffdddd id=ERROR")
			  FONTL=$(echo "$BODY_FONTL bgcolor=#ffdddd")
			  echo $line | awk -v FONTC="$FONTC" -v FONTL="$FONTL" '
			  { printf("<td %s>%s</td>\n<td %s>%s%s%s</td>\n<td %s>%s</td>\n<td %s>%s</td>\n<td %s>%s</td>\n\n",FONTL,$1,FONTL,$3,$4,$5,FONTC,$8,FONTL,$11,FONTL,$13);}
			  '; 
		  fi

	  fi
	  echo "</tr>"
  done
	  )
</table>
<br>
END
fi
#}}}
fi
if [ $IS_PYTHON -eq 1 ]; then
#python{{{
#display table ONLY if installation worked and there has been at leat one FAILURE or ERROR
if [ $IS_INSTALL -eq 1 ] && [ $NUM_PYTHON_TOT -gt 1 ] && [ $NUM_PYTHON_SUC -ne $NUM_PYTHON_TOT ]
then
	cat << END >> report.html
	<table $(echo $SECTION_STYLE)><tr><td $(echo $SECTION_FONT)>List of Python tests</td></tr></table>
	<table $(echo $BODY_STYLE) style="border-collapse:collapse;">
	<tr> 
	<th $BODY_FONT>Result</th> 
	<th $BODY_FONT>Tolerance</th> 
	<th $BODY_FONT>Test id</th>  
	<th $BODY_FONT>Test name</th> 
	<th $BODY_FONT>Field checked</th>
	</tr>
	$(cat python_short.log | while read line
do
	echo "<tr>"
	STATUS=`echo $line | awk '{print $1}'`

	#FAILURE
	if [ "$STATUS" = "FAILURE" ]
	then

		FONTC=$(echo "$BODY_FONTC bgcolor=#ffff00");
		FONTL=$(echo "$BODY_FONTL bgcolor=#ffff00");
		echo $line | awk -v FONTC="$FONTC" -v FONTL="$FONTL" '
		{ printf("<td %s id=FAILURE>%s</td>\n<td %s>%s</td>\n<td %s>%s</td>\n<td %s>%s</td>\n<td %s>%s</td>\n\n",FONTL,$1,FONTC,$3,FONTC,$6,FONTL,$9,FONTL,$11);}
		'; 

	else

		#SUCCESS
		if [ "$STATUS" = "SUCCESS" ]
		then
			FONTC=$(echo "$BODY_FONTC bgcolor=#ddffdd")
			FONTL=$(echo "$BODY_FONTL bgcolor=#ddffdd")
			#do not write anything
			#ERROR
		else
			FONTC=$(echo "$BODY_FONTC bgcolor=#ffdddd id=ERROR")
			FONTL=$(echo "$BODY_FONTL bgcolor=#ffdddd")
			echo $line | awk -v FONTC="$FONTC" -v FONTL="$FONTL" '
			{ printf("<td %s>%s</td>\n<td %s>%s%s%s</td>\n<td %s>%s</td>\n<td %s>%s</td>\n<td %s>%s</td>\n\n",FONTL,$1,FONTL,$3,$4,$5,FONTC,$8,FONTL,$11,FONTL,$13);}
			'; 
		fi

	fi
	echo "</tr>"
done
)
</table>
<br>
END
fi
#}}}
fi

#4. Error report
if [ $IS_MATLAB -eq 1 ] && [ -s matlaberror.log ]; then
#Matlab {{{
cat << END >> report.html
<table $SECTION_STYLE><tr><td $SECTION_FONT>Matlab errors</td></tr></table>
<table $CODE_STYLE><tr><td $CODE_FONT>
<pre style="
white-space: -moz-pre-wrap;
white-space: -pre-wrap;
white-space: -o-pre-wrap;
white-space: pre-wrap;
word-wrap: break-word;
">$(cat matlaberror.log)</pre>
</td></tr></table>
END
#}}}
fi
if [ $IS_PYTHON -eq 1 ] && [ -s pythonerror.log ]; then
	#Python {{{
	cat << END >> report.html
<table $SECTION_STYLE><tr><td $SECTION_FONT>Python errors</td></tr></table>
<table $CODE_STYLE><tr><td $CODE_FONT>
<pre style="
white-space: -moz-pre-wrap;
white-space: -pre-wrap;
white-space: -o-pre-wrap;
white-space: pre-wrap;
word-wrap: break-word;
">$(cat pythonerror.log)</pre>
</td></tr></table>
END
	#}}}
fi
if [ -s petscerror.log ]; then
	#PETSc{{{
cat << END >> report.html
<table $SECTION_STYLE><tr><td $SECTION_FONT>PETSc errors</td></tr></table>
<table $CODE_STYLE><tr><td $CODE_FONT>
<pre style="
white-space: -moz-pre-wrap;
white-space: -pre-wrap;
white-space: -o-pre-wrap;
white-space: pre-wrap;
word-wrap: break-word;
">$(cat petscerror.log)</pre>
</td></tr></table>
END
#}}}
fi

#last: footer
#{{{
cat << END >> report.html
<br>
<table $FOOTER_STYLE><tr><td $FOOTER_FONT><a href="http://issm.jpl.nasa.gov" title="ISSM website" target="_blank">ISSM</a> nightly run report</td></tr></table>
</div>
END
#}}}
