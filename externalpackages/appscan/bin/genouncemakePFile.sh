#!/bin/bash

#----------------------------------------------------------------------
#VARIABLES containing information to generate OunceMake Properties file
#----------------------------------------------------------------------
XMLHEADER="<?xml version=\""1.0\"" encoding=\"UTF-8\"?>"
OUNCEHEADER="<OunceMakeProperties xmlns=\"http://www.ouncelabs.com/namespace\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">"
MAKEOPT="<MakeOptions>-f Makefile</MakeOptions>"

COMPILER_TAG="Compiler"
COMPILER_MACRO="macro"

C_COMPILER_MVALUE="\"CC\""
CPP_COMPILER_MVALUE="\"CXX\""

MAKE_TAG="Make"
MAKE_MACRO="macro"
MAKE_MVALUE="\"MAKE\""


LINKER_TAG="Linker"
LINKER_MACRO="macro"
LINKER_MVALUE="\"LD\""


ARCHIVE_TAG="Executable"
ARCHIVE_MACRO="macro"
ARCHIVE_MVALUE="\"ARCHIVE\""

OUNCEOPT="<Options recursive=\"true\" single_project=\"false\" verbose=\"false\" clean=\"make clean\" no_clean=\"false\" build=\"true\"/>"

OUNCE_MAIN_TAG="OunceMakeProperties"
OUNCE_COMPILER_GPO_TAG="GlobalProjectOptions"
OUNCE_COMPILER_OPT="compiler_options"
OUNCE_COMPILER_OPT_gpp="\"-g++_linux_i386\""
OUNCE_COMPILER_OPT_gcc="\"-gcc_linux_i386\""

OUNCE_INCLUDE_PATH="include_paths"
OUNCE_MACRO_PATH="macros"

OUNCE_COMPILER_FO_TAG="FileOptions"
OUNCE_COMPILER_EXT="extensions"

OUNCE_COMPILER_EXT_CPP="cpp;cxx"
OUNCE_COMPILER_EXT_C="c"

#Did the search for a particular tool pass/fail
FAIL_FILE="failed.status"

#---------------------------------
#COMPILER data extraction function
#---------------------------------

function initialize {
#Create Directory for info
 echo;echo
 echo -e "Initialization....\t\t\t[\c "
 echo -e "\E[33mSTARTED\c "; tput sgr0
 echo "]"
 search_tool which
 FILEPRE=${1}

 case "$1" in
	"gcc"|"cc")
	  TESTFILE="testfoo.c"
	;;
	"g++"|"CC")
	  TESTFILE="testfoo.cpp"
	;;
	*)
	 echo "Error setting prefix"
	;;
 esac

 DATAFOLDER="${FILEPRE}results"
#Declare files to generate and store results
 RESULTFILE="${DATAFOLDER}/${FILEPRE}libinc.result"
 SLIMFILE="${DATAFOLDER}/${FILEPRE}slim.inc"
 SLIM2FILE="${DATAFOLDER}/${FILEPRE}slim2.inc"

#Declare the Tags that enclose the information needed
 START_TAG="#include <...> search starts here:"
 END_TAG="End of search list."

#Declare Variable to track Tags
 MARK_FOUND=0

#For Macros
 MACROSDUMPFILE="${DATAFOLDER}/${FILEPRE}macrosd.result"
 INTFILE1="${DATAFOLDER}/${FILEPRE}interf1.result"
 INTFILE2="${DATAFOLDER}/${FILEPRE}interf2.result"


#The file extensions for (both) options
CONFIG_FOLDER="config"
C_EXT="${CONFIG_FOLDER}/cext.data"
CPP_EXT="${CONFIG_FOLDER}/cppext.data"

if [ -d ${CONFIG_FOLDER} ]; then
	echo
else
	mkdir ${CONFIG_FOLDER}
	echo ${OUNCE_COMPILER_EXT_CPP} > ${C_EXT}
	echo ${OUNCE_COMPILER_EXT_C} > ${CPP_EXT}
fi

 echo -e "Initialization....\t\t\t[\c "
 echo -e "\E[33mDONE\c "; tput sgr0
 echo "]"
}

function purgefiles {
 echo -e "Purging files....\t\t\t[\c "
 echo -e "\E[33mSTARTED\c "; tput sgr0
 echo "]"
  if [ -e ${DATAFOLDER} ]; then
	rm -R ${DATAFOLDER}
	mkdir ${DATAFOLDER}
  else
	mkdir ${DATAFOLDER}
  fi
 echo -e "Purging files....\t\t\t[\c "
 echo -e "\E[33mDONE\c "; tput sgr0
 echo "]"
}


function processGNUInclude {
echo
#Create the dummy file is does not exist
if [ -e ${TESTFILE} ]; then
	echo "Will use an existing ${TESTFILE} file"
else 
	echo -e "int main(){\n return 0;\n}" &>  ${TESTFILE}
	if [ -e ${TESTFILE} ]; then
	  echo -e "Generating new test file\t\t[\c "
	  echo -e "\E[31m${TESTFILE}\c " ; tput sgr0 
	  echo "]"
	else
	  echo -e "Error generating new test file\t\t[\c "
	  echo -e "\E[31m${TESTFILE}\c " ; tput sgr0 
	  echo "]"
   	#TO DO
   	#Will the user like to continue or exit?
	fi
fi


 echo
 echo -e "Generating information for Compiler\t\t[\c "
 echo -e "\E[31m${COMPILER}\c " ; tput sgr0 
 echo "]"

${COMPILER} -v ${TESTFILE} &> ${RESULTFILE}


cat ${RESULTFILE} | while read line
do 

	if [ "${line}" = "${START_TAG}" ] 
	then
		MARK_FOUND=1
	fi

	
	if [ "${line}" = "${END_TAG}" ] 
	then
		MARK_FOUND=0
	fi

	
	if [ "${MARK_FOUND}" = 1 ] && [ "${line}" != "${START_TAG}" ] 
	then
		let "inclib +=  1"
		echo ${line} 1>> ${SLIMFILE}
	fi
done

if [ -e ${SLIMFILE} ]; then
	sed -e :a -e N -e 's/\n/; /' -e ta  ${SLIMFILE}  > ${RESULTFILE}
	echo -e "Inc & Lib results are in file...\t\t[\c "
	echo -e "\033[1m\E[31m${RESULTFILE}\033[0m\c " ; tput sgr0
	echo "]"
fi

#String to contatenate the include path
concatpath="include_paths=\""
#Merge file content
sed -e 's/\n/;/'  ${SLIMFILE} > ${SLIM2FILE}

}


function processGNUMacro {
#*****************************
#Generate macros for system
#*****************************
 echo
 echo -e "Generating Macro data for Compiler\t\t[\c "
 echo -e "\E[31m${COMPILER} \c " ; tput sgr0 
 echo "]"

 ${COMPILER} -dM -E ${TESTFILE} >> ${MACROSDUMPFILE}
 sed 's/#define //' ${MACROSDUMPFILE} > ${INTFILE1}
 if [ -e ${INTFILE1} ]; then
	sed 's/ /=/' ${INTFILE1} > ${INTFILE2}
 else
	echo -e "\033[1m Sorry <${INTFILE1}> was not created.. \033[0m"
 fi

 if [ -e ${INTFILE2} ]; then
	#Remove the Macro's that have no values
	sed '/^_.*=$/d' ${INTFILE2} > ${INTFILE1} 
	#Join the individual line and delimit with a (; )
	sed -e :a -e N  -e 's/\n/; /' -e ta  ${INTFILE1}  > ${MACROSDUMPFILE}
	#All done display the results file
	echo -e "Macro results are in file...\t\t\t[\c "
	echo -e "\033[1m\E[31m${MACROSDUMPFILE}\033[0m\c " ; tput sgr0
	echo "]"
 else
	echo -e "\033[1m Sorry <${MACROSDUMPFILE}> was not created.. \033[0m"
 fi
}

function processSUNMacro {
#*******************************************
#This function is for the SUN Compilers ONLY
#*******************************************
 echo
 echo -e "Generating Macro data for Compiler\t\t[\c "
 echo -e "\E[31m${COMPILER} \c " ; tput sgr0 
 echo "]"

 ${COMPILER} -xdumpmacros ${TESTFILE} &> ${MACROSDUMPFILE}
 sed 's/#define //' ${MACROSDUMPFILE} > ${INTFILE1}
 if [ -e ${INTFILE1} ]; then
	sed 's/ /=/' ${INTFILE1} > ${INTFILE2}
 else
	echo -e "\033[1m Sorry <${INTFILE1}> was not created.. \033[0m"
 fi

 if [ -e ${INTFILE2} ]; then
	#Remove the Macro's that have no values
	sed '/^_.*=$/d' ${INTFILE2} > ${INTFILE1} 
	#Join the individual line and delimit with a (; )
	sed -e :a -e N  -e 's/\n/; /' -e ta  ${INTFILE1}  > ${MACROSDUMPFILE}
	#All done display the results file
	echo -e "Macro results are in file...\t\t\t[\c "
	echo -e "\033[1m\E[31m${MACROSDUMPFILE}\033[0m\c " ; tput sgr0
	echo "]"
 else
	echo -e "\033[1m Sorry <${MACROSDUMPFILE}> was not created.. \033[0m"
 fi
}

#--------------------------
#DATA Processing variables
#--------------------------
PROPFILE="ouncemake_properties.xml"
PROBEFILE="probe.results"
ERRORLOG="error.log"

#---------------------------------------
#FUNCTION TO CREATE FILE
#--------------------------------------

function create_file {
#If a previous copy of file exist, delete it
if [ -e ${PROPFILE} ]; then
	rm ${PROPFILE}
fi

echo;echo
echo "CREATING THE PROPERTIES XML FILE"
#insert file headers
echo ${XMLHEADER} >> ${PROPFILE}
echo ${OUNCEHEADER} >> ${PROPFILE}
echo ${MAKEOPT} >> ${PROPFILE}
case "${1}" in
	'gcc')
	  search_tool gcc
	;;
	'g++')
	  search_tool g++
	;;
	'bothGNU')
	  search_tool gcc
	  search_tool g++
	;;
	'cc')
	  search_tool cc
	;;
	'CC')
	  search_tool CC
	;;
	'bothSUN')
	  search_tool gcc
	  search_tool g++
	;;
	'*')
	;;
esac
search_tool make
search_tool ld
search_tool ar
echo "	${OUNCEOPT}" >> ${PROPFILE}

case "${1}" in
	'gcc')
	  singleGNUC gcc
	;;
	'g++')
	  singleGNUC g++
	;;
	'bothGNU')
	  bothGNUC gcc
	  bothGNUC g++
	;;
	'cc')
	  singleSUNC cc
	;;
	'CC')
	  singleSUNC CC
	;;
	'bothSUN')
	  bothSUNC cc
	  bothSUNC CC
	;;
	'*')
	;;
esac

echo "</${OUNCE_MAIN_TAG}>" >> ${PROPFILE}

echo
echo -e "Successully created Properties file\t\t[\c "
echo -e "\033[1m\E[31m${PROPFILE}\033[0m\c " ; tput sgr0
echo "]"
echo
}

#when the system is unable to locate the tool
#the user is prompted to provide that information
function usersupplypath {
	echo "Please provide the path to ${1}."
	read TOOLPATH
}


function search_tool {
which $1 &> ${PROBEFILE}
#Why this approach? I had to be able to non-interactively  
#know if the pattern was successfully found. 

sed -n -e "/.*no[ ]${1}[ ]in[ ].*/w $FAIL_FILE" ${PROBEFILE} 

#If the file is empty, the process PASSED
#else the Fail criteria was found. Process FAILED
if [ -e ${FAIL_FILE} ]; then
	exec < ${FAIL_FILE}
	read FAILPASS
	if [${FAILPASS} = ""]; then
	  SEARCHTOOL="PASSED"
	  rm ${FAIL_FILE}
	else
	 SEARCHTOOL="FAILED"
	fi
fi


if [ ${SEARCHTOOL} == "FAILED" ] ; then
  echo "Failed to find ${1} in PATH" >> ${ERRORLOG}
  echo -e "check for $1 in ${TLOC}\t\t\t[$SEARCHTOOL]"
elif [ ${SEARCHTOOL} == "PASSED" ]; then
  #Redirects stdin to a file so that
  exec < ${PROBEFILE}
  #I can read line by line
  read TLOC
  echo -e "check for $1 in ${TLOC}\t\t\t[\c "
  echo -e "\E[33m${SEARCHTOOL}\c "; tput sgr0
  echo "]"
else
  echo "Unable to process the request"
  #exit
fi

case "$1" in
	'gcc')
	#set up the sub strings 
	  TAG_L="<${COMPILER_TAG} ${COMPILER_MACRO}=${C_COMPILER_MVALUE}>"
	  MVALUE="$TLOC"
	  TAG_R="</${COMPILER_TAG}>"
	;;
	'g++')
	  TAG_L="<${COMPILER_TAG} ${COMPILER_MACRO}=${CPP_COMPILER_MVALUE}>"
	  MVALUE="$TLOC"
	  TAG_R="</${COMPILER_TAG}>"
	;;
	'make')
	  TAG_L="<${MAKE_TAG} ${MAKE_MACRO}=${MAKE_MVALUE}>"
	  MVALUE="$TLOC"
	  TAG_R="</${MAKE_TAG}>"
	;;
	'ld')
	  TAG_L="<${LINKER_TAG} ${LINKER_MACRO}=${LINKER_MVALUE}>"
	  MVALUE="$TLOC"
	  TAG_R="</${LINKER_TAG}>"
	;;
	'ar')
	  TAG_L="<${ARCHIVE_TAG} ${ARCHIVE_MACRO}=${ARCHIVE_MVALUE}>"
	  MVALUE="$TLOC"
	  TAG_R="</${ARCHIVE_TAG}>"
	;;
	'*')
	;;	
esac

#then concatenate
FULL_STRING="${TAG_L}${MVALUE}${TAG_R}"
#and store in file
echo "	$FULL_STRING" >> ${PROPFILE}
}
function singleGNUC {

echo "	<${OUNCE_COMPILER_GPO_TAG}" >> ${PROPFILE}

case "$1" in 
	'gcc')
	  echo -e "\t  ${OUNCE_COMPILER_OPT}=${OUNCE_COMPILER_OPT_gcc}" >> ${PROPFILE}
	;;
	'g++')
	  echo -e "\t  ${OUNCE_COMPILER_OPT}=${OUNCE_COMPILER_OPT_gpp}" >> ${PROPFILE}
	;;
	'*')
	;;	
esac

exec < ${RESULTFILE}
read INCLUDE_DATA
echo -e "\t  ${OUNCE_INCLUDE_PATH}=\"${INCLUDE_DATA}\"" >> ${PROPFILE}
exec < ${MACROSDUMPFILE}
read MACRO_DATA
#Using single quotes to wrap double quotes in data
echo -e "\t  ${OUNCE_MACRO_PATH}='${MACRO_DATA}' />" >> ${PROPFILE}

}


function bothGNUC {

echo -e "\t<${OUNCE_COMPILER_FO_TAG}" >> ${PROPFILE}

case "$1" in 
	'gcc')
	  FILEPRE="$1"
	  exec < ${C_EXT}
	  read CEXTENSIONS
	  if [ ${CEXTENSIONS} = "" ]; then
	  	CEXTENSIONS=${OUNCE_COMPILER_EXT_C}
	  else
		OUNCE_COMPILER_EXT_C=${CEXTENSIONS}
	  fi
	  echo -e "\t  ${OUNCE_COMPILER_EXT}=\"${OUNCE_COMPILER_EXT_C}\"" >> ${PROPFILE}
	;;
	'g++')
	  FILEPRE="$1"
	  exec < ${CPP_EXT}
	  read CPPEXTENSIONS
	  if [ ${CPPEXTENSIONS} = "" ]; then
	  	CPPEXTENSIONS=${OUNCE_COMPILER_EXT_CPP}
	  else
		OUNCE_COMPILER_EXT_CPP=${CPPEXTENSIONS}
	  fi
	  echo -e "\t  ${OUNCE_COMPILER_EXT}=\"${OUNCE_COMPILER_EXT_CPP}\"" >> ${PROPFILE}
	;;
	'*')
	;;	
esac

#To make sure that the correct folder and file is read
DATAFOLDER="${FILEPRE}results"
#Declare files to generate and store results
RESULTFILE="${DATAFOLDER}/${FILEPRE}libinc.result"

exec < ${RESULTFILE}
read INCLUDE_DATA
echo -e "\t  ${OUNCE_INCLUDE_PATH}=\"${INCLUDE_DATA}\"" >> ${PROPFILE}
exec < ${MACROSDUMPFILE}
read MACRO_DATA
echo -e "\t  ${OUNCE_MACRO_PATH}='${MACRO_DATA}'/>" >> ${PROPFILE}

}


clear
echo
echo -e "\t \E[34m\033[1mOUNCEMAKE PROPERTIES FILE GENERATOR (1.0)b\033[0m";tput sgr0

#Start of the program
case "$1" in
	'-t')
	  case "$2" in 
		"gcc")
		  COMPILER=${2}
		 initialize ${2}
		 purgefiles
		 processGNUInclude
		 processGNUMacro
		 create_file ${2}
		;;
		"g++")
		  COMPILER=${2}
		 initialize ${2}
		 purgefiles
		 processGNUInclude
		 processGNUMacro
		 create_file ${2}
		;;
		"bothGNU")
		  COMPILER="gcc"
		 initialize ${2}
		 purgefiles
		 processGNUInclude
		 processGNUMacro
		  COMPILER="g++"
		 initialize ${2}
		 purgefiles
		 processGNUInclude
		 processGNUMacro
		 create_file ${2}
		;;
		"cc")
		  COMPILER=${2}
		 initialize ${2}
		 purgefiles
		 processSUNInclude
		 processSUNMacro
		 create_file ${2}
		;;
		"CC")
		  COMPILER=${2}
		 initialize ${2}
		 purgefiles
		# processSUNInclude
		 processSUNMacro
		# create_file ${2}
		;;
		"bothSUN")
		  COMPILER="cc"
		 initialize cc
		 purgefiles
		 processSUNInclude
		 processSUNMacro
		  COMPILER="CC"
		 initialize CC
		 purgefiles
		 processSUNInclude
		 processSUNMacro
		 create_file both
		;;
		'*')
		  echo "Usage: -t [gcc|g++|bothGNU]"
		  exit $?
		;;
	  esac  
	;;
	'-h')
	  echo -e "Usage: -t [Options] -v "
	  echo -e "Options:"
	  echo -e "\t[gcc]\tC project only";
	  echo -e "\t[g++]\tCPP project only"
	  echo -e "\t[bothGNU]\tMixed C and CPP project"
	  echo
	  echo -e "\t-v \tDisplay version information"
	  exit $?
	;;
	'-v')
	  echo "Version 1.0 beta"
	  exit $?
	;;
	*)
	  echo -e "Usage: -t [Options] -v "
	  echo -e "Options:"
	  echo -e "\t[gcc]\tC project only";
	  echo -e "\t[g++]\tCPP project only"
	  echo -e "\t[both]\tMixed C and CPP project"
	  echo
	  echo -e "\t-v \tDisplay version information"
	  exit $?
	;;
esac
echo "[ P R O C E S S  C O M P L E T E !!! ]"
echo
