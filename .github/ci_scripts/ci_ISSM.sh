#!/bin/bash

# Get configuration
# Source config file
if [ $# -ne 1 ]; then
	#no config file specified: exit
	echo "No config file specified. Exiting..." >&2 # Error message to stderr.
	exit 1
fi
if [ ! -f "$1" ]; then
	echo "File $1 not found!" >&2 # Error message to stderr.
	exit 1
fi

# Initialize test suite variables (to avoid "-eq: unary operator expected")
MATLAB_TEST=0
PYTHON_TEST=0
JAVASCRIPT_TEST=0
EXAMPLES_TEST=0

# Initialize resource variables (to avoid "i<=: syntax error: operand expected" in for loops)
NUMCPUS_INSTALL=1
NUMCPUS_RUN=1

# Source configuration script
source $1;


## External Packages

# Number of packages
NUMPACKAGES=$(($(echo ${EXTERNALPACKAGES} | wc -w ) / 2))
source ${ISSM_DIR}/etc/environment.sh

for ((i=1;i<=$NUMPACKAGES;i++)); do
	NUM1=$((2*$i-1))
	NUM2=$((2*$i))
	PACKAGENAME=$(echo $EXTERNALPACKAGES | cut -d " " -f $NUM1-$NUM1)
	PACKAGEINST=$(echo $EXTERNALPACKAGES | cut -d " " -f $NUM2-$NUM2)

	# Install if requested or if previous install has not been successful
	cd ${ISSM_DIR}/externalpackages/$PACKAGENAME

	# Do not install package if it already exists
	PREFIX=$(egrep "PREFIX=" ./$PACKAGEINST | cut -d'"' -f 2 | envsubst)
	echo "======================================================";
	echo "       Installing $PACKAGENAME                        ";
	echo "======================================================";

	./$PACKAGEINST $NUMCPUS_INSTALL &> compil.log
	if [ $? -ne 0 ] && [ "${PACKAGENAME}" != "boost" ]; then
		cat compil.log
		echo "======================================================";
		echo "    ERROR: installation of $PACKAGENAME failed        ";
		echo "======================================================";
		exit 1;
	fi
	source ${ISSM_DIR}/etc/environment.sh
done

# ISSM compilation
cd $ISSM_DIR
echo "======================================================";
echo "                    Reconfiguring                     ";
echo "======================================================";
autoreconf -ivf
if [ $? -ne 0 ]; then
	echo "autoreconf failed!"
	exit 1
fi

eval "./configure ${ISSM_CONFIG}"
if [ $? -ne 0 ]; then
	echo "ISSM configuration failed (see options below)"
	echo $ISSM_CONFIG
	echo "ISSM configuration failed!"
	exit 1
fi

# Compile and install ISSM
echo "======================================================"
echo "                    Compiling ISSM                    "
echo "======================================================"
if [ $NUMCPUS_INSTALL -gt 1 ]; then
	echo "Making with ${NUMCPUS_INSTALL} CPUs"
	make -j $NUMCPUS_INSTALL install
else
	make install
fi
if [ $? -ne 0 ]; then
	echo "ISSM_COMPILATION failed!"
	exit 1
fi

#Make nightlylog dir
mkdir $ISSM_DIR/nightlylog/

# MATLAB tests {{{
if [ $MATLAB_TEST -eq 1 ]; then

	# Launch MATLAB and the nightly run script
	cat > ${ISSM_DIR}/nightlylog/matlab_run$i.m << EOF
	warning off %necessary to avoid a log of several Go for parallel runs
	try,
		$(if [ "${MATLAB_NROPTIONS}" = "" ]; then
			echo "runme('output','nightly','rank',${i},'numprocs',${NUMCPUS_RUN},'quionerror',1);"
		else
			echo "runme(${MATLAB_NROPTIONS},'output','nightly','rank',${i},'numprocs',${NUMCPUS_RUN},'quionerror',1);"
		fi)
	catch me,
		%An error occured, get report and exit
		message=getReport(me)
		directory=strsplit(pwd,'/');
		fid=fopen([issmdir '/nightlylog/matlaberror.log'], 'at');
		fprintf(fid,'\nMatlab error occured in: %s\n\n',directory{end});
		fprintf(fid,'%s',message);
		fclose(fid);
		quit(1);
	end
	disp('MATLABEXITEDCORRECTLY');
	quit(0);
EOF
	cd $ISSM_DIR/test/NightlyRun
	$MATLAB_PATH/bin/matlab -nodisplay -nosplash -r "addpath ${ISSM_DIR}/src/m/dev; devpath; addpath ${ISSM_DIR}/nightlylog; matlab_run$i"
	if [ $? -ne 0 ]; then
		echo "MATLAB returned an error message!"
		exit 1
	fi
fi
# }}}

# Python tests # {{{
#
# TODO: Figure out why "Running Python test for Rank $i" is printed twice for each CPU
#
if [ $PYTHON_TEST -eq 1 ]; then
	# Launch all tests on different CPUs
	if [ -z $PYTHONPATH ]; then
		export PYTHONPATH="${ISSM_DIR}/src/m/dev"
	else
		export PYTHONPATH="${ISSM_DIR}/src/m/dev:${PYTHONPATH}"
	fi
	export PYTHONSTARTUP="${PYTHONPATH}/devpath.py"
	export PYTHONUNBUFFERED=1 # We don't want Python to buffer otherwise issm.exe output is not captured
	for (( i=1;i<=$NUMCPUS_RUN;i++ ))
	do
		cd $ISSM_DIR/test/NightlyRun
		echo "--------------Running Python test for Rank $i---------------------"
		./runme.py --output=nightly --rank=$i --numprocs=$NUMCPUS_RUN $PYTHON_NROPTIONS &> $ISSM_DIR/nightlylog/python_log$i.log &
		echo "--------------Running Python test for Rank $i---------------------"
	done

	# Concatenate logs
	cd $ISSM_DIR/nightlylog
	if [ -f python_log.log ]; then
		rm python_log.log
	fi

	for job in `jobs -p`
	do
		echo "Waiting on: ${job}"
		wait $job
	done

	for (( i=1;i<=$NUMCPUS_RUN;i++ ))
	do
		echo "This is the concatenation phase for rank: python_log$i.log"
		cat python_log$i.log >> python_log.log
	done
fi
# }}}
