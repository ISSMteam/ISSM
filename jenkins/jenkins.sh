#!/bin/bash
################################################################################
# This script manages installation of ISSM on a given Jenkins node using a
# configuration file passed as the only argument. This file also contains
# details about which nightly run tests should be executed after the build has
# been completed. Finally, results of the build and tests are emailed to the
# members of the ISSM development team.
#
# NOTE:
# - Variable OS is set in environment by running 
#	`source $ISSM_DIR/etc/environment.sh`.
#
# TODO:
# - Investigate refactoring parsing of list of changed files
################################################################################

# Override certain aliases
alias grep=$(which grep)

echo "Cleaning up execution directory"
rm -rf ${ISSM_DIR}/execution/*
rm -rf ${ISSM_DIR}/nightlylog
mkdir ${ISSM_DIR}/nightlylog

# Get configuration
# Source config file {{{
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
# }}}

if [[ $EXAMPLES_TEST -eq 1 && $MATLAB_TEST+$PYTHON_TEST+$JAVASCRIPT_TEST -ne 0 ]]; then
	echo "When running examples tests, set *only* EXAMPLES_TEST=1"
	exit 1
fi

# Install ISSM
# Determining installation type depending on last changes to repository {{{
echo "======================================================";
echo "             Determining installation type            "
echo "======================================================";
if [ -f ${ISSM_DIR}/.PREV_COMMIT ]; then
	# Fetch main branch from remote origin (this does not affect local files 
	# like `git pull` would)
	git fetch --quiet origin main

	# Retrieve previous commit SHA
	PREV_COMMIT=$(cat ${ISSM_DIR}/.PREV_COMMIT)

	# Get list of changed files
	CHANGES=$(git diff --name-only ${PREV_COMMIT} FETCH_HEAD)

	if [ ! "${CHANGES}" == "" ]; then
		# Print list of changed files
		echo "   "
		echo "List of changed files"
		echo "---------------------"
		echo "${CHANGES}"
		echo "   "
	fi

	# If the contents of the externalpackages directory were modified in any
	# way, check for changed external packages
	if [ ! -z "$(echo ${CHANGES} | grep externalpackages)" ]; then
		echo "-- checking for changed externalpackages... yes"
		ISSM_EXTERNALPACKAGES="yes"
	else
		echo "-- checking for changed externalpackages... no"
		ISSM_EXTERNALPACKAGES="no"
	fi

	# If the Makefile or m4 directory were changed in any way or if certain
	# binary files from a previous compilation do not exist, reconfigure
	if [ ! -z "$(echo ${CHANGES}| grep -e "Makefile.am" -e "m4" )" ] ||
		[ ! -f "${ISSM_DIR}/bin/issm.exe" ] && [ ! -f "${ISSM_DIR}/bin/issm-bin.js" ] ||
		[ "$ISSM_EXTERNALPACKAGES" == "yes" ]; then
		echo "-- checking for reconfiguration... yes"
		ISSM_RECONFIGURE="yes"
	else
		echo "-- checking for reconfiguration... no"
		ISSM_RECONFIGURE="no"
	fi

	# If source files were changed in any way, recompile
	if [ ! -z "$(echo ${CHANGES} | grep -e "\.cpp" -e "\.h" )" ] ||
		[ "$ISSM_RECONFIGURE" == "yes" ]; then
		echo "-- checking for recompilation... yes"
		ISSM_COMPILATION="yes"
	else
		echo "-- checking for recompilation... no"
		ISSM_COMPILATION="no"
	fi
else
	echo "Fresh copy of repository; building everything"
	echo "-- checking for changed externalpackages... yes"
	echo "-- checking for reconfiguration... yes"
	echo "-- checking for recompilation... yes"
	ISSM_EXTERNALPACKAGES="yes"
	ISSM_RECONFIGURE="yes"
	ISSM_COMPILATION="yes"
fi

# Write out hidden file containing this commit's SHA
git rev-parse HEAD > ${ISSM_DIR}/.PREV_COMMIT

# }}}

## External Packages
#

# Number of packages
NUMPACKAGES=$(($(echo ${EXTERNALPACKAGES} | wc -w ) / 2))

# Jenkins XML files for individual packages
EXTERNAL_TEST_FILE="${ISSM_DIR}/nightlylog/results/external.xml"
mkdir -p ${ISSM_DIR}/nightlylog/results
echo "<testsuite tests=\"$NUMPACKAGES\">" > $EXTERNAL_TEST_FILE

# Need a source here for when builds start midway through installation of externalpackages
source ${ISSM_DIR}/etc/environment.sh

EXTERNALPACKAGES_FAILED=0;

for ((i=1;i<=$NUMPACKAGES;i++)); do
	NUM1=$((2*$i-1))
	NUM2=$((2*$i))
	PACKAGENAME=$(echo $EXTERNALPACKAGES | cut -d " " -f $NUM1-$NUM1)
	PACKAGEINST=$(echo $EXTERNALPACKAGES | cut -d " " -f $NUM2-$NUM2)

	# Install if requested or if previous install has not been successful
	if [ "${ISSM_EXTERNALPACKAGES}" == "yes" ]; then # NOTE: Removed check on if 'install' directory exist
		cd ${ISSM_DIR}/externalpackages/$PACKAGENAME

		# Do not install package if it already exists
		#
		# NOTE: Assumes PREFIX is defined in package installation script
		#
		PREFIX=$(egrep "PREFIX=" ./$PACKAGEINST | cut -d'"' -f 2 | envsubst)
		if [ ! -d "${PREFIX}" ]; then
			echo "======================================================";
			echo "       Installing $PACKAGENAME                        ";
			echo "======================================================";

			./$PACKAGEINST $NUMCPUS_INSTALL &> compil.log
			if [ $? -ne 0 ] && [ "${PACKAGENAME}" != "boost" ]; then
				cat compil.log
				echo "======================================================";
				echo "    ERROR: installation of $PACKAGENAME failed        ";
				echo "======================================================";
				echo "<testcase classname=\"externalpackages\" name=\"$PACKAGENAME\">" >> $EXTERNAL_TEST_FILE
				echo '<failure message="failure">External packages did not install right. Check it.' >> $EXTERNAL_TEST_FILE
				cat compil.log >> $EXTERNAL_TEST_FILE
				echo '</failure>' >> $EXTERNAL_TEST_FILE
				echo '</testcase>' >> $EXTERNAL_TEST_FILE
				echo '</testsuite>' >> $EXTERNAL_TEST_FILE
				exit 1;
			else
				echo "<testcase classname=\"externalpackages\" name=\"$PACKAGENAME\"/>" >> $EXTERNAL_TEST_FILE
			fi
			source ${ISSM_DIR}/etc/environment.sh

			#If external package is rebuilt, we also need to recompile
			ISSM_RECONFIGURE="yes"
			ISSM_COMPILATION="yes"
		else
			echo "======================================================";
			echo "       Skipping $PACKAGENAME                          ";
			echo "======================================================";
			echo "<testcase classname=\"externalpackages\" name=\"$PACKAGENAME\"/>" >> $EXTERNAL_TEST_FILE
		fi
	else
		echo "======================================================";
		echo "       Skipping $PACKAGENAME                          ";
		echo "======================================================";
		echo "<testcase classname=\"externalpackages\" name=\"$PACKAGENAME\"/>" >> $EXTERNAL_TEST_FILE
	fi
done

echo '</testsuite>' >> $EXTERNAL_TEST_FILE

# Source here to include any newly-installed external packages on the path
source ${ISSM_DIR}/etc/environment.sh

if [ "${OS}" == CYGWIN* ]; then
	echo " == WINDOWS ENVIRONMENT DETECTED =="
	source ${ISSM_DIR}/externalpackages/windows/windows_environment.sh
fi

# Set CXX/CC flags for JS runs after external packages to avoid conflicts 
# during their compilation
#
# TODO:
# - Check a different boolean variable as compiling for JavaScript should be
#	independent from running JavaScript tests (one should be able to do the
#	former without having to do the latter).
# - Revisit environment variables (especially EMCC_CFLAGS) once support for
#	Fortran has been accomplished.
#
CXX_PREVIOUS=$CXX
CC_PREVIOUS=$CC
if [ $JAVASCRIPT_TEST -eq 1 ]; then
	export CC=emcc
	export CXX=em++
	export AR=emar
	export RANLIB=emranlib
	#export EMCC_DEBUG=1 # Uncomment to enable debugging
	export EMCC_CFLAGS="-sERROR_ON_UNDEFINED_SYMBOLS=0" # Required after v1.38.14 to avoid undefined symbol warnings from our Fortran object files being treated as errors
	source ${EMSCRIPTEN_ROOT}/emsdk_env.sh
fi

# }}}
# ISSM compilation yes/no (ISSM_COMPILATION) {{{
if [ "${ISSM_COMPILATION}" == "yes" ]; then
	cd $ISSM_DIR
	if [ "${ISSM_RECONFIGURE}" == "yes" ]; then
		echo "======================================================";
		echo "             Cleaning up and reconfiguring            "
		echo "======================================================";
		make uninstall
		make distclean
		./scripts/automakererun.sh
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
	fi

	# Compile and install ISSM
	echo "======================================================"
	echo "                    Compiling ISSM                    "
	echo "======================================================"
	if [ $NUMCPUS_INSTALL -gt 1 ]; then
		echo "Making with ${NUMCPUS_INSTALL} CPUs"

		# To debug compilation/linking, add 'V=1' option to the call to make
		#make -j $NUMCPUS_INSTALL V=1
		make -j $NUMCPUS_INSTALL
	else
		#make V=1
		make
	fi
	if [ $? -ne 0 ] && [ $NUMCPUS_INSTALL -gt 1 ]; then
		echo " "
		echo "Compilation failed, trying to compile with only one thread"
		echo " "
		make
	fi
	if [ $? -ne 0 ]; then
		echo "ISSM_COMPILATION failed!"
		exit 1
	fi
	make install
elif [ "${ISSM_COMPILATION}" == "no" ]; then
	echo "Skipping ISSM compilation"
else
	echo "ISSM_COMPILATION supported values are: yes and no. Exiting..." >&2 # Redirect error messages to stderr
	exit 1
fi
# }}}

# Restore CC/CXX to their previous values
export CC=$CC_PREVIOUS
export CXX=$CXX_PREVIOUS

# MATLAB tests {{{
if [ $MATLAB_TEST -eq 1 ]; then
	MINGW=0
	if [[ $(uname -s) == MINGW* ]]; then
		MINGW=1
		if [ -z "${ISSM_DIR_WIN+x}" ]; then
			export ISSM_DIR_WIN=$(cygpath -w "${ISSM_DIR}")
		fi
	fi

	# Launch all tests on different CPUs
	for (( i=1;i<=$NUMCPUS_RUN;i++ )); do
		# Launch MATLAB and the nightly run script
		cat > ${ISSM_DIR}/nightlylog/matlab_run$i.m << EOF
		warning off %necessary to avoid a log of several Go for parallel runs
		try,
			$(if [ "${MATLAB_NROPTIONS}" = "" ]; then
				echo "runme('output','nightly','rank',${i},'numprocs',${NUMCPUS_RUN});"
			else
				echo "runme(${MATLAB_NROPTIONS},'output','nightly','rank',${i},'numprocs',${NUMCPUS_RUN});"
			fi)
		catch me,
			%An error occured, get report and exit
			message=getReport(me)
			directory=strsplit(pwd,'/');
			fid=fopen([issmdir '/nightlylog/matlaberror.log'], 'at');
			fprintf(fid,'\nMatlab error occured in: %s\n\n',directory{end});
			fprintf(fid,'%s',message);
			fclose(fid);
		end
		disp('MATLABEXITEDCORRECTLY');
		exit
EOF
		cd $ISSM_DIR/test/NightlyRun

		# NOTE: We redirect all output to logfile in order to catch certain errors. For some reason, this does not work under Windows: the logfile option must be used and process must be run in background
		if [[ ${MINGW} -eq 1 ]]; then
			$MATLAB_PATH/bin/matlab -nodesktop -nosplash -nojvm -r "addpath ${ISSM_DIR_WIN}/src/m/dev; devpath; addpath ${ISSM_DIR_WIN}/nightlylog; matlab_run$i" -logfile ${ISSM_DIR_WIN}/nightlylog/matlab_log$i.log &
		else
			$MATLAB_PATH/bin/matlab -nodisplay -nosplash -r "addpath ${ISSM_DIR}/src/m/dev; devpath; addpath ${ISSM_DIR}/nightlylog; matlab_run$i" &> ${ISSM_DIR}/nightlylog/matlab_log$i.log &
		fi
	done

	# Wait for MATLAB to exit
	#
	# TODO:
	# - Replace by adding -wait option to above calls to MATLAB?
	#
	if [[ ${MINGW} -eq 1 ]]; then
		sleep 5;
		echo "Waiting for MATLAB to exit"
		pid=$(ps -W | grep MATLAB | awk '{print $1}')
		echo '-----------------------------'
		echo "pid: ${pid}"
		echo '-----------------------------'

		# Time out after $max_time seconds because sometimes multiple MATLAB processes get locked in race condition
		timer=0
		max_time=7200
		while [[ $timer -lt $max_time && -n "${pid}" ]]; do
			pid=$(ps -W | grep MATLAB | awk '{print $1}')
			timer=$((timer + 1))
			sleep 1;
		done

		# Check if timer hit $max_time
		if [ $timer -eq $max_time ]; then
			echo "Testing timed out at ${timer} seconds"
			# Kill MATLAB processes
			pid=$(ps -W | grep MATLAB | awk '{print $1}')
			echo "${pid}" | xargs /bin/kill -f
			exit 1
		fi
	else
		wait
	fi

	# Concatenate logs
	cd $ISSM_DIR/nightlylog

	if [ -f matlab_log.log ]; then
		rm matlab_log.log
	fi

	for job in `jobs -p`; do
		#echo "Waiting on: ${job}" # Commented out because it really has nothing to do with MATLAB processes
		wait $job
	done

	for (( i=1;i<=$NUMCPUS_RUN;i++ )); do
		cat matlab_log$i.log >> matlab_log.log
	done

	# Filter out Windows characters
	cat matlab_log.log | tr -cd '\11\12\40-\176' > matlab_log.log2 && mv matlab_log.log2 matlab_log.log
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
	export PYTHONSTARTUP="${ISSM_DIR}/src/m/dev/devpath.py"
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

# Clean up ADOL-C tape files
rm -f $ISSM_DIR/execution/*/ADOLC-*

# Examples Tests # {{{
if [ $EXAMPLES_TEST -eq 1 ]; then
	export MATLAB_PATH

	# Download examples datasets if they are not already present
	if [[ -z $(ls -A1q $ISSM_DIR/examples/Data | grep -vw -E "download.sh") ]]; then
		$ISSM_DIR/scripts/DownloadExamplesDatasets.sh
	fi

	$ISSM_DIR/jenkins/examples_tests.sh
fi
# }}}

# Process logs to be JUnit compatible # {{{
cd $ISSM_DIR/nightlylog
source $ISSM_EXT_DIR/shell2junit/install/bin/sh2ju.sh
juLogClean

if [ $MATLAB_TEST -eq 1 ]; then
	# Strip special characters
	sed -i \
		-e 's|\[92m||g' \
		-e 's|\[m||g' \
		-e 's|\x1B||g' \
		matlab_log.log

	# Number tests
	numtests=`cat matlab_log.log | grep "\-\-\-\-\-\-\-\-starting" | wc -l`
	testlist=`cat matlab_log.log | grep "\-\-\-\-\-\-\-\-starting" | sed 's/----------------starting://g' | sed 's/-//g'`

	# NOTE: Matching pattern is "----------------starting:<TEST_NUM>-". Note trailing hyphen, which guards against greedy matching of, for example, "2010" when searching for "201".
	for i in `echo $testlist`; do
		juLog -test=MATLAB-$i -name=Error -error=ERROR awk "/starting:${i}-/{flag=1;next}/finished:${i}-/{flag=0} flag{print}" matlab_log.log
		juLog -test=MATLAB-$i -name=Failure -error=FAILURE awk "/starting:${i}-/{flag=1;next}/finished:${i}-/{flag=0} flag{print}" matlab_log.log
	done

	# Check that MATLAB did not exit in error
	matlabExitedInError=`grep -E "Activation cannot proceed|Error in|Illegal|Invalid MEX-file|Warning: Name is nonexistent or not a directory" matlab_log.log | wc -l`

	if [ $matlabExitedInError -ne 0 ]; then
		echo "----------MATLAB exited in error!----------"
		cat matlab_log.log
		echo "-----------End of matlab_log.log-----------"

		# Clean up execution directory
		rm -rf $ISSM_DIR/execution/*

		exit 1
	fi
fi

if [ $PYTHON_TEST -eq 1 ]; then
	# Strip special characters
	sed -i \
		-e 's|\[92m||g' \
		-e 's|\[m||g' \
		-e 's|\x1B||g' \
		python_log.log

	# Number tests
	numtests=`cat python_log.log | grep "\-\-\-\-\-\-\-\-starting" | wc -l`
	testlist=`cat python_log.log | grep "\-\-\-\-\-\-\-\-starting" | sed 's/----------------starting://g' | sed 's/-//g'`

	# NOTE: Matching pattern is "----------------starting:<TEST_NUM>-". Note trailing hyphen, which guards against greedy matching of, for example, "2010" when searching for "201".
	for i in `echo $testlist`; do
		juLog -test=PYTHON-$i -name=Error -error=ERROR awk "/starting:${i}-/{flag=1;next}/finished:${i}-/{flag=0} flag{print}" python_log.log
		juLog -test=PYTHON-$i -name=Failure -error=FAILURE awk "/starting:${i}-/{flag=1;next}/finished:${i}-/{flag=0} flag{print}" python_log.log
	done

	# Check that Python did not exit in error
	pythonExitedInError=`grep -c -E "Error|No such file or directory|Permission denied|Standard exception|Traceback|bad interpreter|syntax error|error:" python_log.log`

	if [ $pythonExitedInError -ne 0 ]; then
		echo "----------Python exited in error!----------"
		cat python_log.log
		echo "-----------End of python_log.log-----------"

		# Clean up execution directory
		rm -rf $ISSM_DIR/execution/*

		exit 1
	fi
fi

if [ $EXAMPLES_TEST -eq 1 ]; then
	# Strip special characters
	sed -i \
		-e 's|\[92m||g' \
		-e 's|\[m||g' \
		-e 's|\x1B||g' \
		matlab_log_examples.log

	# Inexplicably, there are backspace characters in the error output; remove them
	sed -i -e 's|\x08||g' matlab_log_examples.log

	numtests=`cat matlab_log_examples.log | grep "starting: " | wc -l`
	testlist=`cat matlab_log_examples.log | grep "starting: " | sed 's/starting: //'`

	for i in `echo $testlist`; do
		juLog -test=Example-$i -name=Error -error=ERROR awk "/^starting: ${i}$/{flag=1;next}/^finished: ${i}$/{flag=0} flag{print}" matlab_log_examples.log
		juLog -test=Example-$i -name=Failure -error=FAILURE awk "/^starting: ${i}$/{flag=1;next}/^finished: ${i}$/{flag=0} flag{print}" matlab_log_examples.log
	done

	# Check that MATLAB did not exit in error
	matlabExitedInError=`grep -E "Activation cannot proceed|Error in|Illegal|Invalid MEX-file|Warning: Name is nonexistent or not a directory" matlab_log_examples.log | wc -l`

	if [ $matlabExitedInError -ne 0 ]; then
		echo "----------MATLAB exited in error!----------"
		cat matlab_log_examples.log
		echo "-----------End of matlab_log.log-----------"

		# Clean up execution directory
		rm -rf $ISSM_DIR/execution/*

		exit 1
	fi
fi

# Clean up execution directory
rm -rf $ISSM_DIR/execution/*
# }}}
