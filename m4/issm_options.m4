dnl ISSM Options

dnl TODO:
dnl - Check if we need statements such as,
dnl
dnl.  	  AM_CONDITIONAL([JAVASCRIPT], [test "x${HAVE_JAVASCRIPT}" = "xyes"])
dnl
dnl	  when we have already performed a similar check,
dnl
dnl  	  if test "x${JAVASCRIPT}" = "xno"; then
dnl
dnl - Move library dependency checks from end of file to appropriate places
dnl   inline
dnl - Refactor conditionals that test both -d <file> and -f <file>
dnl

AC_DEFUN([ISSM_OPTIONS],[
	AC_MSG_NOTICE(============================================================================)
	AC_MSG_NOTICE(=                      Checking ISSM specific options                      =)
	AC_MSG_NOTICE(============================================================================)

	dnl ISSM's internal options
	dnl Build info{{{

	dnl Build date
	AC_PATH_PROGS(DATE, date)
	AC_MSG_CHECKING([for build date])
	if test "$DATE"; then
		PACKAGE_DATE=`date`
	else
		PACKAGE_DATE="unknown"
	fi
	AC_DEFINE_UNQUOTED([PACKAGE_BUILD_DATE], "${PACKAGE_DATE}", [build date])
	AC_MSG_RESULT([${PACKAGE_DATE}])

	dnl User name
	AC_MSG_CHECKING([user name])
	if test -n "$USER"
	then
		user_name="$USER"
	else
		if test -n "$LOGNAME"; then
			user_name ="$LOGNAME"
		else
			user_name =`(whoami) 2>/dev/null` || user_name=unknown
		fi
	fi
	AC_DEFINE_UNQUOTED([USER_NAME], "${user_name}", [user name])
	AC_MSG_RESULT([${user_name}])

	AC_MSG_CHECKING([host full OS name and version])
	dnl Normalize some host OS names
	case ${host_os} in
		dnl linux is linux is linux, regardless of RMS
		linux-gnu* | lignux* )	host_os=linux ;;
	esac
	AC_DEFINE_UNQUOTED([HOST_OS], "${host_os}", [host full OS name and version])
	AC_MSG_RESULT([${host_os}])

	AC_MSG_CHECKING([host cpu])
	AC_DEFINE_UNQUOTED([HOST_CPU], "${host_cpu}", [host CPU])
	AC_MSG_RESULT([${host_cpu}])

	AC_MSG_CHECKING([vendor])
	AC_DEFINE_UNQUOTED([HOST_VENDOR], "${host_vendor}", [host vendor])
	AC_MSG_RESULT([${host_vendor}])

	AC_MSG_CHECKING([host OS name])
	host_os_name=`echo $host_os | sed 's/\..*//g'`
	dnl Normalize some OS names
	case ${host_os_name} in
		dnl linux is linux is linux, regardless of RMS.
		linux-gnu* | lignux* )	host_os_name=linux ;;
	esac
	AC_DEFINE_UNQUOTED([HOST_OS_NAME], "${host_os_name}", [host OS name])
	AC_MSG_RESULT([${host_os_name}])

	dnl Parse out the OS version of the host
	AC_MSG_CHECKING([host OS version])
	host_os_version=`echo $host_os | sed 's/^[[^0-9]]*//g'`
	if test -z "$host_os_version"; then
		host_os_version=`(uname -r) 2>/dev/null` || host_os_version=unknown
	fi
	AC_DEFINE_UNQUOTED([HOST_OS_VERSION], "${host_os_version}", [host OS version])
	AC_MSG_RESULT([${host_os_version}])

	dnl Determine host architecture (different than CPU)
	AC_MSG_CHECKING([host OS architecture])
	host_arch=`(uname -m) 2>/dev/null` || host_arch=unknown
	dnl Normalize some names
	case ${host_arch} in
		sun4* )	host_arch=sun4 ;;
		sun3x )	host_arch=sun3 ;;
		sun )	host_arch=`(arch) 2>/dev/null` || host_arch=unknown ;;
		i?86 )	host_arch=i386 ;; # all x86 should show up as i386
	esac
	AC_DEFINE_UNQUOTED([HOST_ARCH], "${host_arch}", [host archictecture])
	AC_MSG_RESULT([${host_arch}])

	dnl }}}
	dnl Debugging {{{
	AC_ARG_ENABLE(
		[debugging],													dnl feature
		AS_HELP_STRING([--enable-debugging], [turn debug support on]),	dnl help string
		[enable_debugging=${enableval}],								dnl action if given
		[enable_debugging=no]											dnl action if not given
	)
	AC_MSG_CHECKING(for debugging support)
	if test "x${enable_debugging}" == "xyes"; then
		AC_DEFINE([_ISSM_DEBUG_], [1], [Macro to enable debugging in ISSM])
	fi
	AC_MSG_RESULT([${enable_debugging}])
	dnl }}}
	dnl Development{{{
	AC_ARG_ENABLE(
		[development],													dnl feature
		AS_HELP_STRING([--enable-development], [turn development on]),  dnl help string
		[enable_development=${enableval}],								dnl action if given
		[enable_development=no]											dnl action if not given
	)
	AC_MSG_CHECKING(for development support)
	if test "x${enable_development}" == "xyes"; then
		AC_DEFINE([_DEVELOPMENT_], [1], [enable development support in ISSM])
	fi
	AM_CONDITIONAL([DEVELOPMENT], [test "x${enable_development}" == "xyes"])
	AC_MSG_RESULT([${enable_development}])
	dnl }}}
	dnl Standalone Options {{{
	AC_ARG_ENABLE(
		[standalone-modules],															dnl feature
		AS_HELP_STRING([--enable-standalone-modules], [produce standalone modules]),	dnl help string
		[enable_standalone_modules=${enableval}],										dnl action if given
		[enable_standalone_modules=no]													dnl action if not given
	)
	AC_MSG_CHECKING(for standalone modules build)
	AM_CONDITIONAL([STANDALONE_MODULES], [test "x${enable_standalone_modules}" == "xyes"])
	AC_MSG_RESULT([${enable_standalone_modules}])

	AC_ARG_ENABLE(
		[standalone-executables],																dnl feature
		AS_HELP_STRING([--enable-standalone-executables], [produce standalone executables]),	dnl help string
		[enable_standalone_executables=${enableval}],											dnl action if given
		[enable_standalone_executables=no]														dnl action if not given
	)
	AC_MSG_CHECKING(for standalone executables build)
	AM_CONDITIONAL([STANDALONE_EXECUTABLES], [test "x${enable_standalone_executables}" == "xyes"])
	AC_MSG_RESULT([${enable_standalone_executables}])

	AC_ARG_ENABLE(
		[standalone-libraries],																dnl feature
		AS_HELP_STRING([--enable-standalone-libraries], [produce standalone libraries]),	dnl help string
		[enable_standalone_libraries=${enableval}],											dnl action if given
		[enable_standalone_libraries=no]													dnl action if not given
	)
	AC_MSG_CHECKING(for standalone libraries build)
	AM_CONDITIONAL([STANDALONE_LIBRARIES], [test "x${enable_standalone_libraries}" == "xyes"])
	AC_MSG_RESULT([${enable_standalone_libraries}])
	dnl }}}
	dnl Version{{{
	AC_ARG_ENABLE(
		[version],													dnl feature
		AS_HELP_STRING([--enable-version], [produce libISSM.so.0]),	dnl help string
		[enable_version=${enableval}],								dnl action if given
		[enable_version=no]											dnl action if not given
	)
	AM_CONDITIONAL([VERSION], [test "x${enable_version}" == "xyes"])
	dnl }}}
	dnl Wrappers build {{{
	AC_ARG_WITH(
		[wrappers],															dnl feature
		AS_HELP_STRING([--with-wrappers = value], [wrappers compilation]),	dnl help string
		[WRAPPERS_VALUE=${withval}],										dnl action if given
		[WRAPPERS_VALUE="yes"]												dnl action if not given
	)
	AC_MSG_CHECKING(for wrappers compilation)
	AM_CONDITIONAL([WRAPPERS], [test "x${WRAPPERS_VALUE}" == "xyes"])
	AC_MSG_RESULT([${WRAPPERS_VALUE}])
	dnl }}}
	dnl Extensions{{{
	ISSMEXT=".exe"
	AC_SUBST([ISSMEXT])
	dnl }}}

	dnl OS{{{
	IS_MAC=no
	IS_MSYS2=no
	SYSTEM_FMEMOPEN=1
	AC_MSG_CHECKING([operating system type])
	case "${host_os}" in
		*darwin*)
			AC_MSG_RESULT([macOS])
			IS_MAC=yes
			AC_DEFINE([_IS_MAC_], [1], [is macOS])
			AC_DEFINE([_IS_MSYS2_], [0], [is Windows (MSYS2 MinGW)])
			dnl For some reason, CXXFLAGS is not empty by default under clang
			export CXXFLAGS="-g -O2 -fPIC -std=c++11 -D_DO_NOT_LOAD_GLOBALS_"

			dnl When standard Dakota installation has been updated to new 
			dnl version, remove the following
			DAKOTA_COMPILER_FLAGS="-Wno-deprecated-register -Wno-return-type"
			export CFLAGS="${DAKOTA_COMPILER_FLAGS}"
			export CXXFLAGS="${CXXFLAGS} ${DAKOTA_COMPILER_FLAGS}"

			dnl NOTE: Commenting out the following, for now, as ISSM seems to 
			dnl 	  compile and run fine, but certain errors (e.g. file not 
			dnl 	  found) were not bubbling up, and instead causing MATLAB 
			dnl 	  to crash.
			dnl
# 			if test "${LDFLAGS}" == ""; then
# 				export LDFLAGS="-Wl,-no_compact_unwind"
# 			else
# 				export LDFLAGS="${LDFLAGS} -Wl,-no_compact_unwind"
# 			fi

			dnl Check if system copy of libc has fmemopen
			AC_MSG_CHECKING([if system copy of libc has fmemopen (macOS-only check)])
			sys_ver=$(sw_vers -productVersion)
			if test $(echo ${sys_ver} | cut -d "." -f 1) -eq 10 && test $(echo ${sys_ver} | cut -d "." -f 2) -lt 13; then
				SYSTEM_FMEMOPEN=0
				AC_MSG_RESULT([no])
			else
				AC_MSG_RESULT([yes])
			fi
		;;
		*linux*)
			AC_MSG_RESULT([Linux])
			AC_DEFINE([_IS_MAC_], [0], [is macOS])
			AC_DEFINE([_IS_MSYS2_], [0], [is Windows (MSYS2 MinGW)])
		;;
		*mingw*)
			AC_MSG_RESULT([Windows (MSYS2 MinGW)])
			IS_MSYS2=yes
			AC_DEFINE([_IS_MAC_], [0], [is macOS])
			AC_DEFINE([_IS_MSYS2_], [1], [is Windows (MSYS2 MinGW)])
			export CXXFLAGS="-D_MSYS2_ -std=c++11"
			export LDFLAGS="-no-undefined"
			export OSLIBS="-Wl,-L/c/msys64/mingw64/lib -Wl,-lstdc++ -Wl,-lmingw32 -Wl,-lgcc_s -Wl,-lmoldname -Wl,-lmingwex -Wl,-lmsvcrt -Wl,-lm -Wl,-lpthread -Wl,-lshell32 -Wl,-luser32 -Wl,-lgdi32 -Wl,-luser32 -Wl,-ladvapi32 -Wl,-lkernel32 -Wl,-lgcc"
		;;
		*)
			AC_MSG_ERROR([unsupported operating system type)])
		;;
	esac

	AM_CONDITIONAL([MAC], [test "x${IS_MAC}" == "xyes"])
	AM_CONDITIONAL([MSYS2], [test "x${IS_MSYS2}" == "xyes"])

	AC_DEFINE_UNQUOTED([_SYSTEM_HAS_FMEMOPEN_], ${SYSTEM_FMEMOPEN}, [does system copy of libc have fmemopen])
	AM_CONDITIONAL([SYSTEM_HAS_FMEMOPEN], [test "${SYSTEM_FMEMOPEN}" == "1"])

	dnl Set default environment variables
	if test ! -z "${COPTFLAGS+x}"; then
		AC_MSG_WARN([If you want to use the optimization flags provided by COPTFLAGS (${COPTFLAGS}), please pass them via CFLAGS])
	fi
	if test -z "${CXXFLAGS+x}"; then
		export CXXFLAGS="-g -O2 -fPIC -std=c++11 -D_DO_NOT_LOAD_GLOBALS_"
	fi
	if test ! -z "${CXXOPTFLAGS+x}"; then
		AC_MSG_WARN([If you want to use the optimization flags provided by CXXOPTFLAGS (${CXXOPTFLAGS}), please pass them via CXXFLAGS])
	fi
	dnl }}}
	dnl Xlib (graphics library){{{
	AC_MSG_CHECKING([for Xlib (graphics library)])
	AC_ARG_WITH(
		[graphics-lib],
		AS_HELP_STRING([--with-graphics-lib=options], [Xlib (graphics library) to use]),
		[GRAPHICS_LIB=${withval}],
		[GRAPHICS_LIB=""]
	)
	if test -n "${GRAPHICS_LIB}"; then
		GRAPHICS_DIR=$(echo ${GRAPHICS_LIB} | sed -e "s/-L//g" | awk '{print $[1]}')
		if test -d "${GRAPHICS_DIR}" || test -f "${GRAPHICS_DIR}"; then
			HAVE_GRAPHICS=yes
			GRAPHICSLIB="${GRAPHICS_LIB}"
			AC_DEFINE([_HAVE_GRAPHICS_], [1], [with Xlib (graphics library) in ISSM src])
			AC_SUBST([GRAPHICSLIB])
		else
			if test -f "${PETSC_ROOT}/conf/petscvariables"; then
				PETSC_REC_GRAPHICS_LIB=$(cat ${PETSC_ROOT}/conf/petscvariables | grep X_LIB)
				AC_MSG_ERROR([Xlib (graphics library) provided (${GRAPHICS_LIB}) does not exist! PETSc suggests the following library: ${PETSC_REC_GRAPHICS_LIB}]);
			fi
			AC_MSG_ERROR([Xlib (graphics library) provided (${GRAPHICS_LIB}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([done])
	dnl }}}
	dnl MATLAB{{{
	dnl See if MATLAB has been provided
	AC_MSG_CHECKING([for MATLAB])
	AC_ARG_WITH(
		[matlab-dir],														dnl feature
		AS_HELP_STRING([--with-matlab-dir=DIR], [MATLAB root directory]),	dnl help string
		[MATLAB_ROOT=${withval}],											dnl action if given
		[MATLAB_ROOT="no"]													dnl action if not given
	)
	if test "x${MATLAB_ROOT}" == "xno"; then
		HAVE_MATLAB=no
	else
		if ! test -d "${MATLAB_ROOT}"; then
			AC_MSG_ERROR([MATLAB directory provided (${MATLAB_ROOT}) does not exist!]);
		fi
		if ! test -f "${MATLAB_ROOT}/extern/include/mex.h"; then
			AC_MSG_ERROR([Couldn't find mex.h... check your installation of MATLAB])
		fi
		HAVE_MATLAB=yes
	fi
	AC_MSG_RESULT([${HAVE_MATLAB}])
	AM_CONDITIONAL([MATLAB], [test "x${HAVE_MATLAB}" == "xyes"])

	dnl Set variables
	if test "x${HAVE_MATLAB}" == "xyes"; then
		AC_DEFINE([_HAVE_MATLAB_], [1], [with MATLAB in ISSM src])

		dnl Set MEX* variable
		AC_MSG_CHECKING([MATLAB's mex compilation flags])

		dnl NOTE: We know $VENDOR cannot be empty at this point, so no need to
		dnl		  check again in the following conditionals
		dnl
		case "${host_os}" in
			*mingw*)
				dnl Value to set MEXEXT to can be found on Windows by running $MATLAB_ROOT/bin/mexext.bat
				MEXEXT=".mexw64"
				MATLABINCL="-I${MATLAB_ROOT}/extern/include"
				MEXOPTFLAGS="-O2 -fwrapv -DNDEBUG -g"
				MEXCFLAGS="-fexceptions -fno-omit-frame-pointer -m64 -DMATLAB_MEX_FILE"
				MEXCXXFLAGS="-fexceptions -fno-omit-frame-pointer -std=c++11 -m64 -DMATLAB_MEX_FILE"
				MEXLINKFLAGS="-m64 -Wl,--no-undefined -shared -static -Wl,${MATLAB_ROOT}/extern/lib/win64/mingw64/mexFunction.def"
				MEXLIB_DIR="${MATLAB_ROOT}/extern/lib/win64/mingw64"
				MEXLIB="-L${MEXLIB_DIR} -lmx -lmex -lmat -lm -lmwlapack -lmwblas"
			;;
			*)
				MEXEXT=$(${MATLAB_ROOT}/bin/mex -v 2>&1 < /dev/null | grep LDEXTENSION | sed -e "s/         LDEXTENSION        = //g")
				MATLABINCL="-I${MATLAB_ROOT}/extern/include"
				MEXLINKFLAGS=$(${MATLAB_ROOT}/bin/mex -v 2>&1 < /dev/null | grep LDFLAGS | sed -e "s/         LDFLAGS            = //g")
				MEXLIB=$(${MATLAB_ROOT}/bin/mex -v 2>&1 < /dev/null | grep CXXLIBS | sed -e "s/         CXXLIBS            = //g")
				if test -z "${MEXEXT}"; then
					echo "#include <mex.h>" > conftest.cpp
					echo "void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){}" >> conftest.cpp
					${MATLAB_ROOT}/bin/mex -v -lmex conftest.cpp > conftest.tmp 2>&1
					rm -f conftest.cpp
					MEXLINKFLAGS=$(cat conftest.tmp | grep LDFLAGS | sed -e "s/LDFLAGS ://g")
					MEXLIB=$(cat conftest.tmp | grep LINKLIBS | sed -e "s/LINKLIBS ://g")
					MEXEXT=$(cat conftest.tmp | grep LDEXT | sed -e "s/LDEXT ://g" | awk '{print $[1]}')
					if test -z "${MEXEXT}"; then
					 cat conftest.tmp
					fi
					rm -f conftest.tmp
				fi

				dnl Make sure mexFunction.map is not in MEXLIB to avoid problems with global variables
				dnl MEXLINKFLAGS=$(echo ${MEXLINKFLAGS} | sed -e "s/,-expo.*mexFunction\\.map\"//g" | sed -e "s/-[[^ ]]*mexFunction\\.map//g")
				MEXLINKFLAGS="" dnl We actually don't need MEXLINK????

				dnl on new versions of macOS, MATLAB adds a -weak prefix to its library which makes libtool trip, remove
				MEXLIB=$(echo $MEXLIB | sed -e "s/weak-//g")
			;;
		esac
		AC_MSG_RESULT([done])
		if test -z "${MEXEXT}"; then
			AC_MSG_ERROR([Couldn't use MATLAB's mex... check manual compilation with MATLAB or error message above])
		fi

		AC_SUBST([MEXEXT])
		AC_SUBST([MEXOPTFLAGS])
		AC_SUBST([MEXCFLAGS])
		AC_SUBST([MEXCXXFLAGS])
		AC_SUBST([MATLABINCL])
		AC_SUBST([MEXLINKFLAGS])
		AC_SUBST([MEXLIB])
	fi
	dnl }}}
	dnl JavaScript{{{
	AC_MSG_CHECKING([for JavaScript])
	AC_ARG_WITH(
		[javascript],
		AS_HELP_STRING([--with-javascript], [compile JavaScript wrappers? (default: no)]),
		[JAVASCRIPT=${withval}],
		[JAVASCRIPT="no"]
	)
	if test "x${JAVASCRIPT}" == "xno"; then
		HAVE_JAVASCRIPT=no
	else
		HAVE_JAVASCRIPT=yes
		AC_DEFINE([_HAVE_JAVASCRIPT_], [1], [with JavaScript])
	fi
	AC_MSG_RESULT([${HAVE_JAVASCRIPT}])
	AM_CONDITIONAL([JAVASCRIPT], [test "x${HAVE_JAVASCRIPT}" == "xyes"])
	JAVASCRIPTWRAPPEREXT=.js
	AC_SUBST([JAVASCRIPTWRAPPEREXT])
	dnl }}}
	dnl Triangle {{{
	AC_MSG_CHECKING([for triangle])
	AC_ARG_WITH(
		[triangle-dir],
		AS_HELP_STRING([--with-triangle-dir=DIR], [Triangle root directory]),
		[TRIANGLE_ROOT=${withval}],
		[TRIANGLE_ROOT="no"]
	)
	if test "x${TRIANGLE_ROOT}" == "xno"; then
		HAVE_TRIANGLE=no
	else
		HAVE_TRIANGLE=yes
		if ! test -d "${TRIANGLE_ROOT}"; then
			AC_MSG_ERROR([Triangle directory provided (${TRIANGLE_ROOT}) does not exist!]);
		fi
		if ! test -f "${TRIANGLE_ROOT}/include/triangle.h"; then
			AC_MSG_ERROR([Couldn't find triangle.h... check your installation of triangle])
		fi
	fi
	AC_MSG_RESULT([${HAVE_TRIANGLE}])
	AM_CONDITIONAL([TRIANGLE], [test "x${HAVE_TRIANGLE}" == "xyes"])

	dnl Triangle libraries and header files
	if test "x${HAVE_TRIANGLE}" == "xyes"; then
		TRIANGLEINCL=-I${TRIANGLE_ROOT}/include
		case "${host_os}" in
			*darwin*)
				if test "x${HAVE_JAVASCRIPT}" == "xyes"; then
					dnl Link to the object file, not the library
					TRIANGLELIB=${TRIANGLE_ROOT}/share/triangle.o
				else
					TRIANGLELIB="-L${TRIANGLE_ROOT}/lib -ltriangle"
				fi
			;;
			*linux*)
				if test "x${HAVE_JAVASCRIPT}" == "xyes"; then
					dnl Link to the object file, not the library
					TRIANGLELIB=${TRIANGLE_ROOT}/share/triangle.o
				else
					TRIANGLELIB="-L${TRIANGLE_ROOT}/lib -ltriangle"
				fi
			;;
			*mingw*)
				if test "x${HAVE_JAVASCRIPT}" == "xyes"; then
					dnl Link to the object file, not the library
					TRIANGLELIB=${TRIANGLE_ROOT}/share/triangle.o
				else
					TRIANGLELIB="-Wl,-L${TRIANGLE_ROOT}/lib -Wl,-ltriangle"
				fi
			;;
		esac
		AC_DEFINE([_HAVE_TRIANGLE_], [1], [with Triangle in ISSM src])
		AC_SUBST([TRIANGLEINCL])
		AC_SUBST([TRIANGLELIB])
	fi
	dnl }}}
	dnl Boost{{{
	AC_MSG_CHECKING([for Boost])
	AC_ARG_WITH(
		[boost-dir],
		AS_HELP_STRING([--with-boost-dir=DIR], [Boost root directory]),
		[BOOST_ROOT=${withval}],
		[BOOST_ROOT="no"]
	)
	if test "x${BOOST_ROOT}" == "xno"; then
		HAVE_BOOST=no
	else
		HAVE_BOOST=yes
		if ! test -d "${BOOST_ROOT}"; then
			AC_MSG_ERROR([Boost directory provided (${BOOST_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_BOOST}])
	AM_CONDITIONAL([BOOST], [test "x${HAVE_BOOST}" == "xyes"])

	dnl Boost libraries and header files
	if test "x${HAVE_BOOST}" == "xyes"; then
		if test -z "${CXXFLAGS+x}"; then
			export CXXFLAGS="-Wno-deprecated"
		else
			export CXXFLAGS="${CXXFLAGS} -Wno-deprecated"
		fi
		BOOSTROOT="${BOOST_ROOT}"
		BOOSTINCL="-I${BOOST_ROOT}/include"
		#BOOSTLIB="-L$BOOST_ROOT/lib -lboost_python"
		AC_MSG_CHECKING(for Boost version)
		BOOST_VERSION=`cat ${BOOST_ROOT}/include/boost/version.hpp | grep "#define BOOST_VERSION " | sed 's/.*BOOST_VERSION //'`
		BOOST_VERSION_MAJOR=`expr ${BOOST_VERSION} / 100000`
		BOOST_VERSION_MINOR=`expr ${BOOST_VERSION} / 100 % 1000`
		AC_MSG_RESULT([${BOOST_VERSION_MAJOR}.${BOOST_VERSION_MINOR}])
		AC_DEFINE([_HAVE_BOOST_], [1], [with Boost in ISSM src])
		AC_SUBST([BOOSTROOT])
		AC_SUBST([BOOSTINCL])
		AC_SUBST([BOOSTLIB])
	fi
	dnl }}}
	dnl Dakota{{{
	AC_MSG_CHECKING([for Dakota])
	AC_ARG_WITH(
		[dakota-dir],
		AS_HELP_STRING([--with-dakota-dir=DIR], [Dakota root directory]),
		[DAKOTA_ROOT=${withval}],
		[DAKOTA_ROOT="no"]
	)
	if test "x${DAKOTA_ROOT}" == "xno"; then
		HAVE_DAKOTA=no
	else
		HAVE_DAKOTA=yes
		if ! test -d "${DAKOTA_ROOT}"; then
			AC_MSG_ERROR([Dakota directory provided (${DAKOTA_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_DAKOTA}])
	AM_CONDITIONAL([DAKOTA], [test "x${HAVE_DAKOTA}" == "xyes"])

	dnl Dakota libraries and header files
	if test "x${HAVE_DAKOTA}" == "xyes"; then
		DAKOTAINCL=-I${DAKOTA_ROOT}/include

		AC_MSG_CHECKING(for Dakota version)
		dnl TODO:
		dnl - Check if this method applies to all other versions of Dakota (it
		dnl   should as long as the Dakota binaries have been compiled). If so,
		dnl   we can remove the other methods of getting the version.
		dnl - Modify src/wrappers/IssmConfig/IssmConfig.cpp so that strlen is
		dnl   not called with _DAKOTA_VERSION_ as an argument so that we can
		dnl   do,
		dnl
		dnl   	AC_DEFINE_UNQUOTED([_DAKOTA_VERSION_], ${DAKOTA_VERSION}, [Dakota version number])
		dnl
		if test -f "${DAKOTA_ROOT}/VERSION"; then
			DAKOTA_VERSION=`cat ${DAKOTA_ROOT}/VERSION | grep 'DAKOTA Version' | sed 's/.*DAKOTA Version //' | sed 's/ .*//'`
		else
			DAKOTA_VERSION_OUTPUT=`${DAKOTA_ROOT}/bin/dakota -v`
			if test -n "${DAKOTA_VERSION_OUTPUT}"; then
				DAKOTA_VERSION=`echo ${DAKOTA_VERSION_OUTPUT} grep "Dakota version" | sed 's/Dakota version //' | sed 's/ .*//'`
			elif test -f "${DAKOTA_ROOT}/../src/src/CommandLineHandler.C"; then
				DAKOTA_VERSION=`cat ${DAKOTA_ROOT}/../src/src/CommandLineHandler.C | grep 'DAKOTA version' | grep 'release' | grep -v // | sed 's/.*DAKOTA version //' | sed 's/ .*//' `
			elif test -f "${DAKOTA_ROOT}/../src/src/CommandLineHandler.cpp"; then
				DAKOTA_VERSION=`cat ${DAKOTA_ROOT}/../src/src/CommandLineHandler.cpp | grep 'DAKOTA version' | grep 'release' | grep -v // | sed 's/.*DAKOTA version //' | sed 's/ .*//' `
			else
				AC_MSG_ERROR([Dakota CommandLineHandler.C or CommandLineHandler.cpp file not found to determine DAKOTA_VERSION!]);
			fi
		fi
		AC_MSG_RESULT([${DAKOTA_VERSION}])
		AC_DEFINE_UNQUOTED(_DAKOTA_VERSION_, "${DAKOTA_VERSION}", [Dakota version number])

		DAKOTAFLAGS=""

		dnl NOTE:
		dnl - See,
		dnl
		dnl 	$ISSM_DIR/externalpackages/dakota/build/src/Makefile.export.Dakota
		dnl
		dnl   for the flags needed by your combination of Boost and Dakota
		dnl   versions
		dnl - We know $DAKOTA_ROOT cannot be empty at this point, so no need to
		dnl   check again in the following conditionals
		dnl
		dnl TODO:
		dnl - Should we also be checking if HAVE_BOOST before adding boost libs?
		dnl - Clean up the following conditionals
		dnl
		case "${host_os}" in
			*darwin*)
				if test "${DAKOTA_VERSION}" == "5.1" || test "${DAKOTA_VERSION}" == "5.2"; then
					DAKOTALIB="-L${DAKOTA_ROOT}/lib -ldakota -lteuchos -lpecos -llhs -lsparsegrid -lsurfpack -lconmin -lddace -lfsudace -ljega -lcport -loptpp -lpsuade -lncsuopt -lcolin -linterfaces -lmomh -lscolib -lpebbl -ltinyxml -lutilib -l3po -lhopspack -lnidr -lamplsolver -lboost_signals -lboost_regex -lboost_filesystem -lboost_system -ldl"
				elif test "${DAKOTA_VERSION}" == "5.3" || test "${DAKOTA_VERSION}" == "5.3.1"; then
					DAKOTAFLAGS="-DDISABLE_DAKOTA_CONFIG_H -DBOOST_MULTI_INDEX_DISABLE_SERIALIZATION -DDAKOTA_PLUGIN -DBOOST_DISABLE_ASSERTS -DDAKOTA_HAVE_BOOST_FS -DHAVE_UNISTD_H -DHAVE_SYSTEM -DHAVE_WORKING_FORK -DHAVE_WORKING_VFORK -DHAVE_SYS_WAIT_H -DHAVE_USLEEP -DDAKOTA_F90 -DDAKOTA_HAVE_MPI -DHAVE_PECOS -DHAVE_SURFPACK -DDAKOTA_COLINY -DDAKOTA_UTILIB -DHAVE_ADAPTIVE_SAMPLING -DHAVE_CONMIN -DDAKOTA_DDACE -DHAVE_FSUDACE -DDAKOTA_HOPS -DHAVE_JEGA -DHAVE_NCSU -DHAVE_NL2SOL -DHAVE_OPTPP -DDAKOTA_OPTPP -DHAVE_PSUADE -DHAVE_AMPL"
					DAKOTALIB="-L${DAKOTA_ROOT}/lib -L${BOOST_ROOT}/lib -ldakota_src -lpecos_src -lscolib -ljega_fe -llhs -lpebbl -lcolin -linterfaces -lmods -lmoga -loptpp -lsampling -lsoga -lsurfpack -lutilib -lconmin -ldakota_src_fortran -lmod -lncsuopt -lsurfpack_fortran -lteuchos -l3po -lamplsolver -lanalyzer -lbose -lcport -ldace -ldfftpack -leutils -lfsudace -lhopspack -ljega -lnidr -lpecos -lpsuade -lrandom -ltinyxml -lutilities -lsparsegrid -lboost_signals -lboost_regex -lboost_filesystem -lboost_system"
					AC_DEFINE([DISABLE_DAKOTA_CONFIG_H], [1], [enabling DAKOTA_CONFIG_H])
					AC_DEFINE([DAKOTA_HAVE_MPI], [1], [enabling Dakota with MPI])
				elif test "${DAKOTA_VERSION}" == "6.1" || test "${DAKOTA_VERSION}" == "6.2"; then
					if test "${BOOST_VERSION_MAJOR}" == "1"; then
						DAKOTAFLAGS="-DHAVE_CONFIG_H -DDISABLE_DAKOTA_CONFIG_H -DBOOST_DISABLE_ASSERTS -DHAVE_UNISTD_H -DHAVE_SYSTEM -DHAVE_WORKING_FORK -DHAVE_WORKING_VFORK -DHAVE_SYS_WAIT_H -DHAVE_USLEEP -DDAKOTA_F90 -DDAKOTA_HAVE_MPI -DHAVE_PECOS -DHAVE_SURFPACK -DHAVE_ADAPTIVE_SAMPLING -DHAVE_ESM -DHAVE_QUESO -DHAVE_QUESO_GPMSA -DHAVE_CONMIN -DHAVE_DDACE -DHAVE_DREAM -DHAVE_FSUDACE -DDAKOTA_HOPS -DHAVE_NCSU -DHAVE_NL2SOL -DHAVE_NOMAD -DHAVE_OPTPP -DDAKOTA_OPTPP -DHAVE_PSUADE -DHAVE_AMPL"
						DAKOTALIB="-L${DAKOTA_ROOT}/lib -ldakota_src -ldakota_src_fortran -lnidr -lteuchos -lpecos -lpecos_src -llhs -llhs_mods -llhs_mod -ldfftpack -lsparsegrid -lsurfpack -lsurfpack -lsurfpack_fortran -lconmin -lddace -ldream -lfsudace -lhopspack -lncsuopt -lcport -lnomad -loptpp -lpsuade -lamplsolver"
						DAKOTALIB+=" -L${BOOST_ROOT}/lib -lboost_filesystem -lboost_program_options -lboost_regex -lboost_serialization -lboost_system"
						DAKOTALIB+=" ${BLASLAPACKLIB}"
					fi
					AC_DEFINE([DISABLE_DAKOTA_CONFIG_H], [1], [enabling DAKOTA_CONFIG_H])
					AC_DEFINE([DAKOTA_HAVE_MPI], [1], [enabling Dakota with MPI])
				elif test "${DAKOTA_VERSION}" == "6.11"; then
					if test "${BOOST_VERSION_MAJOR}" == "1"; then
						if test "${BOOST_VERSION_MINOR}" == "55"; then
							DAKOTAFLAGS="-DDISABLE_DAKOTA_CONFIG_H -DBOOST_MULTI_INDEX_DISABLE_SERIALIZATION -DDAKOTA_PLUGIN -DBOOST_DISABLE_ASSERTS -DDAKOTA_HAVE_BOOST_FS -DHAVE_UNISTD_H -DHAVE_SYSTEM -DHAVE_WORKING_FORK -DHAVE_WORKING_VFORK -DHAVE_SYS_WAIT_H -DHAVE_USLEEP -DDAKOTA_F90 -DDAKOTA_HAVE_MPI -DHAVE_PECOS -DHAVE_SURFPACK -DDAKOTA_UTILIB -DHAVE_ADAPTIVE_SAMPLING -DHAVE_CONMIN -DDAKOTA_DDACE -DHAVE_FSUDACE -DDAKOTA_HOPS -DHAVE_NCSU -DHAVE_NL2SOL -DHAVE_OPTPP -DDAKOTA_OPTPP -DHAVE_PSUADE -DHAVE_AMPL"
							DAKOTALIB="-L${DAKOTA_ROOT}/lib -ldakota_src -ldream -lfsudace -lddace -lnomad -lpecos_src -llhs -llhs_mods -loptpp -lsurfpack -lconmin -ldakota_src_fortran -llhs_mod -lncsuopt -lsurfpack_fortran -lteuchos -lamplsolver -lcport -ldfftpack -lfsudace -lhopspack -lnidr -lpecos -lpsuade -lsparsegrid -L$BOOST_ROOT/lib -lboost_serialization -lboost_signals -lboost_regex -lboost_filesystem -lboost_system ${BLASLAPACKLIB}"
						elif test "${BOOST_VERSION_MINOR}" == "72"; then
							DAKOTAFLAGS="-DHAVE_CONFIG_H -DDISABLE_DAKOTA_CONFIG_H -DBOOST_DISABLE_ASSERTS -DHAVE_UNISTD_H -DHAVE_SYSTEM -DHAVE_WORKING_FORK -DHAVE_WORKING_VFORK -DHAVE_SYS_WAIT_H -DHAVE_USLEEP -DDAKOTA_F90 -DDAKOTA_HAVE_MPI -DHAVE_PECOS -DHAVE_SURFPACK -DHAVE_ADAPTIVE_SAMPLING -DHAVE_ESM -DHAVE_CONMIN -DHAVE_DDACE -DHAVE_DREAM -DHAVE_FSUDACE -DDAKOTA_HOPS -DHAVE_NCSU -DHAVE_NL2SOL -DHAVE_NOMAD -DHAVE_OPTPP -DDAKOTA_OPTPP -DHAVE_PSUADE -DHAVE_AMPL"
							DAKOTALIB="-L${DAKOTA_ROOT}/lib -ldakota_src -ldakota_src_fortran -lnidr -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchosparser -lteuchoscore -lpecos_util -lpecos_src -llhs -llhs_mods -llhs_mod -ldfftpack -lsparsegrid -lsurfpack -lsurfpack -lsurfpack_fortran -lapproxnn -lconmin -lddace -ldream -lfsudace -lhopspack -lncsuopt -lcport -lnomad -loptpp -lpsuade -lamplsolver -L${BOOST_ROOT}/lib -lboost_filesystem -lboost_program_options -lboost_regex -lboost_serialization -lboost_system ${BLASLAPACKLIB}"
						fi
					fi
					AC_DEFINE([DISABLE_DAKOTA_CONFIG_H], [1], [enabling DAKOTA_CONFIG_H])
					AC_DEFINE([DAKOTA_HAVE_MPI], [1], [enabling Dakota with MPI])
				else
					AC_MSG_ERROR([Dakota version not found or version (${DAKOTA_VERSION}) not supported!]);
				fi
			;;
			*linux*)
				if test "${DAKOTA_VERSION}" == "5.1" || test "${DAKOTA_VERSION}" == "5.2"; then
					DAKOTALIB="-L${DAKOTA_ROOT}/lib -ldakota -lteuchos -lpecos -llhs -lsparsegrid -lsurfpack -lconmin -lddace -lfsudace -ljega -lcport -loptpp -lpsuade -lncsuopt -lcolin -linterfaces -lmomh -lscolib -lpebbl -ltinyxml -lutilib -l3po -lhopspack -lnidr -lamplsolver -lboost_signals -lboost_regex -lboost_filesystem -lboost_system -ldl"
				elif test "${DAKOTA_VERSION}" == "5.3" || test "${DAKOTA_VERSION}" == "5.3.1"; then
					DAKOTAFLAGS="-DDISABLE_DAKOTA_CONFIG_H -DBOOST_MULTI_INDEX_DISABLE_SERIALIZATION -DDAKOTA_PLUGIN -DBOOST_DISABLE_ASSERTS -DDAKOTA_HAVE_BOOST_FS -DHAVE_UNISTD_H -DHAVE_SYSTEM -DHAVE_WORKING_FORK -DHAVE_WORKING_VFORK -DHAVE_SYS_WAIT_H -DHAVE_USLEEP -DDAKOTA_F90 -DDAKOTA_HAVE_MPI -DHAVE_PECOS -DHAVE_SURFPACK -DDAKOTA_COLINY -DDAKOTA_UTILIB -DHAVE_ADAPTIVE_SAMPLING -DHAVE_CONMIN -DDAKOTA_DDACE -DHAVE_FSUDACE -DDAKOTA_HOPS -DHAVE_JEGA -DHAVE_NCSU -DHAVE_NL2SOL -DHAVE_OPTPP -DDAKOTA_OPTPP -DHAVE_PSUADE -DHAVE_AMPL"
					DAKOTALIB="-L${DAKOTA_ROOT}/lib -L${BOOST_ROOT}/lib -ldakota_src -lpecos_src -lscolib -ljega_fe -llhs -lpebbl -lcolin -linterfaces -lmods -lmoga -loptpp -lsampling -lsoga -lsurfpack -lutilib -lconmin -ldakota_src_fortran -lmod -lncsuopt -lsurfpack_fortran -lteuchos -l3po -lamplsolver -lanalyzer -lbose -lcport -ldace -ldfftpack -leutils -lfsudace -lhopspack -ljega -lnidr -lpecos -lpsuade -lrandom -ltinyxml -lutilities -lsparsegrid -lboost_signals -lboost_regex -lboost_filesystem -lboost_system"
					AC_DEFINE([DISABLE_DAKOTA_CONFIG_H], [1], [enabling DAKOTA_CONFIG_H])
					AC_DEFINE([DAKOTA_HAVE_MPI], [1], [enabling Dakota with MPI])
				elif test "${DAKOTA_VERSION}" == "6.1" || test "${DAKOTA_VERSION}" == "6.2"; then
					if test "${BOOST_VERSION_MAJOR}" == "1"; then
						DAKOTAFLAGS="-DHAVE_CONFIG_H -DDISABLE_DAKOTA_CONFIG_H -DBOOST_DISABLE_ASSERTS -DHAVE_UNISTD_H -DHAVE_SYSTEM -DHAVE_WORKING_FORK -DHAVE_WORKING_VFORK -DHAVE_SYS_WAIT_H -DHAVE_USLEEP -DDAKOTA_F90 -DDAKOTA_HAVE_MPI -DHAVE_PECOS -DHAVE_SURFPACK -DHAVE_ADAPTIVE_SAMPLING -DHAVE_ESM -DHAVE_QUESO -DHAVE_QUESO_GPMSA -DHAVE_CONMIN -DHAVE_DDACE -DHAVE_DREAM -DHAVE_FSUDACE -DDAKOTA_HOPS -DHAVE_NCSU -DHAVE_NL2SOL -DHAVE_NOMAD -DHAVE_OPTPP -DDAKOTA_OPTPP -DHAVE_PSUADE -DHAVE_AMPL"
						DAKOTALIB="-L${DAKOTA_ROOT}/lib -ldakota_src -ldakota_src_fortran -lnidr -lteuchos -lpecos -lpecos_src -llhs -llhs_mods -llhs_mod -ldfftpack -lsparsegrid -lsurfpack -lsurfpack -lsurfpack_fortran -lconmin -lddace -ldream -lfsudace -lhopspack -lncsuopt -lcport -lnomad -loptpp -lpsuade -lamplsolver"
						DAKOTALIB+=" -L${BOOST_ROOT}/lib -lboost_filesystem -lboost_program_options -lboost_regex -lboost_serialization -lboost_system"
						DAKOTALIB+=" ${BLASLAPACKLIB}"
					fi
					AC_DEFINE([DISABLE_DAKOTA_CONFIG_H], [1], [enabling DAKOTA_CONFIG_H])
					AC_DEFINE([DAKOTA_HAVE_MPI], [1], [enabling Dakota with MPI])
				elif test "${DAKOTA_VERSION}" == "6.11"; then
					if test "${BOOST_VERSION_MAJOR}" == "1"; then
						if test "${BOOST_VERSION_MINOR}" == "55"; then
							DAKOTAFLAGS="-DDISABLE_DAKOTA_CONFIG_H -DBOOST_MULTI_INDEX_DISABLE_SERIALIZATION -DDAKOTA_PLUGIN -DBOOST_DISABLE_ASSERTS -DDAKOTA_HAVE_BOOST_FS -DHAVE_UNISTD_H -DHAVE_SYSTEM -DHAVE_WORKING_FORK -DHAVE_WORKING_VFORK -DHAVE_SYS_WAIT_H -DHAVE_USLEEP -DDAKOTA_F90 -DDAKOTA_HAVE_MPI -DHAVE_PECOS -DHAVE_SURFPACK -DDAKOTA_UTILIB -DHAVE_ADAPTIVE_SAMPLING -DHAVE_CONMIN -DDAKOTA_DDACE -DHAVE_FSUDACE -DDAKOTA_HOPS -DHAVE_NCSU -DHAVE_NL2SOL -DHAVE_OPTPP -DDAKOTA_OPTPP -DHAVE_PSUADE -DHAVE_AMPL"
							DAKOTALIB="-L${DAKOTA_ROOT}/lib -ldakota_src -ldream -lfsudace -lddace -lnomad -lpecos_src -llhs -llhs_mods -loptpp -lsurfpack -lconmin -ldakota_src_fortran -llhs_mod -lncsuopt -lsurfpack_fortran -lteuchos -lamplsolver -lcport -ldfftpack -lfsudace -lhopspack -lnidr -lpecos -lpsuade -lsparsegrid -L$BOOST_ROOT/lib -lboost_serialization -lboost_signals -lboost_regex -lboost_filesystem -lboost_system ${BLASLAPACKLIB}"
						elif test "${BOOST_VERSION_MINOR}" == "72"; then
							DAKOTAFLAGS="-DHAVE_CONFIG_H -DHAVE_CONFIG_H -DDISABLE_DAKOTA_CONFIG_H -DBOOST_DISABLE_ASSERTS -DHAVE_UNISTD_H -DHAVE_SYSTEM -DHAVE_WORKING_FORK -DHAVE_WORKING_VFORK -DHAVE_SYS_WAIT_H -DHAVE_USLEEP -DDAKOTA_F90 -DDAKOTA_HAVE_MPI -DHAVE_PECOS -DHAVE_SURFPACK -DHAVE_ADAPTIVE_SAMPLING -DHAVE_ESM -DHAVE_CONMIN -DHAVE_DDACE -DHAVE_DREAM -DHAVE_FSUDACE -DDAKOTA_HOPS -DHAVE_NCSU -DHAVE_NL2SOL -DHAVE_NOMAD -DHAVE_OPTPP -DDAKOTA_OPTPP -DHAVE_PSUADE -DHAVE_AMPL"
							DAKOTALIB="-L${DAKOTA_ROOT}/lib -ldakota_src -ldakota_src_fortran -lnidr -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchosparser -lteuchoscore -lpecos_util -lpecos_src -llhs -llhs_mods -llhs_mod -ldfftpack -lsparsegrid -lsurfpack -lsurfpack -lsurfpack_fortran -lapproxnn -lconmin -lddace -ldream -lfsudace -lhopspack -lncsuopt -lcport -lnomad -loptpp -lpsuade -lamplsolver -L${BOOST_ROOT}/lib -lboost_filesystem -lboost_program_options -lboost_regex -lboost_serialization -lboost_system ${BLASLAPACKLIB}"
						fi
					fi
					AC_DEFINE([DISABLE_DAKOTA_CONFIG_H], [1], [enabling DAKOTA_CONFIG_H])
					AC_DEFINE([DAKOTA_HAVE_MPI], [1], [enabling Dakota with MPI])
				else
					AC_MSG_ERROR([Dakota version not found or version (${DAKOTA_VERSION}) not supported!]);
				fi
			;;
		esac

		case ${DAKOTA_VERSION} in
			@<:@1-9@:>@*.@<:@0-9@:>@*.@<:@0-9@:>@*)
				DAKOTA_MAJOR=`echo ${DAKOTA_VERSION} | sed -e 's/^\(@<:@0-9@:>@*\)\..*/\1/'`
				DAKOTA_MINOR=`echo ${DAKOTA_VERSION} | sed -e 's/^@<:@0-9@:>@*\.\(@<:@0-9@:>@*\)\..*/\1/'`
				DAKOTA_BUILD=`echo ${DAKOTA_VERSION} | sed -e 's/^@<:@0-9@:>@*\.@<:@0-9@:>@*\.\(@<:@0-9@:>@*\).*/\1/'`
			;;
			@<:@1-9@:>@*.@<:@0-9@:>@*)
				DAKOTA_MAJOR=`echo ${DAKOTA_VERSION} | sed -e 's/^\(@<:@0-9@:>@*\)\..*/\1/'`
				DAKOTA_MINOR=`echo ${DAKOTA_VERSION} | sed -e 's/^@<:@0-9@:>@*\.\(@<:@0-9@:>@*\).*/\1/'`
				DAKOTA_BUILD=0
			;;
			*)
				AC_MSG_ERROR([Dakota version (${DAKOTA_VERSION}) not supported!])
			;;
		esac
		AC_MSG_CHECKING(for Dakota major version)
		AC_MSG_RESULT(${DAKOTA_MAJOR})
		AC_DEFINE_UNQUOTED([_DAKOTA_MAJOR_], ${DAKOTA_MAJOR}, [Dakota major version number])
		AC_MSG_CHECKING(for Dakota minor version)
		AC_MSG_RESULT(${DAKOTA_MINOR})
		AC_DEFINE_UNQUOTED([_DAKOTA_MINOR_], ${DAKOTA_MINOR}, [Dakota minor version number])
		AC_MSG_CHECKING(for Dakota build version)
		AC_MSG_RESULT(${DAKOTA_BUILD})
		AC_DEFINE_UNQUOTED([_DAKOTA_BUILD_], ${DAKOTA_BUILD}, [Dakota build version number])

		AC_DEFINE([_HAVE_DAKOTA_], [1], [with Dakota in ISSM src])
		AC_SUBST([DAKOTAINCL])
		AC_SUBST([DAKOTAFLAGS])
		AC_SUBST([DAKOTALIB])
	fi
	AM_CONDITIONAL([ISSM_DAKOTA], [test "x${DAKOTA_MAJOR}" == "x6"])
	dnl }}}
	dnl Python{{{
	AC_MSG_CHECKING([for Python])
	AC_ARG_WITH(
		[python],
		AS_HELP_STRING([--with-python=EXEC], [Python path, e.g., "/usr/bin/python3"]),
		[PYTHON_PATH=${withval}],
		[PYTHON_PATH="no"]
	)

	if test "x${PYTHON_PATH}" == "xno"; then
		HAVE_PYTHON=no
	else
		HAVE_PYTHON=yes
		if ! test -f "${PYTHON_PATH}"; then
			AC_MSG_ERROR([Python provided (${PYTHON_PATH}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_PYTHON}])
	AM_CONDITIONAL([PYTHON], [test "x${HAVE_PYTHON}" == "xyes"])

	dnl Python specifics
	if test "x${HAVE_PYTHON}" == "xyes"; then

		AC_MSG_CHECKING([for Python version])
		PYTHON_VERSION=$(${PYTHON_PATH} -c "import sys; sys.stdout.write(str(sys.version_info.major)+'.'+str(sys.version_info.minor))")
		AC_MSG_RESULT([${PYTHON_VERSION}])
		
		dnl Make sure major version is 3
		PYTHON_MAJOR=${PYTHON_VERSION%.*}
		AC_DEFINE_UNQUOTED([_PYTHON_MAJOR_], ${PYTHON_MAJOR}, [Python version major])
		if test "x${PYTHON_MAJOR}" != "x3"; then
			AC_MSG_ERROR([Only Python 3 is supported]);
		fi

		AC_MSG_CHECKING([for Python include directory])
		PYTHONINCL=$(${PYTHON_PATH} -c "import sys; import sysconfig; sys.stdout.write(sysconfig.get_config_var('INCLUDEPY'))")
		if ! test -f "${PYTHONINCL}/Python.h"; then
			PYTHONINCL=$(${PYTHON_PATH} -c "from sysconfig import get_paths as gp; print(gp()[['include']])")
			if ! test -f "${PYTHONINCL}/Python.h"; then
				AC_MSG_ERROR([Python.h not found! Please locate this file and contact ISSM developers via forum or email.]);
			fi
		fi
		AC_MSG_RESULT([$PYTHONINCL])
		PYTHONINCL="-I${PYTHONINCL}"

		AC_MSG_CHECKING([for libpython])
		PYTHONLIBDIR=$(${PYTHON_PATH} -c "import sys; import sysconfig; sys.stdout.write(sysconfig.get_config_var('LIBDIR'))")
		if ls ${PYTHONLIBDIR}/libpython${PYTHON_VERSION}m.* 1> /dev/null 2>&1; then
			PYTHONLIB="-L${PYTHONLIBDIR} -lpython${PYTHON_VERSION}m"
		elif ls ${PYTHONLIBDIR}/libpython${PYTHON_VERSION}.* 1> /dev/null 2>&1; then
			PYTHONLIB="-L${PYTHONLIBDIR} -lpython${PYTHON_VERSION}"
		else
			PYTHONLIBDIR=$(${PYTHON_PATH} -c "from sysconfig import get_paths as gp; print(gp()[['stdlib']])")
			if ls ${PYTHONLIBDIR}/../libpython${PYTHON_VERSION}m.* 1> /dev/null 2>&1; then
				PYTHONLIB="-L${PYTHONLIBDIR}/.. -lpython${PYTHON_VERSION}m"
			elif ls ${PYTHONLIBDIR}/../libpython${PYTHON_VERSION}.* 1> /dev/null 2>&1; then
				PYTHONLIB="-L${PYTHONLIBDIR}/.. -lpython${PYTHON_VERSION}"
			else
				AC_MSG_ERROR([libpython not found! Please locate this file and contact ISSM developers via forum or email.]);
			fi
		fi
		AC_MSG_RESULT([$PYTHONLIB])

		case "${host_os}" in
			*darwin*) PYTHONLINK="-dynamiclib" ;;
			*linux*)  PYTHONLINK="-shared" ;;
			*mingw*)  PYTHONLINK="-shared" ;;
		esac
		AC_DEFINE([_HAVE_PYTHON_], [1], [with Python in ISSM src])
		AC_SUBST([PYTHONINCL])
		AC_SUBST([PYTHONLIB])
		PYTHONWRAPPEREXT=".so"
		AC_SUBST([PYTHONWRAPPEREXT])
		AC_SUBST([PYTHONLINK])

		dnl NumPy
		AC_MSG_CHECKING([for NumPy version])
		NUMPY_VERSION=$(${PYTHON_PATH} -c "import numpy; import sys; sys.stdout.write(numpy.version.version)")
		AC_MSG_RESULT([$NUMPY_VERSION])

		AC_MSG_CHECKING([for NumPy include directory])
		NUMPYINCL=$(${PYTHON_PATH} -c "import numpy; import sys; sys.stdout.write(numpy.get_include())")
		AC_MSG_RESULT([$NUMPYINCL])
		if ! test -d "${NUMPYINCL}"; then
			AC_MSG_ERROR([NumPy directory provided (${NUMPYINCL}) does not exist!]);
		fi

		dnl NumPy libraries and header files
		PYTHON_NUMPYINCL="-I${NUMPYINCL} -I${NUMPYINCL}/numpy"
		AC_DEFINE([_HAVE_PYTHON_NUMPY_], [1], [with NumPy in ISSM src])
		AC_SUBST([PYTHON_NUMPYINCL])
	fi
	AM_CONDITIONAL([PYTHON3], [test "xyes" == "xyes"])
	dnl }}}
	dnl Python-OLD{{{
	if test "x${HAVE_PYTHON}" != "xyes"; then
	AC_MSG_CHECKING([for Python])
	AC_ARG_WITH(
		[python-dir],
		AS_HELP_STRING([--with-python-dir=DIR], [Python root directory]),
		[PYTHON_ROOT=${withval}],
		[PYTHON_ROOT="no"]
	)

	AC_ARG_WITH(
		[python-version],
		AS_HELP_STRING([--with-python-version=DIR], [Python forced version]),
		[PYTHON_VERSION=${withval}],
		[PYTHON_VERSION="no"]
	)
	if test "x${PYTHON_ROOT}" == "xno"; then
		HAVE_PYTHON=no
		HAVE_PYTHON3=no
	else
		HAVE_PYTHON=yes
		if ! test -d "${PYTHON_ROOT}"; then
			AC_MSG_ERROR([Python directory provided (${PYTHON_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_PYTHON}])
	AM_CONDITIONAL([PYTHON], [test "x${HAVE_PYTHON}" == "xyes"])

	dnl Python specifics
	if test "x${HAVE_PYTHON}" == "xyes"; then
		if test "x${PYTHON_VERSION}" == "xno"; then
			AC_MSG_CHECKING([for Python version])
			dnl Query Python for its version number. Getting [:3] seems to be
			dnl the best way to do this: it's what "site.py" does in the
			dnl standard library.
			if test -f "${PYTHON_ROOT}/bin/python"; then
				PYTHON_VERSION=$(${PYTHON_ROOT}/bin/python -c "import sys; sys.stdout.write(str(sys.version_info.major)+'.'+str(sys.version_info.minor))")
			elif test -f "${PYTHON_ROOT}/bin/python3"; then
				PYTHON_VERSION=$(${PYTHON_ROOT}/bin/python3 -c "import sys; sys.stdout.write(str(sys.version_info.major)+'.'+str(sys.version_info.minor))")
			else
				AC_MSG_ERROR([Python version could not be determined automatically, please provide option --with-python-version]);
			fi
			AC_MSG_RESULT([${PYTHON_VERSION}])
		else
			AC_MSG_RESULT([enforced Python version is ${PYTHON_VERSION}])
		fi
		dnl Determine major version
		PYTHON_MAJOR=${PYTHON_VERSION%.*}
		AC_DEFINE_UNQUOTED([_PYTHON_MAJOR_], ${PYTHON_MAJOR}, [Python version major])
		if test "x${PYTHON_MAJOR}" == "x3"; then
			HAVE_PYTHON3="yes"
		else
			HAVE_PYTHON3="no"
		fi

		AC_MSG_CHECKING([for Python header file Python.h])
		dnl Python.h might be in different locations:
		if test -f "${PYTHON_ROOT}/include/Python.h"; then
			PYTHONINCL=-I${PYTHON_ROOT}/include
		elif test -f "${PYTHON_ROOT}/include/python${PYTHON_VERSION}/Python.h"; then
			PYTHONINCL=-I${PYTHON_ROOT}/include/python${PYTHON_VERSION}
		elif test -f "${PYTHON_ROOT}/include/python${PYTHON_VERSION}m/Python.h"; then
			PYTHONINCL=-I${PYTHON_ROOT}/include/python${PYTHON_VERSION}m
		elif test -f "${PYTHON_ROOT}/Headers/Python.h"; then
			PYTHONINCL=-I${PYTHON_ROOT}/include/python${PYTHON_VERSION}m
		elif test -f "${PYTHON_ROOT}/Frameworks/Python.framework/Versions/${PYTHON_VERSION}/include/python${PYTHON_VERSION}/Python.h"; then
			PYTHONINCL=-I${PYTHON_ROOT}/Frameworks/Python.framework/Versions/${PYTHON_VERSION}/include/python${PYTHON_VERSION}
		else
			AC_MSG_ERROR([Python.h not found! Please locate this file and contact ISSM developers via forum or email.]);
		fi
		AC_MSG_RESULT([found])

		AC_MSG_CHECKING([for Python library libpython])
		if ls ${PYTHON_ROOT}/lib/x86_64-linux-gnu/libpython${PYTHON_VERSION}m.* 1> /dev/null 2>&1; then
			PYTHONLIB="-L${PYTHON_ROOT}/lib/x86_64-linux-gnu -lpython${PYTHON_VERSION}m"
		elif ls ${PYTHON_ROOT}/lib/x86_64-linux-gnu/libpython${PYTHON_VERSION}.* 1> /dev/null 2>&1; then
			PYTHONLIB="-L${PYTHON_ROOT}/lib/x86_64-linux-gnu -lpython${PYTHON_VERSION}"
		elif ls ${PYTHON_ROOT}/lib/libpython${PYTHON_VERSION}m.* 1> /dev/null 2>&1; then
			PYTHONLIB="-L${PYTHON_ROOT}/lib -lpython${PYTHON_VERSION}m"
		elif ls ${PYTHON_ROOT}/lib/libpython${PYTHON_VERSION}.* 1> /dev/null 2>&1; then
			PYTHONLIB="-L${PYTHON_ROOT}/lib -lpython${PYTHON_VERSION}"
		elif ls ${PYTHON_ROOT}/lib64/libpython${PYTHON_VERSION}m.* 1> /dev/null 2>&1; then
			PYTHONLIB="-L${PYTHON_ROOT}/lib64 -lpython${PYTHON_VERSION}m"
		elif ls ${PYTHON_ROOT}/lib64/libpython${PYTHON_VERSION}.* 1> /dev/null 2>&1; then
			PYTHONLIB="-L${PYTHON_ROOT}/lib64 -lpython${PYTHON_VERSION}"
		elif ls ${PYTHON_ROOT}/Frameworks/Python.framework/Versions/${PYTHON_VERSION}/lib/libpython${PYTHON_VERSION}.* 1> /dev/null 2>&1; then
			PYTHONLIB="-L${PYTHON_ROOT}/Frameworks/Python.framework/Versions/${PYTHON_VERSION}/lib -lpython${PYTHON_VERSION}"
		else
			AC_MSG_ERROR([libpython not found! Please locate this file and contact ISSM developers via forum or email.]);
		fi
		AC_MSG_RESULT([found])

		PYWRAPPEREXT=.so
		case "${host_os}" in
			*darwin*)
				PYTHONLINK="-dynamiclib"
			;;
			*linux*)
				PYTHONLINK="-shared"
			;;
			*mingw*)
				PYTHONLINK="-shared"
			;;
		esac
		AC_DEFINE([_HAVE_PYTHON_], [1], [with Python in ISSM src])
		AC_SUBST([PYTHONINCL])
		AC_SUBST([PYTHONLIB])
		PYTHONWRAPPEREXT=${PYWRAPPEREXT}
		AC_SUBST([PYTHONWRAPPEREXT])
		AC_SUBST([PYTHONLINK])
	fi
	AM_CONDITIONAL([PYTHON3], [test "x${HAVE_PYTHON3}" == "xyes"])
	dnl }}}
	dnl NumPy{{{
	dnl NOTE: You can find NumPy by running,
	dnl
	dnl		>>> import numpy
	dnl		>>> numpy.__file__
	dnl
	dnl TODO:
	dnl - Replace references to python-numpy with numpy (and similar terms)
	dnl	  project-wide
	dnl
	AC_MSG_CHECKING(for python-numpy)
	AC_ARG_WITH(
		[python-numpy-dir],
		AS_HELP_STRING([--with-python-numpy-dir=DIR], [python-numpy root directory]),
		[PYTHON_NUMPY_ROOT=${withval}],
		[PYTHON_NUMPY_ROOT="no"]
	)
	if test "x${PYTHON_NUMPY_ROOT}" == "xno"; then
		HAVE_PYTHON_NUMPY=no
	else
		HAVE_PYTHON_NUMPY=yes
		if ! test -d "${PYTHON_NUMPY_ROOT}"; then
			AC_MSG_ERROR([NumPy directory provided (${PYTHON_NUMPY_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_PYTHON_NUMPY}])

	dnl NumPy libraries and header files
	if test "x${HAVE_PYTHON_NUMPY}" == "xyes"; then
		PYTHON_NUMPYINCL="-I${PYTHON_NUMPY_ROOT} -I${PYTHON_NUMPY_ROOT}/core/include/numpy -I${PYTHON_NUMPY_ROOT}/_core/include -I${PYTHON_NUMPY_ROOT}/_core/include/numpy"
		AC_DEFINE([_HAVE_PYTHON_NUMPY_], [1], [with NumPy in ISSM src])
		AC_SUBST([PYTHON_NUMPYINCL])
	fi
	fi #Starts in pythonb-old
	dnl }}}
	dnl Chaco{{{
	AC_MSG_CHECKING([for Chaco])
	AC_ARG_WITH(
		[chaco-dir],
		AS_HELP_STRING([--with-chaco-dir=DIR], [Chaco root directory]),
		[CHACO_ROOT=${withval}],
		[CHACO_ROOT="no"]
	)
	if test "x${CHACO_ROOT}" == "xno"; then
		HAVE_CHACO=no
	else
		HAVE_CHACO=yes
		if ! test -d "${CHACO_ROOT}"; then
			AC_MSG_ERROR([Chaco directory provided (${CHACO_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_CHACO}])
	AM_CONDITIONAL([CHACO], [test "x${HAVE_CHACO}" == "xyes"])

	dnl Chaco libraries and header files
	if test "x${HAVE_CHACO}" == "xyes"; then
		CHACOINCL="-I${CHACO_ROOT}/include"
		if test "x${IS_MSYS2}" == "xyes"; then
			CHACOLIB="-Wl,-L${CHACO_ROOT}/lib -Wl,-lchacominusblas"
		else
			CHACOLIB="-L${CHACO_ROOT}/lib -lchacominusblas"
		fi
		AC_DEFINE([_HAVE_CHACO_], [1], [with Chaco in ISSM src])
		AC_SUBST([CHACOINCL])
		AC_SUBST([CHACOLIB])
	fi
	dnl }}}
	dnl ESMF{{{
	AC_MSG_CHECKING([for ESMF])
	AC_ARG_WITH(
		[esmf-dir],
		AS_HELP_STRING([--with-esmf-dir=DIR], [ESMF root directory]),
		[ESMF_ROOT=${withval}],
		[ESMF_ROOT="no"]
	)
	if test "x${ESMF_ROOT}" == "xno"; then
		HAVE_ESMF=no
	else
		HAVE_ESMF=yes
		if ! test -d "${ESMF_ROOT}"; then
			AC_MSG_ERROR([ESMF directory provided (${ESMF_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_ESMF}])

	dnl ESMF libraries and header files
	if test "x${HAVE_ESMF}" == "xyes"; then
		ESMFINCL="-I${ESMF_ROOT}/include"
		ESMFLIB="-L${ESMF_ROOT}/lib/libO/Linux.gfortran.64.mpich.default/ -lesmf"
		AC_DEFINE([_HAVE_ESMF_], [1], [with ESMF in ISSM src])
		AC_SUBST([ESMFINCL])
		AC_SUBST([ESMFLIB])
	fi
	AM_CONDITIONAL([ESMF], [test "x${HAVE_ESMF}" == "xyes"])
	dnl }}}
	dnl CoDiPack{{{
	AC_MSG_CHECKING([for CoDiPack])
	AC_ARG_WITH(
		[codipack-dir],
		AS_HELP_STRING([--with-codipack-dir=DIR], [CoDiPack root directory]),
		[CODIPACK_ROOT=${withval}],
		[CODIPACK_ROOT="no"]
	)
	if test "x${CODIPACK_ROOT}" == "xno"; then
		HAVE_CODIPACK=no
	else
		HAVE_CODIPACK=yes
		if ! test -d "${CODIPACK_ROOT}"; then
			AC_MSG_ERROR([CoDiPack directory provided (${CODIPACK_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_CODIPACK}])

	dnl CoDiPack libraries and header files
	if test "x${HAVE_CODIPACK}" == "xyes"; then

		AC_MSG_CHECKING(for CoDiPack version)
		CODIPACK_MAJOR=`cat ${CODIPACK_ROOT}/include/codi.hpp | grep "#define CODI_MAJOR_VERSION" | sed 's/#define CODI_MAJOR_VERSION//' | sed 's/ //g'`
		CODIPACK_MINOR=`cat ${CODIPACK_ROOT}/include/codi.hpp | grep "#define CODI_MINOR_VERSION" | sed 's/#define CODI_MINOR_VERSION//' | sed 's/ //g'`
		if test -z "${CODIPACK_MAJOR}"; then
			 CODIPACK_MAJOR=`cat ${CODIPACK_ROOT}/include/codi/configure.h | grep "#define CODI_MAJOR_VERSION" | sed 's/#define CODI_MAJOR_VERSION//' | sed 's/ //g'`
			 CODIPACK_MINOR=`cat ${CODIPACK_ROOT}/include/codi/configure.h | grep "#define CODI_MINOR_VERSION" | sed 's/#define CODI_MINOR_VERSION//' | sed 's/ //g'`
		fi
		if test -z "${CODIPACK_MAJOR}"; then
			AC_MSG_ERROR([Couldn't determine CoDiPack version])
		fi
		AC_DEFINE_UNQUOTED([_CODIPACK_MAJOR_], ${CODIPACK_MAJOR}, [CoDiPack version major])
		AC_DEFINE_UNQUOTED([_CODIPACK_MINOR_], ${CODIPACK_MINOR}, [CoDiPack version minor])
		AC_MSG_RESULT([${CODIPACK_MAJOR}.${CODIPACK_MINOR}])

		CODIPACKINCL="-I${CODIPACK_ROOT}/include"
		AC_DEFINE([_HAVE_CODIPACK_], [1], [with CoDiPack in ISSM src])
		AC_DEFINE([_HAVE_AD_], [1], [with AD in ISSM src])
		AC_SUBST([CODIPACKINCL])
	fi
	AM_CONDITIONAL([CODIPACK], [test "x${HAVE_CODIPACK}" == "xyes"])
	AM_COND_IF(CODIPACK, [CXXFLAGS+=" -std=c++17"])
	dnl }}}
	dnl Tape Allocation {{{
	AC_MSG_CHECKING(for tape allocation)
	AC_ARG_ENABLE(
		[tape-alloc],																dnl feature
		AS_HELP_STRING([--enable-tape-alloc], [turn tape allocation support on]),
		[enable_tape_alloc=${enableval}],
		[enable_tape_alloc=no]
	)
	if test "x${enable_tape_alloc}" == "xyes"; then
		AC_DEFINE([_AD_TAPE_ALLOC_], [1], [enable a priori tape allocation for AD])
	fi
	AC_MSG_RESULT([${enable_tape_alloc}])
	dnl }}}
	dnl ADOL-C {{{
	AC_MSG_CHECKING([for ADOL-C])
	AC_ARG_WITH(
		[adolc-dir],
		AS_HELP_STRING([--with-adolc-dir=DIR], [ADOL-C root directory]),
		[ADOLC_ROOT=${withval}],
		[ADOLC_ROOT="no"]
	)
	if test "x${ADOLC_ROOT}" == "xno"; then
		HAVE_ADOLC=no
	else
		HAVE_ADOLC=yes
		if ! test -d "${ADOLC_ROOT}"; then
			AC_MSG_ERROR([ADOL-C directory provided (${ADOLC_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_ADOLC}])

	dnl ADOL-C libraries and header files
	if test "x${HAVE_ADOLC}" == "xyes"; then
		ADOLCINCL="-I${ADOLC_ROOT}/include"
		dnl ADOLCLIB="-L${ADOLC_ROOT}/lib64 -ladolc" used to be the path
		ADOLCLIB="-L${ADOLC_ROOT}/lib -ladolc"
		AC_DEFINE([_HAVE_ADOLC_], [1], [with ADOL-C in ISSM src])
		AC_DEFINE([_HAVE_AD_], [1], [with AD in ISSM src])
		AC_SUBST([ADOLCINCL])
		AC_SUBST([ADOLCLIB])
	fi
	AM_CONDITIONAL([ADOLC], [test "x${HAVE_ADOLC}" == "xyes"])
	AM_COND_IF(ADOLC, [CXXFLAGS+=" -std=c++11"])
	dnl }}}
	dnl ADOL-C version{{{
	AC_MSG_CHECKING(for ADOL-C version)
	AC_ARG_WITH(
		[adolc-version],
		AS_HELP_STRING([--with-adolc-version=number], [ADOL-C version]),
		[ADOLC_VERSION=${withval}],
		[ADOLC_VERSION=2]
	)
	AC_DEFINE_UNQUOTED([_ADOLC_VERSION_], ${ADOLC_VERSION}, [ADOL-C version])
	AC_MSG_RESULT(${ADOLC_VERSION})
	dnl }}}
	dnl ATLAS {{{
	AC_MSG_CHECKING(for ATLAS and CBLAS libraries)
	AC_ARG_WITH(
		[atlas-dir],
		AS_HELP_STRING([--with-atlas-dir=DIR], [ATLAS root directory]),
		[ATLAS_ROOT=${withval}],
		[ATLAS_ROOT="no"]
	)
	if test "x${ATLAS_ROOT}" == "xno"; then
		HAVE_ATLAS=no
	else
		HAVE_ATLAS=yes
		if ! test -d "${ATLAS_ROOT}"; then
			AC_MSG_ERROR([ATLAS directory provided (${ATLAS_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_ATLAS}])

	dnl ATLAS libraries and header files
	if test "x${HAVE_ATLAS}" == "xyes"; then
		case "${host_os}" in
			*darwin*)
				ATLASLIB="-L${ATLAS_ROOT}/lib -lcblas -latlas -lm"
			;;
			*linux*)
				ATLASLIB="-L${ATLAS_ROOT}/lib -lcblas -latlas -lm"
			;;
			*mingw*)
				ATLASLIB="-L${ATLAS_ROOT}/lib -lcblas -latlas -lm"
			;;
		esac
		AC_DEFINE([_HAVE_ATLAS_], [1], [with ATLAS in ISSM src])
		AC_SUBST([ATLASLIB])
	fi
	dnl }}}
	dnl GSL{{{
	AC_MSG_CHECKING([for GSL])
	AC_ARG_WITH(
		[gsl-dir],
		AS_HELP_STRING([--with-gsl-dir=DIR], [GSL root directory]),
		[GSL_ROOT=${withval}],
		[GSL_ROOT="no"]
	)
	if test "x${GSL_ROOT}" == "xno"; then
		HAVE_GSL=no
	else
		HAVE_GSL=yes
		if ! test -d "${GSL_ROOT}"; then
			AC_MSG_ERROR([GSL directory provided (${GSL_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_GSL}])

	dnl GSL libraries and header files
	if test "x${HAVE_GSL}" == "xyes"; then
		GSLINCL="-I${GSL_ROOT}/include"
		if test "x${HAVE_ATLAS}" == "xyes"; then
			GSLLIB="-dy -L${GSL_ROOT}/lib -lgsl -L${ATLAS_ROOT}/lib -lcblas -latlas -lm"
		else
			GSLLIB="-L${GSL_ROOT}/lib -lgsl -lgslcblas -lm"
		fi
		AC_DEFINE([_HAVE_GSL_], [1], [with GSL in ISSM src])
		AC_SUBST([GSLINCL])
		AC_SUBST([GSLLIB])
	fi
	AM_CONDITIONAL([GSL], [test "x${HAVE_GSL}" == "xyes"])
	dnl }}}
	dnl AMPI (ADOL-C){{{
	AC_MSG_CHECKING([for AMPI])
	AC_ARG_WITH(
		[ampi-dir],
		AS_HELP_STRING([--with-ampi-dir=DIR], [Adjoinable MPI root directory]),
		[AMPI_ROOT=${withval}],
		[AMPI_ROOT="no"]
	)
	if test "x${AMPI_ROOT}" == "xno"; then
		HAVE_AMPI=no
	else
		HAVE_AMPI=yes
		if ! test -d "${AMPI_ROOT}"; then
			AC_MSG_ERROR([AMPI directory provided (${AMPI_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_AMPI}])

	dnl AMPI libraries and header files
	if test "x${HAVE_AMPI}" == "xyes"; then
		AMPIINCL="-I${AMPI_ROOT}/include"
		if test "x${ADOLC_ROOT}" == "xno"; then
			AC_MSG_ERROR([cannot run AMPI without ADOL-C]);
		fi
		dnl AMPILIB="-dy -L${AMPI_ROOT}/lib -lampiCommon -L${ADOLC_ROOT}/lib -ladolc -L${AMPI_ROOT}/lib -lampiCommon -lampiBookkeeping -lampiTape"
		dnl AMPILIB="-dy -L${AMPI_ROOT}/lib  -L${ADOLC_ROOT}/lib -Wl,--start-group,-lampiCommon,-ladolc,-lampiCommon,-lampiBookkeeping,-lampiTape,-lampiPlainC,-lampiADtoolStubsST,--end-group"
		dnl AMPILIB="-L${AMPI_ROOT}/lib  -L${ADOLC_ROOT}/lib -Wl,--start-group -lampiCommon -ladolc -lampiCommon -lampiBookkeeping -lampiTape -lampiPlainC -lampiADtoolStubsST -Wl,--end-group"
		dnl AMPILIB="${AMPI_ROOT}/lib/libampiCommon.so ${ADOLC_ROOT}/lib/libadolc.so  ${AMPI_ROOT}/lib/libampiCommon.so ${AMPI_ROOT}/lib/libampiBookkeeping.so ${AMPI_ROOT}/lib/libampiTape.so ${AMPI_ROOT}/lib/libampiPlainC.so  ${AMPI_ROOT}/lib/libampiADtoolStubsST.so"
		dnl AMPILIB="-dy -L${AMPI_ROOT}/lib  -L${ADOLC_ROOT}/lib -lampiCommon -ladolc -lampiCommon -lampiBookkeeping -lampiTape -lampiPlainC -lampiADtoolStubsST"
		AMPILIB="-dy -L${AMPI_ROOT}/lib  -lampiCommon -lampiBookkeeping -lampiTape"
		AC_DEFINE([_HAVE_AMPI_], [1], [with AMPI in ISSM src])
		AC_SUBST([AMPIINCL])
		AC_SUBST([AMPILIB])
	fi
	AM_CONDITIONAL([AMPI], [test "x${HAVE_AMPI}" == "xyes"])
	dnl }}}
	dnl MeDiPack (CoDiPack, ADOL-C dev){{{
	AC_MSG_CHECKING([for MeDiPack])
	AC_ARG_WITH(
		[medipack-dir],
		AS_HELP_STRING([--with-medipack-dir=DIR], [MeDiPack root directory]),
		[MEDIPACK_ROOT=${withval}],
		[MEDIPACK_ROOT="no"]
	)
	if test "x${MEDIPACK_ROOT}" == "xno"; then
		HAVE_MEDIPACK=no
	else
		HAVE_MEDIPACK=yes
		if ! test -d "${MEDIPACK_ROOT}"; then
			AC_MSG_ERROR([MeDiPack directory provided (${MEDIPACK_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_MEDIPACK}])

	dnl MeDiPack libraries and header files
	if test "x${HAVE_MEDIPACK}" == "xyes"; then
		if test "x${CODIPACK_ROOT}" == "xno"; then
			AC_MSG_ERROR([cannot run MeDiPack without CoDiPack]);
		fi
		MEDIPACKINCL="-I${MEDIPACK_ROOT}/include -I${MEDIPACK_ROOT}/src"
		dnl Also set _HAVE_AMPI_, because the interface is (almost) the same as
		dnl for AMPI
		AC_DEFINE([_HAVE_AMPI_], [1], [with AMPI in ISSM src])
		AC_DEFINE([_HAVE_MEDIPACK_], [1], [with MeDiPack in ISSM src])
		AC_SUBST([MEDIPACKINCL])
	fi
	AM_CONDITIONAL([MEDIPACK], [test "x${HAVE_MEDIPACK}" == "xyes"])
	dnl }}}
	dnl AdjointPETSc{{{
	AC_MSG_CHECKING([for AdjointPETSc])
	AC_ARG_WITH(
		[adjointpetsc-dir],
		AS_HELP_STRING([--with-adjointpetsc-dir=DIR], [AdjointPETSc root directory]),
		[ADJOINTPETSC_ROOT=${withval}],
		[ADJOINTPETSC_ROOT="no"]
	)
	if test "x${ADJOINTPETSC_ROOT}" == "xno"; then
		HAVE_ADJOINTPETSC=no
	else
		HAVE_ADJOINTPETSC=yes
		if ! test -d "${ADJOINTPETSC_ROOT}"; then
			AC_MSG_ERROR([AdjointPETSc directory provided (${ADJOINTPETSC_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_ADJOINTPETSC}])

	dnl AdjointPETSc libraries and header files
	if test "x${HAVE_ADJOINTPETSC}" == "xyes"; then
		if test "x${CODIPACK_ROOT}" == "xno"; then
			AC_MSG_ERROR([cannot run AdjointPETSc without CoDiPack]);
		fi
		if test "x${PETSC_ROOT}" == "xno"; then
			AC_MSG_ERROR([cannot run AdjointPETSc without PETSc]);
		fi
		ADJOINTPETSCINCL="-I${ADJOINTPETSC_ROOT}/include"
		ADJOINTPETSCLIB="-L${ADJOINTPETSC_ROOT}/lib -ladjoint_petsc"
		dnl Also set _HAVE_AMPI_, because the interface is (almost) the same as
		dnl for AMPI
		AC_DEFINE([_HAVE_AMPI_], [1], [with AMPI in ISSM src])
		AC_DEFINE([_HAVE_ADJOINTPETSC_], [1], [with AdjointPETSc in ISSM src])
		AC_SUBST([ADJOINTPETSCINCL])
		AC_SUBST([ADJOINTPETSCLIB])
	fi
	AM_CONDITIONAL([ADJOINTPETSC], [test "x${HAVE_ADJOINTPETSC}" == "xyes"])
	dnl }}}
	dnl HDF5 {{{
	AC_MSG_CHECKING(for HDF5 libraries)
	AC_ARG_WITH(
		[hdf5-dir],
		AS_HELP_STRING([--with-hdf5-dir=DIR], [HDF5 root directory]),
		[HDF5_ROOT=${withval}],
		[HDF5_ROOT="no"]
	)
	if test "x${HDF5_ROOT}" == "xno"; then
		HAVE_HDF5=no
	else
		HAVE_HDF5=yes
		if ! test -d "${HDF5_ROOT}"; then
			AC_MSG_ERROR([HDF5 directory provided (${HDF5_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_HDF5}])

	dnl HDF5 libraries and header files
	dnl
	dnl TODO: Add check for if we need to link to libhdf5_fortran and
	dnl libhdf5hl_fortran, if and when necessary.
	dnl
	if test "x${HAVE_HDF5}" == "xyes"; then
		case "${host_os}" in
			*darwin*)
				HDF5LIB="-L${HDF5_ROOT}/lib -lhdf5 -lhdf5_hl"
			;;
			*linux*)
				HDF5LIB="-L${HDF5_ROOT}/lib -lhdf5 -lhdf5_hl"
			;;
			*mingw*)
				HDF5LIB="-L${HDF5_ROOT}/lib -lhdf5 -lhdf5_hl"
			;;
		esac
		AC_DEFINE([_HAVE_HDF5_], [1], [with HDF5 in ISSM src])
		AC_SUBST([HDF5LIB])
	fi
	dnl }}}
	dnl PETSc{{{
	AC_MSG_CHECKING([for PETSc])
	AC_ARG_WITH(
		[petsc-dir],
		AS_HELP_STRING([--with-petsc-dir=DIR], [PETSc root directory, necessary for parallel build]),
		[PETSC_ROOT=${withval}],
		[PETSC_ROOT="no"]
	)
	if test "x${PETSC_ROOT}" == "xno"; then
		HAVE_PETSC=no
	else
		HAVE_PETSC=yes
		if ! test -d "${PETSC_ROOT}"; then
			AC_MSG_ERROR([PETSc directory provided (${PETSC_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_PETSC}])
	AM_CONDITIONAL([PETSC], [test "x${HAVE_PETSC}" == "xyes"])

	dnl PETSc libraries and header files
	if test "x${HAVE_PETSC}" == "xyes"; then
		if ! test -f "${PETSC_ROOT}/include/petscversion.h"; then
			AC_MSG_ERROR([PETSc not instaled correctly: file (${PETSC_ROOT}/include/petscversion.h) does not exist!]);
		fi

		AC_MSG_CHECKING(for PETSc version)
		PETSC_MAJOR=`cat ${PETSC_ROOT}/include/petscversion.h | grep "#define PETSC_VERSION_MAJOR" | sed 's/#define PETSC_VERSION_MAJOR//' | sed 's/ //g'`
		PETSC_MINOR=`cat ${PETSC_ROOT}/include/petscversion.h | grep "#define PETSC_VERSION_MINOR" | sed 's/#define PETSC_VERSION_MINOR//' | sed 's/ //g'`
		AC_DEFINE_UNQUOTED([_PETSC_MAJOR_], ${PETSC_MAJOR}, [PETSc version major])
		AC_DEFINE_UNQUOTED([_PETSC_MINOR_], ${PETSC_MINOR}, [PETSc version minor])
		AC_MSG_RESULT([${PETSC_MAJOR}.${PETSC_MINOR}])

		AC_MSG_CHECKING(whether PETSc is the development version)
		PETSC_RELEASE=`cat ${PETSC_ROOT}/include/petscversion.h | grep "#define PETSC_VERSION_RELEASE" | sed 's/#define PETSC_VERSION_RELEASE//' | sed 's/ //g'`
		if test "${PETSC_RELEASE}" == "0"; then
			AC_DEFINE([_HAVE_PETSCDEV_], [1], [with PETSc-dev])
			AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
		fi

		AC_ARG_WITH(
			[petsc-arch],
			AS_HELP_STRING([--with-petsc-arch=DIR], [PETSc arch, necessary for PETSc < 3.0]),
			[PETSC_ARCH=${withval}],
			[PETSC_ARCH=""]
		)

		AC_MSG_CHECKING(for PETSc libraries and header files in ${PETSC_ROOT})
		dnl To get PETSc's libraries,
		dnl
		dnl		cd $ISSM_DIR/externalpackages/petsc/src
		dnl		make getlinklibs
		dnl
		PETSCINCL=" -I${PETSC_ROOT}/include"
		dnl Add other location (not needed anymore since at least PETSc 3.0)
		if test -n "${PETSC_ARCH}" && test -d "${PETSC_ROOT}/${PETSC_ARCH}/include"; then
			PETSCINCL+=" ${PETSC_ROOT}/${PETSC_ARCH}/include"
		fi
		if test -n "${PETSC_ARCH}" && test -d "${PETSC_ROOT}/include/${PETSC_ARCH}"; then
			PETSCINCL+=" ${PETSC_ROOT}/include/${PETSC_ARCH}"
		fi

		case "${host_os}" in
			*darwin*)
				if test ${PETSC_MAJOR} -lt 3; then
					PETSCLIB="-L${PETSC_ROOT}/lib -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetscsnes -lpetscts -lpetsc"
				else
					PETSCLIB="-L${PETSC_ROOT}/lib -lpetsc"
					#if test ${PETSC_MAJOR} -gt 3 || test ${PETSC_MINOR} -ge 3; then
					#	PETSCLIB+=" -lmetis"
					#fi
				fi
			;;
			*linux*)
				if test ${PETSC_MAJOR} -lt 3; then
					PETSCLIB="-L${PETSC_ROOT}/lib -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetscsnes -lpetscts -lmpiuni -lpetsc"
				else
					PETSCLIB="-L${PETSC_ROOT}/lib -lpetsc -ldl"
				fi
				if test "x$host_os_version" = "x3.0.101-0.31.1_1.0502.8394-cray_gem_s"; then
					PETSCLIB="-L${PETSC_ROOT}/lib -lcraypetsc_gnu_real"
				fi
				if test x$HOST = "xmaui01"; then
					PETSCLIB="-L${PETSC_ROOT}/lib -lcraypetsc_intel_real"
				fi
			;;
			*mingw*)
				PETSCLIB="-Wl,-L${PETSC_ROOT}/lib -Wl,-lpetsc"
			;;
		esac
		AC_MSG_RESULT([done])
		AC_DEFINE([_HAVE_PETSC_], [1], [with PETSc in ISSM src])
		AC_SUBST([PETSCINCL])
		AC_SUBST([PETSCLIB])
	fi
	dnl }}}
	dnl MPI{{{
	AC_MSG_CHECKING(for MPI)
	AC_ARG_WITH(
		[mpi-include],
		AS_HELP_STRING([--with-mpi-include=DIR], [MPI include directory, necessary for parallel build]),
		[MPI_INCLUDE=${withval}],
		[MPI_INCLUDE=""]
	)
	AC_ARG_WITH(
		[mpi-libdir],
		AS_HELP_STRING([--with-mpi-libdir=DIR], [MPI library directory, necessary for parallel build]),
		[MPI_LIBDIR=${withval}],
		[MPI_LIBDIR=""]
	)
	AC_ARG_WITH(
		[mpi-libflags],
		AS_HELP_STRING([--with-mpi-libflags=LIBS], [MPI libraries to be used, necessary for parallel build]),
		[MPI_LIBFLAGS=${withval}],
		[MPI_LIBFLAGS=""]
	)
	if test -z "${MPI_INCLUDE}"; then
		HAVE_MPI=no
	else
		HAVE_MPI=yes

		dnl Processing for Windows
		dnl
		dnl NOTE: We know $VENDOR cannot be empty at this point, so no need to
		dnl		  check again in the following conditionals
		dnl
		if test "${VENDOR}" == "intel-win7-32"; then
			MPI_LIBDIR=`cygpath -m ${MPI_LIBDIR}`
			MPI_INCLUDE=`cygpath -m ${MPI_INCLUDE}`
		elif test "${VENDOR}" == "intel-win7-64"; then
			MPI_LIBDIR="/I`cygpath -m ${MPI_LIBDIR}`"
			MPI_INCLUDE=`cygpath -m ${MPI_INCLUDE}`
		elif test "${VENDOR}" == "MSVC-Win64" || test "${VENDOR}" == "MSVC-Win64-par"; then
			MPI_LIBDIR=`cygpath -m ${MPI_LIBDIR}`
			MPI_INCLUDE=`cygpath -m ${MPI_INCLUDE}`
		fi

		if test -z "${MPI_LIBDIR}"; then
			MPILIB="${MPI_LIBFLAGS}"
		else
			MPILIB="${MPI_LIBDIR} ${MPI_LIBFLAGS}"
		fi

		if ! test -f ${MPI_INCLUDE}/mpi.h; then
			 AC_MSG_ERROR([Count not find mpi.h in ${MPI_INCLUDE}!]);
		fi
		MPIINCL="-I${MPI_INCLUDE}"

		AC_DEFINE([_HAVE_MPI_], [1], [with MPI in ISSM src])
		AC_DEFINE([HAVE_MPI], [1], [MPI flag for Dakota (DO NOT REMOVE)])
		AC_SUBST([MPIINCL])
		AC_SUBST([MPILIB])
	fi
	AM_CONDITIONAL([MPI], [test "x${HAVE_MPI}" == "xyes"])
	AC_MSG_RESULT([${HAVE_MPI}])
	dnl }}}
	dnl METIS{{{
	AC_MSG_CHECKING([for METIS])
	AC_ARG_WITH(
		[metis-dir],
		AS_HELP_STRING([--with-metis-dir=DIR], [METIS root directory, necessary for serial build]),
		[METIS_ROOT=${withval}],
		[METIS_ROOT="no"]
	)
	HAVE_METIS=no
	if test "x${METIS_ROOT}" == "xno"; then
		dnl Check if METIS was installed via PETSc
		if test -f ${PETSC_ROOT}/install/include/metis.h; then
			HAVE_METIS=yes
			METIS_ROOT="${PETSC_ROOT}"
		fi
	else
		if ! test -d "${METIS_ROOT}"; then
			AC_MSG_ERROR([METIS directory provided (${METIS_ROOT}) does not exist!]);
		fi
		HAVE_METIS=yes
	fi
	if test "${HAVE_METIS}" = "yes"; then
		if test -f ${METIS_ROOT}/include/metis.h; then
			 METIS_H=${METIS_ROOT}/include/metis.h
		elif test -f ${METIS_ROOT}/metis.h; then
			 METIS_H=${METIS_ROOT}/metis.h
		elif test -f ${METIS_ROOT}/../include/metis.h; then
			 METIS_H=${METIS_ROOT}/../include/metis.h
		else
			 AC_MSG_ERROR([Count not find METIS header file!]);
		fi
		METIS_VERSION=$(grep "#define METIS_VER_MAJOR" ${METIS_H} | sed 's|.*METIS_VER_MAJOR[[:space:]]*||')
		dnl METIS libraries and header files
		if test "x${METIS_VERSION}" == "x4"; then
			METISINCL="-I${METIS_ROOT}/Lib"
			case "${host_os}" in
				*darwin*)
					METISLIB="-L${METIS_ROOT} -lmetis"
				;;
				*linux*)
					METISLIB="-L${METIS_ROOT} -lmetis"
				;;
				*mingw*)
					METISLIB="-Wl,-L${METIS_ROOT}/lib -Wl,-lmetis"
				;;
			esac
		elif test "x${METIS_VERSION}" == "x5"; then
			METISINCL="-I${METIS_ROOT}/include"
			case "${host_os}" in
				*darwin*)
					METISLIB="-L${METIS_ROOT}/lib -lmetis"
				;;
				*linux*)
					METISLIB="-L${METIS_ROOT}/lib -lmetis"
				;;
				*mingw*)
					METISLIB="-Wl,-L${METIS_ROOT}/lib -Wl,-lmetis"
				;;
			esac
		else
			AC_MSG_ERROR([METIS version ${METIS_VERSION} not yet supported! Please contact ISSM developers via forum or email.])
		fi
		AC_DEFINE([_HAVE_METIS_], [1], [with METIS in ISSM src])
		AC_DEFINE_UNQUOTED([_METIS_VERSION_], ${METIS_VERSION}, [METIS version number])
		AC_SUBST([METISINCL])
		AC_SUBST([METISLIB])
	fi
	AC_MSG_RESULT([${HAVE_METIS}])
	AM_CONDITIONAL([METIS], [test "x${HAVE_METIS}" = "xyes"])
	dnl }}}
	dnl ParMETIS{{{
	AC_MSG_CHECKING([for ParMETIS])
	AC_ARG_WITH(
		[parmetis-dir],
		AS_HELP_STRING([--with-parmetis-dir=DIR], [ParMETIS root directory, necessary for parallel build]),
		[PARMETIS_ROOT=${withval}],
		[PARMETIS_ROOT="no"]
	)
	HAVE_PARMETIS=no
	if test "x${PARMETIS_ROOT}" == "xno"; then
		dnl Check if ParMETIS was installed via PETSc
		if test -f ${PETSC_ROOT}/install/include/parmetis.h; then
			HAVE_PARMETIS="yes"
			PARMETIS_ROOT="${PETSC_ROOT}"
		fi
	else
		if ! test -d "${PARMETIS_ROOT}"; then
			AC_MSG_ERROR([ParMETIS directory provided (${PARMETIS_ROOT}) does not exist!]);
		fi
		if ! test -d "${METIS_ROOT}"; then
			AC_MSG_ERROR([If supplying path to ParMETIS with option --with-parmetis-dir, must also supply path to METIS with option --with-metis-dir]);
		fi
		HAVE_PARMETIS="yes"
	fi
	if test "${HAVE_PARMETIS}" == "yes"; then
		#PARMETIS_VERSION=$(grep "#define PARMETIS_MAJOR_VERSION" ${PARMETIS_ROOT}/include/parmetis.h | sed 's|.*PARMETIS_MAJOR_VERSION[[:space:]]*||')
		dnl METIS libraries and header files
		#if test "x${PARMETIS_VERSION}" == "x4"; then
			PARMETISINCL="-I${PARMETIS_ROOT}/include"
			case "${host_os}" in
				*darwin*)
					PARMETISLIB="-L${PARMETIS_ROOT}/lib -lparmetis"
				;;
				*linux*)
					PARMETISLIB="-L${PARMETIS_ROOT}/lib -lparmetis"
				;;
				*mingw*)
					PARMETISLIB="-Wl,-L${PARMETIS_ROOT}/lib -Wl,-lparmetis"
				;;
			esac
		#else
		#	AC_MSG_ERROR([ParMETIS version ${PARMETIS_VERSION} not yet supported! Please contact ISSM developers via forum or email.])
		#fi
		AC_DEFINE([_HAVE_PARMETIS_], [1], [with ParMETIS in ISSM src])
		#AC_DEFINE([_PARMETIS_VERSION_], [${PARMETIS_VERSION}], [ParMETIS version number])
		AC_SUBST([PARMETISINCL])
		AC_SUBST([PARMETISLIB])
	fi
	AC_MSG_RESULT([${HAVE_PARMETIS}])
	AM_CONDITIONAL([PARMETIS], [test "x${HAVE_PARMETIS}" = "xyes"])
	dnl }}}
	dnl Toolkit for Advanced Optimization (TAO){{{
	AC_MSG_CHECKING([for TAO])
	AC_ARG_WITH(
		[tao-dir],
		AS_HELP_STRING([--with-tao-dir=DIR], [TAO root directory]),
		[TAO_ROOT=${withval}],
		[TAO_ROOT="no"]
	)
	if test "x${HAVE_PETSC}" == "xyes" && test "x${PETSC_MAJOR}" == "x3" && test ${PETSC_MINOR} -ge 5; then
		dnl In PETSc >= 3.5, TAO is provided
		HAVE_TAO="yes"
		AC_DEFINE([_HAVE_TAO_], [1], [with TAO in ISSM src])
		AC_MSG_RESULT([${HAVE_TAO}])
	else
		if test "x${TAO_ROOT}" == "xno"; then
			HAVE_TAO=no
		else
			HAVE_TAO=yes
			if ! test -d "${TAO_ROOT}"; then
				AC_MSG_ERROR([TAO directory provided (${TAO_ROOT}) does not exist!]);
			fi
		fi
		AC_MSG_RESULT([${HAVE_TAO}])

		dnl TAO libraries and header files
		if test "x${HAVE_TAO}" == "xyes"; then
			TAOINCL="-I${TAO_ROOT} -I${TAO_ROOT}/include -I${TAO_ROOT}/bmake"
			TAOLIB="-L${TAO_ROOT}/lib -ltao -lpetsc"
			AC_DEFINE([_HAVE_TAO_], [1], [with Tao in ISSM src])
			AC_SUBST([TAOINCL])
			AC_SUBST([TAOLIB])
		fi
	fi
	dnl }}}
	dnl M1QN3{{{
	AC_MSG_CHECKING([for M1QN3])
	AC_ARG_WITH(
		[m1qn3-dir],
		AS_HELP_STRING([--with-m1qn3-dir=DIR], [M1QN3 root directory]),
		[M1QN3_ROOT=${withval}],
		[M1QN3_ROOT="no"]
	)
	if test "x${M1QN3_ROOT}" == "xno"; then
		HAVE_M1QN3=no
	else
		HAVE_M1QN3=yes
		if ! test -d "${M1QN3_ROOT}"; then
			AC_MSG_ERROR([M1QN3 directory provided (${M1QN3_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_M1QN3}])

	dnl M1QN3 libraries and header files
	if test "x${HAVE_M1QN3}" == "xyes"; then
		if test "x${IS_MSYS2}" == "xyes"; then
			M1QN3LIB="-Wl,-L${M1QN3_ROOT} -Wl,-lm1qn3 -Wl,-lddot"
		else
			M1QN3LIB="-L${M1QN3_ROOT} -lm1qn3 -lddot"
		fi
		AC_DEFINE([_HAVE_M1QN3_], [1], [with M1QN3 in ISSM src])
		AC_SUBST([M1QN3LIB])
	fi
	dnl }}}
	dnl PROJ{{{
	AC_MSG_CHECKING([for PROJ])
	AC_ARG_WITH(
		[proj-dir],
		AS_HELP_STRING([--with-proj-dir=DIR], [PROJ root directory]),
		[PROJ_ROOT=${withval}],
		[PROJ_ROOT="no"]
	)
	if test "x${PROJ_ROOT}" == "xno"; then
		HAVE_PROJ=no
	else
		HAVE_PROJ=yes
		if ! test -d "${PROJ_ROOT}"; then
			AC_MSG_ERROR([PROJ directory provided (${PROJ_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_PROJ}])
	AM_CONDITIONAL([PROJ], [test "x${HAVE_PROJ}" == "xyes"])

	dnl PROJ libraries and header files
	if test "x${HAVE_PROJ}" == "xyes"; then
		PROJROOT="${PROJ_ROOT}"
		PROJINCL="-I${PROJ_ROOT}/include"
		PROJLIB="-L${PROJ_ROOT}/lib -lproj"
		AC_DEFINE([_HAVE_PROJ_], [1], [with PROJ in ISSM src])
		AC_SUBST([PROJROOT])
		AC_SUBST([PROJINCL])
		AC_SUBST([PROJLIB])
	fi
	dnl }}}
	dnl ScaLAPACK{{{
	dnl NOTE: User should supply path to root directory or libraries, but not both
	dnl
	AC_MSG_CHECKING([for ScaLAPACK])
	AC_ARG_WITH(
		[scalapack-dir],
		AS_HELP_STRING([--with-scalapack-dir=DIR], [ScaLAPACK root directory]),
		[SCALAPACK_ROOT=${withval}],
		[SCALAPACK_ROOT="no"]
	)
	AC_ARG_WITH(
		[scalapack-lib],
		AS_HELP_STRING([--with-scalapack-lib=LIBS], [ScaLAPACK libraries to link to]),
		[SCALAPACKLIB=${withval}],
		[SCALAPACKLIB="no"]
	)
	if test "x${SCALAPACK_ROOT}" == "xno" && test "x${SCALAPACKLIB}" == "xno"; then
		HAVE_SCALAPACK=no
		SCALAPACKLIB=""
	elif test "x${SCALAPACK_ROOT}" != "xno"; then
		if ! test -d "${SCALAPACK_ROOT}"; then
			AC_MSG_ERROR([ScaLAPACK directory provided (${SCALAPACK_ROOT}) does not exist!]);
		fi
		HAVE_SCALAPACK=yes
		if test "x${IS_MSYS2}" == "xyes"; then
			SCALAPACKLIB="-Wl,-L${SCALAPACK_ROOT}/lib -Wl,-lscalapack"
		else
			if test -f ${SCALAPACK_ROOT}/libscalapack-openmpi.so; then
				SCALAPACKLIB="-L${SCALAPACK_ROOT} -lscalapack-openmpi"
			else
			 SCALAPACKLIB="-L${SCALAPACK_ROOT}/lib -lscalapack"
			fi
		fi
	elif test "x${SCALAPACKLIB}" != "xno"; then
		dnl Value of SCALAPACKLIB should be valid here, so no need to set it (as above)
		HAVE_SCALAPACK=yes
	else
		AC_MSG_ERROR([use --with-scalapack-dir or --with-scalapack-lib, but not both])
	fi
	AC_MSG_RESULT([${HAVE_SCALAPACK}])

	dnl ScaLAPACK libraries and header files
	if test "x${HAVE_SCALAPACK}" == "xyes"; then
		AC_DEFINE([_HAVE_SCALAPACK_], [1], [with ScaLAPACK in ISSM src])
		AC_SUBST([SCALAPACKLIB])
	fi
	dnl }}}
	dnl BLAS/LAPACK{{{
	AC_MSG_CHECKING([for BLAS/LAPACK])
	AC_ARG_WITH(
		[blas-dir],
		[AS_HELP_STRING([--with-blas-dir=DIR], [BLAS root directory])],
		[BLAS_ROOT=$withval],
		[BLAS_ROOT="no"]
	)
	AC_ARG_WITH(
		[lapack-dir],
		[AS_HELP_STRING([--with-lapack-dir=DIR], [LAPACK root directory])],
		[LAPACK_ROOT=$withval],
		[LAPACK_ROOT="no"]
	)
	AC_ARG_WITH(
		[blas-lapack-dir],
		AS_HELP_STRING([--with-blas-lapack-dir=DIR], [BLAS/LAPACK root directory]),
		[BLASLAPACK_ROOT=$withval],
		[BLASLAPACK_ROOT="no"]
	)
	AC_ARG_WITH(
		[blas-lapack-lib],
		AS_HELP_STRING([--with-blas-lapack-lib=LIBFLAGS], [BLAS/LAPACK libflags]),
		[BLASLAPACK_LIB=$withval],
		[BLASLAPACK_LIB="no"]
	)
	if (test "x${BLAS_ROOT}" = "xno" || test "x${LAPACK_ROOT}" = "xno") && test "x${BLASLAPACK_ROOT}" = "xno" && test "x${BLASLAPACK_LIB}" = "xno"; then
		HAVE_BLASLAPACK=no
	else
		HAVE_BLASLAPACK=yes
		if test -z "${BLASLAPACK_LIB}"; then
			if ! test -d "${BLAS_ROOT}" || ! test -d "${LAPACK_ROOT}"; then
				if ! test -d "${BLASLAPACK_ROOT}"; then
					AC_MSG_ERROR([Use either --with-blas-dir and --with-lapack-dir *or* --with-blaslapack-dir *or* --with-blaslapack-lib and supply libflags]);
				fi
			fi
		fi
	fi
	AC_MSG_RESULT([${HAVE_BLASLAPACK}])

	dnl BLAS/LAPACK libraries and header files
	if test "x${HAVE_BLASLAPACK}" == "xyes"; then
		if test "x${BLASLAPACK_LIB}" != "xno"; then
			BLASLAPACKLIB="${BLASLAPACK_LIB}"
		else
			case "${host_os}" in
				*darwin*)
					BLASLAPACKLIB="-L${BLASLAPACK_ROOT}/lib"
					if ls ${BLASLAPACK_ROOT}/lib/libopenblas.* 1> /dev/null 2>&1; then
						BLASLAPACKLIB+=" -lopenblas"
					elif ls ${BLASLAPACK_ROOT}/lib/libf2clapack.* 1> /dev/null 2>&1; then
						BLASLAPACKLIB+=" -lf2clapack -lf2cblas"
					elif ls ${BLASLAPACK_ROOT}/lib/libflapack.* 1> /dev/null 2>&1; then
						BLASLAPACKLIB+=" -lflapack -lfblas"
					else
						BLASLAPACKLIB+=" -llapack -lblas"
					fi
				;;
				*linux*)
					BLASLAPACKLIB="-L${BLASLAPACK_ROOT}/lib"
					if ls ${BLASLAPACK_ROOT}/lib/libopenblas.* 1> /dev/null 2>&1; then
						BLASLAPACKLIB+=" -lopenblas"
					elif ls ${BLASLAPACK_ROOT}/lib/libf2clapack.* 1> /dev/null 2>&1; then
						BLASLAPACKLIB+=" -lf2clapack -lf2cblas"
					elif ls ${BLASLAPACK_ROOT}/lib/libflapack.* 1> /dev/null 2>&1; then
						BLASLAPACKLIB+=" -lflapack -lfblas"
					else
						BLASLAPACKLIB+=" -llapack -lblas"
					fi
				;;
				*mingw*)
					if test -d "${BLASLAPACK_ROOT}"; then
						BLASLAPACKLIB="-Wl,-L${BLASLAPACK_ROOT}/lib"
						if ls ${BLASLAPACK_ROOT}/lib/libopenblas.* 1> /dev/null 2>&1; then
							BLASLAPACKLIB+=" -lopenblas"
						elif ls ${BLASLAPACK_ROOT}/lib/libf2clapack.* 1> /dev/null 2>&1; then
							BLASLAPACKLIB+=" -lf2clapack -lf2cblas"
						elif ls ${BLASLAPACK_ROOT}/lib/libflapack.* 1> /dev/null 2>&1; then
							BLASLAPACKLIB="-Wl,-L${BLASLAPACK_ROOT}/lib -Wl,-lflapack -Wl,-lfblas"
						else
							BLASLAPACKLIB+=" -Wl,-llapack -Wl,-lblas"
						fi
					else
						BLASLAPACKLIB="${LAPACK_ROOT}/lib/liblapack.a ${BLAS_ROOT}/lib/libblas.a"
					fi
				;;
			esac
		fi
		AC_DEFINE([_HAVE_BLASLAPACK_], [1], [with BLAS/LAPACK in ISSM src])
		AC_SUBST([BLASLAPACKLIB])
	fi
	dnl }}}
	dnl Math Kernel Library (MKL){{{
	AC_MSG_CHECKING([for MKL])
	AC_ARG_WITH(
		[mkl-libflags],
		AS_HELP_STRING([--with-mkl-libflags=LIBS], [MKL libraries to be used]),
		[MKL_LIBFLAGS=${withval}],
		[MKL_LIBFLAGS="no"]
	)
	if test "x${MKL_LIBFLAGS}" == "xno"; then
		HAVE_MKL=no
	else
		HAVE_MKL=yes
		MKLLIB="${MKL_LIBFLAGS}"
		AC_DEFINE([_HAVE_MKL_], [1], [with MKL in ISSM src])
		AC_SUBST([MKLLIB])
		AC_SUBST([MKLINCL])
	fi
	AC_MSG_RESULT([${HAVE_MKL}])
	dnl }}}
	dnl PlaLAPACK{{{
	dnl TODO: 	Handle user supplying path to root directory *or* individual
	dnl 		arguments (like ScaLAPACK)
	dnl
	AC_MSG_CHECKING(for PlaLAPACK)
	AC_ARG_WITH(
		[plapack-lib],
		AS_HELP_STRING([--with-plapack-lib=<LIB>], [PlaLAPACK library]),
		[PLAPACK_LIB=${withval}],
		[PLAPACK_LIB=""]
	)
	AC_ARG_WITH(
		[plapack-include],
		AS_HELP_STRING([--with-plapack-include=<INC>], [PlaLAPACK include]),
		[PLAPACK_INCLUDE=${withval}],
		[PLAPACK_INCLUDE=""]
	)

	if test -n "${PLAPACK_LIB}"; then
		if test -n "${PLAPACK_INCLUDE}"; then
			HAVE_PLAPACK=yes
			PLAPACKINCL="${PLAPACK_INCLUDE}"
			PLAPACKLIB="${PLAPACK_LIB}"
			AC_DEFINE([_HAVE_PLAPACK_], [1], [with PlaLAPACK in ISSM src])
			AC_SUBST([PLAPACKINCL])
			AC_SUBST([PLAPACKLIB])
		else
			HAVE_PLAPACK=no
		fi
	else
		HAVE_PLAPACK=no
	fi
	AC_MSG_RESULT([${HAVE_PLAPACK}])
	dnl }}}
	dnl MPLAPACK{{{
	AC_MSG_CHECKING([for MPLAPACK])
	AC_ARG_WITH(
		[mplapack-dir],
		AS_HELP_STRING([--with-mplapack-dir=DIR], [MPLAPACK root directory]),
		[MPLAPACK_ROOT=${withval}],
		[MPLAPACK_ROOT="no"]
	)
	if test "x${MPLAPACK_ROOT}" == "xno"; then
		HAVE_MPLAPACK=no
	else
		HAVE_MPLAPACK=yes
		if ! test -d "${MPLAPACK_ROOT}"; then
			AC_MSG_ERROR([MPLAPACK directory provided (${MPLAPACK_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_MPLAPACK}])

	dnl MPLAPACK libraries and header files
	if test "x${HAVE_MPLAPACK}" == "xyes"; then
		MPLAPACKINCL="-I${MPLAPACK_ROOT}/include/mplapack -I${MPLAPACK_ROOT}/include"
		MPLAPACKLIB="-L${MPLAPACK_ROOT}/lib -lmpblas__Float128 -lmplapack__Float128 -lgomp -lquadmath"
		AC_DEFINE([_HAVE_MPLAPACK_], [1], [with MPLAPACK in ISSM src])
		AC_SUBST([MPLAPACKINCL])
		AC_SUBST([MPLAPACKLIB])
	fi
	AM_CONDITIONAL([MPLAPACK], [test "x${HAVE_MPLAPACK}" == "xyes"])
	dnl }}}
	dnl MUMPS{{{
	AC_MSG_CHECKING([for MUMPS])
	AC_ARG_WITH(
		[mumps-dir],
		AS_HELP_STRING([--with-mumps-dir=DIR], [MUMPS root directory]),
		[MUMPS_ROOT=${withval}],
		[MUMPS_ROOT="no"]
	)
	if test "x${MUMPS_ROOT}" == "xno"; then
		HAVE_MUMPS=no
	else
		HAVE_MUMPS=yes
		if ! test -d "${MUMPS_ROOT}"; then
			AC_MSG_ERROR([MUMPS directory provided (${MUMPS_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_MUMPS}])

	dnl MUMPS libraries and header files
	if test "x${HAVE_MUMPS}" == "xyes"; then
		MUMPSINCL="-I${MUMPS_ROOT}/include"
		if test "x${MUMPS_ROOT}" == "x${PETSC_ROOT}"; then
			if test "x${IS_MSYS2}" == "xyes"; then
				MUMPSLIB="-Wl,-L${MUMPS_ROOT}/lib -Wl,-lcmumps -Wl,-ldmumps -Wl,-lsmumps -Wl,-lzmumps -Wl,-lmumps_common -Wl,-lpord"
			else
				MUMPSLIB="-L${MUMPS_ROOT}/lib -ldmumps -lcmumps -lmumps_common -lpord -lzmumps"
			fi
		else
			MUMPSLIB="-L${MUMPS_ROOT}/lib -ldmumps -lmumps_common -lpord"
		fi
		AC_DEFINE([_HAVE_MUMPS_], [1], [with MUMPS in ISSM src])
		AC_SUBST([MUMPSINCL])
		AC_SUBST([MUMPSLIB])
	fi
	AM_CONDITIONAL([MUMPS], [test "x${HAVE_MUMPS}" == "xyes"])
	dnl }}}
	dnl MUMPS2{{{
	if test "x${HAVE_MUMPS}" != "xyes"; then
		AC_MSG_CHECKING(for MUMPS2 (standalone))
		AC_ARG_WITH(
			[mumps2-include],
			AS_HELP_STRING([--with-mumps2-include=DIR], [MUMPS2 include directory, necessary for parallel build]),
			[MUMPS_INCLUDE=${withval}],
			[MUMPS_INCLUDE=""]
		)
		AC_ARG_WITH(
			[mumps2-libflags],
			AS_HELP_STRING([--with-mumps2-libflags=LIBS], [MUMPS2 libraries to be used, necessary for parallel build]),
			[MUMPS_LIBFLAGS=${withval}],
			[MUMPS_LIBFLAGS=""]
		)
		if test -z "${MUMPS_INCLUDE}"; then
			HAVE_MUMPS=no
		else
			HAVE_MUMPS=yes

			if test -z "${MUMPS_LIBDIR}"; then
				MUMPSINCL="-I${MUMPS_INCLUDE}"
				MUMPSLIB="${MUMPS_LIBFLAGS}"
			else
				MUMPSINCL="-I${MUMPS_INCLUDE}"
				MUMPSLIB="-L${MUMPS_LIBDIR} ${MUMPS_LIBFLAGS}"
			fi
			AC_DEFINE([_HAVE_MUMPS_], [1], [with MUMPS])
			AC_SUBST([MUMPSINCL])
			AC_SUBST([MUMPSLIB])
		fi
		AM_CONDITIONAL([MUMPS], [test "x${HAVE_MUMPS}" == "xyes"])
		AC_MSG_RESULT([${HAVE_MUMPS}])
	fi
	dnl }}}
	dnl BLACS{{{
	AC_MSG_CHECKING([for BLACS])
	AC_ARG_WITH(
		[blacs-dir],
		AS_HELP_STRING([--with-blacs-dir=DIR], [BLACS root directory]),
		[BLACS_ROOT=${withval}],
		[BLACS_ROOT="no"]
	)
	if test "x${BLACS_ROOT}" == "xno"; then
		HAVE_BLACS=no
	else
		HAVE_BLACS=yes
		if ! test -d "${BLACS_ROOT}"; then
			AC_MSG_ERROR([BLACS directory provided (${BLACS_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_BLACS}])

	dnl BLACS libraries and header files
	if test "x${HAVE_BLACS}" == "xyes"; then
		BLACSINCL=""
		BLACSLIB="-L${BLACS_ROOT} -lblacs"
		AC_DEFINE([_HAVE_BLACS_], [1], [with BLACS in ISSM src])
		AC_SUBST([BLACSINCL])
		AC_SUBST([BLACSLIB])
	fi
	dnl }}}
	dnl HYPRE{{{
	AC_MSG_CHECKING([for HYPRE])
	AC_ARG_WITH(
		[hypre-dir],
		AS_HELP_STRING([--with-hypre-dir=DIR], [HYPRE root directory]),
		[HYPRE_ROOT=${withval}],
		[HYPRE_ROOT="no"]
	)
	if test "x${HYPRE_ROOT}" == "xno"; then
		HAVE_HYPRE=no
	else
		HAVE_HYPRE=yes
		if ! test -d "${HYPRE_ROOT}"; then
			AC_MSG_ERROR([HYPRE directory provided (${HYPRE_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_HYPRE}])

	dnl HYPRE libraries and header files
	if test "x${HAVE_HYPRE}" == "xyes"; then
		HYPREINCL=""
		HYPRELIB="-L${HYPRE_ROOT}/lib -lHYPRE"
		AC_DEFINE([_HAVE_HYPRE_], [1], [with HYPRE in ISSM src])
		AC_SUBST([HYPREINCL])
		AC_SUBST([HYPRELIB])
	fi
	dnl }}}
	dnl Prometheus{{{
	AC_MSG_CHECKING([for Prometheus])
	AC_ARG_WITH(
		[prometheus-dir],
		AS_HELP_STRING([--with-prometheus-dir=DIR], [Prometheus root directory]),
		[PROMETHEUS_ROOT=${withval}],
		[PROMETHEUS_ROOT="no"]
	)
	if test "x${PROMETHEUS_ROOT}" == "xno"; then
		HAVE_PROMETHEUS=no
	else
		HAVE_PROMETHEUS=yes
		if ! test -d "${PROMETHEUS_ROOT}"; then
			AC_MSG_ERROR([Prometheus directory provided (${PROMETHEUS_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_PROMETHEUS}])

	dnl Prometheus libraries and header files
	if test "x${HAVE_PROMETHEUS}" == "xyes"; then
		PROMETHEUSINCL="-I${PROMETHEUS_ROOT}/include"
		PROMETHEUSLIB="-L${PROMETHEUS_ROOT}/lib -lpromfei -lprometheus"
		AC_DEFINE([_HAVE_PROMETHEUS_], [1], [with Prometheus in ISSM src])
		AC_SUBST([PROMETHEUSINCL])
		AC_SUBST([PROMETHEUSLIB])
	fi
	dnl }}}
	dnl SEMIC{{{
	AC_MSG_CHECKING([for SEMIC])
	AC_ARG_WITH(
		[semic-dir],
		AS_HELP_STRING([--with-semic-dir=DIR], [SEMIC root directory]),
		[SEMIC_ROOT=${withval}],
		[SEMIC_ROOT="no"]
	)
	if test "x${SEMIC_ROOT}" == "xno"; then
		HAVE_SEMIC=no
	else
		HAVE_SEMIC=yes
		if ! test -d "${SEMIC_ROOT}"; then
			AC_MSG_ERROR([SEMIC directory provided (${SEMIC_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_SEMIC}])

	dnl SEMIC libraries and header files
	if test "x${HAVE_SEMIC}" == "xyes"; then
		SEMICINCL="-I${SEMIC_ROOT}"
		if test "x${IS_MSYS2}" == "xyes"; then
			SEMICLIB="-Wl,-L${SEMIC_ROOT}/lib -Wl,-lsurface_physics -Wl,-lutils"
		else
			SEMICLIB="-L${SEMIC_ROOT}/lib -lsurface_physics -lutils"
		fi
		AC_DEFINE([_HAVE_SEMIC_], [1], [with SEMIC in ISSM src])
		AC_SUBST([SEMICLIB])
		AC_SUBST([SEMICINCL])
	fi
	AM_CONDITIONAL([SEMIC], [test "x${HAVE_SEMIC}" == "xyes"])
	dnl }}}
	dnl SPAI{{{
	AC_MSG_CHECKING([for SPAI])
	AC_ARG_WITH(
		[spai-dir],
		AS_HELP_STRING([--with-spai-dir=DIR], [SPAI root directory]),
		[SPAI_ROOT=${withval}],
		[SPAI_ROOT="no"]
	)
	if test "x${SPAI_ROOT}" == "xno"; then
		HAVE_SPAI=no
	else
		HAVE_SPAI=yes
		if ! test -d "${SPAI_ROOT}"; then
			AC_MSG_ERROR([SPAI directory provided (${SPAI_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_SPAI}])

	dnl SPAI libraries and header files
	if test "x${HAVE_SPAI}" == "xyes"; then
		SPAIINCL="-I${SPAI_ROOT}/include"
		SPAILIB="-L${SPAI_ROOT}/lib -lspai"
		AC_DEFINE([_HAVE_SPAI_], [1], [with SPAI in ISSM src])
		AC_SUBST([SPAIINCL])
		AC_SUBST([SPAILIB])
	fi
	dnl }}}
	dnl SuperLU{{{
	AC_MSG_CHECKING([for SuperLU])
	AC_ARG_WITH(
		[superlu-dir],
		AS_HELP_STRING([--with-superlu-dir=DIR], [SuperLU root directory]),
		[SUPERLU_ROOT=${withval}],
		[SUPERLU_ROOT="no"]
	)
	if test "x${SUPERLU_ROOT}" == "xno"; then
		HAVE_SUPERLU=no
	else
		HAVE_SUPERLU=yes
		if ! test -d "${SUPERLU_ROOT}"; then
			AC_MSG_ERROR([SuperLU directory provided (${SUPERLU_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_SUPERLU}])

	dnl SuperLU libraries and header files
	if test "x${HAVE_SUPERLU}" == "xyes"; then
		SUPERLUINCL="-I${SUPERLU_ROOT}/include"
		SUPERLULIB="-L${SUPERLU_ROOT}/lib -lsuperlu_dist"
		AC_DEFINE([_HAVE_SUPERLU_], [1], [with SuperLU in ISSM src])
		AC_SUBST([SUPERLUINCL])
		AC_SUBST([SUPERLULIB])
	fi
	dnl }}}
	dnl SPOOLES{{{
	AC_MSG_CHECKING([for SPOOLES])
	AC_ARG_WITH(
		[spooles-dir],
		AS_HELP_STRING([--with-spooles-dir=DIR], [SPOOLES root directory]),
		[SPOOLES_ROOT=${withval}],
		[SPOOLES_ROOT="no"]
	)
	if test "x${SPOOLES_ROOT}" == "xno"; then
		HAVE_SPOOLES=no
	else
		HAVE_SPOOLES=yes
		if ! test -d "${SPOOLES_ROOT}"; then
			AC_MSG_ERROR([SPOOLES directory provided (${SPOOLES_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_SPOOLES}])

	dnl SPOOLES libraries and header files
	if test "x${HAVE_SPOOLES}" == "xyes"; then
		SPOOLESINCL="-I${SPOOLES_ROOT}/include"
		SPOOLESLIB="-L${SPOOLES_ROOT}/lib -lspooles"
		AC_DEFINE([_HAVE_SPOOLES_], [1], [with SPOOLES in ISSM src])
		AC_SUBST([SPOOLESINCL])
		AC_SUBST([SPOOLESLIB])
	fi
	dnl }}}
	dnl PaStiX{{{
	AC_MSG_CHECKING([for PaStiX])
	AC_ARG_WITH(
		[pastix-dir],
		AS_HELP_STRING([--with-pastix-dir=DIR], [PaStiX root directory]),
		[PASTIX_ROOT=${withval}],
		[PASTIX_ROOT="no"]
	)
	if test "x${PASTIX_ROOT}" == "xno"; then
		HAVE_PASTIX=no
	else
		HAVE_PASTIX=yes
		if ! test -d "${PASTIX_ROOT}"; then
			AC_MSG_ERROR([PaStiX directory provided (${PASTIX_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_PASTIX}])

	dnl PaStiX libraries and header files
	if test "x${HAVE_PASTIX}" == "xyes"; then
		PASTIXINCL="-I${PASTIX_ROOT}/include"
		PASTIXLIB="-L${PASTIX_ROOT}/lib -lpastix_XXbit_mpi_smp_nobubble_int32_simple_real_scotch_i686_pc_linux -lptscotch -lptscotcherr -lpastix"
		AC_DEFINE([_HAVE_PASTIX_], [1], [with PaStiX in ISSM src])
		AC_SUBST([PASTIXINCL])
		AC_SUBST([PASTIXLIB])
	fi
	dnl }}}
	dnl ml{{{
	AC_MSG_CHECKING([for ml])
	AC_ARG_WITH(
		[ml-dir],
		AS_HELP_STRING([--with-ml-dir=DIR],[ml root directory]),
		[ML_ROOT=$withval],
		[ML_ROOT="no"]
	)
	if test "x${ML_ROOT}" == "xno"; then
		HAVE_ML=no
	else
		HAVE_ML=yes
		if ! test -d "${ML_ROOT}"; then
			AC_MSG_ERROR([ml directory provided (${ML_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_ML}])

	dnl ml libraries and header files
	if test "x${HAVE_ML}" == "xyes"; then
		MLINCL=-I"${ML_ROOT}/include"
		MLLIB=-L"${ML_ROOT}/lib -lml"
		AC_DEFINE([_HAVE_ML_], [1], [with ml in ISSM src])
		AC_SUBST([MLINCL])
		AC_SUBST([MLLIB])
	fi
	dnl }}}
	dnl UMFPACK{{{
	AC_MSG_CHECKING([for UMFPACK])
	AC_ARG_WITH(
		[umfpack-dir],
		AS_HELP_STRING([--with-umfpack-dir=DIR], [UMFPACK root directory]),
		[UMFPACK_ROOT=${withval}],
		[UMFPACK_ROOT="no"]
	)
	if test "x${UMFPACK_ROOT}" == "xno"; then
		HAVE_UMFPACK=no
	else
		HAVE_UMFPACK=yes
		if ! test -d "${UMFPACK_ROOT}"; then
			AC_MSG_ERROR([UMFPACK directory provided (${UMFPACK_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_UMFPACK}])

	dnl UMFPACK libraries and header files
	if test "x${HAVE_UMFPACK}" == "xyes"; then
		UMFPACKINCL=""
		UMFPACKLIB="-L${UMFPACK_ROOT}/lib -lumfpack -lumfpack.5.5.1"
		AC_DEFINE([_HAVE_UMFPACK_], [1], [with UMFPACK in ISSM src])
		AC_SUBST([UMFPACKINCL])
		AC_SUBST([UMFPACKLIB])
	fi
	dnl }}}
	dnl libm (GNU math library){{{
	AC_MSG_CHECKING(for libm)
	AC_ARG_WITH(
		[math-lib],
		AS_HELP_STRING([--with-math-lib=LIB], [libm (GNU math library) to use]),
		[MATH_LIB=${withval}],
		[MATH_LIB=""]
	)
	if test -n "${MATH_LIB}"; then
		HAVE_MATH=yes
		MATHLIB="${MATH_LIB}"
		AC_DEFINE([_HAVE_MATH_], [1], [with libm (GNU math library) in ISSM src])
		AC_SUBST([MATHLIB])
	fi
	AC_MSG_RESULT([done])
	dnl }}}
	dnl Fortran{{{
	AC_MSG_CHECKING(for Fortran compilation)
	AC_ARG_WITH(
		[fortran],
		AS_HELP_STRING([--with-fortran=YES], [do we compile Fortran code (default: yes)]),
		[FORTRAN=${withval}],
		[FORTRAN=yes]
	)
	if test "x${FORTRAN}" == "xyes"; then
		HAVE_FORTRAN=yes
		AC_DEFINE([_HAVE_FORTRAN_], [1], [with Fortran capability])
	else
		HAVE_FORTRAN=no
	fi
	AM_CONDITIONAL([FORTRAN], [test "x${FORTRAN}" == "xyes"])
	AC_MSG_RESULT([${FORTRAN}])

	IS_FORTRANDIR_A_DIR=no
	if test "x${FORTRAN}" == "xyes"; then
		dnl Fortran library
		AC_MSG_CHECKING([for Fortran library])
		AC_ARG_WITH(
			[fortran-lib],
			AS_HELP_STRING([--with-fortran-lib=LIB], [Fortran library to use (and, if needed, libraries on which it depends)]),
			[FORTRAN_LIB=${withval}],
			[FORTRAN_LIB=""]
		)
		if test -n "${FORTRAN_LIB}"; then
			FORTRAN_DIR=$(echo ${FORTRAN_LIB} | sed -e "s/-Wl,//g" | sed -e "s/-L//g" | awk '{print $[1]}')
			if test -d "${FORTRAN_DIR}"; then
				FORTRANDIR="${FORTRAN_DIR}"
				IS_FORTRANDIR_A_DIR=yes
				FORTRANLIB="${FORTRAN_LIB}"
				AC_DEFINE([_HAVE_FORTRAN_], [1], [with Fortran library in ISSM src])
				AC_SUBST([FORTRANDIR])
				AC_SUBST([FORTRANLIB])
			elif test -f "${FORTRAN_DIR}"; then
				FORTRANLIB="${FORTRAN_LIB}"
				AC_DEFINE([_HAVE_FORTRAN_], [1], [with Fortran library in ISSM src])
				AC_SUBST([FORTRANDIR])
			else
				AC_MSG_ERROR([Fortran library provided (${FORTRAN_LIB}) does not exist!]);
			fi
		fi
		AC_MSG_RESULT([done])
	fi
	AM_CONDITIONAL([HAVE_FORTRANDIR], [test "x${IS_FORTRANDIR_A_DIR}" == "xyes"])
	dnl }}}
	dnl NeoPZ{{{
	AC_MSG_CHECKING([for NeoPZ])
	AC_ARG_WITH(
		[neopz-dir],
		AS_HELP_STRING([--with-neopz-dir=DIR], [NeoPZ root directory]),
		[NEOPZ_ROOT=${withval}],
		[NEOPZ_ROOT="no"]
	)
	if test "x${NEOPZ_ROOT}" == "xno"; then
		HAVE_NEOPZ=no
	else
		HAVE_NEOPZ=yes
		if ! test -d "${NEOPZ_ROOT}"; then
			AC_MSG_ERROR([NeoPZ directory provided (${NEOPZ_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_NEOPZ}])

	dnl NeoPZ libraries and header files
	if test "x${HAVE_NEOPZ}" == "xyes"; then
		NEOPZLIB="${NEOPZ_ROOT}/lib/libpz.a"
		NEOPZINCL="-I${NEOPZ_ROOT}/include"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Analysis"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Common"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/External"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Frontal"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Geom"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Integral"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/LinearSolvers"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Material"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Matrix"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Mesh"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Multigrid"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/PerfUtil"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Post"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Pre"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Refine"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Save"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Shape"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/SpecialMaps"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/StrMatrix"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/SubStruct"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Topology"
		NEOPZINCL+=" -I${NEOPZ_ROOT}/include/Util"
		AC_DEFINE([_HAVE_NEOPZ_], [1], [with NeoPZ in ISSM src])
		AC_SUBST([NEOPZINCL])
		AC_SUBST([NEOPZLIB])
	fi
	AM_CONDITIONAL([NEOPZ], [test "x${HAVE_NEOPZ}" == "xyes"])
	dnl }}}
	dnl Gmsh{{{
	AC_MSG_CHECKING([for Gmsh])
	AC_ARG_WITH(
		[gmsh-dir],
		AS_HELP_STRING([--with-gmsh-dir=DIR], [Gmsh root directory]),
		[GMSH_ROOT=${withval}],
		[GMSH_ROOT="no"]
	)
	if test "x${GMSH_ROOT}" == "xno"; then
		HAVE_GMSH=no
	else
		HAVE_GMSH=yes
		if ! test -d "${GMSH_ROOT}"; then
			AC_MSG_ERROR([Gmsh directory provided (${GMSH_ROOT}) does not exist!]);
		fi
	fi
	AC_MSG_RESULT([${HAVE_GMSH}])
	AM_CONDITIONAL([GMSH], [test "x${HAVE_GMSH}" == "xyes"])

	if test "x${HAVE_GMSH}" == "xyes"; then
		AC_DEFINE([_HAVE_GMSH_], [1], [with Gmsh in ISSM src])

		AC_MSG_CHECKING(for Gmsh version)
		GMSH_VERSION_MAJOR=`${GMSH_ROOT}/bin/gmsh -info | grep "Version" | sed -e "s/Version@<:@@<:@:blank:@:>@@:>@*:@<:@@<:@:blank:@:>@@:>@//" | cut -d "." -f 1`
		AC_MSG_RESULT([${GMSH_VERSION_MAJOR}])
		AC_DEFINE_UNQUOTED([_GMSH_VERSION_MAJOR_], ${GMSH_VERSION_MAJOR}, [Gmsh major version])
	fi
	dnl }}}
	dnl Capabilities
	dnl with-bamg{{{
	AC_MSG_CHECKING([for BAMG capability compilation])
	AC_ARG_WITH(
		[bamg],
		AS_HELP_STRING([--with-bamg=YES], [compile with BAMG capabilities (default: yes)]),
		[BAMG=${withval}],
		[BAMG=yes]
	)
	HAVE_BAMG=no
	if test "x${BAMG}" == "xyes"; then
		HAVE_BAMG=yes
		AC_DEFINE([_HAVE_BAMG_], [1], [with BAMG meshing capability])
	fi
	AM_CONDITIONAL([BAMG], [test "x${HAVE_BAMG}" == "xyes"])
	AC_MSG_RESULT([${HAVE_BAMG}])
	dnl }}}
	dnl with-ocean{{{
	AC_MSG_CHECKING(for ice/ocean coupling capability compilation)
	AC_ARG_WITH(
		[ocean],
		AS_HELP_STRING([--with-ocean = YES], [compile with ice/ocean coupling capability (default: no)]),
		[OCEAN=${withval}],
		[OCEAN=no]
	)
	HAVE_OCEAN=no
	if test "x${OCEAN}" == "xyes"; then
		HAVE_OCEAN=yes
		AC_DEFINE([_HAVE_OCEAN_], [1], [with ice/ocean coupling capability])
	fi
	AM_CONDITIONAL([OCEAN], [test "x${HAVE_OCEAN}" == "xyes"])
	AC_MSG_RESULT([${HAVE_OCEAN}])
	dnl }}}
	dnl with-kriging{{{
	AC_MSG_CHECKING(for kriging capability compilation)
	AC_ARG_WITH(
		[kriging],
		AS_HELP_STRING([--with-kriging=YES], [compile with kriging capabilities (default: yes)]),
		[KRIGING=${withval}],
		[KRIGING=yes]
	)
	HAVE_KRIGING=no
	if test "x${KRIGING}" == "xyes"; then
		HAVE_KRIGING=yes
		AC_DEFINE([_HAVE_KRIGING_], [1], [with kriging capability])
	fi
	AM_CONDITIONAL([KRIGING], [test "x${HAVE_KRIGING}" == "xyes"])
	AC_MSG_RESULT([${HAVE_KRIGING}])
	dnl }}}
	dnl performancemeasurements{{{
	AC_ARG_ENABLE(
		[performancemeasurements],
		AS_HELP_STRING([--enable-performancemeasurements], [turn performance measurements on]),
		[performancemeasurements=${enableval}],
		[performancemeasurements=no]
	)
	AC_MSG_CHECKING(for performance measurements support)
	HAVE_PERF=no
	if test "x${performancemeasurements}" == "xyes"; then
		HAVE_PERF=yes
		AC_DEFINE([_HAVE_PERFORMANCE_MEASUREMENTS_], [1], [Macro to enable performance measurements in ISSM])
	fi
	AM_CONDITIONAL([PERFORMANCE_MEASUREMENTS], [test "x${HAVE_PERF}" == "xyes"])
	AC_MSG_RESULT([${HAVE_PERF}])
	dnl }}}

	dnl Analyses
	AX_ANALYSES_SELECTION

	dnl Platform specifics
	dnl multithreading{{{
	AC_MSG_CHECKING(for number of threads)
	AC_ARG_WITH(
		[numthreads],
		AS_HELP_STRING([--with-numthreads=NUMTHREADS_VALUE], [number of threads (default: 1)]),
		[NUMTHREADS_VALUE=${withval}],
		[NUMTHREADS_VALUE=1]
	)
	dnl Check that supplied value is an integer
	if test "${NUMTHREADS_VALUE}" != "${NUMTHREADS_VALUE}"; then
		AC_MSG_ERROR([Number of threads provided (${NUMTHREADS_VALUE}) is not an integer!]);
	elif test "${NUMTHREADS_VALUE}" == "0"; then
		AC_MSG_ERROR([Number of threads must be at least 1!]);
	fi
	MULTITHREADING=no
	MULTITHREADINLIB=""
	if test "x${NUMTHREADS_VALUE}" != "x1"; then
		MULTITHREADINGLIB="-lpthread -lrt"
		case "${host_os}" in
			*darwin*)
				MULTITHREADINGLIB="-lpthread"
			;;
			*linux*)
				MULTITHREADINGLIB="-lpthread -lrt"
			;;
			*mingw*)
				MULTITHREADINGLIB=""
			;;
		esac
		AC_DEFINE([_MULTITHREADING_], [1], [with multithreading enabled])
	fi
	AC_DEFINE_UNQUOTED([_NUMTHREADS_], ${NUMTHREADS_VALUE}, [number of threads])
	AC_SUBST([MULTITHREADINGLIB])
	AC_MSG_RESULT([${NUMTHREADS_VALUE}])
	dnl }}}
	dnl 64-bit indices{{{
	AC_MSG_CHECKING([for 64-bit indices])
	AC_ARG_WITH(
		[64bit-indices],
		AS_HELP_STRING([--with-64bit-indices=bool], [use 64-bit indices (default: 0)]),
		[USE_64BIT_INDICES=${withval}],
		[USE_64BIT_INDICES=0]
	)
	if test "x${USE_64BIT_INDICES}" == "x1"; then
		AC_DEFINE([ISSM_USE_64BIT_INDICES], [1], [with 64-bit indices])
	else
		AC_DEFINE([ISSM_USE_64BIT_INDICES], [0], [with 64-bit indices])
	fi
	AC_MSG_RESULT([${USE_64BIT_INDICES}])
	dnl }}}

	dnl Checks {{{
	AC_MSG_CHECKING(consistency between all external packages)

	dnl Check that if PETSc is requested, MPI is specified
	if test "x${HAVE_PETSC}" == "xyes"; then
		if test "x${HAVE_MPI}" == "xno"; then
			AC_MSG_ERROR([PETSc requires MPI!]);
		fi
	fi

	dnl Check that we have MATLAB and/or Python support if we compile the modules
	if test "x${MODULES_VALUE}" == "xyes" && test "${HAVE_MATLAB}" == "xno" && test "${HAVE_PYTHON}" == "xno"; then
		AC_MSG_ERROR([need at least MATLAB and/or Python support to compile modules! (or use --with-modules=no)]);
	fi

	dnl Check that Fortran is provided if Gia is on
	if test "x${HAVE_GIA}" == "xyes" &&  test "${HAVE_FORTRAN}" == "xno"; then
		AC_MSG_ERROR([need Fortran compiler to compile Gia! (or use --without-Gia)]);
	fi

	dnl Check that Fortran is provided if Love is on
	if test "x${HAVE_LOVE}" == "xyes" && test "x${HAVE_FORTRAN}" == "xno"; then
		AC_MSG_ERROR([need Fortran compiler to compile Love! (or use --without-Love)]);
	fi

	dnl Check that if we have MPI, we have METIS
	if test "x${HAVE_METIS}" == "xyes" && test "x${HAVE_MPI}" == "xno"; then
		AC_MSG_ERROR([need MPI if using the METIS partitioner!]);
	fi

	dnl Check that if we run ADOL-C, we don't compile kriging.exe
	if test "x${HAVE_ADOLC}" == "xyes" && test "${HAVE_KRIGING}" == "xyes"; then
		AC_MSG_ERROR([cannot compile kriging.exe under ADOL-C conditions!]);
	fi

	dnl Check that if we run ADOL-C, we don't use PETSc for now
	if test "x${HAVE_ADOLC}" == "xyes" && test "x${HAVE_PETSC}" == "xyes"; then
		AC_MSG_ERROR([cannot compile ISSM with both PETSc and ADOL-C]);
	fi
	if test "x${HAVE_ADOLC}" == "xyes" && test "x${HAVE_CODIPACK}" == "xyes"; then
		AC_MSG_ERROR([cannot compile ISSM with both ADOL-C and CoDiPack]);
	fi
	if test "x${HAVE_ADJOINTMPI}" == "xyes" && test "x${HAVE_MEDIPACK}" == "xyes"; then
		AC_MSG_ERROR([cannot compile ISSM with both MeDiPack and AdjointMPI]);
	fi
	if test "x${HAVE_CODIPACK}" == "xyes" && test "x${HAVE_PETSC}" == "xyes" && test "x${HAVE_ADJOINTPETSC}" == "xno" ; then
		AC_MSG_ERROR([cannot compile ISSM with both CoDiPack and PETSc without adjointpetsc]);
	fi

	AC_MSG_RESULT([done])
	dnl }}}
	dnl optimization{{{
	dnl -- bypass standard optimization -g -O2 -fPIC?
	AC_MSG_CHECKING(for C++ optimization flags)
	AC_ARG_WITH(
		[cxxoptflags],
		AS_HELP_STRING([--with-cxxoptflags=CXXOPTFLAGS], [C++ optimization flags - DEPRECATED - DO NOT USE]),
		[CXXOPTFLAGS=${withval}],
		[CXXOPTFLAGS="DEPRECATED"]
	)
	if test "x${CXXOPTFLAGS}" != "xDEPRECATED"; then
	 AC_MSG_RESULT([DEPRECATED!!])
		AC_MSG_ERROR([--with-cxxoptflags does not exist anymore! you should now add this line to your configuration script: export CXXFLAGS="${CXXOPTFLAGS}"]);
	fi
	AC_MSG_RESULT([DEPRECATED])
	dnl }}}

	dnl Final variable substitution
	AC_SUBST([CFLAGS])
	AC_SUBST([CXXFLAGS])
	AC_SUBST([OSLIBS])
])

dnl =====================================================================
dnl  ISSM_ENABLE_AD – Automatic-Differentiation (CoDiPack + MediPack)
dnl =====================================================================
AC_DEFUN([ISSM_ENABLE_AD], [
  # --- command-line switches ------------------------------------------
  AC_ARG_ENABLE([ad],
    AS_HELP_STRING([--enable-ad],
      [Build ISSM with CoDiPack+MediPack automatic differentiation (disables PETSc)]),
    [enable_ad=$enableval],
    [enable_ad=no])

  AC_ARG_WITH([codipack-dir],
    AS_HELP_STRING([--with-codipack-dir=DIR],
      [Prefix of CoDiPack install]),
    [CODIPACK_ROOT=$withval], [CODIPACK_ROOT=])

  AC_ARG_WITH([medipack-dir],
    AS_HELP_STRING([--with-medipack-dir=DIR],
      [Prefix of MediPack install]),
    [MEDIPACK_ROOT=$withval], [MEDIPACK_ROOT=])

  # --- validation & flag injection ------------------------------------
  if test "x$enable_ad" = "xyes"; then
    if test -z "$CODIPACK_ROOT" || test -z "$MEDIPACK_ROOT"; then
      AC_MSG_ERROR([--enable-ad needs BOTH --with-codipack-dir and --with-medipack-dir])
    fi

    AC_DEFINE([ISSM_USE_AD], [1],
              [Define to 1 if building with automatic differentiation])

    ENABLE_PETSC=no
    AM_CONDITIONAL([USE_AD], [true])

    AM_CPPFLAGS="$AM_CPPFLAGS -I$CODIPACK_ROOT/include -I$MEDIPACK_ROOT/include -DCODI_ForcedInlines"
    AM_LDFLAGS="$AM_LDFLAGS -L$CODIPACK_ROOT/lib -L$MEDIPACK_ROOT/lib"
    LIBS="$LIBS -lcodi -lmedi"
  else
    ENABLE_PETSC=yes
    AM_CONDITIONAL([USE_AD], [false])
  fi

  dnl Export augmented vars once
  AC_SUBST([AM_CPPFLAGS])
  AC_SUBST([AM_LDFLAGS])
  AC_SUBST([LIBS])
])
