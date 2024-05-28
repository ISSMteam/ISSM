dnl ISSM Options

AC_DEFUN([ISSM_OPTIONS],[

	dnl ISSM's internal options
	dnl Debugging {{{
	AC_ARG_ENABLE([debugging],                                        dnl feature
		AS_HELP_STRING([--enable-debugging],[turn debug support on]),  dnl help string
		[enable_debugging=$enableval],                                 dnl action if given
		[enable_debugging=no])                                         dnl action if not given
	if test "x$enable_debugging" = xyes; then
		AC_DEFINE([_ISSM_DEBUG_],[1],[Macro to enable debugging in ISSM])
	fi

	dnl }}}
    dnl Shared build {{{
    AC_ARG_ENABLE([sharedlibs],                                                dnl feature
        AS_HELP_STRING([--enable-sharedlibs], [produce libISSM.so.0]),         dnl help string
        [enable_sharedlibs=$enableval],                                        dnl action if given
        [enable_sharedlibs=no])                                                dnl action if not given
    AM_CONDITIONAL([SHAREDLIBS], [test x$enable_sharedlibs = xyes])
    dnl }}}
    dnl Version{{{
    AC_ARG_ENABLE([version],                                                dnl feature
        AS_HELP_STRING([--enable-version], [produce libISSM.so.0]),         dnl help string
        [enable_version=$enableval],                                        dnl action if given
        [enable_version=no])                                                dnl action if not given
    AM_CONDITIONAL([VERSION], [test x$enable_VERSION = xyes])
    dnl }}}
	dnl Wrappers build {{{
	AC_ARG_WITH([wrappers],
		AS_HELP_STRING([--with-wrappers = value],[wrappers compilation. ]),
		[WRAPPERS_VALUE=$withval],[WRAPPERS_VALUE="yes"])
	AC_MSG_CHECKING(for wrappers compilation)
	AM_CONDITIONAL([WRAPPERS], [test x$WRAPPERS_VALUE = xyes])
	AC_MSG_RESULT($WRAPPERS_VALUE) 
	dnl }}}
	dnl Extensions{{{
	ISSMEXT=".exe"
	AC_SUBST([ISSMEXT])
	dnl }}}

	dnl ISSM's externalpackages
	dnl vendor{{{
	AC_ARG_WITH([vendor],
	  AS_HELP_STRING([--with-vendor = VENDOR], [vendor name, ex: intel]),
	  [VENDOR=$withval],[VENDOR=""]) 
	AC_MSG_CHECKING(for vendor compilers)
	if test -n "$VENDOR"; then

		if  test $VENDOR = intel-win32; then
			export CC=icl
			export CXX=icl
			export CFLAGS="-DWIN32 -D_INTEL_WIN_"
			export CXXFLAGS="-DWIN32 -D_INTEL_WIN_"
		elif  test $VENDOR = intel-win7-32; then
			export CC=cccl
			export CXX=cccl
			export CXXFLAGS="-DWIN32 -D_INTEL_WIN_ -EHsc"
			export CFLAGS="-DWIN32 -D_INTEL_WIN_ -EHsc"
			export AR="ar-lib lib"
			export RANLIB=true
			OSLIBS="kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib"
		elif  test $VENDOR = intel-win7-64; then
			export CC=cccl
			export CXX=cccl
			export CXXFLAGS="-DWIN32 -D_INTEL_WIN_ -EHsc"
			export CFLAGS="-DWIN32 -D_INTEL_WIN_ -EHsc"
			export AR="ar-lib lib"
			export RANLIB=true
			OSLIBS="kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib"
		elif test $VENDOR = intel-linux; then
			export CC=icc
			export CXX=icpc
			export CFLAGS=" -D_INTEL_LINUX_"
			export CXXFLAGS=" -D_INTEL_LINUX_"
		elif test $VENDOR = intel-discover; then
			export CC=icc
			export CXX=icpc
			export CXXFLAGS=" -O3 -D_INTEL_LINUX_ -DMPICH_IGNORE_CXX_SEEK"
			export CFLAGS=" -O3 -D_INTEL_LINUX_ -DMPICH_IGNORE_CXX_SEEK"
		elif test $VENDOR = intel-pleiades; then
			export CC=icc
			export CXX=icpc
			export CXXFLAGS=" -O3 -D_INTEL_LINUX_ "
			export CFLAGS=" -O3 -D_INTEL_LINUX_ "
        else
		AC_MSG_ERROR([unknow compiler vendor!])
		fi
	fi
	AC_SUBST([OSLIBS]) 
	AC_MSG_RESULT(done)
	dnl }}}
	dnl matlab{{{

	dnl 1. See if matlab has been provided
	AC_ARG_WITH([matlab-dir],
		AS_HELP_STRING([--with-matlab-dir=DIR], [matlab root directory. necessary for serial build.]),
		[MATLAB_ROOT=$withval],[MATLAB_ROOT=""]) 

	AC_MSG_CHECKING([whether matlab is enabled])
	if test -d "$MATLAB_ROOT"; then
		HAVE_MATLAB=yes
	else
		HAVE_MATLAB=no
	fi
	if test x$HAVE_MATLAB = xyes; then
		AC_DEFINE([_HAVE_MATLAB_],[1],[with Matlab in ISSM src])
	fi
	AC_MSG_RESULT($HAVE_MATLAB)
	AM_CONDITIONAL([MATLAB], [test x$HAVE_MATLAB = xyes])

	dnl 2. Get Matlab libraries
	if test x$HAVE_MATLAB = xyes; then

		AC_MSG_CHECKING(for matlab headers and libraries in $MATLAB_ROOT)
  		MATLABINCL="-I$MATLAB_ROOT/extern/include"

		dnl 4. get MEXLIB MEXLINK and MEXEXT (experimental)
      dnl OS-dependent variables and checks
  		case "${host_os}" in
  			*linux*)
  				if test "${host_cpu}" = "x86_64";
  				then
  					MEXLIB=-L"$MATLAB_ROOT/bin/glnxa64/ -lmex"
  					MEXLINK="-pthread -shared -W2,--version-script,${MATLAB_ROOT}/extern/lib/glnxa64/mexFunction.map";
  				else
  					MEXLIB=-L"$MATLAB_ROOT/bin/glnx86/ -lmex"
  					MEXLINK="-pthread -shared -W2,--version-script,${MATLAB_ROOT}/extern/lib/glnx86/mexFunction.map";
  				fi
  				MEXEXT=`$MATLAB_ROOT/bin/mexext`
  				MEXEXT=".$MEXEXT"
  			;;
  			*darwin*)
  				dnl mex -v gives all the flags for compilation of mex files
  				dnl if matlab version is 7.9 or more, we must use mexmaci64 (64 bits)
  				MEXLINK="-O -Wl,-flat_namespace -undefined suppress -arch i386 -bundle -Wl,-exported_symbols_list,$MATLAB_ROOT/extern/lib/maci/mexFunction.map"
  				MEXLIB=" -L$MATLAB_ROOT/bin/maci/ -lmx -lmex -lmat -lstdc++ -largeArrayDims"
  				if test $MATLAB_MAJOR -ge 7; then 
  					 if test $MATLAB_MINOR -ge 9; then 
  						  MEXLINK="-O -Wl,-flat_namespace -undefined suppress -bundle -Wl,-exported_symbols_list,$MATLAB_ROOT/extern/lib/maci64/mexFunction.map"
  							 MEXLIB=" -L$MATLAB_ROOT/bin/maci64/ -lmx -lmex -lmat -lstdc++"
  					 fi
  				fi
  				MEXEXT=`$MATLAB_ROOT/bin/mexext`
  				MEXEXT=".$MEXEXT"
  			;;
  			*cygwin*) 
  				if  test $VENDOR = intel-win7-32; then
					MATLABLIB="$MATLAB_ROOT/extern/lib/win32/microsoft"
  					MEXLIB="/link /DLL /export:mexFunction -L$MATLABLIB libmx.lib libmex.lib libmat.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib" 
  				elif  test $VENDOR = intel-win7-64; then
					MATLABLIB="$MATLAB_ROOT/extern/lib/win64/microsoft"
  					MEXLIB="/link /DLL /export:mexFunction -L$MATLABLIB libmx.lib libmex.lib libmat.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib " 
  				fi
  				MEXEXT=".$MEXEXT"
  			;;
      esac
	   AC_MSG_RESULT(done)

		AC_SUBST([MATLABINCL])
		AC_SUBST([MEX])
		MATLABWRAPPEREXT=$MEXEXT
		AC_SUBST([MATLABWRAPPEREXT])
	    AC_SUBST([MEXLIB]) 
		AC_SUBST([MEXLINK])
	fi
	dnl }}}
	dnl triangle {{{
	AC_ARG_WITH([triangle-dir],
			  AS_HELP_STRING([--with-triangle-dir=DIR], [triangle root directory. necessary for serial build]),
			 [TRIANGLE_ROOT=$withval],[TRIANGLE_ROOT=""]) 
	AC_MSG_CHECKING(for triangle headers and libraries)

	if test -d "$TRIANGLE_ROOT"; then

		dnl defaults
		HAVE_TRIANGLE=yes
		TRIANGLEINCL=-I$TRIANGLE_ROOT/

		case "${host_os}" in
				*cygwin*)
				TRIANGLELIB=$TRIANGLE_ROOT/triangle.lib
				;;
				*linux*)
				TRIANGLELIB=$TRIANGLE_ROOT/triangle.a
				;;
				*darwin*)
				TRIANGLELIB=$TRIANGLE_ROOT/triangle.a
				;;
			esac

		AC_DEFINE([_HAVE_TRIANGLE_],[1],[with Triangle in ISSM src])
		AC_SUBST([TRIANGLEINCL])
		AC_SUBST([TRIANGLELIB])

	else
		HAVE_TRIANGLE=no
	fi
	AC_MSG_RESULT($HAVE_TRIANGLE)
	dnl }}}
	dnl dakota{{{
	AC_ARG_WITH([dakota-dir],
	  AS_HELP_STRING([--with-dakota-dir=DIR], [dakota root directory. necessary for serial build]),
	  [DAKOTA_ROOT=$withval],[DAKOTA_ROOT=""]) 
	AC_MSG_CHECKING(for dakota)
	
	if test -d "$DAKOTA_ROOT"; then

		dnl defaults
		HAVE_DAKOTA=yes
		AC_MSG_RESULT($HAVE_DAKOTA)
		DAKOTAINCL=-I$DAKOTA_ROOT/include
		AC_MSG_CHECKING(for dakota version)
		DAKOTA_VERSION=`cat $DAKOTA_ROOT/include/dakota_config.h | grep "#define PACKAGE_VERSION" | sed 's/#define PACKAGE_VERSION//' | sed 's/ //g' | sed -e 's/\"//g' `
		AC_MSG_RESULT($DAKOTA_VERSION)
		AC_DEFINE_UNQUOTED([DAKOTA_VERSION],"$DAKOTA_VERSION",[Dakota version number])
		case "${host_os}" in
			*cygwin*)
				if      test x$DAKOTA_VERSION = x4.2; then
					DAKOTALIB="-L$DAKOTA_ROOT/lib -ldakota -lteuchos -lpecos -lfftw3 -llhs -levidence -lsurfpack -lconmin -lddace -lfsudace -ljega -lcport -lopt -lpsuade -lnewmat -lncsuopt -lgsl -lquadrature -lcoliny -lcolin -lpebbl -lutilib -l3po -lnappspack -lappspack -lconveyor -lshared -lcdd -lamplsolver"
				else if test x$DAKOTA_VERSION = x5.1 || test x$DAKOTA_VERSION = x5.2; then
					DAKOTALIB="-L$DAKOTA_ROOT/lib -ldakota -lteuchos -lpecos -llhs -lsparsegrid -lsurfpack -lconmin -lddace -lfsudace -ljega -lcport -loptpp -lpsuade -lncsuopt -lcolin -linterfaces -lmomh -lscolib -lpebbl -ltinyxml -lutilib -l3po -lhopspack -lnidr -lamplsolver -lboost_signals -lboost_regex -lboost_filesystem -lboost_system"
				else
					AC_MSG_ERROR([Dakota version not found or version ($DAKOTA_VERSION) not supported!]);
				fi
				fi
			;;
			*linux*)
				if      test x$DAKOTA_VERSION = x4.2; then
					DAKOTALIB="-L$DAKOTA_ROOT/lib -ldakota -lteuchos -lpecos -lfftw3 -llhs -levidence -lsurfpack -lconmin -lddace -lfsudace -ljega -lcport -lopt -lpsuade -lnewmat -lncsuopt -lgsl -lquadrature -lcoliny -lcolin -lpebbl -lutilib -l3po -lnappspack -lappspack -lconveyor -lshared -lcdd -lamplsolver"
				else if test x$DAKOTA_VERSION = x5.1 || test x$DAKOTA_VERSION = x5.2; then
					DAKOTALIB="-L$DAKOTA_ROOT/lib -ldakota -lteuchos -lpecos -llhs -lsparsegrid -lsurfpack -lconmin -lddace -lfsudace -ljega -lcport -loptpp -lpsuade -lncsuopt -lcolin -linterfaces -lmomh -lscolib -lpebbl -ltinyxml -lutilib -l3po -lhopspack -lnidr -lamplsolver -lboost_signals -lboost_regex -lboost_filesystem -lboost_system -ldl"
				else
					AC_MSG_ERROR([Dakota version not found or version ($DAKOTA_VERSION) not supported!]);
				fi
				fi
			;;
			*darwin*)
				if      test x$DAKOTA_VERSION = x4.2; then
					DAKOTALIB="-L$DAKOTA_ROOT/lib -ldakota -lteuchos -lpecos -lfftw3 -llhs -levidence -lsurfpack -lconmin -lddace -lfsudace -ljega -lcport -lopt -lpsuade -lnewmat -lncsuopt -lgsl -lquadrature -lcoliny -lcolin -lpebbl -lutilib -l3po -lnappspack -lappspack -lconveyor -lshared -lcdd -lamplsolver" 
					dnl DAKOTALIB+= "-lgslcblas -L/usr/lib -lblas -llapack"
				else if test x$DAKOTA_VERSION = x5.1 || test x$DAKOTA_VERSION = x5.2; then
					DAKOTALIB="-L$DAKOTA_ROOT/lib -ldakota -lteuchos -lpecos -llhs -lsparsegrid -lsurfpack -lconmin -lddace -lfsudace -ljega -lcport -loptpp -lpsuade -lncsuopt -lcolin -linterfaces -lmomh -lscolib -lpebbl -ltinyxml -lutilib -l3po -lhopspack -lnidr -lamplsolver -lboost_signals -lboost_regex -lboost_filesystem -lboost_system"
					dnl DAKOTALIB+= "-lgslcblas -L/usr/lib -lblas -llapack"
				else
					AC_MSG_ERROR([Dakota version not found or version ($DAKOTA_VERSION) not supported!]);
				fi
				fi
			;;
		esac
		AC_DEFINE([_HAVE_DAKOTA_],[1],[with Dakota in ISSM src])
		AC_SUBST([DAKOTAINCL])
		AC_SUBST([DAKOTALIB])

	else
		HAVE_DAKOTA=no
		AC_MSG_RESULT($HAVE_DAKOTA)
	fi
	AM_CONDITIONAL([DAKOTA], [test x$HAVE_DAKOTA = xyes])
	dnl }}}
	dnl boost{{{
	AC_ARG_WITH([boost-dir],
	  AS_HELP_STRING([--with-boost-dir=DIR], [boost root directory.]),
	  [BOOST_ROOT=$withval],[BOOST_ROOT=""]) 
	AC_MSG_CHECKING(for boost)
	
	if test -d "$BOOST_ROOT"; then
		dnl defaults
		HAVE_BOOST=yes
		BOOSTINCL=-I$BOOST_ROOT/include
		BOOSTLIB="-L$BOOST_ROOT/lib -lboost_python"

		AC_DEFINE([_HAVE_BOOST_],[1],[with Boost in ISSM src])
		AC_SUBST([BOOSTINCL])
		AC_SUBST([BOOSTLIB])
	else
		HAVE_BOOST=no
	fi
	AM_CONDITIONAL([BOOST], [test x$HAVE_BOOST = xyes])
	AC_MSG_RESULT($HAVE_BOOST)
	dnl }}}
	dnl python{{{
	AC_ARG_WITH([python-dir],
	  AS_HELP_STRING([--with-python-dir=DIR], [python root directory.]),
	  [PYTHON_ROOT=$withval],[PYTHON_ROOT=""]) 

	AC_MSG_CHECKING(for python)
	if test -d "$PYTHON_ROOT"; then
		HAVE_PYTHON="yes"
		AC_DEFINE([_HAVE_PYTHON_],[1],[with Python in ISSM src])
	else
		HAVE_PYTHON=no
	fi
	AC_MSG_RESULT($HAVE_PYTHON)

	dnl get python version
	if test x$HAVE_PYTHON = xyes; then
		AC_MSG_CHECKING(for python version)
		dnl Query Python for its version number.  Getting [:3] seems to be
		dnl the best way to do this; it's what "site.py" does in the standard
		dnl library.
		PYTHON_VERSION=$($PYTHON_ROOT/bin/python -c "import sys; print sys.version[[:3]]")
		AC_MSG_RESULT($PYTHON_VERSION)

		dnl recover major: 
		PYTHON_MAJOR=${PYTHON_VERSION%.*}
		if test x$PYTHON_MAJOR = x3; then
			dnl are we running python 3?
			HAVE_PYTHON3="yes"
		else
			HAVE_PYTHON3="no"
		fi
		AC_DEFINE_UNQUOTED([_PYTHON_MAJOR_],$PYTHON_MAJOR,[python version major])

		PYTHONINCL=-I$PYTHON_ROOT/include
		PYTHONLIB="-L$PYTHON_ROOT/lib -lpython$PYTHON_VERSION"
		PYTHONEXT=.so

		case "${host_os}" in
			*cygwin*)
			PYTHONLINK="-shared"
			;;
			*linux*)
			PYTHONLINK="-shared"
			;;
			*darwin*)
			PYTHONLINK="-dynamiclib"
			;;
		esac

		AC_SUBST([PYTHONINCL])
		AC_SUBST([PYTHONLIB])
		PYTHONWRAPPEREXT=$PYTHONEXT
		AC_SUBST([PYTHONWRAPPEREXT])
		AC_SUBST([PYTHONLINK])
	fi
	AM_CONDITIONAL([PYTHON], [test x$HAVE_PYTHON = xyes])
	AM_CONDITIONAL([PYTHON3], [test x$HAVE_PYTHON3 = xyes])
	dnl }}}
	dnl python-numpy{{{
	AC_ARG_WITH([python-numpy-dir],
	  AS_HELP_STRING([--with-python-numpy-dir=DIR], [python-numpy root directory.]),
	  [PYTHON_NUMPY_ROOT=$withval],[PYTHON_NUMPY_ROOT=""]) 
	AC_MSG_CHECKING(for python-numpy)
	
	if test -d "$PYTHON_NUMPY_ROOT"; then
		dnl defaults
		HAVE_PYTHON_NUMPY=yes
		PYTHON_NUMPYINCL="-I$PYTHON_NUMPY_ROOT -I$PYTHON_NUMPY_ROOT/core/include/numpy"

		AC_DEFINE([_HAVE_PYTHON_NUMPY_],[1],[with Python-Numpy in ISSM src])
		AC_SUBST([PYTHON_NUMPYINCL])
	else
		HAVE_PYTHON_NUMPY=no
	fi
	AC_MSG_RESULT($HAVE_PYTHON_NUMPY)
	dnl }}}
	dnl chaco{{{
	AC_ARG_WITH([chaco-dir],
	  AS_HELP_STRING([--with-chaco-dir=DIR], [chaco root directory.]),
	  [CHACO_ROOT=$withval],[CHACO_ROOT=""]) 
	AC_MSG_CHECKING(for chaco)
	
	if test -d "$CHACO_ROOT"; then

		dnl defaults
		HAVE_CHACO=yes
		CHACOINCL=-I$CHACO_ROOT/include
		CHACOLIB="-L$CHACO_ROOT/lib -lchacominusblas"

		AC_DEFINE([_HAVE_CHACO_],[1],[with Chaco in ISSM src])
		AC_SUBST([CHACOINCL])
		AC_SUBST([CHACOLIB])

	else
		HAVE_CHACO=no
	fi
	AC_MSG_RESULT($HAVE_CHACO)
	dnl }}}
	dnl scotch{{{
	AC_ARG_WITH([scotch-dir],
	  AS_HELP_STRING([--with-scotch-dir=DIR], [scotch root directory.]),
	  [SCOTCH_ROOT=$withval],[SCOTCH_ROOT=""]) 
	AC_MSG_CHECKING(for scotch)
	
	if test -d "$SCOTCH_ROOT"; then

		dnl defaults
		HAVE_SCOTCH=yes
		SCOTCHINCL="-DNOFILEIO -I$SCOTCH_ROOT/include -DSCOTCH_VERSION=\\\"UNKNOWN\\\""
		SCOTCHLIB="-L$SCOTCH_ROOT/lib -lnfioscotch -lnfioscotcherr -lnfioscotcherrexit -lscotchmetis"

		AC_DEFINE([_HAVE_SCOTCH_],[1],[with Scotch in ISSM src])
		AC_SUBST([SCOTCHINCL])
		AC_SUBST([SCOTCHLIB])

	else
		HAVE_SCOTCH=no
	fi
	AC_MSG_RESULT($HAVE_SCOTCH)
	dnl }}}
	dnl adolc{{{
	AC_ARG_WITH([adolc-dir],
		AS_HELP_STRING([--with-adolc-dir=DIR], [adolc root directory.]),
		[ADOLC_ROOT=$withval],[ADOLC_ROOT="no"]) 
	AC_MSG_CHECKING(for adolc)

	if test "x$ADOLC_ROOT" = "xno"; then
		HAVE_ADOLC=no
	else
		if test -d "$ADOLC_ROOT"; then

			dnl defaults
			HAVE_ADOLC=yes
			ADOLCINCL="-I$ADOLC_ROOT/include"
			ADOLCLIB="-L$ADOLC_ROOT/lib64 -ladolc"

			AC_DEFINE([_HAVE_ADOLC_],[1],[with adolc in ISSM src])
			AC_SUBST([ADOLCINCL])
			AC_SUBST([ADOLCLIB])

		else
			echo  "Specified directory does not exist!"
			exit 1
		fi
	fi
	AM_CONDITIONAL([ADOLC], [test x$HAVE_ADOLC = xyes])
	AC_MSG_RESULT($HAVE_ADOLC)
	dnl }}}
	dnl adolc-version{{{
	AC_ARG_WITH([adolc-version],
		AS_HELP_STRING([--with-adolc-version=number], [adolc version.]),
		[ADOLC_VERSION=$withval],[ADOLC_VERSION=2]) 
	AC_MSG_CHECKING(for adolc-version) 

	AC_DEFINE_UNQUOTED([_ADOLC_VERSION_],$ADOLC_VERSION,[ADOLC version])
	AC_MSG_RESULT($ADOLC_VERSION)
	dnl }}}
	dnl adic2{{{
	AC_ARG_WITH([adic2-dir],
	  AS_HELP_STRING([--with-adic2-dir=DIR], [adic2 root directory.]),
	  [ADIC2_ROOT=$withval],[ADIC2_ROOT="no"]) 
	AC_MSG_CHECKING(for adic2)

	if test "x$ADIC2_ROOT" = "xno"; then
		HAVE_ADIC2=no
	else
		if test -d "$ADIC2_ROOT"; then

			dnl defaults
			HAVE_ADIC2=yes
			ADIC2INCL="-DADIC2_DENSE -I$ADIC2_ROOT/include -I$ADIC2_ROOT/share/runtime_dense/"
			ADIC2LIB=""

			AC_DEFINE([_HAVE_ADIC2_],[1],[with adic2 in ISSM src])
			AC_SUBST([ADIC2INCL])
			AC_SUBST([ADIC2LIB])

		else
			echo  "Specified directory does not exist!"
			exit 1
		fi
	fi
	AM_CONDITIONAL([ADIC2], [test x$HAVE_ADIC2 = xyes])
	AC_MSG_RESULT($HAVE_ADIC2)
	dnl }}}
	dnl gsl{{{
	AC_ARG_WITH([gsl-dir],
	  AS_HELP_STRING([--with-gsl-dir=DIR], [gsl root directory.]),
	  [GSL_ROOT=$withval],[GSL_ROOT=""]) 
	AC_MSG_CHECKING(for gsl)
	
	if test -d "$GSL_ROOT"; then

		dnl defaults
		HAVE_GSL=yes
		GSLINCL="-I$GSL_ROOT/include"
		GSLLIB="-dy -L$GSL_ROOT/lib -lgsl -lgslcblas -lm"

		AC_DEFINE([_HAVE_GSL_],[1],[with gsl in ISSM src])
		AC_SUBST([GSLINCL])
		AC_SUBST([GSLLIB])

	else
		HAVE_GSL=no
	fi
	AM_CONDITIONAL([GSL], [test x$HAVE_GSL = xyes])
	AC_MSG_RESULT($HAVE_GSL)
	dnl }}}
	dnl rose{{{
	AC_ARG_WITH([rose-dir],
	  AS_HELP_STRING([--with-rose-dir=DIR], [rose root directory.]),
	  [ROSE_ROOT=$withval],[ROSE_ROOT=""]) 
	AC_MSG_CHECKING(for rose)
	
	if test -d "$ROSE_ROOT"; then

		dnl defaults
		HAVE_ROSE=yes
		ROSEINCL="-I$ROSE_ROOT/include"
		ROSELIB=""

		AC_DEFINE([_HAVE_ROSE_],[1],[with rose in ISSM src])
		AC_SUBST([ROSEINCL])
		AC_SUBST([ROSELIB])

	else
		HAVE_ROSE=no
	fi
	AM_CONDITIONAL([ROSE], [test x$HAVE_ROSE = xyes])
	AC_MSG_RESULT($HAVE_ROSE)
	dnl }}}
	dnl mpi{{{
	AC_MSG_CHECKING(for mpi)
	AC_ARG_WITH([mpi-lib],
		AS_HELP_STRING([--with-mpi-lib = options],[mpi options, for ex: "-L$MPIROOT -lmpich]),
		[MPILIB=$withval],[MPILIB=""])
	
	AC_ARG_WITH([mpi-include],
	  AS_HELP_STRING([--with-mpi-include=DIR],[mpi include directory, necessary for parallel build]),
	  [MPI_INCLUDE=$withval],[MPI_INCLUDE=""])
	
	if test -z "$MPILIB" ; then
		HAVE_MPI=no
	else
		if test -z "$MPI_INCLUDE" ; then
			HAVE_MPI=no
		else
			HAVE_MPI=yes
			MPIINCL=-I"$MPI_INCLUDE"
			AC_DEFINE([_HAVE_MPI_],[1],[with Mpi in ISSM src])
			AC_DEFINE([HAVE_MPI],[1],[Mpi Flag for Dakota (DO NOT REMOVE)])
			AC_SUBST([MPIINCL])
			AC_SUBST([MPILIB])
		fi
	fi
	AM_CONDITIONAL([MPI], [test x$HAVE_MPI = xyes])
	AC_MSG_RESULT($HAVE_MPI)
	dnl }}}
	dnl petsc{{{
	AC_ARG_WITH([petsc-dir],
	  AS_HELP_STRING([--with-petsc-dir=DIR],[PETSc root directory, necessary for parallel build]),
	  [PETSC_ROOT=$withval],[PETSC_ROOT=""])
		
	if test -d "$PETSC_ROOT"; then
		AC_MSG_CHECKING(for petsc version)
		PETSC_MAJOR=`cat $PETSC_ROOT/include/petscversion.h | grep "#define PETSC_VERSION_MAJOR" | sed 's/#define PETSC_VERSION_MAJOR//' | sed 's/ //g'`
		PETSC_MINOR=`cat $PETSC_ROOT/include/petscversion.h | grep "#define PETSC_VERSION_MINOR" | sed 's/#define PETSC_VERSION_MINOR//' | sed 's/ //g'`
		AC_DEFINE_UNQUOTED([_PETSC_MAJOR_],$PETSC_MAJOR,[PETSc version major])
		AC_DEFINE_UNQUOTED([_PETSC_MINOR_],$PETSC_MINOR,[PETSc version minor])
		AC_MSG_RESULT($PETSC_MAJOR.$PETSC_MINOR)

		PETSC_VERSION_DATE_HG=`cat $PETSC_ROOT/include/petscversion.h | grep "#define PETSC_VERSION_DATE_HG" | sed 's/#define PETSC_VERSION_DATE_HG//' | sed 's/ //g' | sed -e 's/\"//g' `
		PETSC_RELEASE=`cat $PETSC_ROOT/include/petscversion.h | grep "#define PETSC_VERSION_RELEASE" | sed 's/#define PETSC_VERSION_RELEASE//' | sed 's/ //g'`

		AC_MSG_CHECKING(whether petsc is the development version)
		dnl if test x$PETSC_VERSION_DATE_HG = xunknown; then
		if test "$PETSC_RELEASE" = "0"; then
		   AC_DEFINE([_HAVE_PETSCDEV_],[1],[with PETSc-dev])
			AC_MSG_RESULT(yes)
		else
			AC_MSG_RESULT(no)
		fi
	fi
	
	AC_ARG_WITH([petsc-arch],
	  AS_HELP_STRING([--with-petsc-arch=DIR],[PETSc arch , necessary for parallel build]),
	  [PETSC_ARCH=$withval],[PETSC_ARCH=""])

	AC_MSG_CHECKING(for petsc headers and libraries in $PETSC_ROOT for architecture $PETSC_ARCH)
	
	dnl To ge PETSc's libraries:
	dnl cd externalpackages/petsc/src
	dnl make getlinklibs
	if test -d "$PETSC_ROOT"; then

	 PETSCINCL=" -I$PETSC_ROOT/include"
	 dnl Add other location (maybe not needed anymore)
	 if test -d "$PETSC_ROOT/$PETSC_ARCH/include"; then
	  PETSCINCL+=" $PETSC_ROOT/$PETSC_ARCH/include"
	 fi
	 if test -d "$PETSC_ROOT/include/$PETSC_ARCH"; then
	  PETSCINCL+=" $PETSC_ROOT/include/$PETSC_ARCH"
	 fi
	
	 case "${host_os}" in
			*cygwin*)
			if test $PETSC_MAJOR -lt 3 ; then
				PETSCLIB="-L$PETSC_ROOT/lib libpetscksp.lib libpetscdm.lib libpetscmat.lib libpetscvec.lib libpetscsnes.lib libpetscts.lib libmpiuni.lib libpetsc.lib"
			else
				PETSCLIB="-L$PETSC_ROOT/lib libpetsc.lib"
				if test $PETSC_MAJOR -gt 3 || test $PETSC_MINOR -ge 3; then PETSCLIB+=" libmetis.lib"; fi
			fi
			;;
			*linux*)
			if test $PETSC_MAJOR -lt 3 ; then
				PETSCLIB="-L$PETSC_ROOT/lib -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetsc  -lpetscsnes -lpetscts"
			else
				PETSCLIB="-L$PETSC_ROOT/lib -lpetsc -ldl"
				if test $PETSC_MAJOR -gt 3 || test $PETSC_MINOR -ge 3; then PETSCLIB+=" -lmetis"; fi
			fi
			;;
			*darwin*)
			if test $PETSC_MAJOR -lt 3 ; then
				PETSCLIB="-L$PETSC_ROOT/lib -lpetscksp -lpetscdm -lpetscmat -lpetscvec -lpetscsnes -lpetscts -lpetsc"
			else
				PETSCLIB="-L$PETSC_ROOT/lib -lpetsc"
				if test $PETSC_MAJOR -gt 3 || test $PETSC_MINOR -ge 3; then PETSCLIB+=" -lmetis"; fi
			fi
			;;
		esac
		AC_DEFINE([_HAVE_PETSC_],[1],[with PETSc in ISSM src])
		AC_SUBST([PETSCINCL])
		AC_SUBST([PETSCLIB])
		HAVE_PETSC=yes
	else
		HAVE_PETSC=no
	fi
	
	AM_CONDITIONAL([PETSC], [test x$HAVE_PETSC = xyes])
	AC_MSG_RESULT($HAVE_PETSC)
	dnl }}}
	dnl metis{{{
	if test "$HAVE_PETSC" = "yes" && test "x$PETSC_MAJOR" = "x3" && test $PETSC_MINOR -ge 3; then

		dnl in petsc >=3.3, metis is provided
		HAVE_METIS="yes"
		AC_DEFINE([_METIS_VERSION_],[5],[ Metis version number])
		AC_DEFINE([_HAVE_METIS_],[1],[with Metis in ISSM src])
	else
		AC_ARG_WITH([metis-dir],
		  AS_HELP_STRING([--with-metis-dir=DIR],[metis root directory. necessary for serial build]),
		  [METIS_ROOT=$withval],[METIS_ROOT=""])

		AC_MSG_CHECKING(for metis headers and libraries in $METIS_ROOT)

		if test -d "$METIS_ROOT"; then

			dnl first figure out version of metis: does the VERSION file exist?
			if test -e "$METIS_ROOT/VERSION"; then
				METIS_VERSION=4
			else
				METIS_VERSION=5
			fi

			dnl defaults
			HAVE_METIS=yes

			if test "$METIS_VERSION" = "4" ; then
					
				METISINCL=-I"$METIS_ROOT/Lib" 
				case "${host_os}" in
					*cygwin*)
					METISLIB="-L$METIS_ROOT libmetis.lib"
					;;
					*linux*)
					METISLIB=-L"$METIS_ROOT/ -lmetis"
					;;
					*darwin*)
					METISLIB=-L"$METIS_ROOT/ -lmetis"
					;;
				esac

					AC_DEFINE([_METIS_VERSION_],[4],[ Metis version number])
			fi
			if test "$METIS_VERSION" = "5" ; then
		
				case "${host_os}" in
					*cygwin*)
					METISLIB="-L$METIS_ROOT libmetis.lib"
					;;
					*linux*)
					METISLIB=-L"$METIS_ROOT/lib -lmetis"
					;;
					*darwin*)
					METISLIB=-L"$METIS_ROOT/lib -lmetis"
					;;
				esac

				METISINCL=-I"$METIS_ROOT/include" 
				AC_DEFINE([_METIS_VERSION_],[5],[ Metis version number])
			fi

			AC_DEFINE([_HAVE_METIS_],[1],[with Metis in ISSM src])
			AC_SUBST([METISINCL])
			AC_SUBST([METISLIB])
		else
			HAVE_METIS=no
		fi
		AC_MSG_RESULT($HAVE_METIS)
	fi
	AM_CONDITIONAL([METIS], [test x$HAVE_METIS = xyes])
	dnl }}}
	dnl tao{{{
	AC_ARG_WITH([tao-dir],
		AS_HELP_STRING([--with-tao-dir=DIR], [tao root directory.]),
		[TAO_ROOT=$withval],[TAO_ROOT=""]) 
	AC_MSG_CHECKING(for tao)

	if test -d "$TAO_ROOT"; then

	  HAVE_TAO=yes
	  TAOINCL="-I$TAO_ROOT/ -I$TAO_ROOT/include -I$TAO_ROOT/bmake/ "
	  TAOLIB="-L$TAO_ROOT/lib -ltao -lpetsc"

	  AC_DEFINE([_HAVE_TAO_],[1],[with Tao in ISSM src])
	  AC_SUBST([TAOINCL])
	  AC_SUBST([TAOLIB])
	else
		HAVE_TAO=no
	fi
	AC_MSG_RESULT($HAVE_TAO)
	dnl }}}
	dnl slepc{{{
	AC_ARG_WITH([slepc-dir],
	  AS_HELP_STRING([--with-slepc-dir=DIR],[slepc root directory]),
	  [SLEPC_ROOT=$withval],[SLEPC_ROOT=""])
			  
	AC_MSG_CHECKING(for slepc headers and libraries in $SLEPC_ROOT)
	if test -d "$SLEPC_ROOT"; then
		HAVE_SLEPC=yes
		SLEPCINCL=-I"$SLEPC_ROOT/include"
		SLEPCLIB=-L"$SLEPC_ROOT/lib/ -lslepc"

		AC_DEFINE([_HAVE_SLEPC_],[1],[with Slepc in ISSM src])
		AC_SUBST([SLEPCINCL])
		AC_SUBST([SLEPCLIB])
	else
		HAVE_SLEPC=no
	fi
	AC_MSG_RESULT($HAVE_SLEPC)
	dnl }}}
	dnl shapelib{{{
	AC_ARG_WITH([shapelib-dir],
	  AS_HELP_STRING([--with-shapelib-dir=DIR], [shapelib root directory]),
	  [SHAPELIB_ROOT=$withval],[SHAPELIB_ROOT=""])
			  
	AC_MSG_CHECKING(for shapelib headers and libraries in $SHAPELIB_ROOT)
	if test -d "$SHAPELIB_ROOT"; then

		dnl defaults
		HAVE_SHAPELIB=yes
		SHAPELIBINCL=-I"$SHAPELIB_ROOT/include"
		SHAPELIBLIB=-L"$SHAPELIB_ROOT/lib/ -lshape"

		AC_DEFINE([_HAVE_SHAPELIB_],[1],[with Shapelib in ISSM src])
		AC_SUBST([SHAPELIBINCL])
		AC_SUBST([SHAPELIBLIB])
	else
		HAVE_SHAPELIB=no
	fi
	AC_MSG_RESULT($HAVE_SHAPELIB)
	dnl }}}
	dnl scalapack{{{
	AC_ARG_WITH([scalapack-dir],
	  AS_HELP_STRING([--with-scalapack-dir=DIR],[scalapack root directory]),
	  [SCALAPACK_ROOT=$withval],[SCALAPACK_ROOT=""])
			  
	AC_MSG_CHECKING(for scalapack headers and libraries in $SCALAPACK_ROOT)
	if test -d "$SCALAPACK_ROOT"; then

		dnl defaults
		HAVE_SCALAPACK=yes
		if test x$VENDOR = xintel-discover; then
		 SCALAPACKLIB=-L"$SCALAPACK_ROOT/ -lmkl_scalapack_lp64"
		else
		 SCALAPACKLIB=-L"$SCALAPACK_ROOT/ -lscalapack"
		fi

		AC_DEFINE([_HAVE_SCALAPACK_],[1],[with Scalapack in ISSM src])
		AC_SUBST([SCALAPACKLIB])
	else
		HAVE_SCALAPACK=no
	fi
	AC_MSG_RESULT($HAVE_SCALAPACK)
	dnl }}}
	dnl blas{{{
	AC_ARG_WITH([blas-lapack-dir],
	  AS_HELP_STRING([--with-blas-lapack-dir=DIR],[blas-lapack root directory]),
	  [BLASLAPACK_ROOT=$withval],[BLASLAPACK_ROOT=""])
			  
	AC_MSG_CHECKING(for blas and lapack headers and libraries in $BLASLAPACK_ROOT)
	if test -d "$BLASLAPACK_ROOT"; then

		dnl defaults
		HAVE_BLASLAPACK=yes
		BLASLAPACKINCL=""
	
		if test x$VENDOR = xintel-discover; then
		 BLASLAPACKLIB=-L"$BLASLAPACK_ROOT -lmkl_lapack -lmkl -lguide -lpthread"
		else
		dnl: branch on whether we are running on windows or linux.
		case "${host_os}" in
			*cygwin*)
			BLASLAPACKLIB="-L$BLASLAPACK_ROOT libf2cblas.lib  libf2clapack.lib"
			;;
			*linux*)
			BLASLAPACKLIB=-L"$BLASLAPACK_ROOT -lflapack -lfblas " 
			;;
			*darwin*)
			BLASLAPACKLIB=-L"$BLASLAPACK_ROOT -lflapack -lfblas " 
			;;
		esac
		fi

		AC_DEFINE([_HAVE_BLASLAPACK_],[1],[with blas lapack in ISSM src])
		AC_SUBST([BLASLAPACKLIB])
		AC_SUBST([BLASLAPACKINCL])
	else
		HAVE_BLASLAPACK=no
	fi
	AC_MSG_RESULT($HAVE_BLASLAPACK)
	dnl }}}
	dnl mkl{{{
	AC_ARG_WITH([mkl-dir],
	  AS_HELP_STRING([--with-mkl-dir=DIR],[mkl root directory]),
	  [MKL_ROOT=$withval],[MKL_ROOT=""])
			  
	AC_MSG_CHECKING(for mkl headers and libraries in $MKL_ROOT)
	if test -d "$MKL_ROOT"; then

		dnl defaults
		HAVE_MKL=yes
		MKLINCL=""
		MKLLIB=-L"$MKL_ROOT -lmkl -lmkl_lapack -lmkl_scalapack_ilp64   -lmkl_blacs_sgimpt_ilp64 -lguide  -lpthread"
		AC_DEFINE([_HAVE_MKL_],[1],[with mkl in ISSM src])
		AC_SUBST([MKLLIB])
		AC_SUBST([MKLINCL])
	else
		HAVE_MKL=no
	fi
	AC_MSG_RESULT($HAVE_MKL)
	dnl }}}
	dnl plapack{{{
	AC_MSG_CHECKING(for plapack)
	
	AC_ARG_WITH([plapack-lib],
	  AS_HELP_STRING([--with-plapack-lib = lib],[plapack library]),
	  [PLAPACK_LIB=$withval],[PLAPACK_LIB=""])
	
	AC_ARG_WITH([plapack-include],
			  AS_HELP_STRING([--with-plapack-include = include],
							 [plapack include ]),
			  [PLAPACK_INCLUDE=$withval],[PLAPACK_INCLUDE=""])
	  
	if test -n "$PLAPACK_LIB"; then
		if test -n "$PLAPACK_INCLUDE"; then
		
			dnl defaults
			HAVE_PLAPACK=yes
			PLAPACKINCL="$PLAPACK_INCLUDE"
			PLAPACKLIB="$PLAPACK_LIB"

			AC_DEFINE([_HAVE_PLAPACK_],[1],[with Plapack in ISSM src])
			AC_SUBST([PLAPACKINCL])
			AC_SUBST([PLAPACKLIB])
		else
			HAVE_PLAPACK=no
		fi
	else
		HAVE_PLAPACK=no
	fi
	AC_MSG_RESULT($HAVE_PLAPACK)
	dnl }}}
	dnl mumps{{{
	AC_ARG_WITH([mumps-dir],
	  AS_HELP_STRING([--with-mumps-dir=DIR],[mumps root directory]),
	  [MUMPS_ROOT=$withval],[MUMPS_ROOT=""])
			  
	AC_MSG_CHECKING(for mumps headers and libraries in $MUMPS_ROOT)
	if test -d "$MUMPS_ROOT"; then

		dnl defaults
		HAVE_MUMPS=yes
		MUMPSINCL=-I"$MUMPS_ROOT/include"
		if test "$PETSC_MAJOR" = "2" ; then
		MUMPSLIB=-L"$MUMPS_ROOT/lib -ldmumps -lcmumps  -lpord "
		else
		dnl MUMPSLIB=-L"$MUMPS_ROOT/lib -ldmumps -lcmumps  -lmumps_common -lpord -lparmetis -lzmumps"
		MUMPSLIB=-L"$MUMPS_ROOT/lib -ldmumps -lcmumps  -lmumps_common -lpord -lparmetis"
		fi

		AC_DEFINE([_HAVE_MUMPS_],[1],[with Mumps in ISSM src])
		AC_SUBST([MUMPSINCL])
		AC_SUBST([MUMPSLIB])
	else
		HAVE_MUMPS=no
	fi
	AM_CONDITIONAL([MUMPS], [test x$HAVE_MUMPS = xyes])
	AC_MSG_RESULT($HAVE_MUMPS)
	dnl }}}
	dnl blacs{{{
	AC_ARG_WITH([blacs-dir],
		AS_HELP_STRING([--with-blacs-dir=DIR],[blacs root directory]),
			  [BLACS_ROOT=$withval],[BLACS_ROOT=""])
			  
	AC_MSG_CHECKING(for blacs headers and libraries in $BLACS_ROOT)
	if test -d "$BLACS_ROOT"; then

		dnl defaults
		HAVE_BLACS=yes
		BLACSINCL=""
		if test x$VENDOR = xintel-discover; then
		 BLACSLIB=-L"$BLACS_ROOT/ -lmkl_blacs_intelmpi_lp64"
		else
		 BLACSLIB=-L"$BLACS_ROOT/ -lblacs"
		fi
        
		AC_DEFINE([_HAVE_BLACS_],[1],[with Blacs in ISSM src])
		AC_SUBST([BLACSINCL])
		AC_SUBST([BLACSLIB])
	else
		HAVE_BLACS=no
	fi
	AC_MSG_RESULT($HAVE_BLACS)
	dnl }}}
	dnl hypre{{{
	AC_ARG_WITH([hypre-dir],
	  AS_HELP_STRING([--with-hypre-dir=DIR],[hypre root directory]),
			  [HYPRE_ROOT=$withval],[HYPRE_ROOT=""])
			  
	AC_MSG_CHECKING(for hypre headers and libraries in $HYPRE_ROOT)
	if test -d "$HYPRE_ROOT"; then

		dnl defaults
		HAVE_HYPRE=yes
		HYPREINCL=""
		HYPRELIB=-L"$HYPRE_ROOT/lib -lHYPRE"
        
		AC_DEFINE([_HAVE_HYPRE_],[1],[with Blacs in ISSM src])
		AC_SUBST([HYPREINCL])
		AC_SUBST([HYPRELIB])
	else
		HAVE_HYPRE=no
	fi
	AC_MSG_RESULT($HAVE_HYPRE)
	dnl }}}
	dnl prometheus{{{
		AC_ARG_WITH([prometheus-dir],
					AS_HELP_STRING([--with-prometheus-dir=DIR],[prometheus root directory]),
					[PROMETHEUS_ROOT=$withval],[PROMETHEUS_ROOT=""])

		  AC_MSG_CHECKING(for prometheus headers and libraries in $PROMETHEUS_ROOT)
		  if test -d "$PROMETHEUS_ROOT"; then

			dnl defaults
			  HAVE_PROMETHEUS=yes
			  PROMETHEUSINCL=-I"$PROMETHEUS_ROOT/include"
			  PROMETHEUSLIB=-L"$PROMETHEUS_ROOT/lib -lpromfei -lprometheus -lparmetis"

			  AC_DEFINE([_HAVE_PROMETHEUS_],[1],[with Prometheus in ISSM src])
			  AC_SUBST([PROMETHEUSINCL])
			  AC_SUBST([PROMETHEUSLIB])
		  else
				HAVE_PROMETHEUS=no
			fi
			AC_MSG_RESULT($HAVE_PROMETHEUS)
		dnl }}}
dnl spai{{{
	AC_ARG_WITH([spai-dir],
				AS_HELP_STRING([--with-spai-dir=DIR],[spai root directory]),
				[SPAI_ROOT=$withval],[SPAI_ROOT=""])

	  AC_MSG_CHECKING(for spai headers and libraries in $SPAI_ROOT)
	  if test -d "$SPAI_ROOT"; then

		dnl defaults
		  HAVE_SPAI=yes
		  SPAIINCL=-I"$SPAI_ROOT/include"
		  SPAILIB=-L"$SPAI_ROOT/lib -lspai"

		  AC_DEFINE([_HAVE_SPAI_],[1],[with Spai in ISSM src])
		  AC_SUBST([SPAIINCL])
		  AC_SUBST([SPAILIB])
	  else
		HAVE_SPAI=no
		  fi
		  AC_MSG_RESULT($HAVE_SPAI)
		  dnl }}}
dnl superlu{{{ 
	AC_ARG_WITH([superlu-dir],
				AS_HELP_STRING([--with-superlu-dir=DIR],[superlu root directory]),
				[SUPERLU_ROOT=$withval],[SUPERLU_ROOT=""])

	  AC_MSG_CHECKING(for superlu headers and libraries in $SUPERLU_ROOT)
	  if test -d "$SUPERLU_ROOT"; then

		dnl defaults
		  HAVE_SUPERLU=yes
		  SUPERLUINCL=-I"$SUPERLU_ROOT/include"
		  SUPERLULIB=-L"$SUPERLU_ROOT/lib -lsuperlu_4.3"

		  AC_DEFINE([_HAVE_SUPERLU_],[1],[with Superlu in ISSM src])
		  AC_SUBST([SUPERLUINCL])
		  AC_SUBST([SUPERLULIB])
	  else
		HAVE_SUPERLU=no
		  fi
		  AC_MSG_RESULT($HAVE_SUPERLU)
		  dnl }}}
dnl spooles{{{ 
	AC_ARG_WITH([spooles-dir],
				AS_HELP_STRING([--with-spooles-dir=DIR],[spooles root directory]),
				[SPOOLES_ROOT=$withval],[SPOOLES_ROOT=""])

	  AC_MSG_CHECKING(for spooles headers and libraries in $SPOOLES_ROOT)
	  if test -d "$SPOOLES_ROOT"; then

		dnl defaults
		  HAVE_SPOOLES=yes
		  SPOOLESINCL=-I"$SPOOLES_ROOT/include"
		  SPOOLESLIB=-L"$SPOOLES_ROOT/lib -lspooles"

		  AC_DEFINE([_HAVE_SPOOLES_],[1],[with Spooles in ISSM src])
		  AC_SUBST([SPOOLESINCL])
		  AC_SUBST([SPOOLESLIB])
	  else
		HAVE_SPOOLES=no
		  fi
		  AC_MSG_RESULT($HAVE_SPOOLES)
		  dnl }}}
dnl pastix{{{ 
	AC_ARG_WITH([pastix-dir],
				AS_HELP_STRING([--with-pastix-dir=DIR],[pastix root directory]),
				[PASTIX_ROOT=$withval],[PASTIX_ROOT=""])

	  AC_MSG_CHECKING(for pastix headers and libraries in $PASTIX_ROOT)
	  if test -d "$PASTIX_ROOT"; then

		dnl defaults
		  HAVE_PASTIX=yes
		  PASTIXINCL=-I"$PASTIX_ROOT/include"
		  PASTIXLIB=-L"$PASTIX_ROOT/lib -lpastix_XXbit_mpi_smp_nobubble_int32_simple_real_scotch_i686_pc_linux -lptscotch -lptscotcherr -lpastix"

		  AC_DEFINE([_HAVE_PASTIX_],[1],[with Pastix in ISSM src])
		  AC_SUBST([PASTIXINCL])
		  AC_SUBST([PASTIXLIB])
	  else
		HAVE_PASTIX=no
		  fi
		  AC_MSG_RESULT($HAVE_PASTIX)
		  dnl }}}
	dnl ml{{{
	AC_ARG_WITH([ml-dir],
	  AS_HELP_STRING([--with-ml-dir=DIR],[ml root directory]),
			  [ML_ROOT=$withval],[ML_ROOT=""])
			  
	AC_MSG_CHECKING(for ml headers and libraries in $ML_ROOT)
	if test -d "$ML_ROOT"; then

		dnl defaults
		HAVE_ML=yes
		MLINCL=-I"$ML_ROOT/include"
		MLLIB=-L"$ML_ROOT/lib -lml"
        
		AC_DEFINE([_HAVE_ML_],[1],[with Blacs in ISSM src])
		AC_SUBST([MLINCL])
		AC_SUBST([MLLIB])
	else
		HAVE_ML=no
	fi
	AC_MSG_RESULT($HAVE_ML)
	dnl }}}
	dnl umfpack{{{
		AC_ARG_WITH([umfpack-dir],
		  AS_HELP_STRING([--with-umfpack-dir=DIR],[UMFPACK root directory]),
					[UMFPACK_ROOT=$withval],[UMFPACK_ROOT=""])

		AC_MSG_CHECKING(for UMFPACK headers and libraries in $UMFPACK_ROOT)
		if test -d "$UMFPACK_ROOT"; then

			dnl defaults
			HAVE_UMFPACK=yes
			UMFPACKINCL=""
			UMFPACKLIB=-L"$UMFPACK_ROOT/lib -lumfpack -lumfpack.5.5.1"

			AC_DEFINE([_HAVE_UMFPACK_],[1],[with UMFPACK in ISSM src])
			AC_SUBST([UMFPACKINCL])
			AC_SUBST([UMFPACKLIB])
		else
			HAVE_UMFPACK=no
		fi
		AC_MSG_RESULT($HAVE_UMFPACK)
	dnl }}}
dnl math{{{
	AC_MSG_CHECKING(for math library)
	AC_ARG_WITH([math-lib],
	  AS_HELP_STRING([--with-math-lib = otions],[math options, for ex: "/usr/lib/libm.a]),
	  [MATH_LIB=$withval],[MATH_LIB=""])

	dnl check that --with-math-lib may have been provided
	if test -n "$MATH_LIB" ; then
		HAVE_MATH=yes
		MATHLIB="$MATH_LIB"

		AC_DEFINE([_HAVE_MATH_],[1],[with MATH in ISSM src])
		AC_SUBST([MATHLIB])
	fi
	AC_MSG_RESULT(done)
	dnl }}}
	dnl fortran{{{
	AC_ARG_WITH([fortran],
		AS_HELP_STRING([--with-fortran = YES], [do we compile fortran code (default is yes)]),
		[FORTRAN=$withval],[FORTRAN=yes]) 
	AC_MSG_CHECKING(for fortran compilation)
	if test "x$FORTRAN" = "xyes"; then
		dnl defaults
		HAVE_FORTRAN=yes

		AC_DEFINE([_HAVE_FORTRAN_],[1],[with fortran capability])
	else
		HAVE_FORTRAN=no
	fi
	AM_CONDITIONAL([FORTRAN], [test x$FORTRAN = xyes])
	AC_MSG_RESULT($FORTRAN)

	if test "x$FORTRAN" = "xyes"; then
		dnl fortran library  option
		AC_MSG_CHECKING(for fortran library)
		AC_ARG_WITH([fortran-lib],
		  AS_HELP_STRING([--with-fortran-lib = options],[fortran options, for ex: "/usr/lib/gfortran.a]),
			[FORTRAN_LIB=$withval],[FORTRAN_LIB=""])

		dnl check that --with-fortran-lib may have been provided
		if test -n "$FORTRAN_LIB" ; then
			dnl check that library provided EXISTS!
		   FORTRAN_DIR=$(echo $FORTRAN_LIB | sed -e "s/-L//g" | awk '{print $[1]}')
			if test -d "$FORTRAN_DIR" || test -f "$FORTRAN_DIR"; then
				FORTRANLIB="$FORTRAN_LIB"
				AC_DEFINE([_HAVE_FORTRAN_],[1],[with FORTRAN in ISSM src])
				AC_SUBST([FORTRANLIB])
			else
			 if test "x$HAVE_MPI" = "xyes"; then
				FORTRANLIB=$(mpif77 -print-file-name="libgfortran.a")
				if test -f "$FORTRANLIB"; then
					 AC_MSG_ERROR([fortran library provided ($FORTRAN_LIB) does not exist, MPI suggests the following library: $FORTRANLIB]);
				fi
			 fi
				AC_MSG_ERROR([frtran library provided ($FORTRAN_LIB$) does not exist!]);
			fi
		fi
		AC_MSG_RESULT(done)
	fi
	dnl }}}
	dnl graphics{{{
	AC_MSG_CHECKING(for graphics library)
	AC_ARG_WITH([graphics-lib],
	  AS_HELP_STRING([--with-graphics-lib = options],[graphics options, for ex: "/usr/X11/lib/libX11.a]),
	  [GRAPHICS_LIB=$withval],[GRAPHICS_LIB=""])

	dnl check that --with-graphics-lib may have been provided
	
	if test -n "$GRAPHICS_LIB" ; then
		dnl check that library provided EXISTS!
		GRAPHICS_DIR=$(echo $GRAPHICS_LIB | sed -e "s/-L//g" | awk '{print $[1]}')
		if test -d "$GRAPHICS_DIR" || test -f "$GRAPHICS_DIR"; then
			HAVE_GRAPHICS=yes
			GRAPHICSLIB="$GRAPHICS_LIB"
			AC_DEFINE([_HAVE_GRAPHICS_],[1],[with GRAPHICS in ISSM src])
			AC_SUBST([GRAPHICSLIB])
		else
			if test -f "$PETSC_ROOT/conf/petscvariables"; then
				GRAPHICSLIB=$(cat $PETSC_ROOT/conf/petscvariables | grep X_LIB)
				AC_MSG_ERROR([graphics library provided ($GRAPHICS_LIB) does not exist, PETSc suggests the following library: $GRAPHICSLIB]);
			fi
			AC_MSG_ERROR([graphics library provided ($GRAPHICS_LIB$) does not exist!]);
		fi
	fi
	AC_MSG_RESULT(done)
	dnl }}}

	dnl Capabilities
	dnl with-kml{{{
	AC_ARG_WITH([kml],
		AS_HELP_STRING([--with-kml = YES],[compile with kml capabilities (default is yes)]),
		[KML=$withval],[KML=yes]) 
	AC_MSG_CHECKING(for kml capability compilation)

	if test "x$KML" = "xyes"; then
		HAVE_KML=yes
		AC_DEFINE([_HAVE_KML_],[1],[with kml capability])
	else
		HAVE_KML=no
	fi
	AM_CONDITIONAL([KML], [test x$HAVE_KML = xyes])
	AC_MSG_RESULT($HAVE_KML)
	dnl }}}
	dnl with-kriging{{{
	AC_ARG_WITH([kriging],
		AS_HELP_STRING([--with-kriging = YES],[compile with kriging capabilities (default is yes)]),
		[KRIGING=$withval],[KRIGING=yes]) 
	AC_MSG_CHECKING(for kriging capability compilation)

	if test "x$KRIGING" = "xyes"; then
		HAVE_KRIGING=yes
		AC_DEFINE([_HAVE_KRIGING_],[1],[with kriging capability])
	else
		HAVE_KRIGING=no
	fi
	AM_CONDITIONAL([KRIGING], [test x$HAVE_KRIGING = xyes])
	AC_MSG_RESULT($HAVE_KRIGING)
	dnl }}}
	dnl with-steadystate{{{
	AC_ARG_WITH([steadystate],
		AS_HELP_STRING([--with-steadystate = YES],[compile with steadystate capabilities (default is yes)]),
		[STEADYSTATE=$withval],[STEADYSTATE=yes]) 
	AC_MSG_CHECKING(for steadystate capability compilation)

	if test "x$STEADYSTATE" = "xyes"; then

		dnl defaults
		HAVE_STEADYSTATE=yes

		AC_DEFINE([_HAVE_STEADYSTATE_],[1],[with steadystate capability])
	else
		HAVE_STEADYSTATE=no
	fi
	AM_CONDITIONAL([STEADYSTATE], [test x$HAVE_STEADYSTATE = xyes])
	AC_MSG_RESULT($HAVE_STEADYSTATE)
	dnl }}}
	dnl with-transient{{{
	AC_ARG_WITH([transient],
		AS_HELP_STRING([--with-transient = YES], [compile with transient capabilities (default is yes)]),
		[TRANSIENT=$withval],[TRANSIENT=yes]) 
	AC_MSG_CHECKING(for transient capability compilation)

	if test "x$TRANSIENT" = "xyes"; then

		dnl defaults
		HAVE_TRANSIENT=yes

		AC_DEFINE([_HAVE_TRANSIENT_],[1],[with transient capability])
	else
		HAVE_TRANSIENT=no
	fi
	AM_CONDITIONAL([TRANSIENT], [test x$HAVE_TRANSIENT = xyes])
	AC_MSG_RESULT($HAVE_TRANSIENT)
	dnl }}}
	dnl with-thermal{{{
	AC_ARG_WITH([thermal],
		AS_HELP_STRING([--with-thermal = YES], [compile with thermal capabilities (default is yes)]),
		[THERMAL=$withval],[THERMAL=yes]) 
	AC_MSG_CHECKING(for thermal capability compilation)

	if test "x$THERMAL" = "xyes"; then

		dnl defaults
		HAVE_THERMAL=yes

		AC_DEFINE([_HAVE_THERMAL_],[1],[with thermal capability])
	else
		HAVE_THERMAL=no
	fi
	AM_CONDITIONAL([THERMAL], [test x$HAVE_THERMAL = xyes])
	AC_MSG_RESULT($HAVE_THERMAL)
	dnl }}}
	dnl with-prognostic{{{
	AC_ARG_WITH([prognostic],
		AS_HELP_STRING([--with-prognostic = YES], [compile with prognostic capabilities (default is yes)]),
		[PROGNOSTIC=$withval],[PROGNOSTIC=yes]) 
	AC_MSG_CHECKING(for prognostic capability compilation)

	if test "x$PROGNOSTIC" = "xyes"; then

		dnl defaults
		HAVE_PROGNOSTIC=yes

		AC_DEFINE([_HAVE_PROGNOSTIC_],[1],[with prognostic capability])
	else
		HAVE_PROGNOSTIC=no
	fi
	AM_CONDITIONAL([PROGNOSTIC], [test x$HAVE_PROGNOSTIC = xyes])
	AC_MSG_RESULT($HAVE_PROGNOSTIC)
	dnl }}}
	dnl with-control{{{
	AC_ARG_WITH([control],
		AS_HELP_STRING([--with-control = YES], [compile with control capabilities (default is yes)]),
		[CONTROL=$withval],[CONTROL=yes]) 
	AC_MSG_CHECKING(for control capability compilation)

	if test "x$CONTROL" = "xyes"; then

		dnl defaults
		HAVE_CONTROL=yes

		AC_DEFINE([_HAVE_CONTROL_],[1],[with control capability])
	else
		HAVE_CONTROL=no
	fi
	AM_CONDITIONAL([CONTROL], [test x$HAVE_CONTROL = xyes])
	AC_MSG_RESULT($HAVE_CONTROL)
	dnl }}}
	dnl with-hydrology{{{
	AC_ARG_WITH([hydrology],
		AS_HELP_STRING([--with-hydrology = YES], [compile with hydrology capabilities (default is yes)]),
		[HYDROLOGY=$withval],[HYDROLOGY=yes]) 
	AC_MSG_CHECKING(for hydrology capability compilation)

	if test "x$HYDROLOGY" = "xyes"; then

		dnl defaults
		HAVE_HYDROLOGY=yes

		AC_DEFINE([_HAVE_HYDROLOGY_],[1],[with hydrology capability])
	else
		HAVE_HYDROLOGY=no
	fi
	AM_CONDITIONAL([HYDROLOGY], [test x$HAVE_HYDROLOGY = xyes])
	AC_MSG_RESULT($HAVE_HYDROLOGY)
	dnl }}}
	dnl with-diagnostic{{{
	AC_ARG_WITH([diagnostic],
		AS_HELP_STRING([--with-diagnostic = YES], [compile with diagnostic capabilities (default is yes)]),
		[DIAGNOSTIC=$withval],[DIAGNOSTIC=yes]) 
	AC_MSG_CHECKING(for diagnostic capability compilation)

	if test "x$DIAGNOSTIC" = "xyes"; then

		dnl defaults
		HAVE_DIAGNOSTIC=yes

		AC_DEFINE([_HAVE_DIAGNOSTIC_],[1],[with diagnostic capability])
	else
		HAVE_DIAGNOSTIC=no
	fi
	AM_CONDITIONAL([DIAGNOSTIC], [test x$HAVE_DIAGNOSTIC = xyes])
	AC_MSG_RESULT($HAVE_DIAGNOSTIC)
	dnl }}}
	dnl with-balanced{{{
	AC_ARG_WITH([balanced],
		AS_HELP_STRING([--with-balanced = YES], [compile with balanced capabilities (default is yes)]),
		[BALANCED=$withval],[BALANCED=yes]) 
	AC_MSG_CHECKING(for balanced capability compilation)

	if test "x$BALANCED" = "xyes"; then

		dnl defaults
		HAVE_BALANCED=yes

		AC_DEFINE([_HAVE_BALANCED_],[1],[with balanced capability])
	else
		HAVE_BALANCED=no
	fi
	AM_CONDITIONAL([BALANCED], [test x$HAVE_BALANCED = xyes])
	AC_MSG_RESULT($HAVE_BALANCED)
	dnl }}}
	dnl with-responses{{{
	AC_ARG_WITH([responses],
		AS_HELP_STRING([--with-responses = YES], [compile with responses capabilities (default is yes)]),
		[RESPONSES=$withval],[RESPONSES=yes]) 
	AC_MSG_CHECKING(for responses capability compilation)

	if test "x$RESPONSES" = "xyes"; then

		dnl defaults
		HAVE_RESPONSES=yes

		AC_DEFINE([_HAVE_RESPONSES_],[1],[with responses capability])
	else
		HAVE_RESPONSES=no
	fi
	AM_CONDITIONAL([RESPONSES], [test x$HAVE_RESPONSES = xyes])
	AC_MSG_RESULT($HAVE_RESPONSES)
	dnl }}}
	dnl with-slope{{{
	AC_ARG_WITH([slope],
		AS_HELP_STRING([--with-slope = YES], [compile with slope capabilities (default is yes)]),
		[SLOPE=$withval],[SLOPE=yes]) 
	AC_MSG_CHECKING(for slope capability compilation)

	if test "x$SLOPE" = "xyes"; then

		dnl defaults
		HAVE_SLOPE=yes

		AC_DEFINE([_HAVE_SLOPE_],[1],[with slope capability])
	else
		HAVE_SLOPE=no
	fi
	AM_CONDITIONAL([SLOPE], [test x$HAVE_SLOPE = xyes])
	AC_MSG_RESULT($HAVE_SLOPE)
	dnl }}}
	dnl with-groundingline{{{
	AC_ARG_WITH([groundingline],
		AS_HELP_STRING([--with-groundingline = YES], [compile with groundingline capabilities (default is yes)]),
		[GROUNDINGLINE=$withval],[GROUNDINGLINE=yes]) 
	AC_MSG_CHECKING(for groundingline capability compilation)

	if test "x$GROUNDINGLINE" = "xyes"; then

		dnl defaults
		HAVE_GROUNDINGLINE=yes

		AC_DEFINE([_HAVE_GROUNDINGLINE_],[1],[with groundingline capability])
	else
		HAVE_GROUNDINGLINE=no
	fi
	AM_CONDITIONAL([GROUNDINGLINE], [test x$HAVE_GROUNDINGLINE = xyes])
	AC_MSG_RESULT($HAVE_GROUNDINGLINE)
	dnl }}}
	dnl with-rifts{{{
	AC_ARG_WITH([rifts],
		AS_HELP_STRING([--with-rifts = YES], [compile with rifts capabilities (default is yes)]),
		[RIFTS=$withval],[RIFTS=yes]) 
	AC_MSG_CHECKING(for rifts capability compilation)

	if test "x$RIFTS" = "xyes"; then

		dnl defaults
		HAVE_RIFTS=yes

		AC_DEFINE([_HAVE_RIFTS_],[1],[with rifts capability])
	else
		HAVE_RIFTS=no
	fi
	AM_CONDITIONAL([RIFTS], [test x$HAVE_RIFTS = xyes])
	AC_MSG_RESULT($HAVE_RIFTS)
	dnl }}}
	dnl math77{{{
	AC_ARG_WITH([math77-dir],
		AS_HELP_STRING([--with-math77-dir=DIR], [math77 root directory.]),
		[MATH77_ROOT=$withval],[MATH77_ROOT=""]) 
	AC_MSG_CHECKING(for math77)

	if test -d "$MATH77_ROOT"; then

	  HAVE_MATH77=yes
	  MATH77LIB="-L$MATH77_ROOT/ -lmath77"

	  AC_DEFINE([_HAVE_MATH77_],[1],[with math77 in ISSM src])
	  AC_SUBST([MATH77LIB])
	else
		HAVE_MATH77=no
	fi
	AC_MSG_RESULT($HAVE_MATH77)
	dnl }}}
	dnl with-gia{{{
	AC_ARG_WITH([gia],
		AS_HELP_STRING([--with-gia = YES], [compile with gia capabilities (default is yes)]),
		[GIA=$withval],[GIA=no]) 
	AC_MSG_CHECKING(for gia capability compilation)

	if test "x$GIA" = "xyes"; then
	  
	  if test "x$HAVE_MATH77" = "xno"; then
		  HAVE_GIA=no
		  AC_MSG_ERROR([gia requires compilation of math77 library! Reconfigure with --with-math77 option on]);
	  else
		dnl defaults
		HAVE_GIA=yes
		AC_DEFINE([_HAVE_GIA_],[1],[with gia capability])
	  fi

	else
		HAVE_GIA=no
	fi
	AM_CONDITIONAL([GIA], [test x$HAVE_GIA = xyes])
	AC_MSG_RESULT($HAVE_GIA)
	dnl }}}
	dnl with-ios{{{
	AC_ARG_WITH([ios],
		AS_HELP_STRING([--with-ios = YES], [compile with iOS capabilities (default is no, alternatives are yes)]),
		[IOS=$withval],[IOS=no]) 
	AC_MSG_CHECKING(for iOS compilation)

	if test "x$IOS" = "xyes"; then
		dnl defaults
		HAVE_IOS=yes

		AC_DEFINE([_HAVE_IOS_],[1],[with android capability])
	elif test "x$IOS" = "xno"; then
		HAVE_IOS=no
	else
	  AC_MSG_ERROR([--with-ios should be either no or yes])
	fi
	AM_CONDITIONAL([IOS], [test x$HAVE_IOS != xno])
	AC_DEFINE([_HAVE_IOS_],[1],[with ios.])
	AC_MSG_RESULT($HAVE_IOS)
	dnl }}}
	dnl with-android{{{
	AC_ARG_WITH([android],
		AS_HELP_STRING([--with-android = EXE], [compile with android capabilities (default is no, alternatives are exe and jni)]),
		[ANDROID=$withval],[ANDROID=no]) 
	AC_MSG_CHECKING(for android capability compilation)

	if test "x$ANDROID" = "xjni"; then

		dnl defaults
		HAVE_ANDROID=jni
		AC_DEFINE([_HAVE_ANDROID_],[1],[with android capability])
		AC_DEFINE([_HAVE_ANDROID_JNI_],[1],[with android jni])
	elif test "x$ANDROID" = "xexe"; then
		dnl defaults
		HAVE_ANDROID=exe

		AC_DEFINE([_HAVE_ANDROID_],[1],[with android capability])
	elif test "x$ANDROID" = "xno"; then
		HAVE_ANDROID=no
	else
	  AC_MSG_ERROR([--with-android should be either no, exe or jni])
	fi
	AM_CONDITIONAL([ANDROID], [test x$HAVE_ANDROID != xno])
	AM_CONDITIONAL([ANDROIDJNI], [test x$HAVE_ANDROID = xjni])
	AM_CONDITIONAL([ANDROIDEXE], [test x$HAVE_ANDROID = xexe])
	AC_DEFINE([_HAVE_ANDROID_],[1],[with android.])
	AC_MSG_RESULT($HAVE_ANDROID)
	dnl }}}
	dnl with-android-ndk{{{
	AC_ARG_WITH([android-ndk],
	  AS_HELP_STRING([--with-android-ndk=DIR], [android-ndk root directory.]),
	  [ANDROID_NDK_ROOT=$withval],[ANDROID_NDK_ROOT=""]) 
	AC_MSG_CHECKING(with android ndk)
	
	if test -d "$ANDROID_NDK_ROOT"; then
		dnl defaults
		HAVE_ANDROID_NDK=yes
		ANDROID_NDKINCL="-I$ANDROID_NDK_ROOT/arm-linux-android-install/sysroot/usr/include"

		AC_DEFINE([_HAVE_ANDROID_NDK_],[1],[with android ndk in ISSM src])
		AC_SUBST([ANDROID_NDKINCL])
	else
		HAVE_ANDROID_NDK=no
	fi
	AC_MSG_RESULT($HAVE_ANDROID_NDK)
	dnl }}}
	dnl with-3d{{{
	AC_ARG_WITH([3d],
		AS_HELP_STRING([--with-3d = YES], [compile with 3d capabilities (default is yes)]),
		[THREED=$withval],[THREED=yes]) 
	AC_MSG_CHECKING(for 3d capability compilation)

	if test "x$THREED" = "xyes"; then

		dnl defaults
		HAVE_3D=yes

		AC_DEFINE([_HAVE_3D_],[1],[with 3d capability])
	else
		HAVE_3D=no
	fi
	AM_CONDITIONAL([THREED], [test x$HAVE_3D = xyes])
	AC_MSG_RESULT($HAVE_3D)
	dnl }}}
	dnl checks{{{
	AC_MSG_CHECKING(consistency between all libraries)
	
	dnl check that if petsc is requested , mpi should be specified
	if test "$HAVE_PETSC" = "yes" ; then
		if test "$HAVE_MPI" = "NO";  then
			AC_MSG_ERROR([petsc requires mpi!]);
		fi
	fi

	dnl check that we have either python or matlab support if we compile the modules
	if test "$MODULES_VALUE" = "yes"  && test "$HAVE_MATLAB" = "no" && test "$HAVE_PYTHON" = "no"; then
		AC_MSG_ERROR([need at least python or matlab support to compile modules (or use --with-modules=no)]);
	fi

	dnl check that if we have MPI, we have metis
	if test "$HAVE_METIS" = "yes"  && test "$HAVE_MPI" = "no" ; then
	AC_MSG_ERROR([need mpi if using the metis partitioner!]);
	fi

	AC_MSG_RESULT(done)
	dnl }}}

	dnl other options
	dnl optimization{{{
	dnl bypass standard optimization -g -O2 ? 
	AC_ARG_WITH([cxxoptflags],
	  AS_HELP_STRING([--with-cxxoptflags = CXXOPTFLAGS], [optimization using CXX flags, ex: --with-cxxoptflags=-march=opteron -O3]),
	  [CXXOPTFLAGS=$withval],[CXXOPTFLAGS="-g -O2"]) 
	AC_MSG_CHECKING(for c++ optimization flags)
	AC_SUBST([CXXOPTFLAGS])
	AC_MSG_RESULT(done)

	dnl }}}
	dnl multithreading{{{
	AC_ARG_WITH([numthreads],
	  AS_HELP_STRING([--with-numthreads = NUMTHREADS_VALUE],[numthreads, default is 1. ]),
	  [NUMTHREADS_VALUE=$withval],[NUMTHREADS_VALUE=1])
	AC_MSG_CHECKING(for number of threads)
	dnl defaults
	MULTITHREADING=no
	MULTITHREADINLIB=""
	if test "$NUMTHREADS_VALUE" != "1"; then
		
		MULTITHREADINGLIB="-lpthread -lrt"
		case "${host_os}" in
		*cygwin*)
		MULTITHREADINGLIB="-lpthread -lrt"
		;;
		*linux*)
		MULTITHREADINGLIB="-lpthread -lrt"
		;;
		*darwin*)
		MULTITHREADINGLIB="-lpthread"
		;;
		esac

		AC_DEFINE([_MULTITHREADING_],[1],[with numthreads enabled])
		AC_DEFINE_UNQUOTED([_NUMTHREADS_],[$NUMTHREADS_VALUE],[number of threads])
	fi
	AC_SUBST([MULTITHREADINGLIB])
	AC_MSG_RESULT($NUMTHREADS_VALUE) 
	dnl }}}
	dnl 64bit {{{
	AC_ARG_WITH([64bit-indices],
	  AS_HELP_STRING([--with-64bit-indices = bool], [use 64 bit integers, default 0, ex: --with-64bit-indices=1]),
	  [USE_64BIT_INDICES=$withval],[USE_64BIT_INDICES=0]) 
	AC_MSG_CHECKING(for 64 bit indices)

	if test "$USE_64BIT_INDICES" == "1"; then
	AC_DEFINE([ISSM_USE_64BIT_INDICES],[1],[with 64 bits indices])
	else
	AC_DEFINE([ISSM_USE_64BIT_INDICES],[0],[with 64 bits indices])
	fi
	AC_MSG_RESULT($USE_64BIT_INDICES)
	dnl }}}
])
