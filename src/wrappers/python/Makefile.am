AM_CPPFLAGS = @NEOPZINCL@ @DAKOTAINCL@ @BOOSTINCL@ @PYTHONINCL@ @PYTHON_NUMPYINCL@ @PETSCINCL@ @MPIINCL@ @SPOOLESINCL@ @PARMETISINCL@ @METISINCL@ @TRIANGLEINCL@ @CHACOINCL@ @AMPIINCL@ @MEDIPACKINCL@ @ADJOINTPETSCINCL@ @CODIPACKINCL@ @PROJINCL@
AUTOMAKE_OPTIONS = subdir-objects

#define prefix (from http://www.gnu.org/software/autoconf/manual/autoconf-2.67/html_node/Defining-Directories.html)
AM_CPPFLAGS += -DISSM_PREFIX='"$(prefix)"'

#python io{{{
lib_LTLIBRARIES = libISSMPython.la

io_sources = \
	./io/CheckNumPythonArguments.cpp \
	./io/FetchPythonData.cpp \
	./io/WritePythonData.cpp

ALL_CXXFLAGS = -fPIC -D_WRAPPERS_

if MSYS2
ALL_CXXFLAGS += -D_DO_NOT_LOAD_GLOBALS_ # Cannot have undefined symbols under MSYS2
ALL_CXXFLAGS += -D_USE_MATH_DEFINES # Need to open fence to math.h definitions when --std=c++ is used (which is the default)
endif

ALL_CXXFLAGS += $(CXXFLAGS)

libISSMPython_la_SOURCES = $(io_sources)
libISSMPython_la_CXXFLAGS= ${ALL_CXXFLAGS}
#}}}
#api io{{{
lib_LTLIBRARIES += libISSMApi_python.la

if !MSYS2
api_sources= ./io/ApiPrintf.cpp

libISSMApi_python_la_SOURCES = $(api_sources)
libISSMApi_python_la_CXXFLAGS = ${ALL_CXXFLAGS}
endif
#}}}
#Wrappers {{{
lib_LTLIBRARIES += \
	BamgConvertMesh_python.la \
	BamgMesher_python.la \
	BamgTriangulate_python.la \
	ContourToMesh_python.la \
	ContourToNodes_python.la \
	ElementConnectivity_python.la \
	ExpToLevelSet_python.la \
	InterpFromGridToMesh_python.la \
	InterpFromMesh2d_python.la \
	InterpFromMeshToGrid_python.la \
	InterpFromMeshToMesh2d_python.la \
	InterpFromMeshToMesh3d_python.la \
	IssmConfig_python.la \
	MeshPartition_python.la \
	MeshProfileIntersection_python.la \
	NodeConnectivity_python.la \
	Triangle_python.la \
	ProcessRifts_python.la

if CHACO
lib_LTLIBRARIES += Chaco_python.la
endif
#}}}
# Dependencies {{{
deps = $(MATHLIB)

# Triangle
CXXFLAGS_FOR_TRI = -DTRILIBRARY -DANSI_DECLARATORS -DNO_TIMER

# Python
LDFLAGS_FOR_PY = -module $(PYTHONLINK) -shrext $(PYTHONWRAPPEREXT) -no-undefined

LDFLAGS_FOR_PY += --no-warnings
CXXFLAGS_FOR_PY = -D_HAVE_PYTHON_MODULES_

if PYTHON3
CXXFLAGS_FOR_PY += -DNPY_NO_DEPRECATED_API
endif

AM_CXXFLAGS = ${CXXFLAGS_FOR_TRI} ${CXXFLAGS_FOR_PY}
AM_LDFLAGS = ${LDFLAGS_FOR_PY}

if !MSYS2
AM_LDFLAGS += --export-dynamic -rdynamic
endif

if !VERSION
AM_LDFLAGS += -avoid-version
endif

# NOTE:
# - On Linux, We probably do not need the -static flag as long as we only
#	generate static libraries for external packages. Dynamic system libraries
#	will be linked to, whether we like it or not, if no static version is
#	available.
# - On macOS, static linking of binaries is not supported.
#
if STANDALONE_MODULES
if MAC
AM_LDFLAGS += -Wl,-rpath,'@rpath'
AM_LDFLAGS += -Wl,-rpath,/System/Volumes/Data/Applications/Xcode.app/Contents/Developer/Library/Frameworks
AM_LDFLAGS += -Wl,-rpath,/System/Volumes/Data/Library/Developer/CommandLineTools/Library/Frameworks
else
AM_LDFLAGS += -Wl,-lbacktrace
if !MSYS2
AM_LDFLAGS += -Wl,-static -Wl,--disable-new-dtags -Wl,-rpath,'$$ORIGIN'
endif
endif
endif

if !STANDALONE_EXECUTABLES
if MAC
if HAVE_FORTRANDIR
AM_LDFLAGS += -Wl,-rpath,$(FORTRANDIR)
endif
endif
endif

if !STANDALONE_LIBRARIES
if MAC
if BOOST
AM_LDFLAGS += -Wl,-rpath,$(BOOSTROOT)/lib
endif
if PROJ
AM_LDFLAGS += -Wl,-rpath,$(PROJROOT)/lib
endif
endif
endif

AM_CXXFLAGS += -fPIC -D_WRAPPERS_

deps += ./libISSMPython.la ../../c/libISSMModules.la ../../c/libISSMCore.la

if !MSYS2
deps += ./libISSMApi_python.la
endif

if ADOLC
deps += $(ADOLCLIB)
endif

if FORTRAN
deps += $(FLIBS) $(FORTRANLIB)
if !STANDALONE_LIBRARIES
PYTHONLIB += $(FLIBS) $(FORTRANLIB)
endif
endif

deps += $(PYTHONLIB)

# Optimization flags
AM_CXXFLAGS += $(CXXFLAGS)
#}}}
# Module sources and dependencies {{{
libISSMPython_la_LIBADD = ./../../c/libISSMCore.la ./../../c/libISSMModules.la $(PETSCLIB) $(MUMPSLIB) $(SCALAPACKLIB) $(BLASLAPACKLIB) $(PARMETISLIB) $(METISLIB) $(HDF5LIB) $(MPILIB) $(NEOPZLIB) $(MKLLIB) $(GSLLIB) $(MATHLIB)

if STANDALONE_LIBRARIES
if !MSYS2
libISSMPython_la_LDFLAGS = -static
libISSMApi_python_la_LDFLAGS = -static
endif
deps += $(DAKOTALIB) $(PETSCLIB) $(MUMPSLIB) $(SCALAPACKLIB) $(BLASLAPACKLIB) $(PARMETISLIB) $(METISLIB) $(HDF5LIB) $(TAOLIB) $(M1QN3LIB) $(SEMICLIB) $(PLAPACKLIB) $(SUPERLULIB) $(SPOOLESLIB) $(TRIANGLELIB) $(BLACSLIB) $(HYPRELIB) $(SPAILIB) $(PROMETHEUSLIB) $(PASTIXLIB) $(MLLIB) $(MKLLIB) $(MPILIB) $(NEOPZLIB) $(MATHLIB) $(GRAPHICSLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(ADOLCLIB) $(AMPILIB) $(OSLIBS) $(PYTHONLIB)
endif

libISSMApi_python_la_LIBADD = $(PETSCLIB) $(MUMPSLIB) $(SCALAPACKLIB) $(BLASLAPACKLIB) $(PARMETISLIB) $(METISLIB) $(HDF5LIB) $(MPILIB) $(NEOPZLIB) $(GSLLIB) $(MATHLIB)

BamgConvertMesh_python_la_SOURCES = ../BamgConvertMesh/BamgConvertMesh.cpp
BamgConvertMesh_python_la_CXXFLAGS = ${AM_CXXFLAGS}
BamgConvertMesh_python_la_LIBADD = ${deps}

BamgMesher_python_la_SOURCES = ../BamgMesher/BamgMesher.cpp
BamgMesher_python_la_CXXFLAGS = ${AM_CXXFLAGS}
BamgMesher_python_la_LIBADD = ${deps}

BamgTriangulate_python_la_SOURCES = ../BamgTriangulate/BamgTriangulate.cpp
BamgTriangulate_python_la_CXXFLAGS = ${AM_CXXFLAGS}
BamgTriangulate_python_la_LIBADD = ${deps}

if CHACO
Chaco_python_la_SOURCES = ../Chaco/Chaco.cpp
Chaco_python_la_CXXFLAGS = ${AM_CXXFLAGS}
Chaco_python_la_LIBADD = ${deps} $(CHACOLIB)
endif

ContourToMesh_python_la_SOURCES = ../ContourToMesh/ContourToMesh.cpp
ContourToMesh_python_la_CXXFLAGS = ${AM_CXXFLAGS}
ContourToMesh_python_la_LIBADD = ${deps}

ContourToNodes_python_la_SOURCES = ../ContourToNodes/ContourToNodes.cpp
ContourToNodes_python_la_CXXFLAGS = ${AM_CXXFLAGS}
ContourToNodes_python_la_LIBADD = ${deps}

ElementConnectivity_python_la_SOURCES = ../ElementConnectivity/ElementConnectivity.cpp
ElementConnectivity_python_la_CXXFLAGS = ${AM_CXXFLAGS}
ElementConnectivity_python_la_LIBADD = ${deps}

ExpToLevelSet_python_la_SOURCES = ../ExpToLevelSet/ExpToLevelSet.cpp
ExpToLevelSet_python_la_CXXFLAGS = ${AM_CXXFLAGS}
ExpToLevelSet_python_la_LIBADD = ${deps}

InterpFromMesh2d_python_la_SOURCES = ../InterpFromMesh2d/InterpFromMesh2d.cpp
InterpFromMesh2d_python_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMesh2d_python_la_LIBADD = ${deps}

InterpFromGridToMesh_python_la_SOURCES = ../InterpFromGridToMesh/InterpFromGridToMesh.cpp
InterpFromGridToMesh_python_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromGridToMesh_python_la_LIBADD = ${deps}

InterpFromMeshToGrid_python_la_SOURCES = ../InterpFromMeshToGrid/InterpFromMeshToGrid.cpp
InterpFromMeshToGrid_python_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMeshToGrid_python_la_LIBADD = ${deps}

InterpFromMeshToMesh2d_python_la_SOURCES = ../InterpFromMeshToMesh2d/InterpFromMeshToMesh2d.cpp
InterpFromMeshToMesh2d_python_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMeshToMesh2d_python_la_LIBADD = ${deps}

InterpFromMeshToMesh3d_python_la_SOURCES = ../InterpFromMeshToMesh3d/InterpFromMeshToMesh3d.cpp
InterpFromMeshToMesh3d_python_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMeshToMesh3d_python_la_LIBADD = ${deps}

IssmConfig_python_la_SOURCES = ../IssmConfig/IssmConfig.cpp
IssmConfig_python_la_CXXFLAGS = ${AM_CXXFLAGS}
IssmConfig_python_la_LIBADD = ${deps}

MeshPartition_python_la_SOURCES = ../MeshPartition/MeshPartition.cpp
MeshPartition_python_la_CXXFLAGS = ${AM_CXXFLAGS}
MeshPartition_python_la_LIBADD = ${deps}

MeshProfileIntersection_python_la_SOURCES = ../MeshProfileIntersection/MeshProfileIntersection.cpp
MeshProfileIntersection_python_la_CXXFLAGS = ${AM_CXXFLAGS}
MeshProfileIntersection_python_la_LIBADD = ${deps}

NodeConnectivity_python_la_SOURCES = ../NodeConnectivity/NodeConnectivity.cpp
NodeConnectivity_python_la_CXXFLAGS = ${AM_CXXFLAGS}
NodeConnectivity_python_la_LIBADD = ${deps}

Triangle_python_la_SOURCES = ../Triangle/Triangle.cpp
Triangle_python_la_CXXFLAGS = ${AM_CXXFLAGS}
Triangle_python_la_LIBADD = ${deps} $(TRIANGLELIB)

ProcessRifts_python_la_SOURCES = ../ProcessRifts/ProcessRifts.cpp
ProcessRifts_python_la_CXXFLAGS = ${AM_CXXFLAGS}
ProcessRifts_python_la_LIBADD = ${deps}
#}}}
