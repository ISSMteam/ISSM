AM_CPPFLAGS = @NEOPZINCL@ @DAKOTAINCL@ @BOOSTINCL@ @MATLABINCL@ @PETSCINCL@ @MPIINCL@ @SPOOLESINCL@ @PARMETISINCL@ @METISINCL@ @TRIANGLEINCL@ @CHACOINCL@ @AMPIINCL@ @MEDIPACKINCL@ @ADJOINTPETSCINCL@ @CODIPACKINCL@ @PROJINCL@
AUTOMAKE_OPTIONS = subdir-objects

EXEEXT=$(MATLABWRAPPEREXT)

#define prefix (from http://www.gnu.org/software/autoconf/manual/autoconf-2.67/html_node/Defining-Directories.html)
AM_CPPFLAGS += -DISSM_PREFIX='"$(prefix)"'

#matlab io{{{
lib_LTLIBRARIES = libISSMMatlab.la

io_sources = \
	./io/CheckNumMatlabArguments.cpp \
	./io/FetchMatlabData.cpp \
	./io/WriteMatlabData.cpp

ALL_CXXFLAGS = -fPIC -D_WRAPPERS_

if MSYS2
ALL_CXXFLAGS += -D_DO_NOT_LOAD_GLOBALS_ # Cannot have undefined symbols under MSYS2
ALL_CXXFLAGS += -D_USE_MATH_DEFINES # Need to open fence to math.h definitions when --std=c++ is used (which is the default)
endif

ALL_CXXFLAGS += $(CXXFLAGS)

libISSMMatlab_la_SOURCES = $(io_sources)
libISSMMatlab_la_CXXFLAGS = ${ALL_CXXFLAGS}
#}}}
#api io{{{
lib_LTLIBRARIES += libISSMApi_matlab.la

if !MSYS2
api_sources= ./io/ApiPrintf.cpp

libISSMApi_matlab_la_SOURCES = $(api_sources)
libISSMApi_matlab_la_CXXFLAGS = ${ALL_CXXFLAGS}
endif
#}}}
#Wrappers {{{
lib_LTLIBRARIES += \
	BamgConvertMesh_matlab.la \
	BamgMesher_matlab.la \
	BamgTriangulate_matlab.la \
	ContourToMesh_matlab.la \
	ContourToNodes_matlab.la \
	DistanceToMaskBoundary_matlab.la \
	ElementConnectivity_matlab.la \
	ExpSimplify_matlab.la \
	ExpToLevelSet_matlab.la \
	InterpFromGridToMesh_matlab.la \
	InterpFromGrid.la \
	InterpFromMesh2d_matlab.la \
	InterpFromMeshToGrid_matlab.la \
	InterpFromMeshToMesh2d_matlab.la \
	InterpFromMeshToMesh3d_matlab.la \
	IssmConfig_matlab.la \
	MeshPartition_matlab.la \
	MeshProfileIntersection_matlab.la \
	NodeConnectivity_matlab.la \
	PointCloudFindNeighbors_matlab.la \
	ProcessRifts_matlab.la \
	PropagateFlagsFromConnectivity_matlab.la \
	Triangle_matlab.la

if CHACO
lib_LTLIBRARIES += Chaco_matlab.la
endif
if KRIGING
lib_LTLIBRARIES += Kriging_matlab.la
endif
if PROJ
lib_LTLIBRARIES += CoordTransform_matlab.la
endif
#}}}

# Dependencies {{{
deps = $(MATHLIB)

# Triangle
CXXFLAGS_FOR_TRI = -DTRILIBRARY -DANSI_DECLARATORS -DNO_TIMER

# MATLAB
CXXFLAGS_FOR_MEX = -D_HAVE_MATLAB_MODULES_ $(MEXOPTFLAGS) $(MEXCXXFLAGS)
CPPFLAGS_FOR_MEX = $(MATLABINCL)
LDFLAGS_FOR_MEX = -module -shrext $(MEXEXT) -no-undefined
LIBADD_FOR_MEX = $(MEXLINKFLAGS) $(MEXLIB)

AM_CPPFLAGS += ${CXXFLAGS_FOR_MEX}
AM_CXXFLAGS = ${CXXFLAGS_FOR_TRI} ${CPPFLAGS_FOR_MEX}
AM_LDFLAGS = ${LDFLAGS_FOR_MEX}

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

deps += ./libISSMMatlab.la ../../c/libISSMCore.la ../../c/libISSMModules.la

if !MSYS2
deps += ./libISSMApi_matlab.la
endif

if ADOLC
deps += $(ADOLCLIB)
endif

if FORTRAN
deps += $(FLIBS) $(FORTRANLIB)
if !STANDALONE_LIBRARIES
MEXLIB += $(FLIBS) $(FORTRANLIB)
endif
endif

deps += ${LIBADD_FOR_MEX}

# Optimization flags
AM_CXXFLAGS += $(CXXFLAGS)
#}}}
# Module sources and dependencies {{{
libISSMMatlab_la_LIBADD = ./../../c/libISSMCore.la ./../../c/libISSMModules.la $(PETSCLIB) $(MUMPSLIB) $(SCALAPACKLIB) $(BLASLAPACKLIB) $(PARMETISLIB) $(METISLIB) $(HDF5LIB) $(MPILIB) $(NEOPZLIB) $(MKLLIB) $(GSLLIB) $(MATHLIB) $(MEXLIB)

if VERSION
libISSMMatlab_la_LDFLAGS =
libISSMApi_matlab_la_LDFLAGS =
else
libISSMMatlab_la_LDFLAGS = -avoid-version
libISSMApi_matlab_la_LDFLAGS = -avoid-version
endif

if STANDALONE_LIBRARIES
if !MSYS2
libISSMMatlab_la_LDFLAGS += -static
libISSMApi_matlab_la_LDFLAGS += -static
endif
deps += $(DAKOTALIB) $(PETSCLIB) $(MUMPSLIB) $(SCALAPACKLIB) $(BLASLAPACKLIB) $(PARMETISLIB) $(METISLIB) $(HDF5LIB) $(TAOLIB) $(M1QN3LIB) $(SEMICLIB) $(PLAPACKLIB) $(SUPERLULIB) $(SPOOLESLIB) $(TRIANGLELIB) $(BLACSLIB) $(HYPRELIB) $(SPAILIB) $(PROMETHEUSLIB) $(PASTIXLIB) $(MLLIB) $(MKLLIB) $(MPILIB) $(NEOPZLIB) $(MATHLIB) $(GRAPHICSLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(ADOLCLIB) $(AMPILIB) $(OSLIBS) ${LIBADD_FOR_MEX}
endif

libISSMApi_matlab_la_LIBADD = $(PETSCLIB) $(MUMPSLIB) $(SCALAPACKLIB) $(BLASLAPACKLIB) $(PARMETISLIB) $(METISLIB) $(HDF5LIB) $(MPILIB) $(NEOPZLIB) $(GSLLIB) $(MATHLIB) $(MEXLIB)

BamgConvertMesh_matlab_la_SOURCES = ../BamgConvertMesh/BamgConvertMesh.cpp
BamgConvertMesh_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
BamgConvertMesh_matlab_la_LIBADD = ${deps}

BamgMesher_matlab_la_SOURCES = ../BamgMesher/BamgMesher.cpp
BamgMesher_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
BamgMesher_matlab_la_LIBADD = ${deps}

BamgTriangulate_matlab_la_SOURCES = ../BamgTriangulate/BamgTriangulate.cpp
BamgTriangulate_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
BamgTriangulate_matlab_la_LIBADD = ${deps}

if CHACO
Chaco_matlab_la_SOURCES = ../Chaco/Chaco.cpp
Chaco_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
Chaco_matlab_la_LIBADD = ${deps} $(CHACOLIB)
endif

ContourToMesh_matlab_la_SOURCES = ../ContourToMesh/ContourToMesh.cpp
ContourToMesh_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
ContourToMesh_matlab_la_LIBADD = ${deps}

ContourToNodes_matlab_la_SOURCES = ../ContourToNodes/ContourToNodes.cpp
ContourToNodes_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
ContourToNodes_matlab_la_LIBADD = ${deps}

if PROJ
CoordTransform_matlab_la_SOURCES = ../CoordTransform/CoordTransform.cpp
CoordTransform_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
CoordTransform_matlab_la_LIBADD = ${deps} $(PROJLIB)
endif

DistanceToMaskBoundary_matlab_la_SOURCES = ../DistanceToMaskBoundary/DistanceToMaskBoundary.cpp
DistanceToMaskBoundary_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
DistanceToMaskBoundary_matlab_la_LIBADD = ${deps}

ElementConnectivity_matlab_la_SOURCES = ../ElementConnectivity/ElementConnectivity.cpp
ElementConnectivity_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
ElementConnectivity_matlab_la_LIBADD = ${deps}

ExpToLevelSet_matlab_la_SOURCES = ../ExpToLevelSet/ExpToLevelSet.cpp
ExpToLevelSet_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
ExpToLevelSet_matlab_la_LIBADD = ${deps}

ExpSimplify_matlab_la_SOURCES = ../ExpSimplify/ExpSimplify.cpp
ExpSimplify_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
ExpSimplify_matlab_la_LIBADD = ${deps}

#Very specific case for this mex: do not compile like other modules
InterpFromGrid_la_SOURCES = ../InterpFromGrid/InterpFromGrid.cpp
InterpFromGrid_la_CXXFLAGS = ${MEXLIB}
InterpFromGrid_la_LIBADD = 

InterpFromGridToMesh_matlab_la_SOURCES = ../InterpFromGridToMesh/InterpFromGridToMesh.cpp
InterpFromGridToMesh_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromGridToMesh_matlab_la_LIBADD = ${deps}

InterpFromMesh2d_matlab_la_SOURCES = ../InterpFromMesh2d/InterpFromMesh2d.cpp
InterpFromMesh2d_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMesh2d_matlab_la_LIBADD = ${deps}

InterpFromMeshToGrid_matlab_la_SOURCES = ../InterpFromMeshToGrid/InterpFromMeshToGrid.cpp
InterpFromMeshToGrid_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMeshToGrid_matlab_la_LIBADD = ${deps}

InterpFromMeshToMesh2d_matlab_la_SOURCES = ../InterpFromMeshToMesh2d/InterpFromMeshToMesh2d.cpp
InterpFromMeshToMesh2d_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMeshToMesh2d_matlab_la_LIBADD = ${deps}

InterpFromMeshToMesh3d_matlab_la_SOURCES = ../InterpFromMeshToMesh3d/InterpFromMeshToMesh3d.cpp
InterpFromMeshToMesh3d_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
InterpFromMeshToMesh3d_matlab_la_LIBADD = ${deps}

IssmConfig_matlab_la_SOURCES = ../IssmConfig/IssmConfig.cpp
IssmConfig_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
IssmConfig_matlab_la_LIBADD = ${deps}

Kriging_matlab_la_SOURCES = ../Kriging/Kriging.cpp
Kriging_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
Kriging_matlab_la_LIBADD = ${deps}

MeshPartition_matlab_la_SOURCES = ../MeshPartition/MeshPartition.cpp
MeshPartition_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
MeshPartition_matlab_la_LIBADD = ${deps}

MeshProfileIntersection_matlab_la_SOURCES = ../MeshProfileIntersection/MeshProfileIntersection.cpp
MeshProfileIntersection_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
MeshProfileIntersection_matlab_la_LIBADD = ${deps}

NodeConnectivity_matlab_la_SOURCES = ../NodeConnectivity/NodeConnectivity.cpp
NodeConnectivity_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
NodeConnectivity_matlab_la_LIBADD = ${deps}

PointCloudFindNeighbors_matlab_la_SOURCES = ../PointCloudFindNeighbors/PointCloudFindNeighbors.cpp
PointCloudFindNeighbors_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
PointCloudFindNeighbors_matlab_la_LIBADD = ${deps}

PropagateFlagsFromConnectivity_matlab_la_SOURCES = ../PropagateFlagsFromConnectivity/PropagateFlagsFromConnectivity.cpp
PropagateFlagsFromConnectivity_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
PropagateFlagsFromConnectivity_matlab_la_LIBADD = ${deps}

Triangle_matlab_la_SOURCES = ../Triangle/Triangle.cpp
Triangle_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
Triangle_matlab_la_LIBADD = ${deps} $(TRIANGLELIB)

ProcessRifts_matlab_la_SOURCES = ../ProcessRifts/ProcessRifts.cpp
ProcessRifts_matlab_la_CXXFLAGS = ${AM_CXXFLAGS}
ProcessRifts_matlab_la_LIBADD = ${deps}
#}}}
