AM_CPPFLAGS = @NEOPZINCL@ @DAKOTAINCL@ @SHAPELIBINCL@ @TRIANGLEINCL@ @PETSCINCL@ @AMPIINCL@ @MEDIPACKINCL@ @MPIINCL@ @PARMETISINCL@ @METISINCL@ @CHACOINCL@ @SCOTCHINCL@ @PLAPACKINCL@ @MKLINCL@ @MUMPSINCL@ @SPAIINCL@ @HYPREINCL@ @PROMETHEUSINCL@ @SUPERLUINCL@ @SPOOLESINCL@ @PASTIXINCL@ @MLINCL@ @TAOINCL@ @ADOLCINCL@ @CODIPACKINCL@ @ADJOINTPETSCINCL@ @GSLINCL@ @BOOSTINCL@ @ESMFINCL@ @PROJINCL@ @MPLAPACKINCL@
AM_FCFLAGS = @SEMICINCL@

AUTOMAKE_OPTIONS = subdir-objects

EXEEXT=$(ISSMEXT)

# Library declaration {{{
lib_LTLIBRARIES = libISSMCore.la
if !MSYS2
lib_LTLIBRARIES += libISSMOverload.la
endif
if WRAPPERS
lib_LTLIBRARIES += libISSMModules.la
endif
#}}}

# Core sources
issm_sources =

# BAMG sources {{{
if BAMG
issm_sources += \
	./bamg/BamgGeom.cpp \
	./bamg/BamgMesh.cpp \
	./bamg/BamgOpts.cpp \
	./bamg/CrackedEdge.cpp \
	./bamg/Curve.cpp \
	./bamg/Edge.cpp \
	./bamg/GeomEdge.cpp \
	./bamg/GeomSubDomain.cpp \
	./bamg/GeomVertex.cpp \
	./bamg/Geometry.cpp \
	./bamg/ListofIntersectionTriangles.cpp \
	./bamg/EigenMetric.cpp \
	./bamg/Metric.cpp \
	./bamg/BamgQuadtree.cpp \
	./bamg/SetOfE4.cpp \
	./bamg/SubDomain.cpp \
	./bamg/AdjacentTriangle.cpp \
	./bamg/Triangle.cpp \
	./bamg/BamgVertex.cpp \
	./bamg/VertexOnEdge.cpp \
	./bamg/VertexOnGeom.cpp \
	./bamg/VertexOnVertex.cpp \
	./bamg/Mesh.cpp \
	./shared/Bamg/BigPrimeNumber.cpp \
	./modules/Bamgx/Bamgx.cpp \
	./modules/BamgConvertMeshx/BamgConvertMeshx.cpp \
	./modules/BamgTriangulatex/BamgTriangulatex.cpp

# Do not include AmrBamg with AD
if ADOLC
issm_sources += \
	./shared/Numerics/isnan.cpp \
	./shared/MemOps/MemOps.cpp
else
if CODIPACK
issm_sources += ./shared/Numerics/isnan.cpp
else
issm_sources += ./classes/AmrBamg.cpp
endif
endif
endif
#}}}
# Core sources {{{
issm_sources += \
	./datastructures/DataSet.cpp \
	./classes/gauss/GaussSeg.cpp \
	./classes/gauss/GaussTria.cpp \
	./classes/gauss/GaussTetra.cpp \
	./classes/gauss/GaussPenta.cpp \
	./classes/IoModel.cpp \
	./classes/FemModel.cpp \
	./classes/Loads/Friction.cpp \
	./classes/Constraints/SpcTransient.cpp \
	./classes/DependentObject.cpp \
	./classes/Contours.cpp \
	./classes/Vertices.cpp \
	./classes/Nodes.cpp \
	./classes/Numberedcostfunction.cpp \
	./classes/Misfit.cpp \
	./classes/Cfsurfacesquare.cpp \
	./classes/Cfsurfacesquaretransient.cpp \
	./classes/Cfdragcoeffabsgrad.cpp \
	./classes/Cfdragcoeffabsgradtransient.cpp \
	./classes/Cfrheologybbarabsgrad.cpp \
	./classes/Cfrheologybbarabsgradtransient.cpp \
	./classes/Cfsurfacelogvel.cpp \
	./classes/Cflevelsetmisfit.cpp \
	./classes/Regionaloutput.cpp \
	./classes/Nodalvalue.cpp \
	./classes/Node.cpp \
	./classes/Vertex.cpp \
	./classes/Hook.cpp \
	./classes/Radar.cpp \
	./classes/BarystaticContributions.cpp \
	./classes/ExternalResults/Results.cpp \
	./classes/Elements/Element.cpp \
	./classes/Elements/Elements.cpp \
	./classes/Elements/ElementHook.cpp \
	./classes/Elements/Seg.cpp \
	./classes/Elements/SegRef.cpp \
	./classes/Elements/Tria.cpp \
	./classes/Elements/TriaRef.cpp \
	./classes/Elements/Tetra.cpp \
	./classes/Elements/TetraRef.cpp \
	./classes/Elements/Penta.cpp \
	./classes/Elements/PentaRef.cpp \
	./classes/Materials/Materials.cpp \
	./classes/Materials/Matice.cpp \
	./classes/Materials/Matlitho.cpp \
	./classes/Materials/Matestar.cpp \
	./classes/Constraints/Constraints.cpp \
	./classes/Constraints/SpcStatic.cpp \
	./classes/Constraints/SpcDynamic.cpp \
	./classes/Loads/Channel.cpp \
	./classes/Loads/Loads.cpp \
	./classes/Loads/Penpair.cpp \
	./classes/Loads/Pengrid.cpp \
	./classes/Loads/Moulin.cpp \
	./classes/Loads/Numericalflux.cpp \
	./classes/Loads/Neumannflux.cpp \
	./classes/matrix/ElementMatrix.cpp \
	./classes/matrix/ElementVector.cpp \
	./classes/Params/Parameters.cpp \
	./classes/Params/BoolParam.cpp \
	./classes/Params/ControlParam.cpp \
	./classes/Params/IntParam.cpp \
	./classes/Params/IntVecParam.cpp \
	./classes/Params/IntMatParam.cpp \
	./classes/Params/DoubleParam.cpp \
	./classes/Params/FileParam.cpp \
	./classes/Params/StringArrayParam.cpp \
	./classes/Params/DoubleMatParam.cpp \
	./classes/Params/DoubleTransientMatParam.cpp \
	./classes/Params/DoubleMatArrayParam.cpp \
	./classes/Params/DoubleVecParam.cpp \
	./classes/Params/StringParam.cpp \
	./classes/Params/MatrixParam.cpp \
	./classes/Params/VectorParam.cpp \
	./classes/Params/TransientParam.cpp \
	./classes/Params/TransientArrayParam.cpp \
	./classes/Params/TransientGriddedFieldParam.cpp \
	./classes/Params/DataSetParam.cpp \
	./classes/Profiler.cpp \
	./shared/Matrix/MatrixUtils.cpp \
	./shared/io/Disk/pfopen.cpp \
	./shared/io/Disk/pfclose.cpp \
	./shared/io/Disk/WriteLockFile.cpp \
	./shared/io/Print/PrintfFunction.cpp \
	./shared/io/Comm/IssmComm.cpp \
	./shared/io/Marshalling/IoCodeConversions.cpp \
	./shared/io/Marshalling/Marshalling.cpp \
	./shared/LatLong/Ll2xyx.cpp \
	./shared/LatLong/Xy2llx.cpp \
	./shared/FSanalyticals/fsanalyticals.cpp \
	./shared/Enum/EnumToStringx.cpp \
	./shared/Enum/StringToEnumx.cpp \
	./shared/Numerics/Verbosity.cpp \
	./shared/Numerics/GaussPoints.cpp \
	./shared/Numerics/cross.cpp \
	./shared/Numerics/cubic.cpp \
	./shared/Numerics/NewtonSolveDnorm.cpp \
	./shared/Numerics/ODE1.cpp \
	./shared/Numerics/extrema.cpp \
	./shared/Numerics/legendre.cpp \
	./shared/Numerics/XZvectorsToCoordinateSystem.cpp \
	./shared/Numerics/Normals.cpp \
	./shared/Numerics/Interpolation.cpp \
	./shared/Exceptions/Exceptions.cpp \
	./shared/Sorting/binary_search.cpp \
	./shared/Elements/Cuffey.cpp \
	./shared/Elements/BuddJacka.cpp \
	./shared/Elements/CuffeyTemperate.cpp \
	./shared/Elements/StressIntensityIntegralWeight.cpp \
	./shared/Elements/Paterson.cpp \
	./shared/Elements/Arrhenius.cpp \
	./shared/Elements/NyeCO2.cpp \
	./shared/Elements/NyeH2O.cpp \
	./shared/Elements/LliboutryDuval.cpp \
	./shared/Elements/PrintArrays.cpp \
	./shared/Elements/PddSurfaceMassBalance.cpp \
	./shared/Elements/PddSurfaceMassBalanceSicopolis.cpp \
	./shared/Elements/ComputeDelta18oTemperaturePrecipitation.cpp \
	./shared/Elements/ComputeMungsmTemperaturePrecipitation.cpp \
	./shared/Elements/ComputeD18OTemperaturePrecipitationFromPD.cpp \
	./shared/Elements/DrainageFunctionWaterfraction.cpp \
	./shared/Elements/EstarComponents.cpp \
	./shared/Random/random.cpp \
	./shared/Random/randomgenerator.cpp \
	./shared/String/DescriptorIndex.cpp \
	./toolkits/issm/IssmToolkitUtils.cpp \
	./toolkits/issm/IssmSolver.cpp \
	./toolkits/mpi/issmmpi.cpp \
	./toolkits/mpi/commops/DetermineLocalSize.cpp \
	./toolkits/mpi/commops/DetermineGlobalSize.cpp \
	./toolkits/mpi/commops/DetermineRowRankFromLocalSize.cpp \
	./toolkits/mpi/commops/GetOwnershipBoundariesFromRange.cpp \
	./toolkits/ToolkitOptions.cpp \
	./modules/MmeToInputFromIdx/MmeToInputFromIdx.cpp\
	./modules/MmeToInputx/MmeToInputx.cpp\
	./modules/ModelProcessorx/ModelProcessorx.cpp \
	./modules/ModelProcessorx/ElementsAndVerticesPartitioning.cpp \
	./modules/ModelProcessorx/EdgesPartitioning.cpp \
	./modules/ModelProcessorx/FacesPartitioning.cpp \
	./modules/ModelProcessorx/CreateParameters.cpp \
	./modules/ModelProcessorx/Autodiff/CreateParametersAutodiff.cpp \
	./modules/ModelProcessorx/CreateFaces.cpp \
	./modules/ModelProcessorx/CreateEdges.cpp \
	./modules/ModelProcessorx/CreateSingleNodeToElementConnectivity.cpp \
	./modules/ModelProcessorx/CreateNumberNodeToElementConnectivity.cpp \
	./modules/ModelProcessorx/CreateElementsVerticesAndMaterials.cpp \
	./modules/ModelProcessorx/CreateNodes.cpp \
	./modules/ParseToolkitsOptionsx/ParseToolkitsOptionsx.cpp \
	./modules/NodesDofx/NodesDofx.cpp \
	./modules/NodalValuex/NodalValuex.cpp \
	./modules/VertexCoordinatesx/VertexCoordinatesx.cpp \
	./modules/ElementCoordinatesx/ElementCoordinatesx.cpp \
	./modules/OutputResultsx/OutputResultsx.cpp \
	./modules/InputDepthAverageAtBasex/InputDepthAverageAtBasex.cpp \
	./modules/InputDuplicatex/InputDuplicatex.cpp \
	./modules/InputExtrudex/InputExtrudex.cpp \
	./modules/SurfaceAreax/SurfaceAreax.cpp \
	./modules/AllocateSystemMatricesx/AllocateSystemMatricesx.cpp \
	./modules/CreateJacobianMatrixx/CreateJacobianMatrixx.cpp \
	./modules/SystemMatricesx/SystemMatricesx.cpp \
	./modules/CreateNodalConstraintsx/CreateNodalConstraintsx.cpp \
	./modules/UpdateDynamicConstraintsx/UpdateDynamicConstraintsx.cpp \
	./modules/IoModelToConstraintsx/IoModelToConstraintsx.cpp \
	./modules/SetActiveNodesLSMx/SetActiveNodesLSMx.cpp \
	./modules/InputUpdateFromConstantx/InputUpdateFromConstantx.cpp \
	./modules/InputUpdateFromSolutionx/InputUpdateFromSolutionx.cpp \
	./modules/GeothermalFluxx/GeothermalFluxx.cpp \
	./modules/GetSolutionFromInputsx/GetSolutionFromInputsx.cpp \
	./modules/GetVectorFromInputsx/GetVectorFromInputsx.cpp \
	./modules/InputUpdateFromVectorx/InputUpdateFromVectorx.cpp \
	./modules/FloatingiceMeltingRatex/FloatingiceMeltingRatex.cpp \
	./modules/FloatingiceMeltingRatePicox/FloatingiceMeltingRatePicox.cpp \
	./modules/FrontalForcingsx/FrontalForcingsx.cpp \
	./modules/ConfigureObjectsx/ConfigureObjectsx.cpp \
	./modules/SpcNodesx/SpcNodesx.cpp \
	./modules/SurfaceMassBalancex/SurfaceMassBalancex.cpp \
	./modules/SurfaceMassBalancex/Gembx.cpp \
	./modules/Reducevectorgtofx/Reducevectorgtofx.cpp \
	./modules/Reduceloadx/Reduceloadx.cpp \
	./modules/ConstraintsStatex/ConstraintsStatex.cpp \
	./modules/ResetConstraintsx/ResetConstraintsx.cpp \
	./modules/ResetFSBasalBoundaryConditionx/ResetFSBasalBoundaryConditionx.cpp \
	./modules/Solverx/Solverx.cpp \
	./modules/StochasticForcingx/StochasticForcingx.cpp \
	./modules/Mergesolutionfromftogx/Mergesolutionfromftogx.cpp \
	./modules/UpdateMmesx/UpdateMmesx.cpp \
	./cores/ProcessArguments.cpp \
	./cores/ResetBoundaryConditions.cpp \
	./cores/WrapperCorePointerFromSolutionEnum.cpp \
	./cores/WrapperPreCorePointerFromSolutionEnum.cpp \
	./cores/CorePointerFromSolutionEnum.cpp \
	./cores/ad_core.cpp \
	./main/EnvironmentInit.cpp \
	./main/EnvironmentFinalize.cpp \
	./analyses/EnumToAnalysis.cpp \
	./solutionsequences/solutionsequence_la.cpp \
	./solutionsequences/solutionsequence_la_theta.cpp \
	./solutionsequences/solutionsequence_linear.cpp \
	./solutionsequences/solutionsequence_nonlinear.cpp \
	./solutionsequences/solutionsequence_newton.cpp \
	./solutionsequences/solutionsequence_fct.cpp \
	./solutionsequences/solutionsequence_schurcg.cpp \
	./solutionsequences/solutionsequence_sampling.cpp \
	./solutionsequences/convergence.cpp \
	./classes/Options/Options.cpp \
	./classes/RiftStruct.cpp \
	./modules/ModelProcessorx/Transient/UpdateElementsTransient.cpp \
	./modules/ModelProcessorx/Transient/UpdateParametersTransient.cpp \
	./cores/transient_core.cpp \
	./cores/steadystate_core.cpp \
	./cores/masstransport_core.cpp \
	./cores/mmemasstransport_core.cpp \
	./cores/oceantransport_core.cpp \
	./cores/depthaverage_core.cpp \
	./cores/extrudefrombase_core.cpp \
	./cores/extrudefromtop_core.cpp \
	./cores/thermal_core.cpp \
	./cores/smb_core.cpp \
	./cores/bmb_core.cpp \
	./cores/debris_core.cpp \
	./solutionsequences/solutionsequence_thermal_nonlinear.cpp \
	./modules/ControlInputSetGradientx/ControlInputSetGradientx.cpp \
	./modules/GetVectorFromControlInputsx/GetVectorFromControlInputsx.cpp \
	./modules/SetControlInputsFromVectorx/SetControlInputsFromVectorx.cpp \
	./modules/ModelProcessorx/Control/CreateParametersControl.cpp \
	./modules/ModelProcessorx/Control/UpdateElementsAndMaterialsControl.cpp \
	./modules/SurfaceAbsVelMisfitx/SurfaceAbsVelMisfitx.cpp \
	./modules/SurfaceRelVelMisfitx/SurfaceRelVelMisfitx.cpp \
	./modules/SurfaceLogVelMisfitx/SurfaceLogVelMisfitx.cpp \
	./modules/SurfaceLogVxVyMisfitx/SurfaceLogVxVyMisfitx.cpp \
	./modules/SurfaceAverageVelMisfitx/SurfaceAverageVelMisfitx.cpp \
	./modules/ThicknessAbsMisfitx/ThicknessAbsMisfitx.cpp \
	./modules/Gradjx/Gradjx.cpp \
	./modules/DragCoefficientAbsGradientx/DragCoefficientAbsGradientx.cpp \
	./modules/ThicknessAlongGradientx/ThicknessAlongGradientx.cpp \
	./modules/ThicknessAcrossGradientx/ThicknessAcrossGradientx.cpp \
	./modules/RheologyBbarAbsGradientx/RheologyBbarAbsGradientx.cpp \
	./modules/RheologyBAbsGradientx/RheologyBAbsGradientx.cpp \
	./shared/Numerics/BrentSearch.cpp \
	./cores/control_core.cpp \
	./cores/controltao_core.cpp \
	./cores/controlm1qn3_core.cpp \
	./cores/controladm1qn3_core.cpp \
	./cores/controlvalidation_core.cpp \
	./cores/adjointstressbalance_core.cpp \
	./cores/adjointbalancethickness_core.cpp \
	./cores/adjointbalancethickness2_core.cpp \
	./cores/AdjointCorePointerFromSolutionEnum.cpp \
	./solutionsequences/solutionsequence_adjoint_linear.cpp \
	./cores/hydrology_core.cpp \
	./solutionsequences/solutionsequence_hydro_nonlinear.cpp \
	./solutionsequences/solutionsequence_shakti_nonlinear.cpp \
	./solutionsequences/solutionsequence_glads_nonlinear.cpp \
	./cores/stressbalance_core.cpp \
	./solutionsequences/solutionsequence_stokescoupling_nonlinear.cpp \
	./cores/balancethickness_core.cpp \
	./cores/balancethickness2_core.cpp \
	./cores/balancevelocity_core.cpp \
	./cores/dummy_core.cpp \
	./cores/surfaceslope_core.cpp \
	./cores/bedslope_core.cpp \
	./cores/damage_core.cpp \
	./cores/levelsetfunctionslope_core.cpp \
	./cores/movingfront_core.cpp \
	./cores/groundingline_core.cpp \
	./modules/GroundinglineMigrationx/GroundinglineMigrationx.cpp \
	./classes/Loads/Riftfront.cpp \
	./modules/ConstraintsStatex/RiftConstraintsState.cpp \
	./modules/ModelProcessorx/CreateOutputDefinitions.cpp \
	./modules/OutputDefinitionsResponsex/OutputDefinitionsResponsex.cpp \
	./modules/InterpFromMeshToMesh2dx/InterpFromMeshToMesh2dx.cpp \
	./classes/Inputs/Inputs.cpp \
	./classes/Inputs/BoolInput.cpp \
	./classes/Inputs/DoubleInput.cpp \
	./classes/Inputs/IntInput.cpp \
	./classes/Inputs/ElementInput.cpp \
	./classes/Inputs/SegInput.cpp \
	./classes/Inputs/TriaInput.cpp \
	./classes/Inputs/PentaInput.cpp \
	./classes/Inputs/DatasetInput.cpp \
	./classes/Inputs/ControlInput.cpp \
	./classes/Inputs/TransientInput.cpp \
	./classes/Inputs/TransientFileInput.cpp \
	./classes/Inputs/ArrayInput.cpp \
	./classes/Inputs/IntArrayInput.cpp
#}}}
# CoDiPack/MeDiPack sources {{{
if MEDIPACK
issm_sources += ./toolkits/codipack/ampi_interface.cpp
endif
if CODIPACK
issm_sources += \
	./toolkits/codipack/CoDiPackGlobal.cpp \
	./toolkits/codipack/CoDiPackDebug.cpp
endif
#}}}
# DAKOTA sources {{{
if DAKOTA
issm_sources += \
	./classes/Dakota/IssmDirectApplicInterface.h \
	./classes/Dakota/IssmParallelDirectApplicInterface.cpp \
	./modules/InputUpdateFromDakotax/InputUpdateFromDakotax.cpp \
	./modules/InputUpdateFromVectorDakotax/InputUpdateFromVectorDakotax.cpp \
	./modules/InputUpdateFromMatrixDakotax/InputUpdateFromMatrixDakotax.cpp \
	./modules/AverageOntoPartitionx/AverageOntoPartitionx.cpp \
	./modules/ModelProcessorx/Dakota/CreateParametersDakota.cpp \
	./modules/ModelProcessorx/Dakota/UpdateElementsAndMaterialsDakota.cpp \
	./cores/dakota_core.cpp
if SYSTEM_HAS_FMEMOPEN
issm_sources += ./modules/QmuStatisticsx/QmuStatisticsx.cpp
endif
endif
#}}}
# PETSc sources {{{
if PETSC
issm_sources += \
	./toolkits/petsc \
	./toolkits/petsc/patches \
	./toolkits/petsc/patches/VecToMPISerial.cpp \
	./toolkits/petsc/patches/MatToMPISerial.cpp \
	./toolkits/petsc/patches/NewVec.cpp \
	./toolkits/petsc/patches/PetscOptionsDetermineSolverType.cpp \
	./toolkits/petsc/patches/NewMat.cpp \
	./toolkits/petsc/patches/VecFree.cpp \
	./toolkits/petsc/patches/KSPFree.cpp \
	./toolkits/petsc/patches/MatFree.cpp \
	./toolkits/petsc/patches/ISSMToPetscMatrixType.cpp \
	./toolkits/petsc/patches/ISSMToPetscInsertMode.cpp \
	./toolkits/petsc/patches/ISSMToPetscNormMode.cpp \
	./toolkits/petsc/objects/PetscMat.cpp \
	./toolkits/petsc/objects/PetscVec.cpp \
	./toolkits/petsc/objects/PetscSolver.cpp
endif
#}}}
# MUMPS sources {{{
if MUMPS
issm_sources += ./toolkits/mumps/MumpsSolve.cpp
endif
#}}}
# GSL sources {{{
if GSL
issm_sources += ./toolkits/gsl/DenseGslSolve.cpp
endif
#}}}
# PROJ sources {{{
if PROJ
issm_sources += ./modules/CoordinateSystemTransformx/CoordinateSystemTransformx.cpp
endif
#}}}
# Analyses {{{
if ADJOINTBALANCETHICKNESS
issm_sources += ./analyses/AdjointBalancethicknessAnalysis.cpp
endif
if ADJOINTBALANCETHICKNESS2
issm_sources += ./analyses/AdjointBalancethickness2Analysis.cpp
endif
if ADJOINTHORIZ
issm_sources += ./analyses/AdjointHorizAnalysis.cpp
endif
if AGE
issm_sources += ./analyses/AgeAnalysis.cpp
endif
if BALANCETHICKNESS
issm_sources += ./analyses/BalancethicknessAnalysis.cpp
endif
if BALANCETHICKNESS2
issm_sources += ./analyses/Balancethickness2Analysis.cpp
endif
if BALANCETHICKNESSSOFT
issm_sources += ./analyses/BalancethicknessSoftAnalysis.cpp
endif
if BALANCEVELOCITY
issm_sources += ./analyses/BalancevelocityAnalysis.cpp
endif
if L2PROJECTIONBASE
issm_sources += ./analyses/L2ProjectionBaseAnalysis.cpp
endif
if DAMAGEEVOLUTION
issm_sources += ./analyses/DamageEvolutionAnalysis.cpp
issm_sources += ./modules/Damagex/Damagex.cpp
endif
if DEBRIS
issm_sources += ./analyses/DebrisAnalysis.cpp
endif
if STRESSBALANCE
issm_sources += ./analyses/StressbalanceAnalysis.cpp
endif
if UZAWAPRESSURE
issm_sources += ./analyses/UzawaPressureAnalysis.cpp
endif
if STRESSBALANCESIA
issm_sources += ./analyses/StressbalanceSIAAnalysis.cpp
endif
if STRESSBALANCEVERTICAL
issm_sources += ./analyses/StressbalanceVerticalAnalysis.cpp
endif
if ENTHALPY
issm_sources += ./analyses/EnthalpyAnalysis.cpp
endif
if GLHEIGHTADVECTION
issm_sources += ./analyses/GLheightadvectionAnalysis.cpp
endif
if HYDROLOGYSHREVE
issm_sources += ./analyses/HydrologyShreveAnalysis.cpp
endif
if HYDROLOGYTWS
issm_sources += ./analyses/HydrologyTwsAnalysis.cpp
endif

if HYDROLOGYSHAKTI
issm_sources += ./analyses/HydrologyShaktiAnalysis.cpp
endif
if HYDROLOGYPISM
issm_sources += ./analyses/HydrologyPismAnalysis.cpp
endif
if HYDROLOGYGLADS
issm_sources += ./analyses/HydrologyGlaDSAnalysis.cpp
endif
if HYDROLOGYDCINEFFICIENT
issm_sources += ./analyses/HydrologyDCInefficientAnalysis.cpp
endif
if HYDROLOGYDCEFFICIENT
issm_sources += ./analyses/HydrologyDCEfficientAnalysis.cpp
endif
if HYDROLOGYARMAPW
issm_sources += ./analyses/HydrologyArmapwAnalysis.cpp
endif
if L2PROJECTIONEPL
issm_sources += ./analyses/L2ProjectionEPLAnalysis.cpp
endif
if MELTING
issm_sources += ./analyses/MeltingAnalysis.cpp
endif
if MASSTRANSPORT
issm_sources += ./analyses/MasstransportAnalysis.cpp \
				./analyses/MmemasstransportAnalysis.cpp
endif
if OCEANTRANSPORT
issm_sources += ./analyses/OceantransportAnalysis.cpp
endif
if SMB
issm_sources += ./analyses/SmbAnalysis.cpp
endif

if FREESURFACEBASE
issm_sources += ./analyses/FreeSurfaceBaseAnalysis.cpp
endif
if FREESURFACETOP
issm_sources += ./analyses/FreeSurfaceTopAnalysis.cpp
endif
if EXTRUDEFROMBASE
issm_sources += ./analyses/ExtrudeFromBaseAnalysis.cpp
endif
if EXTRUDEFROMTOP
issm_sources += ./analyses/ExtrudeFromTopAnalysis.cpp
endif
if DEPTHAVERAGE
issm_sources += ./analyses/DepthAverageAnalysis.cpp
endif
if THERMAL
issm_sources += ./analyses/ThermalAnalysis.cpp
endif
if SMOOTH
issm_sources += ./analyses/SmoothAnalysis.cpp
endif
if LEVELSET
issm_sources += ./analyses/LevelsetAnalysis.cpp
issm_sources += ./modules/Calvingx/Calvingx.cpp
issm_sources += ./modules/KillIcebergsx/KillIcebergsx.cpp
endif
if EXTRAPOLATION
issm_sources += ./analyses/ExtrapolationAnalysis.cpp
endif
#}}}
# Love sources (only if have Fortran) {{{
if LOVE
if FORTRAN
issm_sources += \
	./cores/love_core.cpp \
	./analyses/LoveAnalysis.cpp \
	./modules/Zgesvx/Zgesvx.cpp \
	./modules/Zgesvx/Zgesvx.h \
	./modules/Zgesvx/dcabs1.f \
	./modules/Zgesvx/dlamch.f \
	./modules/Zgesvx/ieeeck.f \
	./modules/Zgesvx/ilaenv.f \
	./modules/Zgesvx/iparmq.f \
	./modules/Zgesvx/izamax.f \
	./modules/Zgesvx/lsame.f \
	./modules/Zgesvx/xerbla.f \
	./modules/Zgesvx/zgemm.f \
	./modules/Zgesvx/zgeru.f \
	./modules/Zgesvx/zgesv.f \
	./modules/Zgesvx/zgetf2.f \
	./modules/Zgesvx/zgetrf2.f \
	./modules/Zgesvx/zgetrf.f \
	./modules/Zgesvx/zgetrs.f \
	./modules/Zgesvx/zlaswp.f \
	./modules/Zgesvx/zscal.f \
	./modules/Zgesvx/zswap.f \
	./modules/Zgesvx/ztrsm.f
endif
endif
#}}}
# Esa sources {{{
if ESA
issm_sources += \
	./cores/esa_core.cpp \
	./analyses/EsaAnalysis.cpp
endif
#}}}
# Ocean sources {{{
if OCEAN
issm_sources += ./modules/OceanExchangeDatax/OceanExchangeDatax.cpp
endif
#}}}
# Sampling sources {{{
if SAMPLING
issm_sources += \
	./cores/sampling_core.cpp \
	./analyses/SamplingAnalysis.cpp
endif
#}}}
# SLC sources {{{
if SEALEVELCHANGE
issm_sources += \
	./cores/sealevelchange_core.cpp \
	./analyses/SealevelchangeAnalysis.cpp\
	./classes/GrdLoads.cpp\
	./classes/SealevelGeometry.cpp

# GIA Ivins physics (only if have Fortran)
if FORTRAN
issm_sources += \
	./modules/GiaDeflectionCorex/GiaDeflectionCorex.cpp \
	./modules/GiaDeflectionCorex/distme.f \
	./modules/GiaDeflectionCorex/freed.f \
	./modules/GiaDeflectionCorex/ojrule.f \
	./modules/GiaDeflectionCorex/pwise.f \
	./modules/GiaDeflectionCorex/qwise.f \
	./modules/GiaDeflectionCorex/stot.f \
	./modules/GiaDeflectionCorex/what0.f
endif

endif
#}}}
# METIS sources {{{
if METIS
issm_sources += ./modules/MeshPartitionx/MeshPartitionx.cpp
issm_sources += ./toolkits/metis/patches/METIS_PartMeshNodalPatch.cpp
endif
#}}}
# Esmf sources {{{
if ESMF
issm_sources += ./main/esmfbinders.cpp
endif
#}}}
# SEMIC sources {{{
if SEMIC
if FORTRAN
issm_sources += ./modules/SurfaceMassBalancex/run_semic.f90
issm_sources += ./modules/SurfaceMassBalancex/run_semic_transient.f90
endif
endif
#}}}
# Wrapper sources
# NEOPZ sources {{{
neopz_sources = ./classes/AmrNeopz.cpp
#}}}
# Modules sources {{{
modules_sources = \
	./shared/Threads/LaunchThread.cpp \
	./shared/Threads/PartitionRange.cpp \
	./shared/Exp/exp.cpp \
	./shared/Triangle/AssociateSegmentToElement.cpp \
	./shared/Triangle/GridInsideHole.cpp \
	./shared/Triangle/OrderSegments.cpp \
	./shared/Triangle/SplitMeshForRifts.cpp \
	./shared/Triangle/TriangleUtils.cpp \
	./modules/Trianglex/Trianglex.cpp \
	./modules/ProcessRiftsx/ProcessRiftsx.cpp \
	./modules/PointCloudFindNeighborsx/PointCloudFindNeighborsx.cpp \
	./modules/PointCloudFindNeighborsx/PointCloudFindNeighborsxt.cpp \
	./modules/InterpFromGridToMeshx/InterpFromGridToMeshx.cpp \
	./modules/InterpFromMesh2dx/InterpFromMesh2dx.cpp \
	./modules/InterpFromMesh2dx/InterpFromMesh2dxt.cpp \
	./modules/InterpFromMeshToMesh3dx/InterpFromMeshToMesh3dx.cpp \
	./modules/InterpFromMeshToGridx/InterpFromMeshToGridx.cpp \
	./modules/MeshProfileIntersectionx/MeshProfileIntersectionx.cpp \
	./modules/ContourToMeshx/ContourToMeshx.cpp \
	./modules/ContourToMeshx/ContourToMeshxt.cpp \
	./modules/ExpToLevelSetx/ExpToLevelSetx.cpp \
	./modules/ExpToLevelSetx/ExpToLevelSetxt.cpp \
	./modules/ContourToNodesx/ContourToNodesx.cpp \
	./modules/DistanceToMaskBoundaryx/DistanceToMaskBoundaryx.cpp \
	./modules/DistanceToMaskBoundaryx/DistanceToMaskBoundaryxt.cpp \
	./modules/NodeConnectivityx/NodeConnectivityx.cpp \
	./modules/ElementConnectivityx/ElementConnectivityx.cpp \
	./modules/PropagateFlagsFromConnectivityx/PropagateFlagsFromConnectivityx.cpp

if CHACO
modules_sources += \
	./modules/Chacox/Chacox.cpp \
	./modules/Chacox/input_parse.cpp \
	./modules/Chacox/chaco_seconds.cpp \
	./modules/Chacox/user_params.cpp
endif
if SCOTCH
modules_sources += ./modules/Scotchx/Scotchx.cpp
endif
#}}}

# Kriging (wrapper and executable)
# Kriging sources {{{
if KRIGING
issm_sources += \
	./classes/kriging/Observations.cpp \
	./classes/kriging/GaussianVariogram.cpp \
	./classes/kriging/ExponentialVariogram.cpp \
	./classes/kriging/SphericalVariogram.cpp \
	./classes/kriging/PowerVariogram.cpp \
	./classes/kriging/Quadtree.cpp \
	./classes/kriging/Covertree.cpp \
	./classes/kriging/Observation.cpp \
	./modules/Krigingx/pKrigingx.cpp

modules_sources += \
	./modules/Krigingx/Krigingx.cpp \
	./modules/Krigingx/pKrigingx.cpp
endif
#}}}
# Library flags and sources {{{
ALL_CXXFLAGS = $(CXXFLAGS)

if MSYS2
ALL_CXXFLAGS += -D_DO_NOT_LOAD_GLOBALS_ # Cannot have undefined symbols under MSYS2
ALL_CXXFLAGS += -D_USE_MATH_DEFINES # Need to open fence to math.h definitions when --std=c++ is used (which is the default)
endif

libISSMCore_la_SOURCES  = $(issm_sources)
libISSMCore_la_CXXFLAGS = $(ALL_CXXFLAGS)
libISSMCore_la_FFLAGS = $(AM_FFLAGS)
libISSMCore_LIB_ADD = $(CHACOLIB) $(DAKOTALIB) $(PETSCLIB) $(MUMPSLIB) $(SCALAPACKLIB) $(BLASLAPACKLIB) $(PARMETISLIB) $(METISLIB) $(HDF5LIB) $(TAOLIB) $(M1QN3LIB) $(SEMICLIB) $(PLAPACKLIB) $(MPLAPACKLIB) $(SUPERLULIB) $(SPOOLESLIB) $(BLACSLIB) $(HYPRELIB) $(SPAILIB) $(PROMETHEUSLIB) $(PASTIXLIB) $(MLLIB) $(SCOTCHLIB) $(MKLLIB) $(MPILIB) $(MATHLIB) $(GRAPHICSLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(ADOLCLIB) $(AMPILIB) $(PROJLIB) $(ADJOINTPETSCLIB)
if FORTRAN
libISSMCore_LIB_ADD += $(FLIBS) $(FORTRANLIB)
endif
libISSMCore_LIB_ADD += $(OSLIBS)

if MSYS2
libISSMCore_la_LIBADD = ${libISSMCore_LIB_ADD}
else
if !STANDALONE_LIBRARIES
libISSMCore_la_LIBADD = ${libISSMCore_LIB_ADD}
endif
endif

if NEOPZ
libISSMCore_la_SOURCES += $(neopz_sources)
endif

if WRAPPERS
libISSMModules_la_SOURCES = $(modules_sources)
libISSMModules_la_SOURCES += $(bamg_sources)
if KRIGING
libISSMModules_la_SOURCES += $(kriging_sources)
endif

libISSMModules_la_CXXFLAGS = $(ALL_CXXFLAGS)

libISSMModules_LIB_ADD = $(TRIANGLELIB) $(CHACOLIB) $(BLASLAPACKLIB) $(OSLIBS)
libISSMModules_la_LIBADD = ./libISSMCore.la

if MSYS2
libISSMModules_la_LIBADD += ${libISSMModules_LIB_ADD}
else
if !STANDALONE_LIBRARIES
libISSMModules_la_LIBADD += ${libISSMModules_LIB_ADD}
endif
endif
endif

AM_LDFLAGS =

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
if STANDALONE_EXECUTABLES
if MAC
AM_LDFLAGS += -Wl,-rpath,'@loader_path/../lib'
else
if MSYS2
AM_LDFLAGS += -Wl,-static
else
AM_LDFLAGS += -Wl,--disable-new-dtags -Wl,-rpath,'$$ORIGIN/../lib'
endif
endif
endif

if !STANDALONE_EXECUTABLES
if MAC
if HAVE_FORTRANDIR
AM_LDFLAGS += -Wl,-rpath,$(FORTRANDIR)
endif
if BOOST
AM_LDFLAGS += -Wl,-rpath,$(BOOSTROOT)/lib
endif
if PROJ
AM_LDFLAGS += -Wl,-rpath,$(PROJROOT)/lib
endif
endif
endif

if VERSION
libISSMCore_la_LDFLAGS =
if WRAPPERS
libISSMModules_la_LDFLAGS =
endif
else
libISSMCore_la_LDFLAGS = -avoid-version
if WRAPPERS
libISSMModules_la_LDFLAGS = -avoid-version
endif
endif

if STANDALONE_LIBRARIES
if !MSYS2
libISSMCore_la_LDFLAGS += -static
if WRAPPERS
libISSMModules_la_LDFLAGS += -static
endif
endif
endif

if !STANDALONE_LIBRARIES
if MAC
if BOOST
libISSMCore_la_LDFLAGS += -Wl,-rpath,$(BOOSTROOT)/lib
endif
if PROJ
libISSMCore_la_LDFLAGS += -Wl,-rpath,$(PROJROOT)/lib
endif
endif
endif

if !MSYS2
libISSMOverload_la_SOURCES = ./shared/String/ApiPrintf.cpp
libISSMOverload_la_CFLAGS  = -D_C_ -fPIC $(CFLAGS) $(CXXFLAGS)

if VERSION
libISSMOverload_la_LDFLAGS =
else
libISSMOverload_la_LDFLAGS = -avoid-version
endif

if STANDALONE_LIBRARIES
if !MSYS2
libISSMOverload_la_LDFLAGS += -static
endif
endif
endif
#}}}
# Executables {{{
if JAVASCRIPT
bin_PROGRAMS =
else
bin_PROGRAMS = issm issm_slc
endif

# Standard libraries
LDADD = ./libISSMCore.la

if !MSYS2
LDADD += ./libISSMOverload.la
endif

if WRAPPERS
LDADD += ./libISSMModules.la $(TRIANGLELIB)
endif

# External packages
LDADD += $(DAKOTALIB) $(SEMICLIB) $(M1QN3LIB) $(CHACOLIB) $(PETSCLIB) $(MUMPSLIB) $(SCALAPACKLIB) $(BLASLAPACKLIB) $(PARMETISLIB) $(METISLIB) $(NEOPZLIB) $(TAOLIB) $(PLAPACKLIB) $(MPLAPACKLIB) $(SUPERLULIB) $(SPOOLESLIB) $(BLACSLIB) $(HDF5LIB) $(HYPRELIB) $(SPAILIB) $(PROMETHEUSLIB) $(PASTIXLIB) $(MLLIB) $(SCOTCHLIB) $(MKLLIB) $(MATHLIB) $(GRAPHICSLIB) $(MULTITHREADINGLIB) $(GSLLIB) $(AMPILIB) $(ADOLCLIB) $(MPILIB) $(ESMFLIB)

if FORTRAN
LDADD += $(FLIBS) $(FORTRANLIB)
endif

LDADD += $(OSLIBS)

issm_SOURCES = main/issm.cpp
issm_CXXFLAGS= $(CXXFLAGS)

issm_slc_SOURCES = main/issm_slc.cpp
issm_slc_CXXFLAGS= $(CXXFLAGS)

if OCEAN
bin_PROGRAMS += issm_ocean
issm_ocean_SOURCES = main/issm_ocean.cpp
issm_ocean_CXXFLAGS= $(CXXFLAGS)
endif

if KRIGING
bin_PROGRAMS += kriging
kriging_SOURCES = main/kriging.cpp
kriging_CXXFLAGS= $(CXXFLAGS)
endif

if ISSM_DAKOTA
bin_PROGRAMS += issm_dakota
issm_dakota_SOURCES = main/issm_dakota.cpp
issm_dakota_CXXFLAGS= $(CXXFLAGS)
bin_PROGRAMS += issm_post
issm_post_SOURCES = main/issm_post.cpp
issm_post_CXXFLAGS= $(CXXFLAGS)
endif
#}}}
