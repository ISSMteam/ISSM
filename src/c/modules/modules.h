/*! \file modules.h: 
 *  \brief header file for all the ISSM modules
 */

#ifndef _ISSM_MODULES_H_
#define _ISSM_MODULES_H_

#ifdef _HAVE_BAMG_
#include "./Bamgx/Bamgx.h"
#include "./BamgConvertMeshx/BamgConvertMeshx.h"
#include "./BamgTriangulatex/BamgTriangulatex.h"
#endif

/*Modules: */
#include "./AllocateSystemMatricesx/AllocateSystemMatricesx.h"
#include "./AverageOntoPartitionx/AverageOntoPartitionx.h"
#include "./Calvingx/Calvingx.h"
#include "./Chacox/Chacox.h"
#include "./ConfigureObjectsx/ConfigureObjectsx.h"
#include "./ContourToMeshx/ContourToMeshx.h"
#include "./ContourToNodesx/ContourToNodesx.h"
#include "./ControlInputSetGradientx/ControlInputSetGradientx.h"
#include "./CreateNodalConstraintsx/CreateNodalConstraintsx.h"
#include "./CreateJacobianMatrixx/CreateJacobianMatrixx.h"
#include "./Damagex/Damagex.h"
#include "./DragCoefficientAbsGradientx/DragCoefficientAbsGradientx.h"
#include "./DistanceToMaskBoundaryx/DistanceToMaskBoundaryx.h"
#include "./ExpToLevelSetx/ExpToLevelSetx.h"
#include "./ElementConnectivityx/ElementConnectivityx.h"
#include "./GeothermalFluxx/GeothermalFluxx.h"
#include "./GetSolutionFromInputsx/GetSolutionFromInputsx.h"
#include "./GetVectorFromInputsx/GetVectorFromInputsx.h"
#include "./GetVectorFromControlInputsx/GetVectorFromControlInputsx.h"
#include "./GiaDeflectionCorex/GiaDeflectionCorex.h"
#include "./SetControlInputsFromVectorx/SetControlInputsFromVectorx.h"
#include "./SetActiveNodesLSMx/SetActiveNodesLSMx.h"
#include "./Gradjx/Gradjx.h"
#include "./GroundinglineMigrationx/GroundinglineMigrationx.h"
#include "./FrontalForcingsx/FrontalForcingsx.h"
#include "./InputDepthAverageAtBasex/InputDepthAverageAtBasex.h"
#include "./InputDuplicatex/InputDuplicatex.h"
#include "./InputExtrudex/InputExtrudex.h"
#include "./InterpFromMesh2dx/InterpFromMesh2dx.h"
#include "./InterpFromGridToMeshx/InterpFromGridToMeshx.h"
#include "./InterpFromMeshToMesh2dx/InterpFromMeshToMesh2dx.h"
#include "./InterpFromMeshToMesh3dx/InterpFromMeshToMesh3dx.h"
#include "./InterpFromMeshToGridx/InterpFromMeshToGridx.h"
#include "./InputUpdateFromConstantx/InputUpdateFromConstantx.h"
#include "./InputUpdateFromSolutionx/InputUpdateFromSolutionx.h"
#include "./InputUpdateFromDakotax/InputUpdateFromDakotax.h"
#include "./InputUpdateFromVectorx/InputUpdateFromVectorx.h"
#include "./InputUpdateFromVectorDakotax/InputUpdateFromVectorDakotax.h"
#include "./InputUpdateFromMatrixDakotax/InputUpdateFromMatrixDakotax.h"
#include "./IoModelToConstraintsx/IoModelToConstraintsx.h"
#include "./KillIcebergsx/KillIcebergsx.h"
#include "./Krigingx/Krigingx.h"
#include "./FloatingiceMeltingRatex/FloatingiceMeltingRatex.h"
#include "./FloatingiceMeltingRatePicox/FloatingiceMeltingRatePicox.h"
#include "./Mergesolutionfromftogx/Mergesolutionfromftogx.h"
#include "./MeshPartitionx/MeshPartitionx.h"
#include "./MeshProfileIntersectionx/MeshProfileIntersectionx.h"
#include "./SurfaceAbsVelMisfitx/SurfaceAbsVelMisfitx.h"
#include "./SurfaceRelVelMisfitx/SurfaceRelVelMisfitx.h"
#include "./SurfaceLogVelMisfitx/SurfaceLogVelMisfitx.h"
#include "./SurfaceLogVxVyMisfitx/SurfaceLogVxVyMisfitx.h"
#include "./SurfaceAverageVelMisfitx/SurfaceAverageVelMisfitx.h"
#include "./ModelProcessorx/ModelProcessorx.h"
#include "./MmeToInputFromIdx/MmeToInputFromIdx.h"
#include "./MmeToInputx/MmeToInputx.h"
#include "./ParseToolkitsOptionsx/ParseToolkitsOptionsx.h"
#include "./NodalValuex/NodalValuex.h"
#include "./NodeConnectivityx/NodeConnectivityx.h"
#include "./NodesDofx/NodesDofx.h"
#include "./OutputDefinitionsResponsex/OutputDefinitionsResponsex.h"
#include "./OutputResultsx/OutputResultsx.h"
#include "./ConstraintsStatex/ConstraintsStatex.h"
#include "./PointCloudFindNeighborsx/PointCloudFindNeighborsx.h"
#include "./PropagateFlagsFromConnectivityx/PropagateFlagsFromConnectivityx.h"
#include "./QmuStatisticsx/QmuStatisticsx.h"
#include "./Reduceloadx/Reduceloadx.h"
#include "./Reducevectorgtofx/Reducevectorgtofx.h"
#include "./ResetConstraintsx/ResetConstraintsx.h"
#include "./ResetFSBasalBoundaryConditionx/ResetFSBasalBoundaryConditionx.h"
#include "./RheologyBbarAbsGradientx/RheologyBbarAbsGradientx.h"
#include "./RheologyBAbsGradientx/RheologyBAbsGradientx.h"
#include "./Scotchx/Scotchx.h"
#include "./StochasticForcingx/StochasticForcingx.h"
#include "./SurfaceMassBalancex/SurfaceMassBalancex.h"
#include "./Solverx/Solverx.h"
#include "./SystemMatricesx/SystemMatricesx.h"
#include "./SpcNodesx/SpcNodesx.h"
#include "./SurfaceAreax/SurfaceAreax.h"
#include "./Trianglex/Trianglex.h"
#include "./ProcessRiftsx/ProcessRiftsx.h"
#include "./ThicknessAbsMisfitx/ThicknessAbsMisfitx.h"
#include "./ThicknessAlongGradientx/ThicknessAlongGradientx.h"
#include "./ThicknessAcrossGradientx/ThicknessAcrossGradientx.h"
#include "./UpdateDynamicConstraintsx/UpdateDynamicConstraintsx.h"
#include "./UpdateMmesx/UpdateMmesx.h"
#include "./VertexCoordinatesx/VertexCoordinatesx.h"
#include "./ElementCoordinatesx/ElementCoordinatesx.h"
#include "./Zgesvx/Zgesvx.h"

#ifdef _HAVE_OCEAN_
#include "./OceanExchangeDatax/OceanExchangeDatax.h"
#endif

#endif
