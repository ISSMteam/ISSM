/*
 * cores.h:
 */

#ifndef _CORES_H_
#define _CORES_H_

/*forward declarations: */
class FemModel;
class Parameters;
class SealevelGeometry;
template <class doubletype> class Matrix;
template <class doubletype> class Vector;

#include "../shared/io/Comm/IssmComm.h"
#include "../shared/Numerics/types.h"

/*cores: */
void adjointstressbalance_core(FemModel* femmodel);
void adjointbalancethickness_core(FemModel* femmodel);
void adjointbalancethickness2_core(FemModel* femmodel);
void stressbalance_core(FemModel* femmodel);
void hydrology_core(FemModel* femmodel);
void thermal_core(FemModel* femmodel);
void surfaceslope_core(FemModel* femmodel);
void levelsetfunctionslope_core(FemModel* femmodel);
void movingfront_core(FemModel* femmodel);
void groundingline_core(FemModel* femmodel);
void bedslope_core(FemModel* femmodel);
void meshdeformation_core(FemModel* femmodel);
void control_core(FemModel* femmodel);
void controltao_core(FemModel* femmodel);
void controlm1qn3_core(FemModel* femmodel);
void controladm1qn3_core(FemModel* femmodel);
void controlvalidation_core(FemModel* femmodel);
void masstransport_core(FemModel* femmodel);
void oceantransport_core(FemModel* femmodel);
void depthaverage_core(FemModel* femmodel);
void extrudefrombase_core(FemModel* femmodel);
void extrudefromtop_core(FemModel* femmodel);
void balancethickness_core(FemModel* femmodel);
void balancethickness2_core(FemModel* femmodel);
void balancevelocity_core(FemModel* femmodel);
void slopecompute_core(FemModel* femmodel);
void steadystate_core(FemModel* femmodel);
void transient_core(FemModel* femmodel);
void transient_precore(FemModel* femmodel);
void dakota_core(FemModel* femmodel);
void ad_core(FemModel* femmodel);
void adgradient_core(FemModel* femmodel);
void dummy_core(FemModel* femmodel);
void gia_core(FemModel* femmodel);
void love_core(FemModel* femmodel);
void esa_core(FemModel* femmodel);
void smb_core(FemModel* femmodel);
void bmb_core(FemModel* femmodel);
void damage_core(FemModel* femmodel);
void sampling_core(FemModel* femmodel);
void debris_core(FemModel* femmodel);

/*sealevel change cores:*/
#ifdef _HAVE_SEALEVELCHANGE_
void sealevelchange_core(FemModel* femmodel);
void sealevelchange_initialgeometry(FemModel* femmodel);
SealevelGeometry* sealevelchange_geometry(FemModel* femmodel);
#endif
void grd_core(FemModel* femmodel,SealevelGeometry* slgeom);
void solidearthexternal_core(FemModel* femmodel);
void dynstr_core(FemModel* femmodel);
void couplerinput_core(FemModel* femmodel);
void coupleroutput_core(FemModel* femmodel);

//optimization
int GradJSearch(IssmDouble* search_vector,FemModel* femmodel,int step);
IssmDouble objectivefunction(IssmDouble search_scalar,FemModel* femmodel);

//diverse
void ProcessArguments(int* solution,char** pbinname,char** poutbinname,char** ptoolkitsname,char** plockname,char** prestartname, char** prootpath,char** pmodelname, int argc,char **argv);
void WriteLockFile(char* filename);
void ResetBoundaryConditions(FemModel* femmodel, int analysis_type);
void PrintBanner(void);
void EarthMassTransport(FemModel* femmodel);

//solution configuration
void CorePointerFromSolutionEnum(void (**psolutioncore)(FemModel*),Parameters* parameters,int solutiontype);
void WrapperCorePointerFromSolutionEnum(void (**psolutioncore)(FemModel*),Parameters* parameters,int solutiontype,bool nodakotacore=false);
void WrapperPreCorePointerFromSolutionEnum(void (**psolutioncore)(FemModel*),Parameters* parameters,int solutiontype);
void AdjointCorePointerFromSolutionEnum(void (**padjointcore)(FemModel*),int solutiontype);

#endif
