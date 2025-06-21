/*
 * FemModel.h:
 */

#ifndef _FEMMODEL_H_
#define _FEMMODEL_H_

/*Headers:*/
/*{{{*/
#include "../toolkits/toolkits.h"
class DataSet;
class Parameters;
class Inputs;
class Nodes;
class Vertices;
class Results;
class Constraints;
class Loads;
class Materials;
class SealevelMasks;
class Profiler;
class Elements;
#if defined(_HAVE_NEOPZ_) && !defined(_HAVE_AD_)
#include "./AmrNeopz.h"
#endif
#if defined(_HAVE_BAMG_) && !defined(_HAVE_AD_)
#include "./AmrBamg.h"
#endif
/*}}}*/

class FemModel {

	/*no private members, as we need access to these datasets quite often!:*/

	public:

		int          analysis_counter;     //counter into analysis_type_list
		int         *analysis_type_list;   //list of analyses this femmodel is going to carry out
		int          nummodels;
		int          solution_type;

		Profiler*    profiler;             //keep time, cpu and mem statistics while we are running.

		Elements    *elements;             //elements (one set for all analyses)
		Materials   *materials;            //one set of materials, for each element
		Parameters  *parameters;           //one set of parameters, independent of the analysis_type
		Inputs     *inputs;              //one set of inputs, independent of the analysis_type
		Results     *results;              //results that cannot be fit into the elements
		Vertices    *vertices;             //one set of vertices

		/*Analysis dependent datasets*/
		Constraints  *constraints;
		Constraints **constraints_list;
		Loads        *loads;
		Loads       **loads_list;
		Nodes        *nodes;
		Nodes       **nodes_list;

		//FIXME: do we want only one class and have virtual functions? or keep 2 classes?
		#if defined(_HAVE_NEOPZ_) && !defined(_HAVE_AD_)
		AmrNeopz *amr;		  //adaptive mesh refinement object. It keeps coarse mesh and execute refinement process
		#endif

		#if defined(_HAVE_BAMG_) && !defined(_HAVE_AD_)
		AmrBamg *amrbamg; //adaptive mesh refinement object. It keeps coarse mesh and execute refinement process
		#endif

		/*constructors, destructors: */
		FemModel(void);
		FemModel(int argc,char** argv,ISSM_MPI_Comm comm_init,bool trace=false);
		FemModel(char* rootpath, char* inputfilename, char* outputfilename, char* toolkitsfilename, char* lockfilename, char* restartfilename, char* modelname, ISSM_MPI_Comm incomm, int solution_type,IssmPDouble* X);
		~FemModel();

		/*Methods:*/
		int  AnalysisIndex(int);
		void CheckPoint(void);
		void CheckPointAD(int step);
		void CleanUp(void);
		FemModel* copy();
		void Echo();
		int  GetElementsWidth(){return 3;};//just tria elements in this first version
		void InitFromFiles(char* rootpath, char* inputfilename, char* outputfilename, char* petscfilename, char* lockfilename, char* restartfilename, char* modelname, const int solution_type,bool trace,IssmPDouble* X=NULL);
		void InitFromFids(char* rootpath, FILE* IOMODEL, FILE* toolkitsoptionsfid, int in_solution_type, bool trace, IssmPDouble* X=NULL);
		void Marshall(MarshallHandle* marshallhandle);
		void Restart(int verboselevel=2);
		void RestartAD(int step);
		void SetCurrentConfiguration(int configuration_type);
		void SetCurrentConfiguration(int configuration_type,int analysis_type);
		int  Size(void);
		void SolutionAnalysesList(int** panalyses,int* pnumanalyses,IoModel* iomodel,int solutiontype);
		void Solve(void);

		/*Modules*/
		void AverageButtressingx(IssmDouble* ptheta);
		void BalancethicknessMisfitx(IssmDouble* pV);
		void CalvingRateVonmisesx();
		void CalvingRateLevermannx();
		void CalvingFluxLevelsetx();
		void CalvingMeltingFluxLevelsetx();
		void DeviatoricStressx();
		void Divergencex(IssmDouble* pdiv);
		void ElementOperationx(void (Element::*function)(void));
		void ElementResponsex(IssmDouble* presponse,int response_enum);
		void FloatingAreax(IssmDouble* pV, bool scaled);
		void GetInputLocalMinMaxOnNodesx(IssmDouble** pmin,IssmDouble** pmax,IssmDouble* ug);
		void GetLocalVectorWithClonesGset(IssmDouble** plocal_ug,Vector<IssmDouble> *ug);
		void GetLocalVectorWithClonesVertices(IssmDouble** plocal_vector,Vector<IssmDouble> *vector);
		void SyncLocalVectorWithClonesVertices(IssmDouble* local_vector);
		void SyncLocalVectorWithClonesVerticesAdd(IssmDouble* local_vector);
		void GetLocalVectorWithClonesNodes(IssmDouble** plocal_vector,Vector<IssmDouble> *vector);
		void GroundedAreax(IssmDouble* pV, bool scaled);
		void IcefrontAreax();
		void IcefrontMassFluxx(IssmDouble* presponse, bool scaled);
		void IcefrontMassFluxLevelsetx(IssmDouble* presponse, bool scaled);
		void GroundinglineMassFluxx(IssmDouble* presponse, bool scaled);
		void IceMassx(IssmDouble* pV, bool scaled);
		void IceVolumex(IssmDouble* pV, bool scaled);
		void IceVolumeAboveFloatationx(IssmDouble* pV, bool scaled);
		void InputToP0(int inputenum,int outputenum);
		void InputMakeDiscontinuous(int enum_in);
		void MassFluxx(IssmDouble* presponse);
		void MaxAbsVxx(IssmDouble* presponse);
		void MaxAbsVyx(IssmDouble* presponse);
		void MaxAbsVzx(IssmDouble* presponse);
		void MaxDivergencex(IssmDouble* pdiv);
		void MaxVelx(IssmDouble* presponse);
		void MaxVxx(IssmDouble* presponse);
		void MaxVyx(IssmDouble* presponse);
		void MaxVzx(IssmDouble* presponse);
		void MinVelx(IssmDouble* presponse);
		void MinVxx(IssmDouble* presponse);
		void MinVyx(IssmDouble* presponse);
		void MinVzx(IssmDouble* presponse);
		void MmeToInputFromId(int id, int rootenum, int interpolationenum);
		void DistanceToFieldValue(int fieldenum,IssmDouble fieldvalue,int distanceenum);
		void ResetLevelset();
		void StrainRateparallelx();
		void StrainRateperpendicularx();
		void StrainRateeffectivex();
		void StressIntensityFactorx();
		void RignotMeltParameterizationx();
		void TotalCalvingFluxLevelsetx(IssmDouble* pGbmb, bool scaled);
		void TotalCalvingMeltingFluxLevelsetx(IssmDouble* pGbmb, bool scaled);
		void TotalFloatingBmbx(IssmDouble* pFbmb, bool scaled);
		void TotalGroundedBmbx(IssmDouble* pGbmb, bool scaled);
		void TotalSmbx(IssmDouble* pSmb, bool scaled);
		void TotalSmbMeltx(IssmDouble* pSmbMelt, bool scaled);
		void TotalSmbRefreezex(IssmDouble* pSmbRefreeze, bool scaled);
		#ifdef  _HAVE_DAKOTA_
		void DakotaResponsesx(double* d_responses,char** responses_descriptors,int numresponsedescriptors,int d_numresponses);
		#endif
		void CostFunctionx(IssmDouble* pJ,IssmDouble** pJlist,int* pn);
		void OutputControlsx(Results **presults);
		void RequestedDependentsx(void);
		void RequestedOutputsx(Results **presults,char** requested_outputs, int numoutputs,bool save_results=true);
		void RequestedOutputsx(Results **presults,int* requested_outputs, int numoutputs,bool save_results=true);
		void Responsex(IssmDouble* presponse,int response_descriptor_enum);
		void Responsex(IssmDouble* presponse,const char* response_descriptor);
		void SurfaceAbsMisfitx( IssmDouble* pJ);
		void OmegaAbsGradientx( IssmDouble* pJ);
		void EtaDiffx( IssmDouble* pJ);
		void ThicknessAverage();
		void ThicknessAbsGradientx( IssmDouble* pJ);
		void ThicknessPositivex(IssmDouble* pJ);
		#ifdef _HAVE_ESA_
		void EsaGeodetic2D(Vector<IssmDouble>* pUp, Vector<IssmDouble>* pNorth, Vector<IssmDouble>* pEast, Vector<IssmDouble>* pGravity, Vector<IssmDouble>* pX, Vector<IssmDouble>* pY, IssmDouble* xx, IssmDouble* yy);
      void EsaGeodetic3D(Vector<IssmDouble>* pUp, Vector<IssmDouble>* pNorth, Vector<IssmDouble>* pEast, Vector<IssmDouble>* pGravity, IssmDouble* latitude, IssmDouble* longitude, IssmDouble* radius, IssmDouble* xx, IssmDouble* yy, IssmDouble* zz);
		#endif
		void HydrologyEPLupdateDomainx(IssmDouble* pEplcount);
		void HydrologyIDSupdateDomainx(IssmDouble* pIDScount);
		void TimeAdaptx(IssmDouble* pdt);
		void UpdateConstraintsExtrudeFromBasex();
		void UpdateConstraintsExtrudeFromTopx();
		void UpdateConstraintsL2ProjectionEPLx(IssmDouble* pL2count);
		void InitTransientInputx(int* transientinput_enum,int numoutputs);
		void StackTransientInputx(int* input_enum,int* transientinput_enum,IssmDouble hydrotime,int numoutputs);
		void StackTransientInputonBasex(int* input_enum,int* transientinput_enum,IssmDouble hydrotime,int numoutputs);
		void AverageTransientInputx(int* transientinput_enum,int* averagedinput_enum,IssmDouble init_time,IssmDouble end_time,int numoutputs,int averaging_method);
		void AverageTransientInputonBasex(int* transientinput_enum,int* averagedinput_enum,IssmDouble init_time,IssmDouble end_time,int numoutputs,int averaging_method);
		void UpdateConstraintsx(void);
		int  UpdateVertexPositionsx(void);

		#ifdef _HAVE_JAVASCRIPT_
		FemModel(IssmDouble* buffer, int buffersize, char* toolkits, char* solution, char* modelname,ISSM_MPI_Comm incomm, bool trace=false);
		void CleanUpJs(char** poutput, size_t* psize);
		void InitFromBuffers(char* buffer, int buffersize, char* toolkits, int solution_type,bool trace,IssmPDouble* X=NULL);
		#endif

		/*AMR*/
		#if !defined(_HAVE_AD_)
		void ReMesh(void);
		void BedrockFromMismipPlus(void);
		void AdjustBaseThicknessAndMask(void);
		void GetMesh(Vertices* femmodel_vertices,Elements* femmodel_elements,IssmDouble** px, IssmDouble** py, int** pelementslist);
		void GetMesh(int** elementslist, IssmDouble** x, IssmDouble** y, int* numberofvertices, int* numberofelements);
		void SetMesh(int** elementslist, IssmDouble** x, IssmDouble** y, int* numberofvertices, int* numberofelements);
		void GetMeshOnPartition(Vertices* femmodel_vertices,Elements* femmodel_elements,IssmDouble** px, IssmDouble** py, IssmDouble** pz, int** pelementslist,int** psidtoindex);
		void CreateElements(int newnumberofelements,int elementswidth,int* newelementslist,bool* my_elements,Elements* elements);
		void CreateMaterials(int newnumberofelements,bool* my_elements,Materials* materials);
		void CreateConstraints(Vertices* newfemmodel_vertices,int analysis_enum,Constraints* newfemmodel_constraints);
		void GetInputs(int* pnumP0inputs,IssmDouble** pP0inputs,int** pP0input_enums,int** pP0input_interp,int* pnumP1inputs,IssmDouble** pP1inputs,int** pP1input_enums,int** pP1input_interp);
		void InterpolateInputs(Vertices* newfemmodel_vertices,Elements* newfemmodel_elements,Inputs* new_inputs);
		void UpdateElements(int newnumberofelements,int* newelementslist,bool* my_elements,int analysis_counter,Elements* newelements);
		void WriteMeshInResults(void);
		void WriteErrorEstimatorsInResults(void);
		void SmoothedDeviatoricStressTensor(IssmDouble** ptauxx,IssmDouble** ptauyy,IssmDouble** ptauxy); //nodal values, just for SSA-P1: TauXX, TauYY, TauXY
		void ZZErrorEstimator(IssmDouble** pelementerror);
		void SmoothedGradThickness(IssmDouble** pdHdx,IssmDouble** pdHdy);
		void ThicknessZZErrorEstimator(IssmDouble** pelementerror);
		void MeanGroundedIceLevelSet(IssmDouble** pmasklevelset);
		void GetElementCenterCoordinates(IssmDouble** pxc,IssmDouble** pyc);
		void GetZeroLevelSetPoints(IssmDouble** pzerolevelset_points,int &numberofpoints,int levelset_type);
		#endif

		#if defined(_HAVE_BAMG_) && !defined(_HAVE_AD_)
		void ReMeshBamg(int* pnewnumberofvertices,int* pnewnumberofelements,IssmDouble** pnewx,IssmDouble** pnewy,IssmDouble** pnewz,int** pnewelementslist);
		void InitializeAdaptiveRefinementBamg(void);
		void GethmaxVerticesFromZeroLevelSetDistance(IssmDouble* hmaxvertices,int levelset_type);
		void GethmaxVerticesFromEstimators(IssmDouble* hmaxvertices,int errorestimator_type);
		void GetVerticeDistanceToZeroLevelSet(IssmDouble** pverticedistance,int leveset_type);
		#endif

		#if defined(_HAVE_NEOPZ_) && !defined(_HAVE_AD_)
		void ReMeshNeopz(int* pnewnumberofvertices,int* pnewnumberofelements,IssmDouble** pnewx,IssmDouble** pnewy,IssmDouble** pnewz,int** pnewelementslist);
		void InitializeAdaptiveRefinementNeopz(void);
		void GetElementDistanceToZeroLevelSet(IssmDouble** pelementdistance,int levelset_type);
		void SetRefPatterns(void);
		#endif
};

#endif
