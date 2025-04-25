/*!\file:  Element.h
 * \brief abstract class for Element object
 * This class is a place holder for the Tria and the Penta elements.
 * It is derived from Element, so DataSets can contain them.
*/

#ifndef _ELEMENT_H_
#define _ELEMENT_H_

/*Headers:*/
/*{{{*/
#include "../../datastructures/datastructures.h"
#include "../../toolkits/toolkits.h"
class DataSet;
class Parameters;
class Parameter;
class Elements;
class Loads;
class Nodes;
class Node;
class Vertices;
class Vertex;
class Materials;
class Material;
class Matlitho;
class Inputs;
class Inputs;
class Input;
class Input;
class ElementInput;
class DatasetInput;
class IoModel;
class SealevelGeometry;
class Gauss;
class GrdLoads;
class ElementVector;
template <class doublematrix> class Matrix;
template <class doubletype> class Vector;
class ElementMatrix;
class ElementVector;
class BarystaticContributions;
/*}}}*/

class Element: public Object{

	public:
		int          id;
		int          sid;
		int          lid;
		Inputs     *inputs;
		Node       **nodes;
		Vertex     **vertices;
		Material    *material;
		Parameters  *parameters;
		bool         isonsurface;
		bool         isonbase;

		int* element_type_list;
		int  element_type;

	public:
		/*Constructors/Destructores*/
		Element();
		~Element();

		/*Functions*/
		/*bool               AllActive(void);*/
		/*bool               AnyActive(void);*/
		bool               AnyFSet(void);
      void               ArmaProcess_pre18Oct2022(bool isstepforarma,int arorder,int maorder,IssmDouble telapsed_arma,IssmDouble tstep_arma,IssmDouble* termconstant,IssmDouble* trend,IssmDouble* arlagcoefs,IssmDouble* malagcoefs,bool isfieldstochastic,int enum_type);
      void               ArmaProcess(bool isstepforarma,int arorder,int maorder,int numparams,int numbreaks,IssmDouble tstep_arma,IssmDouble* polyparams,IssmDouble* arlagcoefs,IssmDouble* malagcoefs,IssmDouble* datebreaks,bool isfieldstochastic,int enum_type);
		void               BasinLinearFloatingiceMeltingRate(IssmDouble* deepwaterel,IssmDouble* upperwatermelt,IssmDouble* upperwaterel,IssmDouble* perturbation);
		void               CalvingSetZeroRate(void);
		void               CalvingRateToVector();
		void               CalvingRateToVector(bool isvelvector);
		void               ComputeLambdaS(void);
		void               ComputeNewDamage();
		void               ComputeStrainRate();
		void               CoordinateSystemTransform(IssmDouble** ptransform,Node** nodes,int numnodes,int* cs_array);
		void               DeepEcho();
		void               DeleteMaterials(void);
		void               Delta18oParameterization(void);
		void               Delta18opdParameterization(void);
		void               SmbGradCompParameterization(void);
		IssmDouble         Divergence(void);
		void               dViscositydBFS(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input);
		void               dViscositydBHO(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               dViscositydBSSA(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               dViscositydBMOLHO(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vxbase_input,Input* vybase_input, Input* vxshear_input ,Input* vyshear_input,Input* thickness_input,Input* n_input, IssmDouble zeta);
		void               dViscositydDSSA(IssmDouble* pdmudB,int dim,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               Echo();
		void               FindParam(bool* pvalue,int paramenum);
		void               FindParam(int* pvalue,int paramenum);
		void               FindParam(IssmDouble* pvalue,int paramenum);
		IssmDouble         FindParam(int paramenum);
		void               FindParam(int** pvalues,int* psize,int paramenum);
		IssmDouble         FloatingArea(IssmDouble* mask, bool scaled);
		void	             GetDofList(int** pdoflist,int approximation_enum,int setenum);
		void	             GetDofListPressure(int** pdoflist,int setenum);
		void	             GetDofListVelocity(int** pdoflist,int setenum);
		void	             GetDofListLocal(int** pdoflist,int approximation_enum,int setenum);
		void	             GetDofListLocalPressure(int** pdoflist,int setenum);
		void	             GetDofListLocalVelocity(int** pdoflist,int setenum);
		void               GetInputListOnNodes(IssmDouble* pvalue,int enumtype);
		void               GetInputListOnNodes(IssmDouble* pvalue,int enumtype,IssmDouble defaultvalue);
		void               GetInputListOnNodesVelocity(IssmDouble* pvalue,int enumtype);
		void               GetInputListOnVertices(IssmDouble* pvalue,int enumtype);
		void               GetInputListOnVerticesAtTime(IssmDouble* pvalue,int enumtype,IssmDouble time);
		void               GetInputListOnVertices(IssmDouble* pvalue,int enumtype,IssmDouble defaultvalue);
		void               GetInputLocalMinMaxOnNodes(IssmDouble* min,IssmDouble* max,IssmDouble* ug);
		void               GetInputValue(bool* pvalue,int enum_type);
		void               GetInputValue(int* pvalue,int enum_type);
		void               GetInputValue(IssmDouble* pvalue,int enum_type);
		void               GetInputValue(IssmDouble* pvalue,Gauss* gauss,int enum_type);
		Node*              GetNode(int nodeindex);
		int                GetNodeIndex(Node* node);
		void               GetNodesLidList(int* lidlist);
		void               GetNodesSidList(int* sidlist);
		void               GetPhi(IssmDouble* phi, IssmDouble*  epsilon, IssmDouble viscosity);
		void               GetSolutionFromInputsOneDof(Vector<IssmDouble>* solution,int solutionenum);
		void               GetVectorFromInputs(Vector<IssmDouble>* vector, int name_enum, int type);
		void               GetVectorFromInputs(Vector<IssmDouble>* vector, int name_enum, int type,IssmDouble time);
		void	             GetVerticesLidList(int* lidlist);
		void	             GetVerticesPidList(int* pidlist);
		void               GetVerticesConnectivityList(int* connectivitylist);
		void               GetVerticesCoordinates(IssmDouble** xyz_list);
		void               GetVerticesSidList(int* sidlist);
		IssmDouble         GetXcoord(IssmDouble* xyz_list,Gauss* gauss);
		IssmDouble         GetYcoord(IssmDouble* xyz_list,Gauss* gauss);
		IssmDouble         GetZcoord(IssmDouble* xyz_list,Gauss* gauss);
		void               GradientIndexing(int* indexing,int control_index);
		IssmDouble         GroundedArea(IssmDouble* mask, bool scaled);
		bool               HasNodeOnBase();
		bool               HasNodeOnSurface();
		IssmDouble         IceMass(bool scaled);
		IssmDouble         IceMass(IssmDouble* mask, bool scaled);
		IssmDouble         IceVolume(IssmDouble* mask, bool scaled);
		IssmDouble         IceVolumeAboveFloatation(IssmDouble* mask, bool scaled);
		int                Id();
		void               InputCreate(IssmDouble* vector,Inputs* inputs,IoModel* iomodel,int M,int N,int vector_type,int vector_enum,int code);
		void               InputCreateP1FromConstant(Inputs* inputs,IoModel* iomodel,IssmDouble value,int vector_enum);
		void               InputCreateP0FromConstant(Inputs* inputs,IoModel* iomodel,IssmDouble value,int vector_enum);
		void               ControlInputCreate(IssmDouble* doublearray,IssmDouble* independents_min,IssmDouble* independents_max,Inputs*inputs,IoModel* iomodel,int M,int N,IssmDouble scale,int input_enum,int id);
		void					 DatasetInputAdd(int enum_type,IssmDouble* vector,Inputs* inputs,IoModel* iomodel,int M,int N,int vector_type,int vector_enum,int input_enum);
		void               InputUpdateFromConstant(IssmDouble constant, int name);
		void               InputUpdateFromConstant(int constant, int name);
		void               InputUpdateFromConstant(bool constant, int name);
		void               InputUpdateFromConstant(IssmDouble constant, int name, int type);

		bool               IsAllFloating();
		bool               IsAllGrounded();
		bool               IsGrounded();
		bool               IsOnBase();
		bool               IsOnSurface();
		bool               IsIceInElement();
		bool               IsIceOnlyInElement();
		bool               IsOceanInElement();
		bool               IsOceanOnlyInElement();
		bool		   IsAllMinThicknessInElement();
		bool               IsLandInElement();
		void               Ismip6FloatingiceMeltingRate();
		void               LapseRateBasinSMB(int numelevbins, IssmDouble* lapserates, IssmDouble* elevbins,IssmDouble* refelevation);
		void               LinearFloatingiceMeltingRate();
		void               SpatialLinearFloatingiceMeltingRate();
		void               MantlePlumeGeothermalFlux();
		void               MarshallElement2(MarshallHandle* marshallhandle,int numanalyses);
		void               MigrateGroundingLine(IssmDouble* sheet_ungrounding);
		void               MismipFloatingiceMeltingRate();
		void               MonthlyFactorBasin(IssmDouble* monthlyfac,int enum_type); 
		void               MonthlyPiecewiseLinearEffectBasin(int nummonthbreaks,IssmDouble* monthlyintercepts,IssmDouble* monthlytrends,IssmDouble* monthlydatebreaks,int enum_type); 
		void               BeckmannGoosseFloatingiceMeltingRate();
		void               MungsmtpParameterization(void);
		ElementMatrix*     NewElementMatrix(int approximation_enum=NoneApproximationEnum);
		ElementMatrix*     NewElementMatrixCoupling(int number_nodes,int approximation_enum=NoneApproximationEnum);
		ElementVector*     NewElementVector(int approximation_enum=NoneApproximationEnum);
		void               PicoUpdateBoxid(int* pmax_boxid_basin);
		void               PicoUpdateBox(int loopboxid);
		void               PicoComputeBasalMelt();
		void               PositiveDegreeDay(IssmDouble* pdds,IssmDouble* pds,IssmDouble signorm,bool ismungsm,bool issetpddfac);
		void               PositiveDegreeDaySicopolis(bool isfirnwarming);
		void               PositiveDegreeDayGCM();
		void               ProjectGridDataToMesh(IssmDouble* griddata,IssmDouble* x_grid,IssmDouble* y_grid,int Nx,int Ny,int input_enum);
		void               SmbDebrisEvatt();
		void               RignotMeltParameterization();
		void               ResultInterpolation(int* pinterpolation,int*nodesperelement,int* parray_size, int output_enum);
		void               ResultToPatch(IssmDouble* values,int nodesperelement,int output_enum);
		void               ResultToMatrix(IssmDouble* values,int ncols,int output_enum);
		void               ResultToVector(Vector<IssmDouble>* vector,int output_enum);
		void               SetwiseNodeConnectivity(int* d_nz,int* o_nz,Node* node,bool* flags,int* flagsindices,int* flagsindices_counter,int set1_enum,int set2_enum, int analysis_type);
		void               SetBoolInput(Inputs* inputs,int enum_in,bool value);

		void               SetIntInput(Inputs* inputs,int enum_in,int value);
		void               SmbSemic();
		void               SmbSemicTransient();
		int                Sid();
		void               SmbGemb(IssmDouble timeinputs, int count, int steps);
		void               StrainRateESA(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               StrainRateFS(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input);
		void               StrainRateHO(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               StrainRateHO2dvertical(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               StrainRateMOLHO(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vxshear_input,Input* vyshear_input,Input* thickness_input,Input* n_input,IssmDouble zeta);
		void               StrainRateSSA(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input);
		void               StrainRateSSA1d(IssmDouble* epsilon,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input);
		void               StressMaxPrincipalCreateInput(void);
		void               SubglacialWaterPressure(int output_enum);
		IssmDouble         TotalFloatingBmb(IssmDouble* mask, bool scaled);
		IssmDouble         TotalGroundedBmb(IssmDouble* mask, bool scaled);
		IssmDouble         TotalSmb(IssmDouble* mask, bool scaled);
		IssmDouble         TotalSmbMelt(IssmDouble* mask, bool scaled);
		IssmDouble         TotalSmbRefreeze(IssmDouble* mask, bool scaled);
		void               TransformInvStiffnessMatrixCoord(ElementMatrix* Ke,int cs_enum);
		void               TransformInvStiffnessMatrixCoord(ElementMatrix* Ke,Node** nodes,int numnodes,int cs_enum);
		void               TransformInvStiffnessMatrixCoord(ElementMatrix* Ke,Node** nodes,int numnodes,int* cs_array);
		void               TransformLoadVectorCoord(ElementVector* pe,int cs_enum);
		void               TransformLoadVectorCoord(ElementVector* pe,int* cs_array);
		void               TransformLoadVectorCoord(ElementVector* pe,Node** nodes,int numnodes,int cs_enum);
		void               TransformLoadVectorCoord(ElementVector* pe,Node** nodes,int numnodes,int* cs_array);
		void               TransformLoadVectorCoord(ElementVector* pe,int numnodes,int transformenum){_error_("not implemented yet");};/*Tiling only*/
		void               TransformLoadVectorCoord(ElementVector* pe,int numnodes,int* transformenum_list){_error_("not implemented yet");};/*Tiling only*/
		void               TransformSolutionCoord(IssmDouble* solution,int cs_enum);
		void               TransformSolutionCoord(IssmDouble* solution,int* cs_array);
		void               TransformSolutionCoord(IssmDouble* solution,int numnodes,int cs_enum);
		void               TransformSolutionCoord(IssmDouble* solution,int numnodes,int* cs_array);
		void               TransformSolutionCoord(IssmDouble* solution,Node** nodes,int numnodes,int cs_enum);
		void               TransformSolutionCoord(IssmDouble* solution,Node** nodes,int numnodes,int* cs_array);
		void               TransformStiffnessMatrixCoord(ElementMatrix* Ke,int cs_enum);
		void               TransformStiffnessMatrixCoord(ElementMatrix* Ke,int* cs_array);
		void               TransformStiffnessMatrixCoord(ElementMatrix* Ke,Node** nodes,int numnodes,int cs_enum);
		void               TransformStiffnessMatrixCoord(ElementMatrix* Ke,Node** nodes,int numnodes,int* cs_array);
		void               TransformStiffnessMatrixCoord(ElementMatrix* Ke,int numnodes,int* transformenum_list){_error_("not implemented yet");};/*Tiling only*/
		void               FrictionAlpha2CreateInput(void);
		void               ViscousHeatingCreateInput(void);
		void               ThermalToEnthalpy(IssmDouble * penthalpy,IssmDouble temperature,IssmDouble waterfraction,IssmDouble pressure);
		IssmDouble         TMeltingPoint(IssmDouble pressure);
		void               EnthalpyToThermal(IssmDouble* ptemperature,IssmDouble* pwaterfraction,IssmDouble enthalpy,IssmDouble pressure);
		IssmDouble         EnthalpyDiffusionParameter(IssmDouble enthalpy,IssmDouble pressure);
		IssmDouble         EnthalpyDiffusionParameterVolume(int numvertices,IssmDouble* enthalpy,IssmDouble* pressure);
		IssmDouble         PureIceEnthalpy(IssmDouble pressure);

		/*Virtual functions*/
		virtual void       AddBasalInput(int input_enum, IssmDouble* values, int interpolation_enum){_error_("not implemented");};
		virtual void       AddInput(int input_enum, IssmDouble* values, int interpolation_enum){_error_("not implemented");};
		virtual void       AddControlInput(int input_enum,Inputs* inputs,IoModel* iomodel,IssmDouble* values,IssmDouble* values_min,IssmDouble* values_max, int interpolation_enum,int id){_error_("not supported yet");};
		virtual void       DatasetInputCreate(IssmDouble* array,int M,int N,int* individual_enums,int num_inputs,Inputs* inputs,IoModel* iomodel,int input_enum){_error_("not supported");};
		virtual void       AverageOntoPartition(Vector<IssmDouble>* partition_contributions,Vector<IssmDouble>* partition_areas,IssmDouble* vertex_response,IssmDouble* qmu_part)=0;
		virtual void		 BasalNodeIndices(int* pnumindices,int** pindices,int finiteelement){_error_("not implemented yet");};
		virtual void       CalvingRateParameterization(void){_error_("not implemented yet");};
		virtual void       CalvingRateCalvingMIP(void){_error_("not implemented yet");};
		virtual void       CalvingRateVonmises(void){_error_("not implemented yet");};
		virtual void       CalvingRateVonmisesAD(void){_error_("not implemented yet");};
		virtual void       CalvingRateTest(void){_error_("not implemented yet");};
		virtual void       CalvingCrevasseDepth(void){_error_("not implemented yet");};
		virtual void	    CalvingRateLevermann(void)=0;
		virtual void	    CalvingPollard(void){_error_("not implemented yet");};
		virtual void       CalvingFluxLevelset(void){_error_("not implemented yet");};
		virtual void       CalvingMeltingFluxLevelset(void){_error_("not implemented yet");};
		virtual IssmDouble CharacteristicLength(void)=0;
		virtual void       ComputeBasalStress(void){_error_("not implemented yet");};
		virtual void       ComputeDeviatoricStressTensor(void)=0;
		virtual void       ComputeSigmaNN(void)=0;
		virtual void       ComputeSigmaVM(void){_error_("not implemented yet");};
		virtual void       ComputeStressTensor(void)=0;
		virtual void       ComputeEsaStrainAndVorticity(void)=0;
		virtual void       Configure(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters,Inputs* inputsin)=0;
		virtual void       ControlInputSetGradient(IssmDouble* gradient,int enum_type,int control_index,int offset,int M,int N,int interp)=0;
		virtual void       ControlToVectors(Vector<IssmPDouble>* vector_control, Vector<IssmPDouble>* vector_gradient,int control_enum,int control_interp)=0;
		virtual void       CreateDistanceInputFromSegmentlist(IssmDouble* distances,int distanceenum){_error_("not implemented yet");};
		virtual void       CreateInputTimeAverage(int transientinput_enum,int averagedinput_enum,IssmDouble init_time,IssmDouble end_time,int averaging_method){_error_("not implemented yet "<<this->ObjectEnum());};
		virtual void       ElementResponse(IssmDouble* presponse,int response_enum)=0;
		virtual void       ElementSizes(IssmDouble* phx,IssmDouble* phy,IssmDouble* phz)=0;
		virtual void       ElementCoordinates(Vector<IssmDouble>* vxe,Vector<IssmDouble>* vye,Vector<IssmDouble>* vze, Vector<IssmDouble>* vareae, bool spherical=false)=0;
		virtual void       ElementCoordinates(Vector<IssmDouble>* vlonge,Vector<IssmDouble>* vlate,Vector<IssmDouble>* vareae)=0;
		virtual int        FiniteElement(void)=0;
		virtual IssmDouble FloatingArea(bool scaled)=0;
		virtual void       FSContactMigration(Vector<IssmDouble>* vertex_sigmann,Vector<IssmDouble>* vertex_waterpressure)=0;
		virtual Element*   GetBasalElement(void)=0;
		virtual int        GetElementType(void)=0;
		virtual IssmDouble GetHorizontalSurfaceArea(void){_error_("not implemented");};
		virtual void       GetGroundedPart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlyfloating,int distance_enum, IssmDouble intrusion_distance)=0;
		virtual IssmDouble GetGroundedPortion(IssmDouble* xyz_list)=0;
		virtual void        GetFractionGeometry(IssmDouble* weights, IssmDouble* pphi, int* ppoint1,IssmDouble* pfraction1,IssmDouble* pfraction2, bool* ptrapezeisnegative, IssmDouble* gl)=0;
		virtual void       GetNodalWeightsAndAreaAndCentroidsFromLeveset(IssmDouble* loadweights, IssmDouble* ploadarea, IssmDouble* platbar, IssmDouble* plongbar, IssmDouble late, IssmDouble longe, IssmDouble area,  int levelsetenum)=0;
		virtual void       GetNodalWeightsAndAreaAndCentroidsFromLeveset(IssmDouble* loadweights, IssmDouble* ploadarea, IssmDouble* platbar, IssmDouble* plongbar, IssmDouble late, IssmDouble longe, IssmDouble area, int levelset1enum, int levelset2enum)=0;
		virtual void        GetBarycenterFromLevelset(IssmDouble* platbar, IssmDouble* plongbar,IssmDouble phi,IssmDouble fraction1,IssmDouble fraction2,IssmDouble late, IssmDouble longe, int point1,int istrapeze1, IssmDouble planetradius)=0;

		virtual IssmDouble GetIcefrontArea(){_error_("not implemented");};
		virtual void       GetIcefrontCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum)=0;
		virtual DatasetInput* GetDatasetInput(int inputenum){_error_("not implemented");};
		virtual Input*    GetInput(int inputenum)=0;
		virtual Input*    GetInput(int inputenum,IssmDouble time)=0;
		virtual Input*    GetInput(int inputenum,IssmDouble start_time,IssmDouble end_time,int averaging_method)=0;
		virtual void       GetInputValue(IssmDouble* pvalue,Vertex* vertex,int enumtype){_error_("not implemented yet");};
		virtual void       GetInputValue(IssmDouble* pvalue,Node* node,int enumtype){_error_("not implemented yet");};
		virtual void       GetInputListOnVertices(IssmDouble* pvalue,Input* input,IssmDouble default_value)=0;
		virtual void       GetInputListOnNodes(IssmDouble* pvalue,Input* input,IssmDouble default_value)=0;
		virtual void       GetLevelCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum,IssmDouble level)=0;
		virtual void       GetLevelsetPositivePart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlynegative,IssmDouble* levelsetvalues)=0;
		virtual int        GetVertexIndex(Vertex* vertex){_error_("not implemented");};;
		virtual int        GetNumberOfNodes(void)=0;
		virtual int        GetNumberOfNodes(int enum_type)=0;
		virtual int        GetNumberOfVertices(void)=0;
		virtual void       GetVectorFromControlInputs(Vector<IssmDouble>* gradient,int control_enum,int control_index,int N,const char* data,int offset)=0;
		virtual void       GetVerticesCoordinatesBase(IssmDouble** xyz_list)=0;
		virtual void       GetVerticesCoordinatesTop(IssmDouble** xyz_list)=0;
		virtual IssmDouble GroundedArea(bool scaled)=0;
		virtual IssmDouble IceVolume(bool scaled)=0;
		virtual IssmDouble IceVolumeAboveFloatation(bool scaled)=0;
		virtual IssmDouble IcefrontMassFlux(bool scaled){_error_("not implemented");};
		virtual IssmDouble IcefrontMassFluxLevelset(bool scaled){_error_("not implemented");};
		virtual IssmDouble GroundinglineMassFlux(bool scaled){_error_("not implemented");};
		virtual void       InputDepthAverageAtBase(int enum_type,int average_enum_type)=0;
		virtual void       DatasetInputExtrude(int input_enum,int start){_error_("not implemented yet");};
		virtual void       InputExtrude(int input_enum,int start)=0;
		virtual void       InputUpdateFromSolutionOneDofCollapsed(IssmDouble* solution,int inputenum)=0;
		virtual void       InputUpdateFromSolutionOneDof(IssmDouble* solution,int inputenum)=0;
		#ifdef _HAVE_DAKOTA_
		virtual void  InputUpdateFromMatrixDakota(IssmDouble* matrix, int rows, int ncols, int name, int type)=0;
		virtual void  InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type)=0;
		virtual void  InputScaleFromDakota(IssmDouble* distributed_values, IssmDouble* partition, int npart, int nt, int name)=0;
		#endif
		virtual void  InputUpdateFromIoModel(int index, IoModel* iomodel)=0;
		virtual void  InputUpdateFromVector(IssmDouble* vector, int name, int type)=0;
		virtual bool       IsFaceOnBoundary(void)=0;
		virtual bool       IsIcefront(void)=0;
		virtual bool       IsNodeOnShelfFromFlags(IssmDouble* flags)=0;

		virtual bool       IsZeroLevelset(int levelset_enum)=0;
		virtual void       JacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       JacobianDeterminantBase(IssmDouble* Jdet,IssmDouble* xyz_list_base,Gauss* gauss)=0;
		virtual void       JacobianDeterminantLine(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       JacobianDeterminantSurface(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       JacobianDeterminantTop(IssmDouble* Jdet,IssmDouble* xyz_list_base,Gauss* gauss)=0;
		virtual void       Marshall(MarshallHandle* marshallhandle)=0;
		virtual IssmDouble Masscon(IssmDouble* levelset)=0;
		virtual IssmDouble MassFlux(IssmDouble* segment)=0;
		virtual IssmDouble MassFlux(IssmDouble x1,IssmDouble y1, IssmDouble x2, IssmDouble y2,int segment_id)=0;
		virtual IssmDouble MinEdgeLength(IssmDouble* xyz_list)=0;
		virtual IssmDouble Misfit(int modelenum,int observationenum,int weightsenum)=0;
		virtual IssmDouble MisfitArea(int weightsenum)=0;
		virtual void	   MovingFrontalVelocity(void){_error_("not implemented yet");};
		virtual Gauss*     NewGauss(void)=0;
		virtual Gauss*     NewGauss(int order)=0;
		virtual Gauss*     NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order)=0;
		virtual Gauss*     NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz,int order_vert)=0;
		virtual Gauss*     NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,bool mainlyfloating,int order)=0;
		virtual Gauss*     NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,int order)=0;
      virtual Gauss*     NewGauss(IssmDouble fraction1,IssmDouble fraction2,int order)=0;
		virtual Gauss*     NewGaussBase(int order)=0;
		virtual Gauss*     NewGaussLine(int vertex1,int vertex2,int order)=0;
		virtual Gauss*     NewGaussTop(int order)=0;
		virtual void       NodalFunctions(IssmDouble* basis,Gauss* gauss)=0;
		virtual void       NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       NodalFunctionsDerivativesVelocity(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       NodalFunctionsMINIDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       NodalFunctionsPressure(IssmDouble* basis, Gauss* gauss)=0;
		virtual void       NodalFunctionsP1(IssmDouble* basis,Gauss* gauss)=0;
		virtual void       NodalFunctionsP1Derivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       NodalFunctionsP2(IssmDouble* basis,Gauss* gauss)=0;
		virtual void       NodalFunctionsVelocity(IssmDouble* basis, Gauss* gauss)=0;
		virtual void       NodalFunctionsTensor(IssmDouble* basis, Gauss* gauss)=0;
		virtual int        NodalValue(IssmDouble* pvalue, int index, int natureofdataenum)=0;
		virtual void       NormalBase(IssmDouble* normal,IssmDouble* xyz_list)=0;
		virtual void       NormalSection(IssmDouble* normal,IssmDouble* xyz_list)=0;
		virtual void       NormalTop(IssmDouble* normal,IssmDouble* xyz_list)=0;
		virtual int        NumberofNodesPressure(void)=0;
		virtual int        NumberofNodesVelocity(void)=0;
		virtual void       PotentialUngrounding(Vector<IssmDouble>* potential_sheet_ungrounding)=0;
		virtual int        PressureInterpolation()=0;
      virtual void       Recover3DMOLHOInput(int targetVel_enum, int numnodes, IssmDouble* vb,  IssmDouble* vsh, IssmDouble* n, IssmDouble* H, IssmDouble* s){_error_("not implemented yet");};
		virtual void       ReduceMatrices(ElementMatrix* Ke,ElementVector* pe)=0;
		virtual void       ResetFSBasalBoundaryCondition()=0;
		virtual void       ResetHooks()=0;
		virtual void       SetElementInput(int enum_in,IssmDouble value){_error_("not implemented yet");};
		virtual void       SetElementInput(int enum_in,IssmDouble value, int type)=0;
		virtual void       SetElementInput(Inputs* inputs,int enum_in,IssmDouble values){_error_("not implemented yet");};
		virtual void       SetElementInput(Inputs* inputs,int numindices,int* indices,IssmDouble* values,int enum_in){_error_("not implemented yet");};
		virtual void       SetControlInputsFromVector(IssmDouble* vector,int control_enum,int control_index,int offset,int M,int N)=0;
		virtual void       SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Materials* materials,Parameters* parameters)=0;
		virtual void       SetTemporaryElementType(int element_type_in)=0;
	   virtual Element*   SpawnBasalElement(bool depthaverage_materials=false)=0;
		virtual Element*   SpawnTopElement(void)=0;
		virtual bool       IsSpawnedElement(void)=0;
		virtual IssmDouble StabilizationParameter(IssmDouble u, IssmDouble v, IssmDouble w, IssmDouble diameter, IssmDouble kappa)=0;
		virtual void       StabilizationParameterAnisotropic(IssmDouble* tau_parameter_anisotropic, IssmDouble u, IssmDouble v, IssmDouble w, IssmDouble hx, IssmDouble hy, IssmDouble hz, IssmDouble kappa)=0;
		virtual void	    StrainRateparallel(void)=0;
		virtual void	    StrainRateperpendicular(void)=0;
		virtual void	    StressIntensityFactor(void)=0;
		virtual IssmDouble SurfaceArea(void)=0;
		virtual void       TangentBase(IssmDouble* bed_tangent,IssmDouble* bed_normal){_error_("not implemented yet");};
		virtual int        TensorInterpolation()=0;
		virtual IssmDouble TimeAdapt()=0;
		virtual IssmDouble TotalCalvingFluxLevelset(bool scaled){_error_("not implemented");};
		virtual IssmDouble TotalCalvingMeltingFluxLevelset(bool scaled){_error_("not implemented");};
		virtual IssmDouble TotalFloatingBmb(bool scaled)=0;
		virtual IssmDouble TotalGroundedBmb(bool scaled)=0;
		virtual IssmDouble TotalSmb(bool scaled)=0;
		virtual IssmDouble TotalSmbMelt(bool scaled)=0;
		virtual IssmDouble TotalSmbRefreeze(bool scaled)=0;
		virtual void       Update(Inputs* inputs,int index, IoModel* iomodel,int analysis_counter,int analysis_type,int finite_element)=0;
		virtual void       UpdateConstraintsExtrudeFromBase(void)=0;
		virtual void       UpdateConstraintsExtrudeFromTop(void)=0;
		virtual int        UpdatePotentialUngrounding(IssmDouble* potential_sheet_ungrounding,Vector<IssmDouble>* vec_nodes_on_iceshelf,IssmDouble* nodes_on_iceshelf)=0;
		virtual void       ValueP1DerivativesOnGauss(IssmDouble* dvalue,IssmDouble* values,IssmDouble* xyz_list,Gauss* gauss)=0;
		virtual void       ValueP1OnGauss(IssmDouble* pvalue,IssmDouble* values,Gauss* gauss)=0;
		virtual int        VelocityInterpolation()=0;
		virtual int        VertexConnectivity(int vertexindex)=0;
		virtual void       VerticalSegmentIndices(int** pindices,int* pnumseg)=0;
		virtual void       VerticalSegmentIndicesBase(int** pindices,int* pnumseg)=0;
		virtual void       ViscousHeating(IssmDouble* pphi,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){_error_("not implemented yet");};
		virtual void       WriteFieldIsovalueSegment(DataSet* segments,int fieldenum,IssmDouble fieldvalue){_error_("not implemented yet");};

		#ifdef _HAVE_ESA_
		virtual void          EsaGeodetic2D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast, Vector<IssmDouble>* pGravity, Vector<IssmDouble>* pX, Vector<IssmDouble>* pY,IssmDouble* xx,IssmDouble* yy)=0;
      virtual void          EsaGeodetic3D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,Vector<IssmDouble>* pGravity, IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz)=0;
		#endif
		#ifdef _HAVE_SEALEVELCHANGE_
		virtual IssmDouble    GetArea3D(void)=0;
		virtual IssmDouble    GetAreaSpherical(void)=0;
		virtual IssmDouble    GetTriangleAreaSpherical(IssmDouble xyz_list[3][3])=0;
		virtual void          GiaDeflection(Vector<IssmDouble>* wg,Vector<IssmDouble>* dwgdt,Matlitho* litho, IssmDouble* x,IssmDouble* y)=0;

		virtual void       SealevelchangeGeometrySubElementKernel(SealevelGeometry* slgeom)=0;
		virtual void       SealevelchangeShift(GrdLoads* loads, IssmDouble offset, SealevelGeometry* slgeom)=0;
		virtual void       SealevelchangeGeometryInitial(IssmDouble* xxe, IssmDouble* yye, IssmDouble* zze, IssmDouble* areae, int* lids, int* vcount)=0;
		virtual void       SealevelchangeGeometryCentroidLoads(SealevelGeometry* slgeom, IssmDouble* xxe, IssmDouble* yye, IssmDouble* zze, IssmDouble* areae)=0;
		virtual void       SealevelchangeGeometrySubElementLoads(SealevelGeometry* slgeom, IssmDouble* areae)=0;
		virtual void       SealevelchangeBarystaticLoads(GrdLoads* loads, BarystaticContributions* barycontrib, SealevelGeometry* slgeom)=0;
		virtual void       SealevelchangeConvolution(IssmDouble* sealevelpercpu, GrdLoads* loads, IssmDouble* rotationvector,SealevelGeometry* slgeom)=0;
		virtual void       SealevelchangeOceanAverage(GrdLoads* loads, Vector<IssmDouble>* oceanareas, Vector<IssmDouble>* subelementoceanareas, IssmDouble* sealevelpercpu, SealevelGeometry* slgeom)=0;
		virtual void       SealevelchangeDeformationConvolution(IssmDouble* sealevelpercpu, GrdLoads* loads, IssmDouble* rotationvector,SealevelGeometry* slgeom)=0;
		virtual void       SealevelchangeUpdateViscousFields(IssmDouble lincoeff, int newindex, int offset)=0;
		#endif

};
#endif
