/*! \file Penta.h
 *  \brief: header file for penta object
 */

#ifndef _PENTA_H_
#define _PENTA_H_

/*Headers:*/
/*{{{*/
#include "./Element.h"
#include "./ElementHook.h"
#include "./PentaRef.h"
class Object;
class Parameters;
class Results;
class Inputs;
class Input;
class IoModel;
class Node;
class Material;
class Matlitho;
class Tria;
class ElementMatrix;
class ElementVector;
class GaussPenta;
#include "../../shared/Exceptions/exceptions.h"
#include "../../shared/Enum/Enum.h"
/*}}}*/

class Penta: public Element,public ElementHook,public PentaRef{

	public:

		Penta      **verticalneighbors;           // 2 neighbors: first one under, second one above

		/*Penta constructors and destructor: {{{*/
		Penta(){};
		Penta(int penta_id,int penta_sid,int penta_lid,IoModel* iomodel,int nummodels);
		~Penta();
		/*}}}*/
		/*Object virtual functions definitions: {{{*/
		Object *copy();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*Penta routines:{{{*/
		void           AddBasalInput(int input_enum, IssmDouble* values, int interpolation_enum);
		void           AddInput(int input_enum, IssmDouble* values, int interpolation_enum);
		void           AddControlInput(int input_enum,Inputs* inputs,IoModel* iomodel,IssmDouble* values,IssmDouble* values_min,IssmDouble* values_max, int interpolation_enum,int id);
		void           ControlInputExtrude(int enum_type,int start);
		void           ComputeSigmaVM(void);
		void           DatasetInputExtrude(int enum_type,int start);
		void           DatasetInputCreate(IssmDouble* array,int M,int N,int* individual_enums,int num_inputs,Inputs* inputs,IoModel* iomodel,int input_enum);
		void           AverageOntoPartition(Vector<IssmDouble>* partition_contributions,Vector<IssmDouble>* partition_areas,IssmDouble* vertex_response,IssmDouble* qmu_part);
		void           BasalNodeIndices(int* pnumindices,int** pindices,int finiteelement);
		void           CalvingRateVonmises();
		void           CalvingRateLevermann();
		void           CalvingFluxLevelset();
		void           CalvingMeltingFluxLevelset();
		void           CalvingRateCalvingMIP();
		IssmDouble     CharacteristicLength(void){_error_("not implemented yet");};
		void           ComputeBasalStress(void);
		void           ComputeDeviatoricStressTensor();
		void           ComputeEsaStrainAndVorticity(){_error_("not implemented yet!");};
		void           ComputeSigmaNN(){_error_("not implemented yet");};
		void           ComputeStressTensor();
		//void           ComputeMeanEla(IssmDouble* paltitude, int* pcounter);
		void           Configure(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters,Inputs* inputsin);
		void           ControlInputSetGradient(IssmDouble* gradient,int enum_type,int control_index,int offset,int M,int N,int interp);
		void           ControlToVectors(Vector<IssmPDouble>* vector_control, Vector<IssmPDouble>* vector_gradient,int control_enum,int control_interp);
		void				CreateDistanceInputFromSegmentlist(IssmDouble* distances,int distanceenum);
		ElementMatrix* CreateBasalMassMatrix(void);
		void           CreateInputTimeAverage(int transientinput_enum,int averagedinput_enum,IssmDouble start_time,IssmDouble end_time,int averaging_method);
		void           ElementResponse(IssmDouble* presponse,int response_enum);
		void           ElementCoordinates(Vector<IssmDouble>* vxe,Vector<IssmDouble>* vye,Vector<IssmDouble>* vze,Vector<IssmDouble>* vareae,bool spherical=false){_error_("not implemented yet");};
		void           ElementCoordinates(Vector<IssmDouble>* vlonge,Vector<IssmDouble>* vlate,Vector<IssmDouble>* vareae){_error_("not implemented yet");};
		void           ElementSizes(IssmDouble* hx,IssmDouble* hy,IssmDouble* hz);
		int            FiniteElement(void);
		IssmDouble     FloatingArea(bool scaled);
		void           FSContactMigration(Vector<IssmDouble>* vertex_sigmann,Vector<IssmDouble>* vertex_waterpressure);
		IssmDouble     GetArea3D(void){_error_("not implemented yet!");};
		IssmDouble     GetAreaSpherical(void){_error_("not implemented yet!");};
		void           GetAreaCoordinates(IssmDouble *area_coordinates,IssmDouble* xyz_zero,IssmDouble* xyz_list,int numpoints);
		IssmDouble     GetTriangleAreaSpherical(IssmDouble xyz_list[3][3]){_error_("not implemented yet");};
		Element*       GetBasalElement(void);
		Penta*         GetBasalPenta(void);
		int            GetElementType(void);
		void           GetGroundedPart(int* point1, IssmDouble* fraction1, IssmDouble* fraction2, bool* mainlyfloating, int distance_enum, IssmDouble intrusion_distance);
		IssmDouble     GetGroundedPortion(IssmDouble* xyz_list);
		void           GetFractionGeometry(IssmDouble* weights, IssmDouble* pphi, int* ppoint1,IssmDouble* pfraction1,IssmDouble* pfraction2, bool* ptrapezeisnegative, IssmDouble* gl){_error_("not implemented yet");};
		void           GetFractionGeometry2D(IssmDouble* weights, IssmDouble* pphi, int* ppoint1,IssmDouble* pfraction1,IssmDouble* pfraction2, bool* ptrapezeisnegative, IssmDouble* gl);
		void       GetNodalWeightsAndAreaAndCentroidsFromLeveset(IssmDouble* loadweights, IssmDouble* ploadarea, IssmDouble* platbar, IssmDouble* plongbar, IssmDouble late, IssmDouble longe, IssmDouble area,  int levelsetenum){_error_("not implemented yet");};
		void       GetNodalWeightsAndAreaAndCentroidsFromLeveset(IssmDouble* loadweights, IssmDouble* ploadarea, IssmDouble* platbar, IssmDouble* plongbar, IssmDouble late, IssmDouble longe, IssmDouble area, int levelset1enum, int levelset2enum){_error_("not implemented yet");};
		void        GetBarycenterFromLevelset(IssmDouble* platbar, IssmDouble* plongbar,IssmDouble phi,IssmDouble fraction1,IssmDouble fraction2,IssmDouble late, IssmDouble longe, int point1,int istrapeze1, IssmDouble planetradius){_error_("not implemented yet");};
		IssmDouble		GetIcefrontArea();
		void           GetIcefrontCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum);
		Input*        GetInput(int enumtype);
		Input*        GetInput(int enumtype,IssmDouble time);
		Input*        GetInput(int inputenum,IssmDouble start_time,IssmDouble end_time,int averaging_method){_error_("not implemented yet!");};
		void        GetInputListOnVertices(IssmDouble* pvalue,Input* input,IssmDouble default_value);
		void        GetInputListOnNodes(IssmDouble* pvalue,Input* input,IssmDouble default_value);
		DatasetInput* GetDatasetInput(int inputenum);
		void           GetInputValue(IssmDouble* pvalue,Vertex* vertex,int enumtype);
		void           GetInputValue(IssmDouble* pvalue,Node* node,int enumtype);
		void           GetLevelCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum,IssmDouble level){_error_("not implemented yet");};
		void           GetLevelsetPositivePart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlynegative,IssmDouble* levelsetvalues){_error_("not implemented yet");};
		void				GetLevelsetIntersectionBase(int** pindices, int* pnumiceverts, IssmDouble* fraction, int levelset_enum, IssmDouble level);
		int            GetNumberOfNodes(void);
		int            GetNumberOfNodes(int enum_type);
		int            GetNumberOfVertices(void);
		Penta*         GetLowerPenta(void);
		Penta*         GetUpperPenta(void);
		void           GetVectorFromControlInputs(Vector<IssmDouble>* gradient,int control_enum,int control_index,int N,const char* data,int offset);
		int            GetVertexIndex(Vertex* vertex);
		void           GetVerticesCoordinatesBase(IssmDouble** pxyz_list);
		void           GetVerticesCoordinatesTop(IssmDouble** pxyz_list);
		IssmDouble     GroundedArea(bool scaled);
		IssmDouble     GroundinglineMassFlux(bool scaled);
		IssmDouble		IcefrontMassFluxLevelset(bool scaled);
		IssmDouble     IceVolume(bool scaled);
		IssmDouble     IceVolumeAboveFloatation(bool scaled);
		void           InputDepthAverageAtBase(int enum_type,int average_enum_type);
		void	         InputExtrude(int enum_type,int start);
		void           InputUpdateFromIoModel(int index, IoModel* iomodel);
		void           InputUpdateFromSolutionOneDof(IssmDouble* solutiong,int enum_type);
		void           InputUpdateFromSolutionOneDofCollapsed(IssmDouble* solutiong,int enum_type);
		void           InputUpdateFromVector(IssmDouble* vector, int name, int type);
		bool           IsFaceOnBoundary(void){_error_("not implemented yet");};
		bool           IsIcefront(void);
		bool           IsNodeOnShelfFromFlags(IssmDouble* flags);
		bool           IsZeroLevelset(int levelset_enum);
		void           JacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,Gauss* gauss);
		void           JacobianDeterminantBase(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss);
		void           JacobianDeterminantLine(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss);
		void           JacobianDeterminantSurface(IssmDouble*  pJdet, IssmDouble* xyz_list,Gauss* gauss);
		void           JacobianDeterminantTop(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss);
		IssmDouble     Masscon(IssmDouble* levelset){_error_("not implemented yet");};
		IssmDouble     MassFlux(IssmDouble* segment);
		IssmDouble     MassFlux(IssmDouble x1,IssmDouble y1, IssmDouble x2, IssmDouble y2,int segment_id);
		IssmDouble     MinEdgeLength(IssmDouble* xyz_list);
		IssmDouble     Misfit(int modelenum,int observationenum,int weightsenum){_error_("not implemented yet");};
		IssmDouble     MisfitArea(int weightsenum){_error_("not implemented yet");};
		void           MovingFrontalVelocity(void);
		Gauss*         NewGauss(void);
		Gauss*         NewGauss(int order);
		Gauss*         NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order){_error_("not implemented yet");};
		Gauss*         NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz,int order_vert);
		Gauss*         NewGaussBase(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz);
		Gauss*         NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,bool mainlyfloating,int order);
		Gauss*         NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,int order){_error_("not implemented yet");};
		Gauss*         NewGauss(IssmDouble fraction1,IssmDouble fraction2,int order){_error_("not implemented yet");};
		Gauss*         NewGaussBase(int order);
		Gauss*         NewGaussLine(int vertex1,int vertex2,int order);
		Gauss*         NewGaussTop(int order);
		void           NodalFunctions(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsDerivativesVelocity(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsMINIDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsPressure(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsP1(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsP1Derivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsP2(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsTensor(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsVelocity(IssmDouble* basis,Gauss* gauss);
		void	         NormalBase(IssmDouble* bed_normal, IssmDouble* xyz_list);
		void           NormalSection(IssmDouble* normal,IssmDouble* xyz_list);
		void           NormalSectionBase(IssmDouble* normal,IssmDouble* xyz_list);
		void	         NormalTop(IssmDouble* bed_normal, IssmDouble* xyz_list);
		int            NodalValue(IssmDouble* pvalue, int index, int natureofdataenum);
		int            NumberofNodesPressure(void);
		int            NumberofNodesVelocity(void);
		void           PotentialUngrounding(Vector<IssmDouble>* potential_sheet_ungrounding);
		int            PressureInterpolation();
      void           Recover3DMOLHOInput(int targetVel_enum, int numnodes, IssmDouble* vb,  IssmDouble* vsh, IssmDouble* n, IssmDouble* H, IssmDouble* s);
		void           ReduceMatrices(ElementMatrix* Ke,ElementVector* pe);
		void           ResetFSBasalBoundaryCondition(void);
		void           ResetHooks();
		void           SetElementInput(int enum_in,IssmDouble value);
		void           SetElementInput(int enum_in,IssmDouble value,int type);
		void           SetElementInput(Inputs* inputs,int enum_in,IssmDouble value);
		void           SetElementInput(Inputs* inputs,int enum_in,IssmDouble value,int type){_error_("not implemented yet");};
		void           SetElementInput(Inputs* inputs,int numindices,int* indices,IssmDouble* values,int enum_in);
		void           SetControlInputsFromVector(IssmDouble* vector,int control_enum,int control_index,int offset, int M,int N);
		void           SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Materials* materials,Parameters* parameters);
		void           SetTemporaryElementType(int element_type_in);
	   Element*       SpawnBasalElement(bool depthaverage_materials);
		Element*       SpawnTopElement(void);
		bool           IsSpawnedElement(void);
		Tria*	         SpawnTria(int index1,int index2,int index3);
		IssmDouble     StabilizationParameter(IssmDouble u, IssmDouble v, IssmDouble w, IssmDouble diameter, IssmDouble kappa);
		void           StabilizationParameterAnisotropic(IssmDouble* tau_parameter_anisotropic, IssmDouble u, IssmDouble v, IssmDouble w, IssmDouble hx, IssmDouble hy, IssmDouble hz, IssmDouble kappa);
		void           StressIntensityFactor();
		void           StrainRateparallel();
		void           StrainRateperpendicular();
		IssmDouble     SurfaceArea(void);
		void       	   TangentBase(IssmDouble* bed_tangent,IssmDouble* bed_normal);
		int            TensorInterpolation(){_error_("not implemented yet");};
		IssmDouble     TimeAdapt();
		IssmDouble		TotalCalvingFluxLevelset(bool scaled);
		IssmDouble		TotalCalvingMeltingFluxLevelset(bool scaled);
		IssmDouble     TotalFloatingBmb(bool scaled);
		IssmDouble     TotalGroundedBmb(bool scaled);
		IssmDouble     TotalSmb(bool scaled);
		IssmDouble     TotalSmbMelt(bool scaled);
		IssmDouble     TotalSmbRefreeze(bool scaled);
		void           Update(Inputs* inputs,int index, IoModel* iomodel,int analysis_counter,int analysis_type,int finitelement);
		void           UpdateConstraintsExtrudeFromBase(void);
		void           UpdateConstraintsExtrudeFromTop(void);
		int            UpdatePotentialUngrounding(IssmDouble* potential_sheet_ungrounding,Vector<IssmDouble>* vec_nodes_on_iceshelf,IssmDouble* nodes_on_iceshelf);
		void           ValueP1DerivativesOnGauss(IssmDouble* dvalue,IssmDouble* values,IssmDouble* xyz_list,Gauss* gauss);
		void           ValueP1OnGauss(IssmDouble* pvalue,IssmDouble* values,Gauss* gauss);
		int            VelocityInterpolation();
		int            VertexConnectivity(int vertexindex);
		void           VerticalSegmentIndices(int** pindices,int* pnumseg);
		void           VerticalSegmentIndicesBase(int** pindices,int* pnumseg);
		void           ViscousHeating(IssmDouble* pphi,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input);

		#ifdef _HAVE_DAKOTA_
		void  InputUpdateFromMatrixDakota(IssmDouble* matrix, int nows, int ncols, int name, int type);
		void  InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type);
		void  InputScaleFromDakota(IssmDouble* distributed_values, IssmDouble* partition, int npart, int nt, int name);

		#endif

		#ifdef _HAVE_ESA_
		void    EsaGeodetic2D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,Vector<IssmDouble>* pGravity,Vector<IssmDouble>* pX,Vector<IssmDouble>* pY,IssmDouble* xx,IssmDouble* yy){_error_("not implemented yet!");};
      void    EsaGeodetic3D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,Vector<IssmDouble>* pGravity,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz){_error_("not implemented yet!");};
		#endif
		#ifdef _HAVE_SEALEVELCHANGE_
		void       GiaDeflection(Vector<IssmDouble>* wg,Vector<IssmDouble>* dwgdt,Matlitho* litho, IssmDouble* x,IssmDouble* y){_error_("not implemented yet");};
		void       SealevelchangeGeometrySubElementKernel(SealevelGeometry* slgeom){_error_("not implemented yet");};
		void       SealevelchangeShift(GrdLoads* loads,  IssmDouble offset, SealevelGeometry* slgeom){_error_("not implemented yet");};
		void       SealevelchangeGeometryInitial(IssmDouble* xxe, IssmDouble* yye, IssmDouble* zze, IssmDouble* areae, int* lids,int* vcount){_error_("not implemented yet");};
		void       SealevelchangeGeometryCentroidLoads(SealevelGeometry* slgeom, IssmDouble* xxe, IssmDouble* yye, IssmDouble* zze, IssmDouble* areae){_error_("not implemented yet");};
		void       SealevelchangeGeometrySubElementLoads(SealevelGeometry* slgeom, IssmDouble* areae){_error_("not implemented yet");};
		void       SealevelchangeBarystaticLoads(GrdLoads* loads,  BarystaticContributions* barycontrib, SealevelGeometry* slgeom){_error_("not implemented yet");};
		void       SealevelchangeConvolution(IssmDouble* sealevelpercpu, GrdLoads* loads, IssmDouble* rotationvector,SealevelGeometry* slgeom){_error_("not implemented yet");};
		void       SealevelchangeOceanAverage(GrdLoads* loads, Vector<IssmDouble>* oceanareas, Vector<IssmDouble>* subelementoceanareas, IssmDouble* sealevelpercpu, SealevelGeometry* slgeom){_error_("not implemented yet");};
		void       SealevelchangeDeformationConvolution(IssmDouble* sealevelpercpu, GrdLoads* loads,  IssmDouble* rotationvector,SealevelGeometry* slgeom){_error_("not implemented yet");};
		void       SealevelchangeUpdateViscousFields(IssmDouble lincoeff, int newindex, int offset){_error_("not implemented yet");};
		#endif

		/*}}}*/
};
#endif  /* _PENTA_H */
