/*! \file Tria.h
 *  \brief: header file for tria object
 */

#ifndef _TRIA_H_
#define _TRIA_H_

/*Headers:*/
/*{{{*/
#include "./Element.h"
#include "./ElementHook.h"
#include "./TriaRef.h"
class Parameters;
class Inputs;
class IoModel;
class Results;
class Node;
class Material;
class Matlitho;
class Seg;
class ElementMatrix;
class ElementVector;
class Vertex;
class GaussTria;

#include "../../shared/Exceptions/exceptions.h"
#include "../../shared/Enum/Enum.h"
/*}}}*/

class Tria: public Element,public ElementHook,public TriaRef{

	public:
		int iscollapsed;

		/*Tria constructors, destructors {{{*/
		Tria(){};
		Tria(int tria_id,int tria_sid,int tria_lid,IoModel* iomodel,int nummodels);
		~Tria();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object *copy();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*Update virtual functions resolution: {{{*/
		#ifdef _HAVE_DAKOTA_
		void  InputUpdateFromMatrixDakota(IssmDouble* matrix, int nows, int ncols, int name, int type);
		void  InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type);
		void  InputScaleFromDakota(IssmDouble* distributed_values, IssmDouble* partition, int npart, int nt, int name);
		#endif
		void  InputUpdateFromIoModel(int index, IoModel* iomodel);
		void  InputUpdateFromVector(IssmDouble* vector, int name, int type);
		/*}}}*/
		/*Element virtual functions definitions: {{{*/
		void        AverageOntoPartition(Vector<IssmDouble>* partition_contributions,Vector<IssmDouble>* partition_areas,IssmDouble* vertex_response,IssmDouble* qmu_part);
		void			CalvingRateVonmises();
		void			CalvingRateVonmisesAD();
		void			CalvingRateTest();
		void        CalvingCrevasseDepth();
		void			CalvingRateLevermann();
		void			CalvingPollard();
		void			CalvingFluxLevelset();
		void			CalvingMeltingFluxLevelset();
		void			CalvingRateParameterization();
		void			CalvingRateCalvingMIP();
		IssmDouble  CharacteristicLength(void);
		void        ComputeBasalStress(void);
		void        ComputeDeviatoricStressTensor();
		void        ComputeEsaStrainAndVorticity();
		void        ComputeSigmaNN();
		void        ComputeSigmaVM();
		void        ComputeStressTensor();
		void        ComputeSurfaceNormalVelocity();
		void        Configure(Elements* elements,Loads* loads,Nodes* nodesin,Vertices* verticesin,Materials* materials,Parameters* parameters,Inputs* inputsin);
		void        ControlInputSetGradient(IssmDouble* gradient,int enum_type,int control_index,int offset,int M,int N,int interp);
		void        ControlToVectors(Vector<IssmPDouble>* vector_control, Vector<IssmPDouble>* vector_gradient,int control_enum,int control_interp);
		void        CreateDistanceInputFromSegmentlist(IssmDouble* distances,int distanceenum);
		void        ElementCoordinates(Vector<IssmDouble>* vxe,Vector<IssmDouble>* vye,Vector<IssmDouble>* vze, Vector<IssmDouble>* vareae, bool spherical=false);
		void        ElementCoordinates(Vector<IssmDouble>* vlonge,Vector<IssmDouble>* vlate,Vector<IssmDouble>* vareae);
		int         EdgeOnBaseIndex();
		void        EdgeOnBaseIndices(int* pindex1,int* pindex);
		int         EdgeOnSurfaceIndex();
		void        EdgeOnSurfaceIndices(int* pindex1,int* pindex);
		void        ElementResponse(IssmDouble* presponse,int response_enum);
		void        ElementSizes(IssmDouble* hx,IssmDouble* hy,IssmDouble* hz);
		int         FiniteElement(void);
		IssmDouble  FloatingArea(bool scaled);
		void        FSContactMigration(Vector<IssmDouble>* vertex_sigmann,Vector<IssmDouble>* vertex_waterpressure);
		Element*    GetBasalElement(void){_error_("not implemented yet");};
		void        GetLevelsetPositivePart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlynegative,IssmDouble* levelsetvalues);
		void        GetGroundedPart(int* point1,IssmDouble* fraction1, IssmDouble* fraction2,bool* mainlyfloating, int distance_enum, IssmDouble intrusion_distance);
		IssmDouble  GetGroundedPortion(IssmDouble* xyz_list);
		void        GetFractionGeometry(IssmDouble* weights, IssmDouble* pphi, int* ppoint1,IssmDouble* pfraction1,IssmDouble* pfraction2, bool* ptrapezeisnegative, IssmDouble* gl);
		void        GetNodalWeightsAndAreaAndCentroidsFromLeveset(IssmDouble* loadweights, IssmDouble* ploadarea, IssmDouble* platbar, IssmDouble* plongbar, IssmDouble late, IssmDouble longe, IssmDouble area,  int levelsetenum);
		void        GetNodalWeightsAndAreaAndCentroidsFromLeveset(IssmDouble* loadweights, IssmDouble* ploadarea, IssmDouble* platbar, IssmDouble* plongbar, IssmDouble late, IssmDouble longe, IssmDouble area, int levelset1enum, int levelset2enum);
		void        GetBarycenterFromLevelset(IssmDouble* platbar, IssmDouble* plongbar,IssmDouble phi,IssmDouble fraction1,IssmDouble fraction2,IssmDouble late, IssmDouble longe, int point1,int istrapeze1, IssmDouble planetradius);
		IssmDouble  GetTriangleAreaSpherical(IssmDouble xyz_list[3][3]);
		IssmDouble  GetIcefrontArea();
		void	      GetIcefrontCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum);
		void        GetInputListOnVertices(IssmDouble* pvalue,Input* input,IssmDouble default_value);
		void        GetInputListOnNodes(IssmDouble* pvalue,Input* input,IssmDouble default_value);
		void	      GetLevelCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum,IssmDouble level);
		int         GetVertexIndex(Vertex* vertex);
		int         GetNumberOfNodes(void);
		int         GetNumberOfNodes(int enum_type);
		int         GetNumberOfVertices(void);
		void        GetVectorFromControlInputs(Vector<IssmDouble>* gradient,int control_enum,int control_index,int N,const char* data,int offset);
		void        GetVerticesCoordinatesBase(IssmDouble** pxyz_list);
		void        GetVerticesCoordinatesTop(IssmDouble** pxyz_list);
		IssmDouble  GroundedArea(bool scaled);
		bool        HasEdgeOnBase();
		bool        HasEdgeOnSurface();
		IssmDouble  IceVolume(bool scaled);
		IssmDouble  IceVolumeAboveFloatation(bool scaled);
		IssmDouble  IcefrontMassFlux(bool scaled);
		IssmDouble  IcefrontMassFluxLevelset(bool scaled);
		IssmDouble  GroundinglineMassFlux(bool scaled);
		void        InputDepthAverageAtBase(int enum_type,int average_enum_type);
		void        InputExtrude(int enum_type,int start){_error_("not implemented"); /*For penta only*/};
		void        ControlInputExtrude(int enum_type,int start){/*For penta only*/};
		bool	   	IsFaceOnBoundary(void);
		bool	   	IsIcefront(void);
		bool        IsNodeOnShelfFromFlags(IssmDouble* flags);
		bool        IsZeroLevelset(int levelset_enum);
		IssmDouble  Masscon(IssmDouble* levelset);
		IssmDouble  MassFlux(IssmDouble* segment);
		IssmDouble  MassFlux(IssmDouble x1,IssmDouble y1, IssmDouble x2, IssmDouble y2,int segment_id);
		void        MaterialUpdateFromTemperature(void){_error_("not implemented yet");};
		IssmDouble  Misfit(int modelenum,int observationenum,int weightsenum);
		IssmDouble  MisfitArea(int weightsenum);
		int         NodalValue(IssmDouble* pvalue, int index, int natureofdataenum);
		int         NumberofNodesPressure(void);
		int         NumberofNodesVelocity(void);
		void        PotentialUngrounding(Vector<IssmDouble>* potential_sheet_ungrounding);
		int         PressureInterpolation();
		void        ReduceMatrices(ElementMatrix* Ke,ElementVector* pe);
		void        ResetFSBasalBoundaryCondition(void);
		void        ResetHooks();
		void        ResetLevelsetFromSegmentlist(IssmDouble* segments,int numsegments);
		void        SetElementInput(int enum_in,IssmDouble value);
		void        SetElementInput(int enum_in,IssmDouble value,int type);
		void        SetElementInput(Inputs* inputs,int enum_in,IssmDouble value);
		void        SetElementInput(Inputs* inputs,int numindices,int* indices,IssmDouble* values,int enum_in);
		void        SetControlInputsFromVector(IssmDouble* vector,int control_enum,int control_index,int offset,int M,int N);
		void        SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Materials* materials,Parameters* parameters);
	   Element*    SpawnBasalElement(bool depthaverage_materials);
		Element*    SpawnTopElement(void);
		bool        IsSpawnedElement(void);
		void			StrainRateparallel();
		void			StrainRateperpendicular();
		void        StressIntensityFactor(void){_error_("not implemented yet");};
		IssmDouble  SurfaceArea(void);
		int         TensorInterpolation();
		IssmDouble  TimeAdapt();
		IssmDouble  TotalCalvingFluxLevelset(bool scaled);
		IssmDouble  TotalCalvingMeltingFluxLevelset(bool scaled);
		IssmDouble  TotalFloatingBmb(bool scaled);
		IssmDouble  TotalGroundedBmb(bool scaled);
		IssmDouble  TotalSmb(bool scaled);
		IssmDouble  TotalSmbMelt(bool scaled);
		IssmDouble  TotalSmbRefreeze(bool scaled);
		void        Update(Inputs* inputs,int index, IoModel* iomodel,int analysis_counter,int analysis_type,int finitelement);
		int         UpdatePotentialUngrounding(IssmDouble* vertices_potentially_ungrounding,Vector<IssmDouble>* vec_nodes_on_iceshelf,IssmDouble* nodes_on_iceshelf);
		void        ValueP1DerivativesOnGauss(IssmDouble* dvalue,IssmDouble* values,IssmDouble* xyz_list,Gauss* gauss);
		void        ValueP1OnGauss(IssmDouble* pvalue,IssmDouble* values,Gauss* gauss);
		int         VelocityInterpolation();
		int         VertexConnectivity(int vertexindex);
		void        VerticalSegmentIndices(int** pindices,int* pnumseg){_error_("not implemented yet");};
		void        VerticalSegmentIndicesBase(int** pindices,int* pnumseg){_error_("not implemented yet");};
		void			WriteFieldIsovalueSegment(DataSet* segments,int fieldenum,IssmDouble fieldvalue);

		#ifdef _HAVE_ESA_
		void    EsaGeodetic2D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,Vector<IssmDouble>* pGravity,Vector<IssmDouble>* pX,Vector<IssmDouble>* pY,IssmDouble* xx,IssmDouble* yy);
      void    EsaGeodetic3D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,Vector<IssmDouble>* pGravity,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz);
		#endif
		#ifdef _HAVE_SEALEVELCHANGE_
		void       GiaDeflection(Vector<IssmDouble>* wg,Vector<IssmDouble>* dwgdt,Matlitho* litho, IssmDouble* x,IssmDouble* y);
		void       SealevelchangeGeometryInitial(IssmDouble* xxe, IssmDouble* yye, IssmDouble* zze, IssmDouble* areae, int* lids, int* vcount);
		void       SealevelchangeGeometrySubElementKernel(SealevelGeometry* slgeom);
		void       SealevelchangeGeometrySubElementLoads(SealevelGeometry* slgeom, IssmDouble* areae);
		void       SealevelchangeGeometryCentroidLoads(SealevelGeometry* slgeom, IssmDouble* xxe, IssmDouble* yye, IssmDouble* zze, IssmDouble* areae);
		void       SealevelchangeBarystaticLoads(GrdLoads* loads, BarystaticContributions* barycontrib, SealevelGeometry* slgeom);
		void       SealevelchangeConvolution(IssmDouble* sealevelpercpu, GrdLoads* loads, IssmDouble* rotationvector,SealevelGeometry* slgeom);
		void       SealevelchangeOceanAverage(GrdLoads* loads, Vector<IssmDouble>* oceanareas, Vector<IssmDouble>* subelementoceanareas, IssmDouble* sealevelpercpu, SealevelGeometry* slgeom);
		void       SealevelchangeDeformationConvolution(IssmDouble* sealevelpercpu, GrdLoads* loads, IssmDouble* rotationvector,SealevelGeometry* slgeom);
		void       SealevelchangeShift(GrdLoads* loads,  IssmDouble offset, SealevelGeometry* slgeom);
		void       SealevelchangeUpdateViscousFields(IssmDouble lincoeff, int newindex, int offset);
		#endif
		/*}}}*/
		/*Tria specific routines:{{{*/
		void           AddBasalInput(int input_enum, IssmDouble* values, int interpolation_enum);
		void           AddInput(int input_enum, IssmDouble* values, int interpolation_enum);
		void           AddControlInput(int input_enum,Inputs* inputs,IoModel* iomodel,IssmDouble* values,IssmDouble* values_min,IssmDouble* values_max, int interpolation_enum,int id);
		void           DatasetInputCreate(IssmDouble* array,int M,int N,int* individual_enums,int num_inputs,Inputs* inputs,IoModel* iomodel,int input_enum);
		void           CreateInputTimeAverage(int transientinput_enum,int averagedinput_enum,IssmDouble init_time,IssmDouble end_time,int averaging_method);
		void           GetInputAveragesUpToCurrentTime(int input_enum,IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime);
		IssmDouble     GetArea(void);
		IssmDouble     GetHorizontalSurfaceArea(void);
		IssmDouble 	   GetArea3D(void);
		IssmDouble 	   GetAreaIce(void);
		IssmDouble 	   GetAreaSpherical(void);
		void           GetAreaCoordinates(IssmDouble *area_coordinates,IssmDouble* xyz_zero,IssmDouble* xyz_list,int numpoints);
		int            GetElementType(void);
		Input*        GetInput(int enumtype);
		Input*        GetInput(int enumtype,IssmDouble time);
		Input*        GetInput(int inputenum,IssmDouble start_time,IssmDouble end_time,int averaging_method);
		DatasetInput* GetDatasetInput(int inputenum);
		void           GetInputValue(IssmDouble* pvalue,Node* node,int enumtype);
		void           GetInputValue(IssmDouble* pvalue,Vertex* vertex,int enumtype);
		void		      GetLevelsetIntersection(int** pindices, int* pnumiceverts, IssmDouble* fraction, int levelset_enum, IssmDouble level);
		void           GetMaterialInputValue(IssmDouble* pvalue,Node* node,int enumtype);
		void	         InputUpdateFromSolutionOneDof(IssmDouble* solution,int enum_type);
		void	         InputUpdateFromSolutionOneDofCollapsed(IssmDouble* solution,int enum_type){_error_("not implemented yet");};
		void           JacobianDeterminant(IssmDouble*  pJdet, IssmDouble* xyz_list,Gauss* gauss);
		void           JacobianDeterminantBase(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss);
		void           JacobianDeterminantLine(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss){_error_("not implemented yet");};
		void           JacobianDeterminantSurface(IssmDouble*  pJdet, IssmDouble* xyz_list,Gauss* gauss);
		void           JacobianDeterminantTop(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss);
		IssmDouble     MinEdgeLength(IssmDouble* xyz_list){_error_("not implemented yet");};
		void	       MovingFrontalVelocity(void);
		Gauss*         NewGauss(void);
		Gauss*         NewGauss(int order);
        Gauss*         NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order);
        Gauss*         NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,bool mainlyfloating,int order);
        Gauss*         NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,int order);
        Gauss*         NewGauss(IssmDouble fraction1,IssmDouble fraction2,int order);
        Gauss*         NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz,int order_vert);
		Gauss*         NewGaussBase(int order);
		Gauss*         NewGaussLine(int vertex1,int vertex2,int order){_error_("not implemented yet");};
		Gauss*         NewGaussTop(int order);
		void           NodalFunctions(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsDerivativesVelocity(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsMINIDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){_error_("not implemented yet");};
		void           NodalFunctionsPressure(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsP1(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsP1Derivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss);
		void           NodalFunctionsP2(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsTensor(IssmDouble* basis,Gauss* gauss);
		void           NodalFunctionsVelocity(IssmDouble* basis,Gauss* gauss);
		void           NormalBase(IssmDouble* normal,IssmDouble* xyz_list);
		void           NormalSection(IssmDouble* normal,IssmDouble* xyz_list);
		void           NormalTop(IssmDouble* normal,IssmDouble* xyz_list);
		void           SetTemporaryElementType(int element_type_in){_error_("not implemented yet");};
		void           InputServe(Input* input_in);
		Seg*	       SpawnSeg(int index1,int index2);
		IssmDouble     StabilizationParameter(IssmDouble u, IssmDouble v, IssmDouble w, IssmDouble diameter, IssmDouble kappa){_error_("not implemented yet");};
		void           StabilizationParameterAnisotropic(IssmDouble* tau_parameter_ansiotropic, IssmDouble u, IssmDouble v, IssmDouble w, IssmDouble hx, IssmDouble hy, IssmDouble hz, IssmDouble kappa){_error_("not implemented yet");};
		void           UpdateConstraintsExtrudeFromBase(void);
		void           UpdateConstraintsExtrudeFromTop(void);
		IssmDouble*    SealevelchangeGxL(IssmDouble* G, IssmDouble* Grot, GrdLoads* loads, IssmDouble* polarmotionvector, SealevelGeometry* slgeom, int nel, bool computefuture);
		IssmDouble*    SealevelchangeHorizGxL(int spatial_component, IssmDouble* G, IssmDouble* Grot, GrdLoads* loads, IssmDouble* polarmotionvector, SealevelGeometry* slgeom, int nel, bool computefuture);
		void	       SealevelchangeCollectGrdfield(IssmDouble* grdfieldout, IssmDouble* grdfield, SealevelGeometry* slgeom, int nel, bool percpu, int viscousenum, bool computefuture);
		/*}}}*/

};
#endif  /* _TRIA_H */
