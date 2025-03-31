/*! \file StressbalanceAnalysis.h 
 *  \brief: header file for generic external result object
 */

#ifndef _StressbalanceAnalysis_
#define _StressbalanceAnalysis_

/*Headers*/
#include "./Analysis.h"

class StressbalanceAnalysis: public Analysis{

  public:
		/*Model processing*/
		void CreateConstraints(Constraints* constraints,IoModel* iomodel);
		void CreateLoads(Loads* loads, IoModel* iomodel);
		void CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr=false);
		int  DofsPerNode(int** doflist,int domaintype,int approximation);
		void UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type);
		void UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum);

		/*Finite element Analysis*/
		void           Core(FemModel* femmodel);
		void           PreCore(FemModel* femmodel);
		ElementVector* CreateDVector(Element* element);
		ElementMatrix* CreateJacobianMatrix(Element* element);
		ElementMatrix* CreateKMatrix(Element* element);
		ElementVector* CreatePVector(Element* element);
		void           GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element);
		void           GetSolutionFromInputsHoriz(Vector<IssmDouble>* solution,Element* element);
		void           GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index);
		void           InputUpdateFromSolution(IssmDouble* solution,Element* element);
		void           UpdateConstraints(FemModel* femmodel);

		/*SSA*/
		ElementMatrix* CreateJacobianMatrixSSA(Element* element);
		ElementMatrix* CreateKMatrixSSA(Element* element);
		ElementMatrix* CreateKMatrixSSAFriction(Element* element);
		ElementMatrix* CreateKMatrixSSALateralFriction(Element* element);
		ElementMatrix* CreateKMatrixSSAViscous(Element* element);
		ElementVector* CreatePVectorSSA(Element* element);
		ElementVector* CreatePVectorSSAFront(Element* element);
		ElementVector* CreatePVectorSSADrivingStress(Element* element);
		void           GetBSSA(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss);
		void           GetBSSAFriction(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss);
		void           GetBSSAprime(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss);
		void           InputUpdateFromSolutionSSA(IssmDouble* solution,Element* element);
		void				ComputeSurfaceSlope(IssmDouble* slope,Gauss* gauss_DG,Gauss* gauss_CG,int point1,IssmDouble fraction1,IssmDouble fraction2,int ig,int dim,Element* element);
		void				ComputeHydrologySlope(IssmDouble* hydrologyslope,Gauss* gauss_DG,Gauss* gauss_CG,int point1,IssmDouble fraction1,IssmDouble fraction2,int ig,int dim,Element* element);
		void				NodalFunctionsDerivativesRGB(IssmDouble* dbasis_subelem,Gauss* gauss_DG,Gauss* gauss_CG,int point1,IssmDouble fraction1,IssmDouble fraction2,int ig,int dim,Element* element);

		/*L1L2*/
		ElementMatrix* CreateKMatrixL1L2(Element* element);
		ElementMatrix* CreateKMatrixL1L2Friction(Element* element);
		ElementMatrix* CreateKMatrixL1L2Viscous(Element* element);
		ElementVector* CreatePVectorL1L2(Element* element);
		ElementVector* CreatePVectorL1L2Front(Element* element);
		ElementVector* CreatePVectorL1L2DrivingStress(Element* element);
		void           InputUpdateFromSolutionL1L2(IssmDouble* solution,Element* element);
		/*MOLHO*/
		ElementMatrix* CreateKMatrixMOLHO(Element* element);
		ElementMatrix* CreateKMatrixMOLHOFriction(Element* element);
		ElementMatrix* CreateKMatrixMOLHOViscous(Element* element);
		ElementVector* CreatePVectorMOLHO(Element* element);
		ElementVector* CreatePVectorMOLHOFront(Element* element);
		ElementVector* CreatePVectorMOLHODrivingStress(Element* element);
		void           InputUpdateFromSolutionMOLHO(IssmDouble* solution,Element* element);
		void           GetSolutionFromInputsMOLHO(Vector<IssmDouble>* solution,Element* element);
		/*HO*/
		ElementMatrix* CreateJacobianMatrixHO(Element* element);
		ElementMatrix* CreateKMatrixHO(Element* element);
		ElementMatrix* CreateKMatrixHOFriction(Element* element);
		ElementMatrix* CreateKMatrixHOViscous(Element* element);
		ElementVector* CreatePVectorHO(Element* element);
		ElementVector* CreatePVectorHOFront(Element* element);
		ElementVector* CreatePVectorHODrivingStress(Element* element);
		void           GetBHOFriction(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss);
		void           InputUpdateFromSolutionHO(IssmDouble* solution,Element* element);
		/*FS*/
		ElementVector* CreateDVectorFS(Element* element);
		ElementMatrix* CreateJacobianMatrixFS(Element* element);
		ElementMatrix* CreateKMatrixFS(Element* element);
		ElementMatrix* CreateKMatrixFSFriction(Element* element);
		ElementMatrix* CreateKMatrixFSFrictionNitsche(Element* element);
		ElementMatrix* CreateKMatrixFSShelf(Element* element);
		ElementMatrix* CreateKMatrixFSViscous(Element* element);
		ElementMatrix* CreateKMatrixFSViscousGLS(Element* element);
		ElementMatrix* CreateKMatrixFSViscousLA(Element* element);
		ElementMatrix* CreateKMatrixFSViscousXTH(Element* element);
		ElementMatrix* CreatePressureMassMatrix(Element* element);
		ElementMatrix* CreateSchurPrecondMatrix(Element* element);
		ElementVector* CreatePVectorFS(Element* element);
		ElementVector* CreatePVectorFSFriction(Element* element);
		ElementVector* CreatePVectorFSFront(Element* element);
		ElementVector* CreatePVectorFSShelf(Element* element);
		ElementVector* CreatePVectorFSStress(Element* element);
		ElementVector* CreatePVectorFSViscous(Element* element);
		ElementVector* CreatePVectorFSViscousGLS(Element* element);
		ElementVector* CreatePVectorFSViscousLA(Element* element);
		ElementVector* CreatePVectorFSViscousXTH(Element* element);
		void           GetBFS(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss);
		void           GetBFSFriction(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss);
		void           GetBFSprime(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss);
		void           GetBFSprimeUzawa(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss);
		void           GetBFSprimevel(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss);
		void           GetBFSUzawa(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss);
		void           GetBFSvel(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss);
		void           GetCFS(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss);
		void           GetCFSprime(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss);
		void           GetSolutionFromInputsFS(Vector<IssmDouble>* solution,Element* element);
		void           InitializeXTH(Elements* elements,Parameters* parameters);
		void           InputUpdateFromSolutionFS(IssmDouble* solution,Element* element);
		void           InputUpdateFromSolutionFSXTH_d(Elements* elements,Parameters* parameters);
		void           InputUpdateFromSolutionFSXTH_tau(Elements* elements,Parameters* parameters);

		/*Coupling*/
		ElementMatrix* CreateKMatrixCouplingHOFS(Element* element);
		ElementMatrix* CreateKMatrixCouplingSSAFS(Element* element);
		ElementMatrix* CreateKMatrixCouplingSSAFSFriction(Element* element);
		ElementMatrix* CreateKMatrixCouplingSSAFSViscous(Element* element);
		ElementMatrix* CreateKMatrixCouplingSSAHO(Element* element);
		ElementMatrix* CreateKMatrixCouplingSSAHOFriction(Element* element);
		ElementMatrix* CreateKMatrixCouplingSSAHOViscous(Element* element);
		ElementMatrix* CreateKMatrixHOFS(Element* element);
		ElementMatrix* CreateKMatrixSSAFS(Element* element);
		ElementMatrix* CreateKMatrixSSAHO(Element* element);
		ElementMatrix* CreateKMatrixSSA3d(Element* element);
		ElementMatrix* CreateKMatrixSSA3dFriction(Element* element);
		ElementMatrix* CreateKMatrixSSA3dViscous(Element* element);
		ElementVector* CreatePVectorSSAFS(Element* element);
		ElementVector* CreatePVectorSSAHO(Element* element);
		ElementVector* CreatePVectorCouplingSSAFS(Element* element);
		ElementVector* CreatePVectorCouplingSSAFSFriction(Element* element);
		ElementVector* CreatePVectorCouplingSSAFSViscous(Element* element);
		ElementVector* CreatePVectorHOFS(Element* element);
		ElementVector* CreatePVectorCouplingHOFS(Element* element);
		ElementVector* CreatePVectorCouplingHOFSFriction(Element* element);
		ElementVector* CreatePVectorCouplingHOFSViscous(Element* element);
		void           GetBprimeSSAFS(IssmDouble* Bprime,Element* element,IssmDouble* xyz_list,Gauss* gauss);
		void           GetBprimeSSAFSTria(IssmDouble* Bprime,Element* element,IssmDouble* xyz_list,Gauss* gauss);
		void           GetBSSAFS(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss);
		void           GetBSSAFSTria(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss);
		void           GetBSSAHO(IssmDouble* B,Element* element,IssmDouble* xyz_list,Gauss* gauss);
		void           GetLFSSSA(IssmDouble* L,Element* element,Gauss* gauss);
		void           GetLprimeFSSSA(IssmDouble* Lprime,Element* element,IssmDouble* xyz_list,Gauss* gauss);
		void           GetLprimeSSAFS(IssmDouble* Lprime,Element* element,IssmDouble* xyz_list,Gauss* gauss);
		void           GetLSSAFS(IssmDouble* L,Element* element,Gauss* gauss);
		void           InputUpdateFromSolutionHOFS(IssmDouble* solution,Element* element);
		void           InputUpdateFromSolutionSSAFS(IssmDouble* solution,Element* element);
		void           InputUpdateFromSolutionSSAHO(IssmDouble* solution,Element* element);
};
#endif
