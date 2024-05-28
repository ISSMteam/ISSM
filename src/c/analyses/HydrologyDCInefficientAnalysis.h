/*! \file HydrologyDCInefficientAnalysis.h
 *  \brief: header file for generic external result object
 */

#ifndef _HydrologyDCInefficientAnalysis_
#define _HydrologyDCInefficientAnalysis_

/*Headers*/
#include "./Analysis.h"
class Node;
class Input;
class HydrologyDCInefficientAnalysis: public Analysis{

	public:
		/*Model processing*/
		int  DofsPerNode(int** doflist,int domaintype,int approximation);
		void UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum);
		void UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type);
		void CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr=false);
		void CreateConstraints(Constraints* constraints,IoModel* iomodel);
		void CreateLoads(Loads* loads, IoModel* iomodel);

		/*Finite element Analysis*/
		void           Core(FemModel* femmodel);
		void           PreCore(FemModel* femmodel);
		ElementVector* CreateDVector(Element* element);
		ElementMatrix* CreateJacobianMatrix(Element* element);
		ElementMatrix* CreateKMatrix(Element* element);
		ElementVector* CreatePVector(Element* element);
		void GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element);
		void GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_interp,int control_index);
		void InputUpdateFromSolution(IssmDouble* solution,Element* element);
		void UpdateConstraints(FemModel* femmodel);

		/*Intermediaries*/
		IssmDouble SedimentStoring(Element* element, Gauss* gauss, Input* sed_head_input, Input* base_input);
		IssmDouble SedimentTransmitivity(Element* element,Gauss* gauss,Input* sed_head_input, Input* base_input,Input* SedTrans_input);
		void GetHydrologyDCInefficientHmax(IssmDouble* ph_max,Element* element, Node* innode);
		IssmDouble GetHydrologyKMatrixTransfer(Element* element);
		IssmDouble GetHydrologyPVectorTransfer(Element* element, Gauss* gauss, Input* epl_head_input);
		void ElementizeEplMask(FemModel* femmodel);
		void HydrologyIDSGetMask(Vector<IssmDouble>* vec_mask, Element* element);
		void ElementizeIdsMask(FemModel* femmodel);
		void HydrologyIdsGetActive(Vector<IssmDouble>* active_vec, Element* element);
};
#endif
