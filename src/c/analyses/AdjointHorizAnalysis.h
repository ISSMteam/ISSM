/*! \file AdjointHorizAnalysis.h 
 *  \brief: header file for generic external result object
 */

#ifndef _AdjointHorizAnalysis_
#define _AdjointHorizAnalysis_

/*Headers*/
#include "./Analysis.h"

class AdjointHorizAnalysis: public Analysis{

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
		ElementMatrix* CreateKMatrixFS(Element* element);
		ElementMatrix* CreateKMatrixHO(Element* element);
		ElementMatrix* CreateKMatrixMOLHO(Element* element);
		ElementMatrix* CreateKMatrixMOLHOVerticalIntergrated(Element* element);
		ElementMatrix* CreateKMatrixL1L2(Element* element);
		ElementMatrix* CreateKMatrixSSA(Element* element);
		ElementVector* CreatePVector(Element* element);
		ElementVector* CreatePVectorFS(Element* element);
		ElementVector* CreatePVectorL1L2(Element* element);
		ElementVector* CreatePVectorHO(Element* element);
		ElementVector* CreatePVectorMOLHO(Element* element);
		ElementVector* CreatePVectorSSA(Element* element);
		void           GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element);
		void           GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index);
		void           GradientJBbarFS(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJBbarGradient(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJBinitial(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJBbarL1L2(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJBbarHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJBbarMOLHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJBbarSSA(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJBFS(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJBGradient(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJBHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJBMOLHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJBSSA(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJDragFS(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJDragGradient(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJDragL1L2(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJDragHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJDragMOLHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJDragSSA(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJDragHydroFS(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJDragHydroL1L2(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJDragHydroHO(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJDragHydroSSA(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           GradientJDSSA(Element* element,Vector<IssmDouble>* gradient,int control_interp,int control_index);
		void           InputUpdateFromSolution(IssmDouble* solution,Element* element);
		void           InputUpdateFromSolutionFS(IssmDouble* solution,Element* element);
		void           InputUpdateFromSolutionHoriz(IssmDouble* solution,Element* element);
		void           InputUpdateFromSolutionMOLHO(IssmDouble* solution,Element* element);
		void           UpdateConstraints(FemModel* femmodel);
};
#endif
