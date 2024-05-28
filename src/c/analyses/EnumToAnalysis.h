#ifndef _ENUMTOANALYSIS_
#define _ENUMTOANALYSIS_

class Analysis;

Analysis* EnumToAnalysis(int analysis_enum);

#endif
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
		void GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index);
		void InputUpdateFromSolution(IssmDouble* solution,Element* element);
		void UpdateConstraints(FemModel* femmodel);
