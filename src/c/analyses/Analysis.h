/*!\file:  Analysis.h
 * \brief abstract class for Analysis objects
 */ 

#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

#include "../toolkits/objects/toolkitobjects.h"

// Looks like AD runs without AMPI are missing commmpi.h
// Conditionally including the header

#if !defined(_HAVE_MPI_) && !defined(_HAVE_PETSC_MPI_)
#include "../toolkits/mpi/commops/commops.h"
#endif

class Parameters;
class Inputs;
class IoModel;
class Elements;
class Nodes;
class Constraints;
class Loads;
class Element;
class ElementVector;
class ElementMatrix;
class Gauss;
class FemModel;

class Analysis{

	public: 
		/*Constructor/Destructor*/
		virtual      ~Analysis(){};

		/*Model processing*/
		virtual void CreateConstraints(Constraints* constraints,IoModel* iomodel)=0;
		virtual void CreateLoads(Loads* loads, IoModel* iomodel)=0;
		virtual void CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr=false)=0;
		virtual int  DofsPerNode(int** doflist,int domaintype,int approximation)=0;
		virtual void UpdateElements(Elements* elements,Inputs* inputs, IoModel* iomodel,int analysis_counter,int analysis_type)=0;
		virtual void UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum)=0;

		/*Finite element Analysis*/
		virtual void           Core(FemModel* femmodel)=0;
		virtual void           PreCore(FemModel* femmodel)=0;
		virtual ElementVector* CreateDVector(Element* element)=0;
		virtual ElementMatrix* CreateJacobianMatrix(Element* element)=0;
		virtual ElementMatrix* CreateKMatrix(Element* element)=0;
		virtual ElementVector* CreatePVector(Element* element)=0;
		virtual void           GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element)=0;
		virtual void           GradientJ(Vector<IssmDouble>* gradient,Element* element,int control_type,int control_interp,int control_index)=0;
		virtual void           InputUpdateFromSolution(IssmDouble* solution,Element* element)=0;
		virtual void           UpdateConstraints(FemModel* femmodel)=0;
};
#endif
