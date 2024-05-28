#include "./BalancethicknessSoftAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
int  BalancethicknessSoftAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void BalancethicknessSoftAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void BalancethicknessSoftAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void BalancethicknessSoftAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void BalancethicknessSoftAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void BalancethicknessSoftAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/

/*Finite Element Analysis*/
void           BalancethicknessSoftAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           BalancethicknessSoftAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* BalancethicknessSoftAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* BalancethicknessSoftAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* BalancethicknessSoftAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
ElementVector* BalancethicknessSoftAnalysis::CreatePVector(Element* element){/*{{{*/
_error_("not implemented yet");
}/*}}}*/
void BalancethicknessSoftAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void BalancethicknessSoftAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void BalancethicknessSoftAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void BalancethicknessSoftAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
