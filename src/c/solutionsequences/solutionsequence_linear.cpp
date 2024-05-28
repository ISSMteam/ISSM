/*!\file: solutionsequence_linear.cpp
 * \brief: numerical core of linear solutions
 */ 

#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

void solutionsequence_linear(FemModel* femmodel){

	/*intermediary: */
	Matrix<IssmDouble>*  Kff = NULL;
	Matrix<IssmDouble>*  Kfs = NULL;
	Vector<IssmDouble>*  ug  = NULL;
	Vector<IssmDouble>*  uf  = NULL;
	Vector<IssmDouble>*  pf  = NULL;
	Vector<IssmDouble>*  df  = NULL;
	Vector<IssmDouble>*  ys  = NULL;

	/*Recover parameters: */
	femmodel->UpdateConstraintsx();
	SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel);
	CreateNodalConstraintsx(&ys,femmodel->nodes);
	Reduceloadx(pf, Kfs, ys); delete Kfs;

	femmodel->profiler->Start(SOLVER);
	Solverx(&uf, Kff, pf, NULL, df, femmodel->parameters); 
	femmodel->profiler->Stop(SOLVER);

	/*Clean up*/
	delete Kff; delete pf; delete df;

	//#ifdef  _HAVE_ADOLC_
	//        for (int i =0; i<uf->svector->M; ++i) {
	//          ADOLC_DUMP_MACRO(uf->svector->vector[i]);
	//        }
	//#endif

	Mergesolutionfromftogx(&ug, uf,ys,femmodel->nodes,femmodel->parameters);delete uf; delete ys;
	InputUpdateFromSolutionx(femmodel,ug); 
	delete ug;  
}
