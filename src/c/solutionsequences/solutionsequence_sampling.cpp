/*!\file: solutionsequence_sampling.cpp
 * \brief: numerical core of linear solutions
 */

#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../analyses/analyses.h"
#include "../shared/Random/randomgenerator.h"
//#include <random>

void GaussianVector(Vector<IssmDouble>* ppf,int seed){/*{{{*/

	/*Define seed*/
	rnd::linear_congruential_engine random_engine;
	if(seed>=0){
		int my_rank;
		ISSM_MPI_Comm_rank(ISSM_MPI_COMM_WORLD,&my_rank);
		seed = seed + 783728*my_rank; // change default seed for parallel simulations (by considering an arbitrary shif based on the rank number)
	}
	random_engine.seed(seed);

	/* Define univariate distribution */
	rnd::normal_distribution distribution(0.0,1.0);

	int        *local_indices = NULL;
	IssmDouble *local_vector = NULL;
	ppf->GetLocalVector(&local_vector,&local_indices);

  int M;
	ppf->GetLocalSize(&M);
	for(int i=0;i<M;i++){
		double rdnumber = distribution.generator(random_engine);
		ppf->SetValue(local_indices[i],rdnumber,INS_VAL);
	}
	ppf->Assemble();

	/*Cleanup*/
	random_engine.free_resources();
	xDelete<int>(local_indices);
	xDelete<IssmDouble>(local_vector);
}/*}}}*/
void solutionsequence_sampling(FemModel* femmodel){

	/*intermediary: */
  Matrix<IssmDouble>*  Kff = NULL;
	Matrix<IssmDouble>*  Kfs = NULL;
  Vector<IssmDouble>*  ug  = NULL;
  Vector<IssmDouble>*  uf  = NULL;
	Vector<IssmDouble>*  pf  = NULL;
  Vector<IssmDouble>*  df  = NULL;
  Vector<IssmDouble>*  ys=NULL;

  SamplingAnalysis*    analysis = NULL;

  Vector<IssmDouble>*  Ml = NULL;      // diagonal lumped mass matrix
  Vector<IssmDouble>*  Mscale = NULL;  // square root of diagonal lumped factor matrix

  /*parameters:*/
  int alpha, seed, nsize;

  /*Recover parameters: */
  femmodel->parameters->FindParam(&seed,SamplingSeedEnum);
  femmodel->parameters->FindParam(&alpha,SamplingAlphaEnum);

  /*CreateAnalysis*/

  analysis = new SamplingAnalysis();

  analysis->LumpedMassMatrix(&Ml,femmodel); //Create lumped mass matrix

  if(alpha%2!=0)  analysis->LumpedKMatrix(&Mscale,femmodel);   /* Compute square root of lump mass matrix (for alpha even) or stiffness matrix (for alpha odd) */
  else{
    Mscale=Ml->Duplicate();
    Ml->Copy(Mscale);
  }
  Mscale->Pow(0.5);

  femmodel->UpdateConstraintsx();
  SystemMatricesx(&Kff, &Kfs, &pf, &df, NULL,femmodel);
  CreateNodalConstraintsx(&ys,femmodel->nodes);
  Reduceloadx(pf, Kfs, ys);
  delete Kfs;

  /* Generate random RHS */
  GaussianVector(pf,seed);

  /* Multiply random RHS by square root of mass matrix (for alpha even) or stiffness matrix (for alpha odd) */
  pf->PointwiseMult(pf,Mscale);

  /*Go solve SPDE */
  femmodel->profiler->Start(SOLVER);
  Solverx(&uf, Kff, pf, NULL, df, femmodel->parameters);
  femmodel->profiler->Stop(SOLVER);

  /* Iterate to compute a realization for alpha>2 */
  for(int i=3;i<=alpha;i+=2){

    /*Create RHS */
    uf->Copy(pf);
    pf->PointwiseMult(pf,Ml);

    /*Go solve SPDE */
    femmodel->profiler->Start(SOLVER);
    Solverx(&uf, Kff, pf, NULL, df, femmodel->parameters);
    femmodel->profiler->Stop(SOLVER);

  }

  /* Update input */
  Mergesolutionfromftogx(&ug, uf,ys,femmodel->nodes,femmodel->parameters);
  InputUpdateFromSolutionx(femmodel,ug);

  /*clean-up*/
  delete Kff; delete pf; delete df; delete uf; delete ys;
  delete Ml; delete Mscale;
  delete analysis;
  delete ug;

}
