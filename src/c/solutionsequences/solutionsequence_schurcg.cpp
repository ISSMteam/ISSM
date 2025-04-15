/*!\file: solutionsequence_schurcg.cpp
 * \brief: numerical core of 
 */ 

#include "./solutionsequences.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../analyses/analyses.h"
#include <iostream>
#include <fstream>

#if defined(_HAVE_PETSC_) && !defined(_HAVE_CODIPACK_)
#include <petscversion.h>
#include "../toolkits/petsc/patches/petscpatches.h"

void SchurCGSolver(Vector<IssmDouble>** puf,PMat Kff, PVec pf, PVec uf0,IS isv,IS isp,Parameters* parameters){/*{{{*/

	PMat                  A, B, BT;				/* Saddle point block matrices */
	PMat						IP;						/* Preconditioner or mass matrix */
	int                  nu, np;					/* No of. free nodes in velocity / pressure space */
  PVec                  p,uold,unew;			/* Solution vectors for pressure / vel. */
	PVec						tmpu,tmpu2,resu,resp,tmpp,tmpp2,rhsu,rhsp; /* temp. vectors, arbitrary RHS in vel. / pressure space */
	PVec						gold,gnew,wold,wnew,chi; /* CG intermediaries */
	PVec						f1,f2;					/* RHS of the global system */
	IssmDouble			rho,gamma,tmpScalar,tmpScalar2; /* Step sizes, arbitrary double */
	PKSP						  kspu,kspip;		/* KSP contexts for vel. / pressure systems*/
	KSPConvergedReason	reason;					/* Convergence reason for troubleshooting */
	int						its;						/* No. of iterations for troubleshooting */
	IssmDouble			initRnorm, rnorm, TOL,ELLTOL; /* residual norms, STOP tolerance */
	PC							pcu,pcp;					/* Preconditioner contexts pertaining the KSP contexts*/
	PetscViewer				viewer;					/* Viewer for troubleshooting */
	IssmDouble				t1,t2;					/* Time measurement for bottleneck analysis */

	IssmDouble tmp1,tmp2,tmp3;
	int tmpi;
	IssmDouble tmp4,tmp5,tmp6,tmp7;

	int noIt;

	int precond = 0;

	#if PETSC_VERSION_LT(3,2,0)
	PetscTruth flag,flg;
	#else
	PetscBool flag,flg;
	#endif

	char ksp_type[50];
	char pc_type[50];
	int maxiter;

	#if PETSC_VERSION_LT(3,7,0)
	PetscOptionsGetString(PETSC_NULL,"-ksp_type",ksp_type,49,&flg);
	PetscOptionsGetString(PETSC_NULL,"-pc_type",pc_type,49,&flg);
	PetscOptionsGetReal(PETSC_NULL,"-tol",&TOL,NULL);
	PetscOptionsGetReal(PETSC_NULL,"-elltol",&ELLTOL,NULL);
	PetscOptionsGetInt(PETSC_NULL,"-schur_pc",&precond,NULL);
	PetscOptionsGetInt(PETSC_NULL,"-max_iter",&maxiter,NULL);
	#elif PETSC_VERSION_LT(3,19,0)
	PetscOptionsGetString(NULL,PETSC_NULL,"-ksp_type",ksp_type,49,&flg);
	PetscOptionsGetString(NULL,PETSC_NULL,"-pc_type",pc_type,49,&flg);
	PetscOptionsGetReal(NULL,PETSC_NULL,"-tol",&TOL,NULL);
	PetscOptionsGetReal(NULL,PETSC_NULL,"-elltol",&ELLTOL,NULL);
	PetscOptionsGetInt(NULL,PETSC_NULL,"-schur_pc",&precond,NULL);
	PetscOptionsGetInt(NULL,PETSC_NULL,"-max_iter",&maxiter,NULL);
	#else
	PetscOptionsGetString(NULL,PETSC_NULLPTR,"-ksp_type",ksp_type,49,&flg);
	PetscOptionsGetString(NULL,PETSC_NULLPTR,"-pc_type",pc_type,49,&flg);
	PetscOptionsGetReal(NULL,PETSC_NULLPTR,"-tol",&TOL,NULL);
	PetscOptionsGetReal(NULL,PETSC_NULLPTR,"-elltol",&ELLTOL,NULL);
	PetscOptionsGetInt(NULL,PETSC_NULLPTR,"-schur_pc",&precond,NULL);
	PetscOptionsGetInt(NULL,PETSC_NULLPTR,"-max_iter",&maxiter,NULL);
	#endif

	if(precond){
		_printf0_("Running WITH preconditioner\n");
	}else{
		_printf0_("Running WITHOUT preconditioner\n");
	}

	/*Initialize output*/
	Vector<IssmDouble>* out_uf=new Vector<IssmDouble>(uf0);

	/* Extract block matrices from the saddle point matrix */
	/* [ A   B ] = Kff
    * [ B^T I ]
	 * where A is the elliptic submatrix, B^T represents the incompressibility, 
	 * and I the Schur preconditioner (stored here, because the space was allocated either way) 
	 *         */
	#if PETSC_VERSION_GT(3,8,0)
	MatCreateSubMatrix(Kff,isv,isv,MAT_INITIAL_MATRIX,&A);
	MatCreateSubMatrix(Kff,isv,isp,MAT_INITIAL_MATRIX,&B);
	MatCreateSubMatrix(Kff,isp,isv,MAT_INITIAL_MATRIX,&BT);
	#else
	MatGetSubMatrix(Kff,isv,isv,MAT_INITIAL_MATRIX,&A);
	MatGetSubMatrix(Kff,isv,isp,MAT_INITIAL_MATRIX,&B);
	MatGetSubMatrix(Kff,isp,isv,MAT_INITIAL_MATRIX,&BT);
	#endif

	/* Extract preconditioner matrix on the pressure space*/
	#if PETSC_VERSION_GT(3,8,0)
	MatCreateSubMatrix(Kff,isp,isp,MAT_INITIAL_MATRIX,&IP);
	#else
	MatGetSubMatrix(Kff,isp,isp,MAT_INITIAL_MATRIX,&IP);
	#endif

	/* Get number of velocity / pressure nodes */
	MatGetSize(B,&nu,&np);

	/* Extract initial guesses for uold and pold */
	VecCreate(IssmComm::GetComm(),&p);VecSetSizes(p,PETSC_DECIDE,np);VecSetFromOptions(p);
	VecAssemblyBegin(p);VecAssemblyEnd(p);
	VecCreate(IssmComm::GetComm(),&uold);VecSetSizes(uold,PETSC_DECIDE,nu);VecSetFromOptions(uold);
	VecAssemblyBegin(uold);VecAssemblyEnd(uold);

	VecGetSubVector(out_uf->pvector->vector,isv,&uold);
	VecGetSubVector(out_uf->pvector->vector,isp,&p);

	/* Set up intermediaries */
	VecDuplicate(uold,&f1);VecSet(f1,0.0);
	VecDuplicate(p,&f2);VecSet(f2,0.0);
	VecDuplicate(uold,&tmpu);VecSet(tmpu,0.0);
	VecDuplicate(uold,&tmpu2);VecSet(tmpu2,0.0);
	VecDuplicate(uold,&resu);VecSet(resu,0.0);
	VecDuplicate(p,&tmpp);VecSet(tmpp,0.0);
	VecDuplicate(p,&tmpp2);VecSet(tmpp2,0.0);
	VecDuplicate(p,&rhsp);VecSet(rhsp,0.0);
	VecDuplicate(p,&resp);VecSet(resp,0.0);
	VecDuplicate(uold,&rhsu);VecSet(rhsu,0.0);
	VecDuplicate(p,&gold);VecSet(gold,0.0);
	VecDuplicate(p,&wnew);VecSet(wnew,0.0);
	VecDuplicate(uold,&chi);VecSet(chi,0.0);

	/* Get global RHS (for each block sub-problem respectively)*/
	VecGetSubVector(pf,isv,&f1);
	VecGetSubVector(pf,isp,&f2);

   /* ------------------------------------------------------------ */

	/* Generate initial value for the velocity from the pressure */
	/* a(u0,v) = f1(v)-b(p0,v)  i.e.  Au0 = F1-Bp0 */
	/* u0 = u_DIR on \Gamma_DIR */

	/* Create KSP context */
	KSPCreate(IssmComm::GetComm(),&kspu);
	#if PETSC_VERSION_GE(3,5,0)
	KSPSetOperators(kspu,A,A);
	#else
	KSPSetOperators(kspu,A,A,DIFFERENT_NONZERO_PATTERN);
	#endif
	if (strcmp(ksp_type,"gmres")==0){
		KSPSetType(kspu,KSPGMRES);
	}else if(strcmp(ksp_type,"pipegmres")==0){
		KSPSetType(kspu,KSPPGMRES);
	}else if(strcmp(ksp_type,"cg")==0){
		KSPSetType(kspu,KSPCG);
	}else if(strcmp(ksp_type,"pipecg")==0){
		KSPSetType(kspu,KSPPIPECG);
	}else if(strcmp(ksp_type,"bicg")==0){
		KSPSetType(kspu,KSPBICG);
	}else if(strcmp(ksp_type,"bicgstab")==0){
		KSPSetType(kspu,KSPBCGS);
	}else if(strcmp(ksp_type,"ibicgstab")==0){
		KSPSetType(kspu,KSPIBCGS);
	}else if(strcmp(ksp_type,"minres")==0){
		KSPSetType(kspu,KSPMINRES);
	}else if(strcmp(ksp_type,"cr")==0){
		KSPSetType(kspu,KSPCR);
	}else if(strcmp(ksp_type,"pipecr")==0){
		KSPSetType(kspu,KSPPIPECR);
	}else{
		_error_("Suggested KSP method not implemented yet!\n");
	}

	KSPSetInitialGuessNonzero(kspu,PETSC_TRUE);

	/*Strong rel. residual tolerance needed when solving for the velocity update. 
	 * This is because ISSM uses the dimensional equations, so the initial guess chi = 0
	 * results in a very high initial residual (~ 1e19)*/
	KSPSetTolerances(kspu,ELLTOL,PETSC_DEFAULT,PETSC_DEFAULT,maxiter); //maxiter

	KSPGetPC(kspu,&pcu);
	if (strcmp(pc_type,"bjacobi")==0){
		PCSetType(pcu,PCBJACOBI);
	}else if(strcmp(pc_type,"asm")==0){
		PCSetType(pcu,PCASM);
	}else if(strcmp(pc_type,"gamg")==0){
		PCSetType(pcu,PCGAMG);
	}else if(strcmp(pc_type,"gasm")==0){
		PCSetType(pcu,PCGASM);
	}else if(strcmp(pc_type,"jacobi")==0){
		PCSetType(pcu,PCJACOBI);
	}else if(strcmp(pc_type,"icc")==0){
		PCSetType(pcu,PCICC);
	}else if(strcmp(pc_type,"ilu")==0){
		PCSetType(pcu,PCILU);
	}else if(strcmp(pc_type,"sor")==0){
		PCSetType(pcu,PCSOR);
	}else if(strcmp(pc_type,"eisenstat")==0){
		PCSetType(pcu,PCEISENSTAT);
	}else if(strcmp(pc_type,"none")==0){
		PCSetType(pcu,PCNONE);
	}else{
	_error_("Suggested preconditioner not implemented yet!\n");
	}
	KSPSetUp(kspu);

	/* Create RHS */
	/* RHS = F1-B * pold */
	VecScale(p,-1.);MatMultAdd(B,p,f1,rhsu);VecScale(p,-1.);

	if (VerboseConvergence())
	{

		VecScale(rhsu,-1.);MatMultAdd(A,uold,rhsu,tmpu);VecScale(rhsu,-1.);
		VecNorm(tmpu,NORM_2,&tmp4);

		VecScale(f2,-1.);MatMultAdd(BT,uold,f2,tmpp);VecScale(f2,-1.);
		VecNorm(tmpp,NORM_2,&tmp6);

		KSPInitialResidual(kspu,uold,tmpu,tmpu2,resu,rhsu);
		VecNorm(resu,NORM_2,&tmp5);

	}

	/* Go solve Au0 = F1-Bp0*/
	KSPSolve(kspu,rhsu,uold);
	KSPGetIterationNumber(kspu,&noIt);

	if (VerboseConvergence())
	{

	KSPGetIterationNumber(kspu,&tmpi);
	VecScale(rhsu,-1.);MatMultAdd(A,uold,rhsu,tmpu);VecScale(rhsu,-1.);
	VecNorm(tmpu,NORM_2,&tmp2);
	KSPGetResidualNorm(kspu,&tmp1);

	_printf0_("||Au_00 - rhsu||_euc: " << tmp4 <<"\n||P(-1)(Au_00 - rhsu)||_euc: " << tmp5 << "\n ||Au_n0 - rhsu||_euc: " << tmp2<< "\n||P(-1)(Au_n0 - rhsu)||_euc: " << tmp1 << "\nIteration number: "<<tmpi<<"\n"); 	
	_printf0_("||BTu_00 - f2||_euc: " << tmp6 <<"\n");
	}

	/* Set up u_new */
	VecDuplicate(uold,&unew);VecCopy(uold,unew);
	VecAssemblyBegin(unew);VecAssemblyEnd(unew);

	/* ------------------------------------------------------------- */

	/*Get initial residual*/
	/*(1/mu(x) * g0, q) = b(q,u0) - (f2,q)  i.e.  IP * g0 = BT * u0 - F2*/

	/* Create KSP context */
	KSPCreate(IssmComm::GetComm(),&kspip);
	#if PETSC_VERSION_GE(3,5,0)
	KSPSetOperators(kspip,IP,IP);
	#else
	KSPSetOperators(kspip,IP,IP,DIFFERENT_NONZERO_PATTERN);
	#endif

	/* Create RHS */
	/* RHS = BT * uold - F2 */
	VecScale(uold,-1.);MatMultAdd(BT,uold,f2,rhsp);VecScale(uold,-1.);

	/* Set KSP & PC options */
	KSPSetType(kspip,KSPGMRES);
	KSPSetInitialGuessNonzero(kspip,PETSC_TRUE);

	KSPGetPC(kspip,&pcp);
	PCSetType(pcp,PCJACOBI);
	/* Note: Systems in the pressure space are cheap, so we can afford a good tolerance */
	KSPSetTolerances(kspip,1e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
	KSPSetUp(kspip);

	if (VerboseConvergence())
	{
		KSPInitialResidual(kspip,gold,tmpp,tmpp2,resp,rhsp);
		VecScale(rhsp,-1.);MatMultAdd(IP,gold,rhsp,tmpp);VecScale(rhsp,-1.);
		VecNorm(resp,NORM_2,&tmp5);
		VecNorm(tmpp,NORM_2,&tmp4);
	}

	/* Go solve */
	KSPSolve(kspip,rhsp,gold);

	if (VerboseConvergence())
	{

	KSPGetResidualNorm(kspip,&tmp1);
	VecScale(rhsp,-1.);MatMultAdd(IP,gold,rhsp,tmpp);VecScale(rhsp,-1.);
	VecNorm(tmpp,NORM_2,&tmp2);

	KSPGetIterationNumber(kspip,&tmpi);
	_printf0_("||IP*g00 - rhsp||_euc: " << tmp4 <<"\n||Jac(-1)*(IP*g00-rhsp)||_euc: " << tmp5 << "\n||IP*gn0-rhsp||_euc: " << tmp2<< "\n||Jac(-1)*(IP*gn0-rhsp)||_euc: " << tmp1 << "\nIteration number: "<<tmpi<<"\n"); 	
	}

	/*Initial residual*/
	MatMult(IP,gold,tmpp);VecDot(gold,tmpp,&initRnorm);
	initRnorm = sqrt(initRnorm);
	_printf0_("inner product norm g0: "<<initRnorm<<"\n");

	/*Iterate only if inital residual large enough*/
	if(initRnorm > 1e-5)
	{

	/* Further setup */
	VecDuplicate(gold,&gnew);VecCopy(gold,gnew);
	VecAssemblyBegin(gnew);VecAssemblyEnd(gnew);

	/* ------------------------------------------------------------ */

	/*Set initial search direction*/
	/*w0 = g0*/
	VecDuplicate(gold,&wold);VecCopy(gold,wold);
	VecAssemblyBegin(wold);VecAssemblyEnd(wold);

	/* Count number of iterations */
	int count = 0;

	/* CG iteration*/
	for(;;){
		if(VerboseConvergence())
		{
		 _printf0_("\n###### Step " << count << " ######\n");
		}

		/*Realizing the step size part 2: chim */
		/*a(chim,v) = b(wm,v)  i.e.  A * chim = B * wm */
		/*chim_DIR = 0*/
		VecScale(wold,1.);MatMult(B,wold,rhsu);VecScale(wold,1.);
		VecSet(chi,0.0);

		if(VerboseConvergence())
		{
		VecScale(rhsu,-1.);MatMultAdd(A,chi,rhsu,tmpu);VecScale(rhsu,-1.);
		VecNorm(tmpu,NORM_2,&tmp4);

		KSPInitialResidual(kspu,chi,tmpu,tmpu2,resu,rhsu);
		VecNorm(resu,NORM_2,&tmp5);

		}

			KSPSolve(kspu,rhsu,chi);

		if (VerboseConvergence())
		{
		VecNorm(chi,NORM_2,&tmp1);
		KSPGetResidualNorm(kspu,&tmp2);

		VecScale(rhsu,-1.);MatMultAdd(A,chi,rhsu,tmpu);VecScale(rhsu,-1.);
		VecNorm(tmpu,NORM_2,&tmp3);

		KSPGetIterationNumber(kspu,&tmpi);
		_printf0_("||chi_nk||_euc: "<< tmp1 << "\n||A*chi_0k - rhsu||_euc: "<<tmp4<< "\n||P(-1)*(A*chi_0k-rhsu)||_euc: " << tmp5 << "\n||A*chi_nk - rhsu||_euc: "<< tmp3 <<"\n||P(-1)*(A*chi_nk - rhsu)||_euc: " << tmp2 <<"\nIteration Number: " << tmpi <<"\n");
		}

		/* ---------------------------------------------------------- */

		/*Set step size*/
		/*rhom = [(wm)^T * IP^-1 * (BT * um - F2)]/[(wm)^T * IP^-1 * BT * chim]*/
		MatMult(IP,gold,tmpp);
		VecDot(tmpp,wold,&rho);

		MatMult(BT,chi,tmpp);
		VecDot(tmpp,wold,&tmpScalar);
		rho = rho/tmpScalar;

		/* ---------------------------------------------------------- */

		/*Pressure update*/
		/*p(m+1) = pm + rhom * wm*/
		VecAXPY(p,-1.*rho,wold);

		/*Velocity update*/
		/*u(m+1) = um - rhom * chim*/
		VecWAXPY(unew,rho,chi,uold);
		VecNorm(unew,NORM_2,&tmpScalar);

		if (VerboseConvergence())
		{
		VecScale(p,-1.);MatMultAdd(B,p,f1,rhsu);VecScale(p,-1.);
		VecScale(rhsu,-1.);MatMultAdd(A,unew,rhsu,tmpu);VecScale(rhsu,-1.);
		VecNorm(tmpu,NORM_2,&tmp6);

		VecScale(f2,-1);MatMultAdd(BT,unew,f2,tmpp);VecScale(f2,-1);
		VecNorm(tmpp,NORM_2,&tmp7);
		_printf0_("Momentum balance norm: " << tmp6 <<"\n");
		_printf0_("Incompressibility norm: " << tmp7 <<"\n");
		}

		/* ---------------------------------------------------------- */

		/*Residual update*/
		/*g(m+1) = gm - rhom * BT * chim*/
		MatMult(BT,chi,tmpp);
		MatMult(IP,gold,tmpp2);
		VecWAXPY(rhsp,-1.*rho,tmpp,tmpp2);
		KSPSolve(kspip,rhsp,gnew);

		/* ---------------------------------------------------------- */

		MatMult(IP,gnew,tmpp);

		VecDot(tmpp,gnew,&gamma);
		rnorm = sqrt(gamma);

		/*BREAK if norm(g(m+0),2) < TOL or pressure space has been fully searched*/
		if(rnorm < TOL*initRnorm) 
		{
			break;
		}else if(rnorm > 1000*initRnorm)
		{
		 _printf0_("L2-Norm g: "<< rnorm << "\n");
		 _printf0_("Solver diverged and ends prematurely.\n");
		 break;
		}else{
			_printf0_("L2-Norm g: "<< rnorm << "\n");
		}

		/*Break prematurely if solver doesn't reach desired tolerance in reasonable time frame*/
		if(count > 10./TOL)   
		{ 
		 _printf0_("Ended after " << ceil(5./TOL) << " CG iterations\n");
		 break;
		}

		/* ---------------------------------------------------------- */

		/*Directional update*/
		MatMult(IP,gold,tmpp);
		VecDot(tmpp,gold,&tmpScalar);
		gamma = gamma/tmpScalar;

		VecWAXPY(wnew,gamma,wold,gnew);

		/* Assign new to old iterates */
		VecCopy(wnew,wold);VecCopy(gnew,gold);VecCopy(unew,uold);

		count++;
	}
	}

	/* Restore pressure and velocity sol. vectors to its global form */
	VecRestoreSubVector(out_uf->pvector->vector,isv,&unew);
	VecRestoreSubVector(out_uf->pvector->vector,isp,&p);

	/*return output pointer*/
	*puf=out_uf;

	/* Cleanup */
	KSPDestroy(&kspu);KSPDestroy(&kspip);

	MatDestroy(&A);MatDestroy(&B);MatDestroy(&BT);MatDestroy(&IP);

	VecDestroy(&p);VecDestroy(&uold);VecDestroy(&unew);VecDestroy(&rhsu);VecDestroy(&rhsp);
	VecDestroy(&gold);VecDestroy(&gnew);VecDestroy(&wold);VecDestroy(&wnew);VecDestroy(&chi);
	VecDestroy(&tmpp);VecDestroy(&tmpu);VecDestroy(&f1);VecDestroy(&f2);
	VecDestroy(&tmpu2);VecDestroy(&tmpp2);VecDestroy(&resu);VecDestroy(&resp);

}/*}}}*/
void convergence_schurcg(bool* pconverged, Matrix<IssmDouble>* Kff,Vector<IssmDouble>* pf,Vector<IssmDouble>* uf,Vector<IssmDouble>* old_uf,IssmDouble eps_res,IssmDouble eps_rel,IssmDouble eps_abs,IS isv,IS isp){/*{{{*/

	/*output*/
	bool converged=false;

	/*intermediary*/
	//Vector<IssmDouble>* KU=NULL;
	//Vector<IssmDouble>* KUF=NULL;
	//Vector<IssmDouble>* KUold=NULL;
	//Vector<IssmDouble>* KUoldF=NULL;
	Vector<IssmDouble>* duf=NULL;
	IssmDouble ndu,nduinf,nu;
	IssmDouble nKUF;
	IssmDouble nKUoldF;
	IssmDouble nF;
	IssmDouble solver_residue,res;
	int analysis_type;

	PMat A, B, BT;
	PVec u,p,uold,pold,f1,f2,tmp,res1,res2;
	int n_u,n_p;
	IssmDouble rnorm1, rnorm2;

	if(VerboseModule()) _printf0_("   checking convergence\n");

	/*If uf is NULL in input, f-set is nil, model is fully constrained, therefore converged from 
	 * the get go: */
	if(uf->IsEmpty()){
		*pconverged=true;
		return;
	}

  /* Note: SchurCG also constructs the Schur preconditioner and stores it in the free block of Kff 
   *			[A    B]
	* Kff =  |      |
	*			[B^T IP]
   * To calculate the residual, only the necessary blocks need to be extracted */

	/*Extract A, B, B^T */
	#if PETSC_VERSION_GT(3,8,0)
	MatCreateSubMatrix(Kff->pmatrix->matrix,isv,isv,MAT_INITIAL_MATRIX,&A);
	MatCreateSubMatrix(Kff->pmatrix->matrix,isv,isp,MAT_INITIAL_MATRIX,&B);
	MatCreateSubMatrix(Kff->pmatrix->matrix,isp,isv,MAT_INITIAL_MATRIX,&BT);
	#else
	MatGetSubMatrix(Kff->pmatrix->matrix,isv,isv,MAT_INITIAL_MATRIX,&A);
	MatGetSubMatrix(Kff->pmatrix->matrix,isv,isp,MAT_INITIAL_MATRIX,&B);
	MatGetSubMatrix(Kff->pmatrix->matrix,isp,isv,MAT_INITIAL_MATRIX,&BT);
	#endif

		/*no. of free nodes in velocity/pressure space*/
		MatGetSize(B,&n_u,&n_p);

		/*Extract values corresponding to the free velocity/pressure nodes*/
		VecCreate(IssmComm::GetComm(),&p);VecSetSizes(p,PETSC_DECIDE,n_p);VecSetFromOptions(p);
		VecAssemblyBegin(p);VecAssemblyEnd(p);
		VecCreate(IssmComm::GetComm(),&u);VecSetSizes(u,PETSC_DECIDE,n_u);VecSetFromOptions(u);
		VecAssemblyBegin(u);VecAssemblyEnd(u);

		VecGetSubVector(uf->pvector->vector,isv,&u);
		VecGetSubVector(uf->pvector->vector,isp,&p);

		/*Extract values of the RHS corresponding to the first/second block*/
		VecDuplicate(u,&f1);VecSet(f1,1.0);
		VecDuplicate(p,&f2);VecSet(f2,1.0);
		VecGetSubVector(pf->pvector->vector,isv,&f1);
		VecGetSubVector(pf->pvector->vector,isp,&f2);

		/*Allocate intermediaries*/
		VecDuplicate(u,&res1);VecSet(res1,1.0);
		VecDuplicate(u,&tmp);VecSet(tmp,1.0);
		VecDuplicate(p,&res2);VecSet(res2,1.0);

	/*Display solver caracteristics*/
	if (VerboseConvergence()){

		/*Calculate res1 = A*u + B*p - f1*/
		VecScale(f1,-1.);MatMultAdd(A,u,f1,tmp);MatMultAdd(B,p,tmp,res1);VecScale(f1,-1.);
		/*Calculate res2 = B^T * u - f2*/
		VecScale(f2,-1.);MatMultAdd(BT,u,f2,res2);VecScale(f2,-1.);

		/*compute norm(res1), norm(res2), norm(F) and residue*/
		VecNorm(res1,NORM_2,&rnorm1);VecNorm(res2,NORM_2,&rnorm2);
		nKUF=sqrt(rnorm1*rnorm1 + rnorm2*rnorm2);
		nF=pf->Norm(NORM_TWO);
		solver_residue=nKUF/nF;
		_printf0_("\n" << "   solver residue: norm(KU-F)/norm(F)=" << solver_residue << "\n");
		if(xIsNan<IssmDouble>(solver_residue)){
			//Kff->Echo();
		}

	}
	/*clean up*/
	VecRestoreSubVector(uf->pvector->vector,isv,&u);
	VecRestoreSubVector(uf->pvector->vector,isp,&p);

	/*Extract values corresponding to velocity/pressure on the old solution*/
	VecGetSubVector(old_uf->pvector->vector,isv,&uold);
	VecGetSubVector(old_uf->pvector->vector,isp,&pold);

	/*Force equilibrium (Mandatory)*/

	/*Calculate res1 = A*uold + B*pold - f1*/
	VecScale(f1,-1.);MatMultAdd(A,uold,f1,tmp);MatMultAdd(B,pold,tmp,res1);VecScale(f1,-1.);
	/*Calculate res2 = B^T * uold - f2*/
	VecScale(f2,-1.);MatMultAdd(BT,uold,f2,res2);VecScale(f2,-1.);

	/*compute norm(res1), norm(res2), norm(F) and residue*/
	VecNorm(res1,NORM_2,&rnorm1);VecNorm(res2,NORM_2,&rnorm2);
	nKUoldF=sqrt(rnorm1*rnorm1 + rnorm2*rnorm2);
	nF=pf->Norm(NORM_TWO);
	res=nKUoldF/nF;
	if (xIsNan<IssmDouble>(res)){
		_printf0_("norm nf = " << nF << "f and norm kuold = " << nKUoldF << "f\n");
		_error_("mechanical equilibrium convergence criterion is NaN!");
	}

	MatDestroy(&A);MatDestroy(&B);MatDestroy(&BT);
	VecRestoreSubVector(pf->pvector->vector,isv,&f1);
	VecRestoreSubVector(pf->pvector->vector,isp,&f2);
	VecDestroy(&res1);VecDestroy(&res2);VecDestroy(&tmp);
	VecRestoreSubVector(old_uf->pvector->vector,isv,&uold);
	VecRestoreSubVector(old_uf->pvector->vector,isp,&pold);

	//print
	if(res<eps_res){
		if(VerboseConvergence()) _printf0_(setw(50)<<left<<"   mechanical equilibrium convergence criterion"<<res*100<< " < "<<eps_res*100<<" %\n");
		converged=true;
	}
	else{ 
		if(VerboseConvergence()) _printf0_(setw(50)<<left<<"   mechanical equilibrium convergence criterion"<<res*100<<" > "<<eps_res*100<<" %\n");
		converged=false;
	}

	/*Relative criterion (optional)*/
	if (!xIsNan<IssmDouble>(eps_rel) || (VerboseConvergence())){

		//compute norm(du)/norm(u)
		duf=old_uf->Duplicate(); old_uf->Copy(duf); duf->AYPX(uf,-1.0);
		ndu=duf->Norm(NORM_TWO); nu=old_uf->Norm(NORM_TWO);

		if (xIsNan<IssmDouble>(ndu) || xIsNan<IssmDouble>(nu)) _error_("convergence criterion is NaN!");

		//clean up
		delete duf;

		//print
		if (!xIsNan<IssmDouble>(eps_rel)){
			if((ndu/nu)<eps_rel){
				if(VerboseConvergence()) _printf0_(setw(50) << left << "   Convergence criterion: norm(du)/norm(u)" << ndu/nu*100 << " < " << eps_rel*100 << " %\n");
			}
			else{ 
				if(VerboseConvergence()) _printf0_(setw(50) << left << "   Convergence criterion: norm(du)/norm(u)" << ndu/nu*100 << " > " << eps_rel*100 << " %\n");
				converged=false;
			}
		}
		else _printf0_(setw(50) << left << "   Convergence criterion: norm(du)/norm(u)" << ndu/nu*100 << " %\n");

	}

	/*Absolute criterion (Optional) = max(du)*/
	if (!xIsNan<IssmDouble>(eps_abs) || (VerboseConvergence())){

		//compute max(du)
		duf=old_uf->Duplicate(); old_uf->Copy(duf); duf->AYPX(uf,-1.0);
		ndu=duf->Norm(NORM_TWO); nduinf=duf->Norm(NORM_INF);
		if (xIsNan<IssmDouble>(ndu) || xIsNan<IssmDouble>(nu)) _error_("convergence criterion is NaN!");

		//clean up
		delete duf;

		//print
		if (!xIsNan<IssmDouble>(eps_abs)){
			if ((nduinf)<eps_abs){
				if(VerboseConvergence()) _printf0_(setw(50) << left << "   Convergence criterion: max(du)" << nduinf << " < " << eps_abs << "\n");
			}
			else{
				if(VerboseConvergence()) _printf0_(setw(50) << left << "   Convergence criterion: max(du)" << nduinf << " > " << eps_abs << "\n");
				converged=false;
			}
		}
		else  _printf0_(setw(50) << left << "   Convergence criterion: max(du)" << nduinf << "\n");

	}

	/*assign output*/
	*pconverged=converged;
}/*}}}*/
void solutionsequence_schurcg(FemModel* femmodel){/*{{{*/

	/*intermediary: */
	Matrix<IssmDouble>* Kff = NULL;
	Matrix<IssmDouble>* Kfs = NULL;
	Vector<IssmDouble>* ug  = NULL;
	Vector<IssmDouble>* uf  = NULL;
	Vector<IssmDouble>* old_uf = NULL;
	Vector<IssmDouble>* pf  = NULL;
	Vector<IssmDouble>* df  = NULL;
	Vector<IssmDouble>* ys  = NULL;

	/*parameters:*/
	int max_nonlinear_iterations;
	int configuration_type;
	int precond;
	IssmDouble eps_res,eps_rel,eps_abs;

	/*Recover parameters: */
	femmodel->parameters->FindParam(&max_nonlinear_iterations,StressbalanceMaxiterEnum);
	femmodel->parameters->FindParam(&eps_res,StressbalanceRestolEnum);
	femmodel->parameters->FindParam(&eps_rel,StressbalanceReltolEnum);
	femmodel->parameters->FindParam(&eps_abs,StressbalanceAbstolEnum);
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	femmodel->UpdateConstraintsx();
	int size;
	int  count=0;
	bool converged=false;

	/*Start non-linear iteration using input velocity: */
	GetSolutionFromInputsx(&ug,femmodel);
	Reducevectorgtofx(&uf, ug, femmodel->nodes,femmodel->parameters);

	/*Update once again the solution to make sure that vx and vxold are similar*/
	InputUpdateFromConstantx(femmodel,converged,ConvergedEnum);
	InputUpdateFromSolutionx(femmodel,ug);

	for(;;){

		/*save pointer to old velocity*/
		delete old_uf; old_uf=uf;
		delete ug;

		/*Get stiffness matrix and Load vector*/
		SystemMatricesx(&Kff,&Kfs,&pf,&df,NULL,femmodel);
		CreateNodalConstraintsx(&ys,femmodel->nodes);
		Reduceloadx(pf, Kfs, ys); delete Kfs;

		/*Create pressure matrix of choice*/
		#if PETSC_VERSION_LT(3,7,0)
		PetscOptionsGetInt(PETSC_NULL,"-schur_pc",&precond,NULL);
		#elif PETSC_VERSION_LT(3,19,0)
		PetscOptionsGetInt(NULL,PETSC_NULL,"-schur_pc",&precond,NULL);
		#else
		PetscOptionsGetInt(NULL,PETSC_NULLPTR,"-schur_pc",&precond,NULL);
		#endif

		StressbalanceAnalysis* analysis = new StressbalanceAnalysis();	

		/* Construct Schur preconditioner matrix or mass matrix depending on input */
		if(precond)
		{
			for(Object* & object : femmodel->elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				ElementMatrix* Ie = analysis->CreateSchurPrecondMatrix(element);
				if(Ie) Ie->AddToGlobal(Kff,NULL);
				delete Ie;
			}
		}else{

			for(Object* & object : femmodel->elements->objects){
				Element* element2=xDynamicCast<Element*>(object);
				ElementMatrix* Ie2 = analysis->CreatePressureMassMatrix(element2);
				if(Ie2) Ie2->AddToGlobal(Kff,NULL);
				delete Ie2;
			}
		}

		delete analysis;

		/*Obtain index sets for velocity and pressure components */
		IS isv = NULL;
		IS isp = NULL;
		#if PETSC_VERSION_MAJOR==3

		/*Make indices out of doftypes: */
		if(!(df->pvector->vector))_error_("need doftypes for FS solver!\n");
			DofTypesToIndexSet(&isv,&isp,df->pvector->vector,FSSolverEnum);
		#else
			_error_("Petsc 3.X required");
		#endif
		Kff->Assemble();

		/*Solve*/
		femmodel->profiler->Start(SOLVER);
		_assert_(Kff->type==PetscMatType); 

		SchurCGSolver(&uf,
					Kff->pmatrix->matrix,
					pf->pvector->vector,
					old_uf->pvector->vector,
					isv,
					isp,
					femmodel->parameters);
		femmodel->profiler->Stop(SOLVER);

		/*Merge solution from f set to g set*/
		Mergesolutionfromftogx(&ug, uf,ys,femmodel->nodes,femmodel->parameters);delete ys;

		/*Check for convergence and update inputs accordingly*/
		convergence_schurcg(&converged,Kff,pf,uf,old_uf,eps_res,eps_rel,eps_abs,isv,isp); delete Kff; delete pf; delete df;
		count++;

		if(count>=max_nonlinear_iterations){
			_printf0_("   maximum number of nonlinear iterations (" << max_nonlinear_iterations << ") exceeded\n"); 
			converged=true;
		}

		InputUpdateFromConstantx(femmodel,converged,ConvergedEnum);
		InputUpdateFromSolutionx(femmodel,ug);

		/*Increase count: */
		if(converged==true){
			femmodel->results->AddResult(new GenericExternalResult<IssmDouble>(femmodel->results->Size()+1,StressbalanceConvergenceNumStepsEnum,count));
			break;
		}
	}

	if(VerboseConvergence()) _printf0_("\n   total number of iterations: " << count << "\n");

	/*clean-up*/
	delete uf;
	delete ug;
	delete old_uf;

}/*}}}*/

#else
void solutionsequence_schurcg(FemModel* femmodel){_error_("PETSc needs to be installed");}
#endif
