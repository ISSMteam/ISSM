/*!\file: solutionsequence_fct.cpp
 * \brief: numerical core of flux corrected transport solution
 */

#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../analyses/analyses.h"
#define USEPENALTYMETHOD false

#if defined(_HAVE_PETSC_) & !defined(_HAVE_CODIPACK_)
void CreateDMatrix(Mat* pD,Mat K){/*{{{*/
	/*Create D matrix such that:
	 *
	 * d_ij = max( -k_ij,0,-k_ji) off diagonal
	 *
	 * d_ii = - sum_{i!=j} d_ij for the diagonal
	 *
	 */

	/*Intermediaries*/
	int        ncols,ncols2,rstart,rend;
	double     d,diagD;
	Mat        D        = NULL;
	Mat        K_transp = NULL;
	int*       cols  = NULL;
	int*       cols2 = NULL;
	double*    vals  = NULL;
	double*    vals2 = NULL;

	/*First, we need to transpose K so that we access both k_ij and k_ji*/
	MatTranspose(K,MAT_INITIAL_MATRIX,&K_transp);

	/*Initialize output (D has the same non zero pattern as K)*/
	MatDuplicate(K,MAT_SHARE_NONZERO_PATTERN,&D);

	/*Go through the rows of K an K' and build D*/
	MatGetOwnershipRange(K,&rstart,&rend);
	for(int row=rstart; row<rend; row++){
		diagD = 0.;
		MatGetRow(K       ,row,&ncols, (const int**)&cols, (const double**)&vals);
		MatGetRow(K_transp,row,&ncols2,(const int**)&cols2,(const double**)&vals2);
		_assert_(ncols==ncols2);
		for(int j=0; j<ncols; j++) {
			_assert_(cols[j]==cols2[j]);
			d = max(max(-vals[j],-vals2[j]),0.);
			MatSetValue(D,row,cols[j],d,INSERT_VALUES);
			if(cols[j]!=row) diagD -= d;
		}
		MatSetValue(D,row,row,diagD,INSERT_VALUES);
		MatRestoreRow(K       ,row,&ncols, (const int**)&cols, (const double**)&vals);
		MatRestoreRow(K_transp,row,&ncols2,(const int**)&cols2,(const double**)&vals2);
	}
	MatAssemblyBegin(D,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(  D,MAT_FINAL_ASSEMBLY);

	/*Clean up and assign output pointer*/
	MatFree(&K_transp);
	*pD = D;
}/*}}}*/
void CreateLHS(Mat* pLHS,IssmDouble* pdmax,Mat K,Mat D,Vec Ml,IssmDouble theta,IssmDouble deltat,FemModel* femmodel,int configuration_type){/*{{{*/
	/*Create Left Hand side of Lower order solution
	 *
	 * LHS = [ML - theta*detlat *(K+D)^n+1]
	 *
	 */

	/*Intermediaries*/
	int        dof,ncols,ncols2,rstart,rend;
	double     d,mi,dmax = 0.;
	Mat        LHS   = NULL;
	int*       cols  = NULL;
	int*       cols2 = NULL;
	double*    vals  = NULL;
	double*    vals2 = NULL;

	MatDuplicate(K,MAT_SHARE_NONZERO_PATTERN,&LHS);
	MatGetOwnershipRange(K,&rstart,&rend);
	for(int row=rstart; row<rend; row++){
		MatGetRow(K,row,&ncols, (const int**)&cols, (const double**)&vals);
		MatGetRow(D,row,&ncols2,(const int**)&cols2,(const double**)&vals2);
		_assert_(ncols==ncols2);
		for(int j=0; j<ncols; j++) {
			_assert_(cols[j]==cols2[j]);
			d = -theta*deltat*(vals[j] + vals2[j]);
			if(cols[j]==row){
				VecGetValues(Ml,1,(const int*)&cols[j],&mi);
				d += mi;
			}
			if(fabs(d)>dmax) dmax = fabs(d);
			MatSetValue(LHS,row,cols[j],d,INSERT_VALUES);
		}
		MatRestoreRow(K,row,&ncols, (const int**)&cols, (const double**)&vals);
		MatRestoreRow(D,row,&ncols2,(const int**)&cols2,(const double**)&vals2);
	}
	MatAssemblyBegin(LHS,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(  LHS,MAT_FINAL_ASSEMBLY);

	/*Deal with Dirichlet conditions: 2 options, penalties or zeros in K matrix*/
	if(USEPENALTYMETHOD){
		/*Option 1: Penalty method*/

		/*Broadcast max(dmax)*/
		IssmDouble dmax_all;
		ISSM_MPI_Reduce(&dmax,&dmax_all,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
		ISSM_MPI_Bcast(&dmax_all,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
		dmax = dmax_all;

		dmax = dmax * 1.e+3;
		for(Object* & object: femmodel->constraints->objects){
			Constraint* constraint=xDynamicCast<Constraint*>(object);
			constraint->PenaltyDofAndValue(&dof,&d,femmodel->nodes,femmodel->parameters);
			if(dof!=-1){
				MatSetValue(LHS,dof,dof,dmax,INSERT_VALUES);
			}
		}
	}
	else{
		/*Option 2: zero stiffness matrix, 1 one diagonal*/
		int numrows = femmodel->constraints->Size();
		int* rows = xNew<int>(numrows);
		IssmDouble* rows_spc = xNew<IssmDouble>(numrows);
		numrows = 0;

		dmax = dmax * 1.e+3;
		for(Object* & object: femmodel->constraints->objects){
			Constraint* constraint=xDynamicCast<Constraint*>(object);
			constraint->PenaltyDofAndValue(&dof,&d,femmodel->nodes,femmodel->parameters);
			if(dof!=-1){
				rows[numrows] = dof;
				numrows++;
			}
		}
		MatZeroRows(LHS,numrows,rows,1.,NULL,NULL);
		xDelete<int>(rows);
		xDelete<IssmDouble>(rows_spc);
	}
	MatAssemblyBegin(LHS,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(  LHS,MAT_FINAL_ASSEMBLY);

	/*Clean up and assign output pointer*/
	*pdmax = dmax;
	*pLHS  = LHS;
}/*}}}*/
void CreateRHS(Vec* pRHS,Mat K,Mat D,Vec Ml,Vec p,Vec u,IssmDouble theta,IssmDouble deltat,IssmDouble dmax,FemModel* femmodel,int configuration_type){/*{{{*/
	/*Create Right Hand side of Lower order solution
	 *
	 * RHS = [ML + (1 - theta) deltaT L^n] u^n
	 *
	 * where L = K + D
	 *
	 */

	/*Intermediaries*/
	Vec         Ku  = NULL;
	Vec         Du  = NULL;
	Vec         RHS = NULL;
	int         dof;
	IssmDouble  d;

	/*Initialize vectors*/
	VecDuplicate(u,&Ku);
	VecDuplicate(u,&Du);
	VecDuplicate(u,&RHS);

	/*Create RHS = M*u + (1-theta)*deltat*K*u + (1-theta)*deltat*D*u*/
	MatMult(K,u,Ku);
	MatMult(D,u,Du);
	VecPointwiseMult(RHS,Ml,u);
	VecAXPBYPCZ(RHS,(1-theta)*deltat,(1-theta)*deltat,1,Ku,Du);
	VecAXPBY(RHS,deltat,1,p);
	VecFree(&Ku);
	VecFree(&Du);

	VecAssemblyBegin(RHS);
	VecAssemblyEnd(  RHS);

	/*Deal with Dirichlet conditions: 2 options, penalties or zeros in K matrix*/
	if(USEPENALTYMETHOD){
		/*Option 1: Penalty method*/
		for(Object* & object: femmodel->constraints->objects){
			Constraint* constraint=xDynamicCast<Constraint*>(object);
			constraint->PenaltyDofAndValue(&dof,&d,femmodel->nodes,femmodel->parameters);
			d = d*dmax;
			if(dof!=-1){
				VecSetValues(RHS,1,&dof,(const double*)&d,INSERT_VALUES);
			}
		}
	}
	else{
		/*Option 2: zero stiffness matrix, 1 one diagonal*/
		int  numrows = femmodel->constraints->Size();
		int* rows = xNew<int>(numrows);
		IssmDouble* rows_spc = xNew<IssmDouble>(numrows);
		numrows = 0;
		for(Object* & object: femmodel->constraints->objects){
			Constraint* constraint=xDynamicCast<Constraint*>(object);
			constraint->PenaltyDofAndValue(&dof,&d,femmodel->nodes,femmodel->parameters);
			if(dof!=-1){
				rows[numrows] = dof;
				rows_spc[numrows] = d;
				numrows++;
			}
		}
		VecSetValues(RHS,numrows,rows,rows_spc,INSERT_VALUES);
		xDelete<int>(rows);
		xDelete<IssmDouble>(rows_spc);
	}
	VecAssemblyBegin(RHS);
	VecAssemblyEnd(  RHS);

	/*Assign output pointer*/
	*pRHS = RHS;

}/*}}}*/
void RichardsonUdot(Vec* pudot,Vec u,Vec Ml,Mat K,Mat Mc){/*{{{*/
	/*Use Richardson's formula to get udot using 5 steps and an initial guess of 0
	 *
	 * udot_new = udot_old + Ml^-1 (K^(n+1) u - Mc udot_old)
	 *
	 * */

	/*Intermediaries*/
	Vec udot  = NULL;
	Vec temp1 = NULL;
	Vec temp2 = NULL;

	/*Initialize vectors*/
	VecDuplicate(u,&udot);
	VecDuplicate(u,&temp1);
	VecDuplicate(u,&temp2);

	/*Initial guess*/
	VecZeroEntries(udot);

	/*Richardson's iterations*/
	for(int i=0;i<5;i++){
		MatMult(Mc,udot,temp1);
		MatMult(K, u,   temp2);
		VecAXPBY(temp2,-1.,1.,temp1);       // temp2 = (K^(n+1) u -- Mc udot_old)
		VecPointwiseDivide(temp1,temp2,Ml); // temp1 = Ml^-1 temp2
		VecAXPBY(udot,1.,1.,temp1);
	}

	/*Clean up and assign output pointer*/
	VecFree(&temp1);
	VecFree(&temp2);
	*pudot=udot;

}/*}}}*/
void CreateRis(IssmDouble** pRi_plus_serial,IssmDouble** pRi_minus_serial,Mat Mc,Mat D,IssmDouble* ml_serial,Vec uvec,IssmDouble* u,IssmDouble* udot,IssmDouble* ulmin,IssmDouble* ulmax,IssmDouble deltat){/*{{{*/

	/*Intermediaries*/
	Vec         Ri_plus  = NULL;
	Vec         Ri_minus = NULL;
	int         ncols,ncols2,rstart,rend;
	double      d;
	int        *cols     = NULL;
	int        *cols2    = NULL;
	double     *vals     = NULL;
	double     *vals2    = NULL;

	/*Initialize vectors*/
	VecDuplicate(uvec,&Ri_plus);
	VecDuplicate(uvec,&Ri_minus);

	/*Get global extremas*/
	IssmDouble uLmin_global =  ulmin[0];
	IssmDouble uLmax_global =  ulmax[0];
	VecGetSize(uvec,&ncols);
	for(int i=1;i<ncols;i++) if(ulmin[i]<uLmin_global) uLmin_global = ulmin[i];
	for(int i=1;i<ncols;i++) if(ulmax[i]>uLmax_global) uLmax_global = ulmax[i];
	//printf("%g %g\n",uLmin_global,uLmax_global);

	/*Go through D and calculate Ris*/
	MatGetOwnershipRange(D,&rstart,&rend);
	for(int row=rstart; row<rend; row++){
		double Pi_plus  = 0.;
		double Pi_minus = 0.;
		MatGetRow(Mc,row,&ncols, (const int**)&cols, (const double**)&vals);
		MatGetRow(D ,row,&ncols2,(const int**)&cols2,(const double**)&vals2);
		_assert_(ncols==ncols2);
		for(int j=0; j<ncols; j++) {
			_assert_(cols[j]==cols2[j]);
			d = vals[j]*(udot[row] - udot[cols[j]]) + vals2[j]*(u[row] - u[cols[j]]);
			if(row!=cols[j]){
				if(d>0.){
					Pi_plus  += d;
				}
				else{
					Pi_minus += d;
				}
			}
		}

		/*Compute Qis and Ris*/
		//double Qi_plus  = ml_serial[row]/deltat*(3. - u[row]);
		//double Qi_minus = ml_serial[row]/deltat*(2. - u[row]);
		double Qi_plus  = ml_serial[row]/deltat*(ulmax[row] - u[row]);
		double Qi_minus = ml_serial[row]/deltat*(ulmin[row] - u[row]);
		//double Qi_plus  = ml_serial[row]/deltat*(uLmax_global - u[row]);
		//double Qi_minus = ml_serial[row]/deltat*(uLmin_global - u[row]);
		_assert_(Qi_plus  >= 0.);
		_assert_(Qi_minus <= 0.);
		d = 1.;
		if(Pi_plus!=0.) d = min(1.,Qi_plus/Pi_plus);
		VecSetValue(Ri_plus,row,d,INSERT_VALUES);
		d = 1.;
		if(Pi_minus!=0.) d = min(1.,Qi_minus/Pi_minus);
		VecSetValue(Ri_minus,row,d,INSERT_VALUES);

		MatRestoreRow(Mc,row,&ncols, (const int**)&cols, (const double**)&vals);
		MatRestoreRow(D ,row,&ncols2,(const int**)&cols2,(const double**)&vals2);
	}
	VecAssemblyBegin(Ri_plus);
	VecAssemblyEnd(  Ri_plus);
	VecAssemblyBegin(Ri_minus);
	VecAssemblyEnd(  Ri_minus);

	/*Serialize Ris*/
	IssmDouble* Ri_plus_serial  = NULL;
	IssmDouble* Ri_minus_serial = NULL;
	VecToMPISerial(&Ri_plus_serial,Ri_plus,IssmComm::GetComm());
	VecToMPISerial(&Ri_minus_serial,Ri_minus,IssmComm::GetComm());
	VecFree(&Ri_plus);
	VecFree(&Ri_minus);

	/*Assign output pointer*/
	*pRi_plus_serial  = Ri_plus_serial;
	*pRi_minus_serial = Ri_minus_serial;
}/*}}}*/
void CreateFbar(Vec* pFbar,IssmDouble* Ri_plus,IssmDouble* Ri_minus,Mat Mc,Mat D,IssmDouble* udot,IssmDouble* u,Vec uvec){/*{{{*/

	/*Intermediaries*/
	Vec         Fbar = NULL;
	int         ncols,ncols2,rstart,rend;
	double      d,f_ij;
	int        *cols     = NULL;
	int        *cols2    = NULL;
	double     *vals     = NULL;
	double     *vals2    = NULL;

	/*Build fbar*/
	VecDuplicate(uvec,&Fbar);
	MatGetOwnershipRange(D,&rstart,&rend);
	for(int row=rstart; row<rend; row++){
		MatGetRow(Mc,row,&ncols, (const int**)&cols, (const double**)&vals);
		MatGetRow(D ,row,&ncols2,(const int**)&cols2,(const double**)&vals2);
		_assert_(ncols==ncols2);
		d = 0.;
		for(int j=0; j<ncols; j++) {
			_assert_(cols[j]==cols2[j]);
			if(row==cols[j]) continue;
			f_ij = vals[j]*(udot[row] - udot[cols[j]]) + vals2[j]*(u[row] - u[cols[j]]);
			if(f_ij>0){
				d += min(Ri_plus[row],Ri_minus[cols[j]]) * f_ij;
			}
			else{
				d += min(Ri_minus[row],Ri_plus[cols[j]]) * f_ij;
			}
		}
		VecSetValue(Fbar,row,d,INSERT_VALUES);
		MatRestoreRow(Mc,row,&ncols, (const int**)&cols, (const double**)&vals);
		MatRestoreRow(D ,row,&ncols2,(const int**)&cols2,(const double**)&vals2);
	}
	VecAssemblyBegin(Fbar);
	VecAssemblyEnd(  Fbar);

	/*Assign output pointer*/
	*pFbar = Fbar;
}/*}}}*/
void UpdateSolution(Vec u,Vec udot,Vec Ml,Vec Fbar,IssmDouble deltat){/*{{{*/

	/*Intermediary*/
	Vec temp1 = NULL;

	/*Compute solution u^n+1 = u_L + deltat Ml^-1 fbar*/
	VecDuplicate(u,&temp1);
	VecPointwiseDivide(temp1,Fbar,Ml); //temp1 = Ml^-1 temp2
	VecAXPBY(udot,1.,1.,temp1);
	VecAXPY(u,deltat,temp1);

	/*CLean up and return*/
	VecFree(&temp1);
}/*}}}*/
#endif
void solutionsequence_fct(FemModel* femmodel){/*{{{*/

	/*intermediary: */
	IssmDouble           theta,deltat,dmax;
	int                  configuration_type,analysis_type;
	Matrix<IssmDouble>*  K  = NULL;
	Matrix<IssmDouble>*  Mc = NULL;
	Vector<IssmDouble>*  p  = NULL;
	Vector<IssmDouble>*  ug = NULL;
	Vector<IssmDouble>*  uf = NULL;
	MasstransportAnalysis* manalysis = NULL;
	DamageEvolutionAnalysis* danalysis = NULL;

	/*Recover parameters: */
	femmodel->parameters->FindParam(&deltat,TimesteppingTimeStepEnum);
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	femmodel->UpdateConstraintsx();
	theta = 0.5; //0.5==Crank-Nicolson, 1==Backward Euler

	/*Create analysis*/
	femmodel->parameters->FindParam(&analysis_type,AnalysisTypeEnum);
	switch(analysis_type){
		case MasstransportAnalysisEnum:
			manalysis = new MasstransportAnalysis();
			manalysis->MassMatrix(&Mc,femmodel);
			manalysis->FctKMatrix(&K,NULL,femmodel);
			manalysis->FctPVector(&p,femmodel);
			break;
		case DamageEvolutionAnalysisEnum:
			danalysis = new DamageEvolutionAnalysis();
			danalysis->MassMatrix(&Mc,femmodel);
			danalysis->FctKMatrix(&K,NULL,femmodel);
			danalysis->FctPVector(&p,femmodel);
			break;
		default: _error_("analysis type " << EnumToStringx(analysis_type) << " not supported for FCT\n");
	}
	delete manalysis;
	delete danalysis;

	#ifdef _HAVE_PETSC_
	#ifdef _HAVE_CODIPACK_
		_error_("No CoDiPack handling for PETSc and fct");
	#else
	/*Convert matrices to PETSc matrices*/
	Mat D_petsc  = NULL;
	Mat LHS      = NULL;
	Vec RHS      = NULL;
	Vec u        = NULL;
	Vec udot     = NULL;
	Mat K_petsc  = K->pmatrix->matrix;
	Mat Mc_petsc = Mc->pmatrix->matrix;
	Vec p_petsc	 = p->pvector->vector;

	/*Get previous solution u^n*/
	GetSolutionFromInputsx(&ug,femmodel);
	Reducevectorgtofx(&uf, ug, femmodel->nodes,femmodel->parameters);
	delete ug;

	/*Compute lumped mass matrix*/
	Vec Ml_petsc = NULL;
	VecDuplicate(uf->pvector->vector,&Ml_petsc);
	MatGetRowSum(Mc_petsc,Ml_petsc);

	/*Create D Matrix*/
	CreateDMatrix(&D_petsc,K_petsc);

	/*Create LHS: [ML - theta*detlat *(K+D)^n+1]*/
	CreateLHS(&LHS,&dmax,K_petsc,D_petsc,Ml_petsc,theta,deltat,femmodel,configuration_type);

	/*Create RHS: [ML + (1 - theta) deltaT L^n] u^n */
	CreateRHS(&RHS,K_petsc,D_petsc,Ml_petsc,p_petsc,uf->pvector->vector,theta,deltat,dmax,femmodel,configuration_type);
	delete uf;
	delete p;

	/*Go solve lower order solution*/
	femmodel->profiler->Start(SOLVER);
	SolverxPetsc(&u,LHS,RHS,NULL,NULL, femmodel->parameters);
	femmodel->profiler->Stop(SOLVER);
	MatFree(&LHS);
	VecFree(&RHS);

	/*Richardson to calculate udot*/
	RichardsonUdot(&udot,u,Ml_petsc,K_petsc,Mc_petsc);
	delete K;

	/*Serialize u and udot*/
	IssmDouble* udot_serial = NULL;
	IssmDouble* u_serial    = NULL;
	IssmDouble* ml_serial   = NULL;
	VecToMPISerial(&udot_serial,udot    ,IssmComm::GetComm());
	VecToMPISerial(&u_serial   ,u       ,IssmComm::GetComm());
	VecToMPISerial(&ml_serial  ,Ml_petsc,IssmComm::GetComm());

	/*Anti diffusive fluxes*/
	Vec         Fbar            = NULL;
	IssmDouble *Ri_minus_serial = NULL;
	IssmDouble *Ri_plus_serial  = NULL;
	IssmDouble *ulmin           = NULL;
	IssmDouble *ulmax           = NULL;
	femmodel->GetInputLocalMinMaxOnNodesx(&ulmin,&ulmax,u_serial);
	CreateRis(&Ri_plus_serial,&Ri_minus_serial,Mc_petsc,D_petsc,ml_serial,u,u_serial,udot_serial,ulmin,ulmax,deltat);
	CreateFbar(&Fbar,Ri_plus_serial,Ri_minus_serial,Mc_petsc,D_petsc,udot_serial,u_serial,u);
	xDelete<IssmDouble>(Ri_plus_serial);
	xDelete<IssmDouble>(Ri_minus_serial);
	xDelete<IssmDouble>(ulmin);
	xDelete<IssmDouble>(ulmax);

	/*Clean up*/
	MatFree(&D_petsc);
	delete Mc;
	xDelete<IssmDouble>(udot_serial);
	xDelete<IssmDouble>(u_serial);
	xDelete<IssmDouble>(ml_serial);

	/*Compute solution u^n+1 = u_L + deltat Ml^-1 fbar*/
	UpdateSolution(u,udot,Ml_petsc,Fbar,deltat);
	uf =new Vector<IssmDouble>(u);
	VecFree(&u);
	VecFree(&Fbar);
	VecFree(&udot);
	VecFree(&Ml_petsc);

	/*Update Element inputs*/
	InputUpdateFromSolutionx(femmodel,uf);
	delete uf;

	#endif
	#else
	_error_("PETSc needs to be installed");
	#endif
}/*}}}*/
