/*!\file: convergence.cpp
 * \brief: figure out if convergence has been reached
 */

#include "../classes/classes.h"
#include "../modules/modules.h"
#include "../shared/shared.h"

void convergence(bool* pconverged, Matrix<IssmDouble>* Kff,Vector<IssmDouble>* pf,Vector<IssmDouble>* uf,Vector<IssmDouble>* old_uf,IssmDouble eps_res,IssmDouble eps_rel,IssmDouble eps_abs){

	/*output*/
	bool converged=false;

	/*intermediary*/
	Vector<IssmDouble>* KUold=NULL;
	Vector<IssmDouble>* KUoldF=NULL;
	Vector<IssmDouble>* duf=NULL;
	IssmDouble ndu,nduinf,nu;
	IssmDouble nKUoldF;
	IssmDouble nF;
	IssmDouble res;
	int analysis_type;
	int ufsize;

	if(VerboseModule()) _printf0_("   checking convergence\n");

	/*If uf is NULL in input, f-set is nil, model is fully constrained, therefore converged from
	 * the get go: */
	uf->GetSize(&ufsize);
	if(ufsize==0){
		*pconverged=true;
		return;
	}

	/* if(uf->IsEmpty()){ */
	/* 	*pconverged=true; */
	/* 	return; */
	/* } */

	/*Force equilibrium (Mandatory)*/

	/*compute K[n]U[n-1]F = K[n]U[n-1] - F*/
	_assert_(uf); _assert_(Kff);
	KUold=uf->Duplicate(); Kff->MatMult(old_uf,KUold);
	KUoldF=KUold->Duplicate();KUold->Copy(KUoldF); KUoldF->AYPX(pf,-1.0);
	nKUoldF=KUoldF->Norm(NORM_TWO);
	nF=pf->Norm(NORM_TWO)+1.e-14;
	res=nKUoldF/nF;
	if (xIsNan<IssmDouble>(res)){
		_printf0_("norm nf = " << nF << "f and norm kuold = " << nKUoldF << "f\n");
		_error_("mechanical equilibrium convergence criterion is NaN!");
	}

	//clean up
	delete KUold;
	delete KUoldF;

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
		ndu=duf->Norm(NORM_TWO); nu=old_uf->Norm(NORM_TWO)+1.e-14;

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
		nduinf=duf->Norm(NORM_INF);
		if (xIsNan<IssmDouble>(nduinf)) _error_("convergence criterion is NaN!");

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
}
