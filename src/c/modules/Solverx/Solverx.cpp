/*!\file Solverx
 * \brief solver
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Solverx.h"
#include "../../shared/shared.h"
#include "../../classes/Params/Parameters.h"

void	Solverx(Vector<IssmDouble>** puf, Matrix<IssmDouble>* Kff, Vector<IssmDouble>* pf, Vector<IssmDouble>* uf0,Vector<IssmDouble>* df, Parameters* parameters){

	/*Create Solver Object*/
	Solver<IssmDouble>* solver=new Solver<IssmDouble>(Kff,pf,uf0,df,parameters);

	#if defined(_HAVE_CODIPACK_)
	auto& tape = IssmDouble::getTape();
	auto pos = tape.getPosition();
	#endif

	/*Solve:*/
	if(VerboseModule()) _printf0_("   Solving matrix system\n");
	Vector<IssmDouble>* uf=solver->Solve();

	/*Check convergence, if failed, try recovery model*/
	if(!checkconvergence(Kff,pf,uf,parameters)){
		_printf0_("WARNING: Solver failed, Trying Recovery Mode\n");
		delete uf;

		#if defined(_HAVE_CODIPACK_)
		tape.resetTo(pos);
		#endif

		ToolkitsOptionsFromAnalysis(parameters,RecoveryAnalysisEnum);
		uf=solver->Solve();

		if(!checkconvergence(Kff,pf,uf,parameters)) _error_("Recovery solver failed...");
	}

	/*clean up and assign output pointers:*/
	_assert_(puf);
	delete solver;
	*puf=uf;
}
bool checkconvergence(Matrix<IssmDouble>* Kff,Vector<IssmDouble>* pf,Vector<IssmDouble>* uf,Parameters* parameters){

	/*Recover parameters: */
	IssmDouble solver_residue_threshold;
	parameters->FindParam(&solver_residue_threshold,SettingsSolverResidueThresholdEnum);

	/*don't check convergence if NaN*/
	if(xIsNan<IssmDouble>(solver_residue_threshold)) return true;

	/*compute KUF = KU - F = K*U - F*/
	Vector<IssmDouble>* KU  = uf->Duplicate(); Kff->MatMult(uf,KU);
	Vector<IssmDouble>* KUF = KU->Duplicate(); KU->Copy(KUF); KUF->AYPX(pf,-1.);
	delete KU;

	/*compute norm(KUF), norm(F)*/
	IssmDouble nKUF=KUF->Norm(NORM_TWO);
	IssmDouble nF=pf->Norm(NORM_TWO);
	delete KUF;

	/*Check solver residue*/
	IssmDouble solver_residue = 0.;
	if(nF>0.)  solver_residue = nKUF/(nF);
	if(VerboseConvergence()) _printf0_("\n   solver residue: norm(KU-F)/norm(F)=" << solver_residue << "\n");
	if(xIsNan<IssmDouble>(solver_residue)) _error_("Solver residue is NaN");

	/*Check convergence*/
	if(solver_residue>solver_residue_threshold){
		_printf0_("solver residue too high!: norm(KU-F)/norm(F)=" << solver_residue << " > "<<solver_residue_threshold<<" (md.settings.solver_residue_threshold)\n");
		return false;
	}
	else{
		return true;
	}

}
