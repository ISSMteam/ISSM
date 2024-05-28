/*!\file Reduceloadx
 * \brief reduce loads (wring out boundary conditions)
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Reduceloadx.h"
#include "../../shared/io/io.h"

void	Reduceloadx( Vector<IssmDouble>* pf, Matrix<IssmDouble>* Kfs, Vector<IssmDouble>* y_s,bool flag_ys0){

	/*intermediary*/
	Vector<IssmDouble>*     y_s0   = NULL;
	Vector<IssmDouble>*     Kfsy_s = NULL;
	int         Kfsm,Kfsn;
	int         global_m,global_n;
	bool        fromlocalsize = true;
	bool        oldalloc  = false;

	if(VerboseModule()) _printf0_("   Dirichlet lifting applied to load vector\n");

	Kfs->GetSize(&global_m,&global_n);
	if(pf && global_m*global_n){

		/*Some checks in debugging mode*/
		_assert_(y_s);

		/*pf = pf - Kfs * y_s;*/
		Kfs->GetLocalSize(&Kfsm,&Kfsn);
		if(oldalloc)
		 Kfsy_s=new Vector<IssmDouble>(Kfsm,fromlocalsize);
		else
		 Kfsy_s=new Vector<IssmDouble>(Kfsm,global_m);

		if (flag_ys0){

			/*Create y_s0, full of 0: */
			y_s0=y_s->Duplicate();
			y_s0->Set(0.0);
			y_s0->Assemble();

			Kfs->MatMult(y_s0,Kfsy_s);
		}
		else{
			Kfs->MatMult(y_s,Kfsy_s);
		}

		pf->AXPY(Kfsy_s,-1.);
	}

	/*Free resources: and return*/
	delete y_s0;
	delete Kfsy_s;

}
