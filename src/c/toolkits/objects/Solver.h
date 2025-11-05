/*!\file:  Solver.h
 */ 

#ifndef _SOLVER_CLASS_H_
#define _SOLVER_CLASS_H_

/*Headers:*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "./Matrix.h"
#include "./Vector.h"
#include "../issm/issmtoolkit.h"
#include "../petsc/petscincludes.h"
class Parameters;

template <class doubletype> 
class Solver{

	private:
		Matrix<doubletype>* Kff;
		Vector<doubletype>* pf;
		Vector<doubletype>* uf0;
		Vector<doubletype>* df;
		Parameters* parameters;

	public:
		/*Constructors, destructors:*/
		Solver(){/*{{{*/
		}
		/*}}}*/
		Solver(Matrix<doubletype>* Kff_in, Vector<doubletype>* pf_in, Vector<doubletype>* uf0_in,Vector<doubletype>* df_in, Parameters* parameters_in){/*{{{*/

			/*In debugging mode, check that stiffness matrix and load vectors are not NULL (they can be empty)*/
			_assert_(Kff_in);
			_assert_(pf_in);

			/*initialize fields: */
			this->Kff=Kff_in;
			this->pf=pf_in;
			this->uf0=uf0_in;
			this->df=df_in;
			this->parameters=parameters_in;
		}
		/*}}}*/
		~Solver(){/*{{{*/
		}
		/*}}}*/

		/*Methods: */
		Vector<doubletype>* Solve(void){ /*{{{*/

			/*output: */
			Vector<doubletype>* uf=NULL;

			/*Initialize vector: */
			uf=new Vector<doubletype>();

			/*According to matrix type, use specific solvers: */
			switch(Kff->type){
				#ifdef _HAVE_PETSC_
				case PetscMatType:{
					PetscVec<doubletype>* uf0_vector = NULL;
					PetscVec<doubletype>* df_vector  = NULL;
					if(uf0) uf0_vector = uf0->pvector;
					if(df)  df_vector  = df->pvector;
					PetscSolve(&uf->pvector,Kff->pmatrix,pf->pvector,uf0_vector,df_vector,parameters);
					break;
								  }
				#endif
				case IssmMatType:{
					IssmSolve(&uf->ivector,Kff->imatrix,pf->ivector,parameters);
					break;
								 }
				default:
					_error_("Matrix type: " << Kff->type << " not supported yet!");
			}

			/*allocate output pointer: */
			return uf;
		}
		/*}}}*/
};
#endif //#ifndef _SOLVER_H_
