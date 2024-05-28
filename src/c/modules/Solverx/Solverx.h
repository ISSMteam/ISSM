/*!\file:  Solverx.h
 * \brief solver
 */ 

#ifndef _SOLVERX_H
#define _SOLVERX_H

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../toolkits/toolkits.h"

/* local prototypes: */
void	Solverx(Vector<IssmDouble>** puf, Matrix<IssmDouble>* Kff, Vector<IssmDouble>* pf, Vector<IssmDouble>* uf0,Vector<IssmDouble>* df, Parameters* parameters);
bool checkconvergence(Matrix<IssmDouble>* Kff,Vector<IssmDouble>* pf,Vector<IssmDouble>* uf,Parameters* parameters);

#endif  /* _SOLVERX_H */
