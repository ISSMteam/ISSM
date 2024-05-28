/*!\file:  Mergesolutionfromftogx.h
 * \brief merge solution back from f set into g set
 */ 

#ifndef _MERGESOLUTIONFROMFTOGX_H
#define _MERGESOLUTIONFROMFTOGX_H

#include "../../classes/classes.h"

/* local prototypes: */
void	Mergesolutionfromftogx( Vector<IssmDouble>** pug, Vector<IssmDouble>* uf, Vector<IssmDouble>* ys, Nodes* nodes, Parameters* parameters, bool flag_ys0=false);

#endif  /* _MERGESOLUTIONFROMFTOGX_H */
