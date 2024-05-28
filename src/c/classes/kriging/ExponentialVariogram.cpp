/*!\file ExponentialVariogram.c
 * \brief: implementation of the ExponentialVariogram object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*ExponentialVariogram constructors and destructor*/
ExponentialVariogram::ExponentialVariogram(){/*{{{*/
	this->nugget = 0.2;
	this->sill   = 1;
	this->range  = SQRT3;
	return;
}
/*}}}*/
ExponentialVariogram::ExponentialVariogram(Options* options){/*{{{*/

	/*Defaults*/
	this->nugget = 0.2;
	this->sill   = 1;
	this->range  = SQRT3;

	/*Overwrite from options*/
	if(options->GetOption("nugget")) options->Get(&this->nugget,"nugget");
	if(options->GetOption("sill"))   options->Get(&this->sill,"sill");
	if(options->GetOption("range"))  options->Get(&this->range,"range");

	/*Checks*/
	if(nugget==sill) _error_("nugget and sill cannot be equal (constant semivariogram not allowed)");
}
/*}}}*/
ExponentialVariogram::~ExponentialVariogram(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* ExponentialVariogram::copy(void){/*{{{*/
	   return new ExponentialVariogram(*this);
}
/*}}}*/
void ExponentialVariogram::Echo(void){/*{{{*/
	_printf_("ExponentialVariogram\n");
	_printf_("   nugget: " << this->nugget << "\n");
	_printf_("   sill  : " << this->sill << "\n");
	_printf_("   range : " << this->range << "\n");
}
/*}}}*/

/*Variogram function*/
double ExponentialVariogram::Covariance(double deltax,double deltay){/*{{{*/
	/*The covariance can be deduced from the variogram from the following
	 * relationship:
	 *    2 gamma = C(x,x) + C(y,y) -2 C(x,y)
	 * so
	 *    C(h) = sill - gamma                                            */
	double h,a,cova;

	/*Calculate length*/
	h=sqrt(deltax*deltax + deltay*deltay);

	/*If h is too small, return sill*/
	if(h<0.0000001) return sill;

	/*compute covariance*/
	a     = 1./3.;
	cova = (sill-nugget)*exp(-h/(a*range));
	return cova;
}
/*}}}*/
double ExponentialVariogram::SemiVariogram(double deltax,double deltay){/*{{{*/
	/*http://en.wikipedia.org/wiki/Variogram*/
	double h,a,gamma;

	/*Calculate length*/
	h=sqrt(deltax*deltax + deltay*deltay);

	/*return semi-variogram*/
	a     = 1./3.;
	gamma = (sill-nugget)*(1-exp(-h/(a*range))) + nugget;
	return gamma;
}
/*}}}*/
