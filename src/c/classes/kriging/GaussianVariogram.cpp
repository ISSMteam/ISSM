/*!\file GaussianVariogram.c
 * \brief: implementation of the GaussianVariogram object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*GaussianVariogram constructors and destructor*/
GaussianVariogram::GaussianVariogram(){/*{{{*/
	this->nugget = 0.2;
	this->sill   = 1;
	this->range  = SQRT3;
	return;
}
/*}}}*/
GaussianVariogram::GaussianVariogram(Options* options){/*{{{*/

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
GaussianVariogram::~GaussianVariogram(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* GaussianVariogram::copy(void){/*{{{*/
	   return new GaussianVariogram(*this);
}
/*}}}*/
void GaussianVariogram::Echo(void){/*{{{*/
	_printf_("GaussianVariogram\n");
	_printf_("   nugget: " << this->nugget << "\n");
	_printf_("   sill  : " << this->sill << "\n");
	_printf_("   range : " << this->range << "\n");
}
/*}}}*/

/*Variogram function*/
double GaussianVariogram::Covariance(double deltax,double deltay){/*{{{*/
	/*The covariance can be deduced from the variogram from the following
	 * relationship:
	 *    2 gamma = C(x,x) + C(y,y) -2 C(x,y)
	 * so
	 *    C(h) = sill - gamma                                            */
	double h2,a,cova;

	/*Calculate length square*/
	h2=deltax*deltax + deltay*deltay;

	/*If h is too small, return sill*/
	if(h2<0.0000001) return sill;

	/*compute covariance*/
	a     = 1./3.;
	cova = (sill-nugget)*exp(-h2/(a*range*range));

	return cova;
}
/*}}}*/
double GaussianVariogram::SemiVariogram(double deltax,double deltay){/*{{{*/
	/*http://en.wikipedia.org/wiki/Variogram*/
	double h2,a,gamma;

	/*Calculate length square*/
	h2=deltax*deltax + deltay*deltay;

	/*return semi-variogram*/
	a     = 1./3.;
	gamma = (sill-nugget)*(1.-exp(-h2/(a*range*range))) + nugget;

	//if(h2>1000*1000) _printf_("gamma = " << gamma << " h= " << sqrt(h2) << "\n");
	_printf_("h = " << sqrt(h2) << " gamma = " << gamma << "\n");
	return gamma;
}
/*}}}*/
