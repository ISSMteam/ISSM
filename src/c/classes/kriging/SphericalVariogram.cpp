/*!\file SphericalVariogram.c
 * \brief: implementation of the SphericalVariogram object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*SphericalVariogram constructors and destructor*/
SphericalVariogram::SphericalVariogram(){/*{{{*/
	this->nugget = 0.2;
	this->sill   = 1;
	this->range  = SQRT3;
	return;
}
/*}}}*/
SphericalVariogram::SphericalVariogram(Options* options){/*{{{*/

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
SphericalVariogram::~SphericalVariogram(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* SphericalVariogram::copy(void){/*{{{*/
	return new SphericalVariogram(*this);
}
/*}}}*/
void SphericalVariogram::Echo(void){/*{{{*/
	_printf_("SphericalVariogram\n");
	_printf_("   nugget: " << this->nugget << "\n");
	_printf_("   sill  : " << this->sill << "\n");
	_printf_("   range : " << this->range << "\n");
}
/*}}}*/

/*Variogram function*/
double SphericalVariogram::Covariance(double deltax,double deltay){/*{{{*/
	/*The covariance can be deduced from the variogram from the following
	 * relationship:
	 *    2 gamma = C(x,x) + C(y,y) -2 C(x,y)
	 * so
	 *    C(h) = sill - gamma                                            */
	double h,cova;

	/*Calculate length square*/
	h=sqrt(deltax*deltax + deltay*deltay);

	/*return covariance*/
	if(h<=range)
	 cova = (sill-nugget)*(1 - (3*h)/(2*range) + pow(h,3)/(2*pow(range,3)) );
	else
	 cova = 0.;

	return cova;
}
/*}}}*/
double SphericalVariogram::SemiVariogram(double deltax,double deltay){/*{{{*/
	/*http://en.wikipedia.org/wiki/Variogram*/
	double h,gamma;

	/*Calculate length square*/
	h=sqrt(deltax*deltax + deltay*deltay);

	/*return semi-variogram*/
	if(h<=range)
	 gamma = (sill-nugget)*( (3*h)/(2*range) - pow(h,3)/(2*pow(range,3)) ) + nugget;
	else
	 gamma = sill;

	return gamma;
}
/*}}}*/
