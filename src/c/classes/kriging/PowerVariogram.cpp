/*!\file PowerVariogram.c
 * \brief: implementation of the PowerVariogram object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"

/*PowerVariogram constructors and destructor*/
PowerVariogram::PowerVariogram(){/*{{{*/
	this->nugget = 0.2;
	this->slope  = 1.;
	this->power  = 1.;
	return;
}
/*}}}*/
PowerVariogram::PowerVariogram(Options* options){/*{{{*/

	/*Defaults*/
	this->nugget = 0.2;
	this->slope  = 1.;
	this->power  = 1.;

	/*Overwrite from options*/
	if(options->GetOption("nugget")) options->Get(&this->nugget,"nugget");
	if(options->GetOption("slope"))  options->Get(&this->slope,"slope");
	if(options->GetOption("power"))  options->Get(&this->power,"power");

	/*Checks*/
	if(power<=0 || power>=2) _error_("power must be betwwen 0 and 2 (0 < power < 2)");
	if(slope<=0) _error_("slope must be positive");
}
/*}}}*/
PowerVariogram::~PowerVariogram(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* PowerVariogram::copy(void){/*{{{*/
	   return new PowerVariogram(*this);
}
/*}}}*/
void PowerVariogram::Echo(void){/*{{{*/
	_printf_("PowerVariogram\n");
	_printf_("   nugget: " << this->nugget << "\n");
	_printf_("   slope : " << this->slope << "\n");
	_printf_("   power : " << this->power << "\n");
}
/*}}}*/

/*Variogram function*/
double PowerVariogram::Covariance(double deltax,double deltay){/*{{{*/
	/*The covariance can be deduced from the variogram from the following
	 * relationship:
	 *    2 gamma = C(x,x) + C(y,y) -2 C(x,y)
	 * so
	 *    C(h) = sill - gamma                                            */
	double h,cova;

	/*Calculate length square*/
	h=sqrt(deltax*deltax + deltay*deltay);

	/*return covariance*/
	cova = 9999. - this->slope*pow(h,this->power);

	return cova;
}
/*}}}*/
double PowerVariogram::SemiVariogram(double deltax,double deltay){/*{{{*/
	/*http://en.wikipedia.org/wiki/Variogram*/
	double h,gamma;

	/*Calculate length square*/
	h=sqrt(deltax*deltax + deltay*deltay);

	/*return semi-variogram*/
	gamma = this->nugget + this->slope*pow(h,this->power);

	//if(h>1000) _printf_("gamma = " << gamma << " h=" << h << "\n");
	return gamma;
}
/*}}}*/
