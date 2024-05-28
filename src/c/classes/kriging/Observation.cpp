/*!\file Observation.c
 * \brief: implementation of the Observation object
 */

#include <stdlib.h>
#include <cmath>
#include <utility>
#include "../classes.h"

/*Observation constructors and destructor*/
Observation::Observation(){/*{{{*/
	return;
}
/*}}}*/
Observation::Observation(double x_in,double y_in,int xi_in,int yi_in,int index_in,double value_in){/*{{{*/

	this->x      = x_in;
	this->y      = y_in;
	this->xi     = xi_in;
	this->yi     = yi_in;
	this->index  = index_in;
	this->value  = value_in;
	this->weight = 1.;

}
/*}}}*/
Observation::Observation(double x_in, double y_in,double value_in){
	this->x = x_in;
	this->y = y_in;
	this->value = value_in;

	this->xi     = 0;
	this->yi     = 0;
	this->index  = 0;
	this->weight = 0.;
}
Observation::~Observation(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* Observation::copy(void){/*{{{*/

	Observation* observation = new Observation(this->x,this->y,this->xi,this->yi,this->index,this->value);

	observation->weight = this->weight;

	return (Object*) observation;

}
/*}}}*/
void Observation::Echo(void){/*{{{*/

	_printf_("Observation\n");
	_printf_("   index : " << this->index << "\n");
	_printf_("   x     : " << this->x << "\n");
	_printf_("   y     : " << this->y << "\n");
	_printf_("   xi    : \n"); printbinary(this->xi); _printf_("\n");
	_printf_("   yi    : \n"); printbinary(this->yi); _printf_("\n");
	_printf_("   weight: " << this->weight << "\n");
	_printf_("   value : " << this->value << "\n");
}
/*}}}*/

/*Observations functions*/
void Observation::WriteXYObs(double* px,double* py,double* pobs){/*{{{*/
	*px   = this->x;
	*py   = this->y;
	*pobs = this->value;
}
/*}}}*/

/*Covertree*/
bool Observation::operator==(const Observation& ob) const{/*{{{*/
	return (ob.x == this->x && ob.y == this->y && ob.value == this->value);
}/*}}}*/
double Observation::distance(const Observation& ob) const{/*{{{*/
	return std::sqrt( (std::pow( (ob.x - this->x), 2 ) + std::pow((ob.y - this->y), 2) ));
}
/*}}}*/
void Observation::print(void) const{/*{{{*/

	_printf_("Observation\n");
	_printf_("   x     : " << this->x << "\n");
	_printf_("   y     : " << this->y << "\n");
	_printf_("   value : " << this->value << "\n");
}
/*}}}*/
void Observation::WriteXYObs(const Observation& ob, double* px, double* py, double* pobs){/*{{{*/
    *px   = ob.x;
    *py   = ob.y;
    *pobs = ob.value;
}/*}}}*/
