/*
 * GenericParam.h
 *
 *  Created on: Aug 29, 2012
 *      Author: utke
 */

#ifndef GENERICPARAM_H_
#define GENERICPARAM_H_

/*Headers:*/
#ifdef HAVE_CONFIG_H
        #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "./Param.h"
#include "../../shared/shared.h"

/**
 * here we have a class that holds an instance of P
 * but because it should live side by side with
 * the other instances derived from Param it - unfortunately -
 * inherits all the accessors that are useless in this context
 */
template <class P> class GenericParam: public Param{

	private:
		P myP;

	public:
		/*GenericParam constructors, destructors:*/
		GenericParam(int enumVal) { this->enum_type = enumVal; };
		GenericParam(P Pin, int enumVal) : myP(Pin) { this->enum_type = enumVal; };
		~GenericParam(){};

		/*Object virtual functions definitions:*/
		// unfortunately,  having to implement such a printer method implies that
		// any structured P must provide the friend << operator
		Param* copy() { return new GenericParam<P>(*this); };
		void   DeepEcho() {
			_printf_("GenericParam:\n");
			_printf_("   enum:  " << enum_type << " (" << EnumToStringx(enum_type) << ")\n");
			_printf_("   value: " << myP << "\n");;
		}
		void  Echo() {DeepEcho();};

		// the "copy"  has to implement the base class abstract function but I
		// would prefer to drop this not to hide a "new" in here because it does
		// not clarify  ownership of the newed up instance...  use the default
		// copy constructor instead
		void Marshall(MarshallHandle* marshallhandle){
			if(this->InstanceEnum()!=FemModelCommEnum){
				_printf_("   WARNING: parameter "<<EnumToStringx(this->enum_type)<<" is a GenericParam and cannot be marshalled\n");
			}
			/*Nothing for now*/
		}
		int  ObjectEnum(){return GenericParamEnum;};

		/*Param virtual function definitions:*/
		P&       GetParameterValue() { return myP;}
		const P& GetParameterValue()const { return myP;};
};

#endif /* GENERICPARAM_H_ */
