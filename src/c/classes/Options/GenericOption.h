/*! \file GenericOption.h 
 *  \brief: header file for generic option object
 */

#ifndef _GENERIC_OPTION_
#define _GENERIC_OPTION_

/*Headers:{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <cstring>
#include "../../shared/shared.h"
#include "../../datastructures/datastructures.h"
#include "./Option.h"
/*}}}*/

template <class OptionType> 
class GenericOption: public Option {

	public:
		char       *name;
		OptionType  value;
		int         size[2];

		/*GenericOption constructors, destructors*/
		GenericOption(){ /*{{{*/
			name = NULL;
			size[0] = 0;
			size[1] = 0;
		} /*}}}*/
		~GenericOption(){ /*{{{*/
			if(name) xDelete<char>(name);
		} /*}}}*/

		/*Object virtual functions definitions:*/
		Object* copy(){/*{{{*/
			_error_("Not implemented yet");
		};/*}}}*/
		void DeepEcho(){ /*{{{*/
			char  indent[81]="";
			this->DeepEcho(indent);
		} /*}}}*/
		void DeepEcho(char* indent){ /*{{{*/

			char  cstr[81];

			_printf_(indent << "          name: \"" << name << "\"\n");
			_printf_(indent << "          size: " << size[0] <<"x"<<size[1]<< "\n");
			_printf_(indent << "         value: " << value << "\n");
		} /*}}}*/
		void Echo(){ /*{{{*/
			this->DeepEcho();
		} /*}}}*/
		int  Id(){/*{{{*/
			_error_("Not implemented yet");
		};/*}}}*/
		int  ObjectEnum(){/*{{{*/
			return GenericOptionEnum;
		};/*}}}*/

		/*GenericOption functions: */
		void  Get(OptionType* pvalue){/*{{{*/
			*pvalue=value; 
		};/*}}}*/
		char* Name(){/*{{{*/
			return name;
		};/*}}}*/
};

#if defined(_HAVE_AD_) && !defined(_WRAPPERS_) 
/*We hook off this specific specialization when not running ADOLC, otherwise we get a redeclaration with the next specialization*/
template <> inline void GenericOption<IssmPDouble*>::Get(IssmPDouble** pvalue){ /*{{{*/

	/*Copy vector*/
	int numel = this->size[0]*this->size[1];
	IssmPDouble* outvalue=xNew<IssmPDouble>(numel);
	for(int i=0;i<numel;i++) outvalue[i]=this->value[i];

	/*Assign output pointer*/
	*pvalue=outvalue;
} /*}}}*/
#endif
template <> inline void GenericOption<IssmDouble*>::Get(IssmDouble** pvalue){ /*{{{*/

	/*Copy vector*/
	int numel = this->size[0]*this->size[1];
	IssmDouble* outvalue=xNew<IssmDouble>(numel);
	for(int i=0;i<numel;i++) outvalue[i]=this->value[i];

	/*Assign output pointer*/
	*pvalue=outvalue;
} /*}}}*/
template <> inline void GenericOption<char*>::Get(char** pvalue){ /*{{{*/

	int   stringsize=strlen(this->value)+1;
	char* outstring=xNew<char>(stringsize);
	xMemCpy<char>(outstring,this->value,stringsize);

	*pvalue=outstring;
} 
/*}}}*/

/*Special destructors when there is a pointer*/
template <> inline GenericOption<char*>::~GenericOption(){ /*{{{*/
	if(name)  xDelete<char>(name);
	if(value) xDelete<char>(value);
} 
/*}}}*/

#endif  /* _OPTIONOBJECT_H */
