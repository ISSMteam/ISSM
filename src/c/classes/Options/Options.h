#ifndef _CONTAINER_OPTIONS_H_
#define _CONTAINER_OPTIONS_H_

/*forward declarations */
class Option;
#include "../../datastructures/datastructures.h"
#include "./GenericOption.h"

class Options: public DataSet{

	public:

		/*constructors, destructors*/
		Options();
		~Options();

		/*numerics*/
		int     AddOption(Option* in_oobject);
		Option* GetOption(const char* name);

		template <class OptionType> void Get(OptionType* pvalue,const char* name){ /*{{{*/

			/*Get option*/
			GenericOption<OptionType>* genericoption=xDynamicCast<GenericOption<OptionType>*>(GetOption(name));

			/*If the pointer is not NULL, the option has been found*/
			if(genericoption){
				genericoption->Get(pvalue);
			}
			/*Else, the Option does not exist, no default provided*/
			else{
				_error_("option of name \"" << name << "\" not found, and no default value has been provided");
			}
		}
		/*}}}*/
		template <class OptionType> void Get(OptionType* pvalue,const char* name,OptionType default_value){ /*{{{*/

			/*Get option*/
			GenericOption<OptionType>* genericoption=xDynamicCast<GenericOption<OptionType>*>(GetOption(name));

			/*If the pointer is not NULL, the option has been found*/
			if(genericoption){
				genericoption->Get(pvalue);
			}
			else{
				if(GetOption(name)) _printf_("WARNING: option "<<name<<" found but fetched format not consistent, defaulting...\n");
				*pvalue=default_value;
			}
		}
		/*}}}*/

};

#endif //ifndef _INPUTS_H_

template <> inline void Options::Get(char** pvalue,const char* name,char* default_value){ /*{{{*/

	/*Get option*/
	GenericOption<char*>* genericoption=xDynamicCast<GenericOption<char*>*>(GetOption(name));

	/*If the pointer is not NULL, the option has been found*/
	if(genericoption){
		genericoption->Get(pvalue);
	}
	else{
		/*Make a copy*/
		int   stringsize=strlen(default_value)+1;
		char* outstring=xNew<char>(stringsize);
		xMemCpy<char>(outstring,default_value,stringsize);
		*pvalue=outstring;
	}
}
/*}}}*/
