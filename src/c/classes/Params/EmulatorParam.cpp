/*!\file EmulatorParam.c
 * \brief: implementation of the EmulatorParam object
 */

/*header files: */
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "shared/shared.h"
/*}}}*/
#include <pybind11/numpy.h>
namespace py = pybind11;

/*EmulatorParam constructors and destructor*/
EmulatorParam::EmulatorParam(){/*{{{*/
	this->module_dir = NULL;
	this->pt_name    = NULL;
	this->py_name    = NULL;
	return;
}
/*}}}*/
EmulatorParam::EmulatorParam(int in_enum_type, char* module_dir_in, char* pt_name_in, char* py_name_in){/*{{{*/

	this->enum_type=in_enum_type;

	/*Copy path to emulator*/
	this->module_dir = xNew<char>(strlen(module_dir_in)+1);
	xMemCpy<char>(this->module_dir, module_dir_in,(strlen(module_dir_in)+1));
	this->pt_name = xNew<char>(strlen(pt_name_in)+1);
	xMemCpy<char>(this->pt_name, pt_name_in,(strlen(pt_name_in)+1));
	this->py_name = xNew<char>(strlen(py_name_in)+1);
	xMemCpy<char>(this->py_name, py_name_in,(strlen(py_name_in)+1));

	/*Activate interpretor*/
   this->guard = NULL;
	try{
		/*What if multi-emulator are activated?*/
		this->guard = new py::scoped_interpreter();

		py::module_ sys = py::module_::import("sys");
		sys.attr("path").attr("append")(this->module_dir);
		std::string pt_path(this->module_dir);
      if(!pt_path.empty() && pt_path.back() != '/'){
      	pt_path += "/";
      }
      pt_path += this->pt_name;
		std::string py_module_name(this->py_name);
      std::size_t dot = py_module_name.rfind('.');
      if(dot != std::string::npos){
	       py_module_name = py_module_name.substr(0,dot);
      }

		this->mod = py::module_::import(py_module_name.c_str());
		this->mod.attr("init_model")(pt_path.c_str(), "auto");
	}
	catch(const py::error_already_set& e){
		_printf_("EmulatorParam: Python exception in constructor\n");
		_printf_("   " << e.what() << "\n");
		this->mod = py::module_();
		delete this->guard;
		this->guard = NULL;
		throw;
	}	
}
/*}}}*/
EmulatorParam::~EmulatorParam(){/*{{{*/
	xDelete<char>(this->module_dir);
	xDelete<char>(this->pt_name);
	xDelete<char>(this->py_name);

   this->mod = py::module_();
   delete this->guard;
	this->guard = NULL;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* EmulatorParam::copy() {/*{{{*/

	_error_("not implemented");

}
/*}}}*/
void EmulatorParam::DeepEcho(void){/*{{{*/

	_error_("not implemented");

}
/*}}}*/
void EmulatorParam::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int  EmulatorParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void EmulatorParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	_error_("Not implemented yet");

}
/*}}}*/
int EmulatorParam::ObjectEnum(void){/*{{{*/

	return EmulatorParamEnum;

}
/*}}}*/

/*EmulatorParam virtual functions definitions: */
