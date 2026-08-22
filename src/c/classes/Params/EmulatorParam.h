/*! \file EmulatorParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _EMULATORPARAM_H_
#define _EMULATORPARAM_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "./Param.h"
#include "../../shared/shared.h"
/*}}}*/
#include <pybind11/embed.h>
namespace py = pybind11;
class EmulatorParam: public Param{

	public:
		char*                module_dir;
		char*                   pt_name;
		char*                   py_name;
		py::scoped_interpreter*   guard;
		py::module_                 mod;

		/*EmulatorParam constructors, destructors: {{{*/
		EmulatorParam();
		EmulatorParam(int enum_type, char* module_dir_in, char* pt_name_in, char* py_name_in);
		EmulatorParam(int enum_type, char* module_dir_in, char* pt_name_in, char* py_name_in, int* edge_src, int nsrc, int* edge_dst, int ndst, int num_nodes);
		~EmulatorParam();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Param* copy();
		void   DeepEcho();
		void   Echo();
		int    Id(); 
		void   Marshall(MarshallHandle* marshallhandle);
		int    ObjectEnum();
		/*}}}*/
		/*Param virtual function definitions: {{{*/
		/*}}}*/
};
#endif
