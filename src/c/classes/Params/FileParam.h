/*! \file FileParam.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _FILEPARAM_H_
#define _FILEPARAM_H_

/*Headers:*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Param.h"
#include "../../shared/shared.h"

class FileParam: public Param{

	private:
		FILE* value;

	public:
		/*FileParam constructors, destructors:*/
		FileParam();
		FileParam(int enum_type,FILE* fid);
		~FileParam();

		/*Object virtual functions definitions:*/
		Param* copy();
		void   DeepEcho();
		void   Echo();
		void   Marshall(MarshallHandle* marshallhandle);
		int    ObjectEnum(){return FileParamEnum;}

		/*Param virtual function definitions:*/
		void  GetParameterValue(FILE** pfid){*pfid=value;};
};
#endif  /* _INTPARAM_H */
