/*!\file:  ExternalResult.h
 * \brief abstract class for ExternalResult object
 */ 

#ifndef _EXTERNALRESULT_H_
#define _EXTERNALRESULT_H_

/*Headers:*/
/*{{{*/

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../datastructures/datastructures.h"
#include "../Node.h"
/*}}}*/

class ExternalResult: public Object{

	public: 

		virtual         ~ExternalResult(){};
		virtual int    GetResultEnum(void)=0;
		virtual char*  GetResultName(void)=0;
		virtual int    GetStep(void)=0;
		virtual double GetValue(void)=0;
		virtual void   WriteData(FILE* fid,bool io_gather)=0;
};
#endif
