/*\file PrintfFunction.c
 *\brief: this function is used by the _printf_ macro, to take into account the 
 *fact we may be running on a cluster. 
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdarg.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <iomanip>

#include "./Print.h"
#include "../Comm/IssmComm.h"
#include "../../String/sharedstring.h"
#include "../../MemOps/MemOps.h"

int PrintfFunctionOnCpu0(const string & message){

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	if(my_rank==0){
		#if _IS_MSYS2_
		printf("%s",message.c_str());
		#else
		ApiPrintf(message.c_str());
		#endif
	}
	return 1;
}
int PrintfFunctionOnAllCpus(const string & message){

	#if _IS_MSYS2_
	printf("%s",message.c_str());
	#else
	ApiPrintf(message.c_str());
	#endif

	return 1;
}
