/*\file Print.h
 *\brief: print I/O for ISSM
 */

#ifndef _ISSM_PRINT_H_
#define _ISSM_PRINT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif 

/*Only include forward declaration to save compile time*/
#include <iosfwd>
#include <sstream>

using namespace std;
/*macros:*/
/* _printf_{{{*/
/* macro to print some string on all cpus */
#define _printf_(StreamArgs)\
  do{std::ostringstream aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy; \
	  aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy << StreamArgs; \
	  PrintfFunctionOnAllCpus(aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy.str());}while(0)
/*}}}*/
/* _printf0_ {{{*/
/* macro to print some string only on cpu 0 */
#define _printf0_(StreamArgs)\
  do{std::ostringstream aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy; \
	  aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy << StreamArgs; \
	  PrintfFunctionOnCpu0(aLoNgAnDwEiRdLoCaLnAmeFoRtHiSmAcRoOnLy.str());}while(0)
/*}}}*/

/*functions: */
int PrintfFunctionOnCpu0(const string & message);
int PrintfFunctionOnAllCpus(const string & message);

#endif	
