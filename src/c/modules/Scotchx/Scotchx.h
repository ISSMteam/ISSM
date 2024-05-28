/*!\file:  Scotchxx.h
 * \brief header file for Scotch partitioner
 */

#ifndef _SCOTCHX_H
#define _SCOTCHX_H

#undef __FUNCT__
#define __FUNCT__  "Scotchx"

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include <string.h>
#include "../../classes/classes.h"

#ifdef _HAVE_SCOTCH_ //only works if scotch library has been compiled in.

	#define GMAP

	#ifdef _PETSC_SCOTCH_
		#include "scotch_module.h"
		#include "scotch_common.h"
		#include "scotch_gmap.h"
	#endif

	#ifdef _HAVE_MPI_
		#include "ptscotch.h"
	#else
		#include "scotch.h"
	#endif

	/*
	**  The static variables.
	*/

	static int                  C_partNbr = 2;        /* Default number of parts     */
	static int                  C_paraNum = 0;        /* Number of parameters        */
	static int                  C_paraNbr = 0;        /* No parameters for mapping   */
	static int                  C_fileNum = 0;        /* Number of file in arg list  */
	static int                  C_fileNbr = 4;        /* Number of files for mapping */
	static File                 C_fileTab[C_FILENBR] = { /* File array               */
								  { "-", NULL, "r" },
								  { "-", NULL, "r" },
								  { "-", NULL, "w" },
								  { "-", NULL, "w" } };

	static const char *         C_usageList[] = {     /* Usage */
	  "gmap [<input source file> [<input target file> [<output mapping file> [<output log file>]]]] <options>",
	  "gpart [<nparts>] [<input source file> [<output mapping file> [<output log file>]]] <options>",
	  "  -h         : Display this help",
	  "  -m<strat>  : Set mapping strategy (see user's manual)",
	  "  -s<obj>    : Force unity weights on <obj>:",
	  "                 e  : edges",
	  "                 v  : vertices",
	  "  -V         : Print program version and copyright",
	  "  -v<verb>   : Set verbose mode to <verb>:",
	  "                 m  : mapping information",
	  "                 s  : strategy information",
	  "                 t  : timing information",
	  "",
	  "See default strategy with option '-vs'",
	  NULL };

#endif

/* local prototypes: */
int gmapx ( int (**pmaptabi)[2], int argcm, char *argvm[], int nvi, int ne2i, int *ir, int *jc, int *vli, int *vwi, int *ewi, char archtyp[], int nai, int *api);

#endif  /* _SCOTCHX_H */
