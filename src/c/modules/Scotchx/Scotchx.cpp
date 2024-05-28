/* Copyright 2004,2007,2008 ENSEIRB, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : gmap.c                                  **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a graph mapping software.       **/
/**                This module contains the main function. **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 05 jan 1993     **/
/**                                 to     12 may 1993     **/
/**                # Version 1.1  : from : 15 oct 1993     **/
/**                                 to     15 oct 1993     **/
/**                # Version 1.3  : from : 06 apr 1994     **/
/**                                 to     18 may 1994     **/
/**                # Version 2.0  : from : 06 jun 1994     **/
/**                                 to     17 nov 1994     **/
/**                # Version 2.1  : from : 07 apr 1995     **/
/**                                 to     18 jun 1995     **/
/**                # Version 3.0  : from : 01 jul 1995     **/
/**                                 to     02 oct 1995     **/
/**                # Version 3.1  : from : 07 nov 1995     **/
/**                                 to     25 apr 1996     **/
/**                # Version 3.2  : from : 24 sep 1996     **/
/**                                 to     26 may 1998     **/
/**                # Version 3.3  : from : 19 oct 1998     **/
/**                                 to   : 30 mar 1999     **/
/**                # Version 3.4  : from : 03 feb 2000     **/
/**                                 to   : 03 feb 2000     **/
/**                # Version 4.0  : from : 16 jan 2004     **/
/**                                 to   : 27 dec 2004     **/
/**                # Version 5.0  : from : 23 dec 2007     **/
/**                                 to   : 18 jun 2008     **/
/**                                                        **/
/************************************************************/

#include "./Scotchx.h"

int
gmapx (
  int                 (**pmaptabi)[2],
  int                 argcm,
  char                *argvm[],
  int                 nvi,
  int                 ne2i,
  int                 *ir,
  int                 *jc,
  int                 *vli,
  int                 *vwi,
  int                 *ewi,
  char                archtyp[],
  int                 nai,
  int                 *api)
{ 
#ifdef _HAVE_SCOTCH_ //only works if Scotch library has been compiled in.

  SCOTCH_Graph        grafdat;                    /* Source graph            */
  SCOTCH_Num          grafflag;                   /* Source graph properties */
  SCOTCH_Arch         archdat;                    /* Target architecture     */
  SCOTCH_Strat        stradat;                    /* Mapping strategy        */
  SCOTCH_Mapping      mapdat;                     /* Mapping data            */
  Clock               runtime[2];                 /* Timing variables        */
  SCOTCH_Num          nvert =0;
  SCOTCH_Num          nedge2=0;
  SCOTCH_Num*         adjir  =NULL;
  SCOTCH_Num*         adjjc  =NULL;
  SCOTCH_Num*         vertlab=NULL;
  SCOTCH_Num*         vertwgt=NULL;
  SCOTCH_Num*         edgewgt=NULL;
  SCOTCH_Num          napar =0;
  SCOTCH_Num*         archpar=NULL;
  SCOTCH_Num          (*maptab)[2]=NULL;
  int                 (*maptabi)[2]=NULL;
  int                 flagval;
  int                 i,j,k;

/*  reset static variables from previous runs (jes, 4/27/10)  */

  C_partNbr = 2;        /* Default number of parts     */
  C_paraNum = 0;        /* Number of parameters        */
  C_paraNbr = 0;        /* No parameters for mapping   */
  C_fileNum = 0;        /* Number of file in arg list  */
  C_fileNbr = 4;        /* Number of files for mapping */
  for (i=0; i<C_FILENBR; i++) {
    C_fileTab[i].name = "-";
    C_fileTab[i].pntr = NULL;
    if (i < 2)
      C_fileTab[i].mode = "r";
    else
      C_fileTab[i].mode = "w";
  }

/*  convert input arguments to scotch data types  */

  nvert =(SCOTCH_Num)nvi;
  nedge2=(SCOTCH_Num)ne2i;

  if (ir && jc) {
    adjir = (SCOTCH_Num *) malloc(nedge2*sizeof(SCOTCH_Num));
    for (i=0; i<nedge2; i++)
      adjir[i]=(SCOTCH_Num)ir[i];
    adjjc = (SCOTCH_Num *) malloc((nvert+1)*sizeof(SCOTCH_Num));
    for (i=0; i<(nvert+1); i++)
      adjjc[i]=(SCOTCH_Num)jc[i];
  }

  if (vli) {
    vertlab = (SCOTCH_Num *) malloc(nvert*sizeof(SCOTCH_Num));
    for (i=0; i<nvert; i++)
      vertlab[i]=(SCOTCH_Num)vli[i];
  }

  if (vwi) {
    vertwgt = (SCOTCH_Num *) malloc(nvert*sizeof(SCOTCH_Num));
    for (i=0; i<nvert; i++)
      vertwgt[i]=(SCOTCH_Num)vwi[i];
  }

  if (ewi) {
    edgewgt = (SCOTCH_Num *) malloc(nedge2*sizeof(SCOTCH_Num));
    for (i=0; i<nedge2; i++)
      edgewgt[i]=(SCOTCH_Num)ewi[i];
  }

  napar =(SCOTCH_Num)nai;

  if (api) {
    archpar = (SCOTCH_Num *) malloc(nai*sizeof(SCOTCH_Num));
    for (i=0; i<nai; i++)
      archpar[i]=(SCOTCH_Num)api[i];
  }

/*  start scotch processing  */

  flagval = C_FLAGNONE;                           /* Default behavior */
  i = strlen (argvm[0]);
  if ((i >= 5) && (strncmp (argvm[0] + i - 5, "gpart", 5) == 0)) {
    flagval |= C_FLAGPART;
    C_paraNbr = 1;                                /* One more parameter       */
    C_fileNbr = 3;                                /* One less file to provide */
    errorProg ("gpart");
  }
  else
    errorProg ("gmap");

  intRandResetStatic ();
  intRandInit ();

  if ((argcm >= 2) && (argvm[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (0);
  }

  grafflag = 0;                                   /* Use vertex and edge weights  */
  SCOTCH_stratInit (&stradat);                    /* Set default mapping strategy */

  for (i = 0; i < C_FILENBR; i ++)                /* Set default stream pointers */
    C_fileTab[i].pntr = (C_fileTab[i].mode[0] == 'r') ? stdin : stdout;
  for (i = 1; i < argcm; i ++) {                   /* Loop for all option codes                        */
    if ((argvm[i][0] != '-') || (argvm[i][1] == '\0') || (argvm[i][1] == '.')) { /* If found a file name */
      if (C_paraNum < C_paraNbr) {                /* If number of parameters not reached              */
        if ((C_partNbr = atoi (argvm[i])) < 1)     /* Get the number of parts                          */
          errorPrint ("main: invalid number of parts (\"%s\")", argvm[i]);
        C_paraNum ++;
        continue;                                 /* Process the other parameters */
      }
      if (C_fileNum < C_fileNbr)                  /* A file name has been given */
        C_fileTab[C_fileNum ++].name = argvm[i];
      else
        errorPrint ("main: too many file names given");
    }
    else {                                        /* If found an option name */
      switch (argvm[i][1]) {
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (0);
        case 'M' :
        case 'm' :
          SCOTCH_stratExit (&stradat);
          SCOTCH_stratInit (&stradat);
          SCOTCH_stratGraphMap (&stradat, &argvm[i][2]);
          break;
        case 'S' :
        case 's' :                                /* Source graph parameters */
          for (j = 2; argvm[i][j] != '\0'; j ++) {
            switch (argvm[i][j]) {
              case 'E' :
              case 'e' :
                grafflag |= 2;                    /* Do not load edge weights */
                break;
              case 'V' :
              case 'v' :
                grafflag |= 1;                    /* Do not load vertex weights */
                break;
              default :
                errorPrint ("main: invalid source graph option (\"%c\")", argvm[i][j]);
            }
          }
          break;
        case 'V' :
          fprintf (stderr, "gmap/gpart, version %s - F. Pellegrini\n", SCOTCH_VERSION);
          fprintf (stderr, "Copyright 2004,2007,2008 ENSEIRB, INRIA & CNRS, France\n");
          fprintf (stderr, "This software is libre/free software under CeCILL-C -- see the user's manual for more information\n");
          return  (0);
        case 'v' :                                /* Output control info */
          for (j = 2; argvm[i][j] != '\0'; j ++) {
            switch (argvm[i][j]) {
              case 'M' :
              case 'm' :
                flagval |= C_FLAGVERBMAP;
                break;
              case 'S' :
              case 's' :
                flagval |= C_FLAGVERBSTR;
                break;
              case 'T' :
              case 't' :
                flagval |= C_FLAGVERBTIM;
                break;
              default :
                errorPrint ("main: unprocessed parameter \"%c\" in \"%s\"", argvm[i][j], argvm[i]);
            }
          }
          break;
        default :
          errorPrint ("main: unprocessed option (\"%s\")", argvm[i]);
      }
    }
  }
  if ((flagval && C_FLAGPART) != 0) {              /* If program run as the partitioner            */
    C_fileTab[3].name = C_fileTab[2].name;        /* Put provided file names at their right place */
    C_fileTab[2].name = C_fileTab[1].name;
    C_fileTab[1].name = "-";
  }

  fileBlockOpen (C_fileTab, C_FILENBR);           /* Open all files */

  clockInit  (&runtime[0]);
  clockStart (&runtime[0]);

  SCOTCH_graphInit (&grafdat);                    /* Create graph structure         */
  SCOTCH_graphLoad (&grafdat, C_filepntrsrcinp, -1, grafflag, nvert, nedge2, adjir, adjjc, vertlab, vertwgt, edgewgt); /* Read source graph */

  SCOTCH_archInit (&archdat);                     /* Create architecture structure          */
  if ((flagval & C_FLAGPART) != 0)                /* If program run as the partitioner      */
    SCOTCH_archCmplt (&archdat, C_partNbr);       /* Create a complete graph of proper size */
  else
    SCOTCH_archLoad (&archdat, C_filepntrtgtinp, archtyp, napar, archpar); /* Read target architecture */

  clockStop  (&runtime[0]);                       /* Get input time */
  clockInit  (&runtime[1]);
  clockStart (&runtime[1]);

  SCOTCH_graphMapInit    (&grafdat, &mapdat, &archdat, NULL);
  SCOTCH_graphMapCompute (&grafdat, &mapdat, &stradat); /* Perform mapping */

  clockStop  (&runtime[1]);                       /* Get computation time */
  clockStart (&runtime[0]);

  SCOTCH_graphMapSave (&nvert, &maptab, &grafdat, &mapdat, C_filepntrmapout); /* Write mapping */

/*  convert output arguments from scotch data types  */

  if (maptab) {
    *pmaptabi = (int (*)[2]) malloc(nvert*2*sizeof(int));
    maptabi  = *pmaptabi;
    for (j=0; j<2; j++)
      for (i=0; i<nvert; i++)
          maptabi[i][j]=(int)maptab[i][j];
    free(maptab);
  }

  clockStop (&runtime[0]);                        /* Get output time */

  if (flagval && C_FLAGVERBSTR) {
    fprintf (C_filepntrlogout, "S\tStrat=");
    SCOTCH_stratSave (&stradat, C_filepntrlogout);
    putc ('\n', C_filepntrlogout);
  }
  if (flagval && C_FLAGVERBTIM) {
    fprintf (C_filepntrlogout, "T\tMapping\t\t%g\nT\tI/O\t\t%g\nT\tTotal\t\t%g\n",
             (double) clockVal (&runtime[1]),
             (double) clockVal (&runtime[0]),
             (double) clockVal (&runtime[0]) +
             (double) clockVal (&runtime[1]));
  }
  if (flagval && C_FLAGVERBMAP)
    SCOTCH_graphMapView (&grafdat, &mapdat, C_filepntrlogout);

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  SCOTCH_graphMapExit (&grafdat, &mapdat);
  SCOTCH_graphExit    (&grafdat);
  SCOTCH_stratExit    (&stradat);
  SCOTCH_archExit     (&archdat);

  if (archpar) free(archpar);
  if (edgewgt) free(edgewgt);
  if (vertwgt) free(vertwgt);
  if (vertlab) free(vertlab);
  if (adjjc) free(adjjc);
  if (adjir) free(adjir);

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (un)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  return (0);

#else //#ifdef _HAVE_SCOTCH_ 
  return(0);
#endif
}
