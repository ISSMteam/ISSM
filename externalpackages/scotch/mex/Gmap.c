
#define THISFUNCTION "Gmap"

/*  Gmap structures and prototypes  */

#include "mat.h"
#include "mex.h"
#include "matrix.h"

void GmapUsage( void );


int
Gmapx (
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
	int                 *api);

/******************************/
/*                            */
/* This is the main function. */
/*                            */
/******************************/

void mexFunction( int nlhs,
				  mxArray *plhs[],
				  int nrhs,
				  const mxArray *prhs[] )
{
	int     argcm;
	char    **argvm=NULL;
	int     nvert =0,nedge2=0,napar =0;
	mwIndex *ir=NULL,*jc=NULL;
	int     *adjir=NULL,*adjjc=NULL;
	double  *vld=NULL,*vwd=NULL,*ewd=NULL,*apd=NULL;
	int     *vli=NULL,*vwi=NULL,*ewi=NULL,*api=NULL;
	char    *archtyp=NULL;
	int     (*maptabi)[2]=NULL;
	double* maptabd=NULL;
	int     i,j,k,imi=0,imo=0,isi=0,ierr;

	/* Check for proper number of arguments */
   
	if      (nrhs == 0 && nlhs == 0) {
		GmapUsage();
		return;
	}
	else if (nrhs <  6 || nlhs >  1) {
		GmapUsage();
		mexErrMsgTxt(" ");
	}

/*  load matlab argument list and convert to integer (note that converting here
	and in the x-layer is inefficient, but it makes the x-layer more general)  */

	argvm = (char **) calloc(nrhs,sizeof(char *));

	if (!(mxIsNumeric(prhs[imi]) &&
		  (mxGetM(prhs[imi]) == 1 && mxGetN(prhs[imi]) == 1))) {
		argvm[isi] = (char *) calloc(4+1,sizeof(char));
		strcpy(argvm[isi],"Gmap");
		mexPrintf("%s -- Using \"%s\" entry point.\n",
				  THISFUNCTION,argvm[isi]);
		isi++;
	}
	else {
		argvm[isi] = (char *) calloc(5+1,sizeof(char));
		strcpy(argvm[isi],"Gpart");
		mexPrintf("%s -- Using \"%s\" entry point.\n",
				  THISFUNCTION,argvm[isi]);
		isi++;

		argvm[isi] = (char *) calloc(17,sizeof(char));
		sprintf(argvm[isi],"%d",(int)mxGetScalar(prhs[imi]));
		mexPrintf("%s -- Number of parts is %s.\n",
				  THISFUNCTION,argvm[isi]);
		isi++;
		imi++;
	}

	if (!mxIsNumeric(prhs[imi]) || (!mxIsEmpty(prhs[imi]) && !mxIsSparse(prhs[imi]))) {
		mexPrintf("%s -- Adjacency matrix must be numeric and sparse.\n",THISFUNCTION);
		mexErrMsgTxt(" ");
	}
	else {
		nvert =mxGetM(prhs[imi]);
		nedge2=mxGetNzmax(prhs[imi]);
		if (mxGetNzmax(prhs[imi])) {
			ir    =mxGetIr(prhs[imi]);
			adjir = (int *) malloc(mxGetNzmax(prhs[imi])*sizeof(int));
			for (i=0; i<mxGetNzmax(prhs[imi]); i++)
				adjir[i]=(int)ir[i];
		}
		if (mxGetN(prhs[imi])) {
			jc    =mxGetJc(prhs[imi]);
			adjjc = (int *) malloc((mxGetN(prhs[imi])+1)*sizeof(int));
			for (i=0; i<(mxGetN(prhs[imi])+1); i++)
				adjjc[i]=(int)jc[i];
		}
		mexPrintf("%s -- Adjacency matrix is of size %d by %d with %d non-zeroes.\n",
				  THISFUNCTION,mxGetM(prhs[imi]),mxGetN(prhs[imi]),mxGetNzmax(prhs[imi]));
	}
	imi++;

	if (!mxIsNumeric(prhs[imi])) {
		mexPrintf("%s -- Vertex label vector must be numeric.\n",THISFUNCTION);
		mexErrMsgTxt(" ");
	}
	else {
		if (mxGetM(prhs[imi])*mxGetN(prhs[imi])) {
			vld=mxGetPr(prhs[imi]);
			vli = (int *) malloc(mxGetM(prhs[imi])*mxGetN(prhs[imi])*sizeof(int));
			for (i=0; i<mxGetM(prhs[imi])*mxGetN(prhs[imi]); i++)
				vli[i]=(int)vld[i];
		}
		mexPrintf("%s -- Vertex label vector is of size %d by %d.\n",
				  THISFUNCTION,mxGetM(prhs[imi]),mxGetN(prhs[imi]));
	}
	imi++;

	if (!mxIsNumeric(prhs[imi])) {
		mexPrintf("%s -- Vertex weight vector must be numeric.\n",THISFUNCTION);
		mexErrMsgTxt(" ");
	}
	else {
		if (mxGetM(prhs[imi])*mxGetN(prhs[imi])) {
			vwd=mxGetPr(prhs[imi]);
			vwi = (int *) malloc(mxGetM(prhs[imi])*mxGetN(prhs[imi])*sizeof(int));
			for (i=0; i<mxGetM(prhs[imi])*mxGetN(prhs[imi]); i++)
				vwi[i]=(int)vwd[i];
		}
		mexPrintf("%s -- Vertex weight vector is of size %d by %d.\n",
				  THISFUNCTION,mxGetM(prhs[imi]),mxGetN(prhs[imi]));
	}
	imi++;

	if (!mxIsNumeric(prhs[imi]) || (!mxIsEmpty(prhs[imi]) && !mxIsSparse(prhs[imi]))) {
		mexPrintf("%s -- Edge weight matrix must be numeric and sparse.\n",THISFUNCTION);
		mexErrMsgTxt(" ");
	}
	else {
		if (mxGetM(prhs[imi])) {
			ewd=mxGetPr(prhs[imi]);
			ewi = (int *) malloc(mxGetM(prhs[imi])*sizeof(int));
			for (i=0; i<mxGetNzmax(prhs[imi]); i++)
				ewi[i]=(int)ewd[i];
		}
		mexPrintf("%s -- Edge weight matrix is of size %d by %d with %d non-zeroes.\n",
				  THISFUNCTION,mxGetM(prhs[imi]),mxGetN(prhs[imi]),mxGetNzmax(prhs[imi]));
	}
	imi++;

	if (!((strlen (argvm[0]) >= 5) &&
		  (strncmp (argvm[0] + strlen (argvm[0]) - 5, "Gpart", 5) == 0))) {
		if (!mxIsChar(prhs[imi])) {
			mexPrintf("%s -- Architecture type must be character.\n",THISFUNCTION);
			mexErrMsgTxt(" ");
		}
		else {
			if (mxGetM(prhs[imi])*mxGetN(prhs[imi])) {
				archtyp = (char *) calloc(mxGetM(prhs[imi])*mxGetN(prhs[imi])+1,sizeof(char));
				mxGetString(prhs[imi],archtyp,mxGetM(prhs[imi])*mxGetN(prhs[imi])+1);
			}
			mexPrintf("%s -- Architecture type is \"%s\".\n",
					  THISFUNCTION,archtyp);
		}
		imi++;

		if (!mxIsNumeric(prhs[imi])) {
			mexPrintf("%s -- Architecture parameter vector must be numeric.\n",THISFUNCTION);
			mexErrMsgTxt(" ");
		}
		else {
			napar =mxGetM(prhs[imi])*mxGetN(prhs[imi]);
			if (mxGetM(prhs[imi])*mxGetN(prhs[imi])) {
				apd=mxGetPr(prhs[imi]);
				api = (int *) malloc(mxGetM(prhs[imi])*mxGetN(prhs[imi])*sizeof(int));
				for (i=0; i<mxGetM(prhs[imi])*mxGetN(prhs[imi]); i++)
					api[i]=(int)apd[i];
			}
			mexPrintf("%s -- Architecture parameter vector is of size %d by %d.\n",
					  THISFUNCTION,mxGetM(prhs[imi]),mxGetN(prhs[imi]));
		}
		imi++;
	}

	while (imi < nrhs) {
		if (!mxIsChar(prhs[imi])) {
			mexPrintf("%s -- prhs[%d] must be character.\n",THISFUNCTION,imi);
			mexErrMsgTxt(" ");
		}
		else {
			argvm[isi] = (char *) calloc(mxGetM(prhs[imi])*mxGetN(prhs[imi])+1,sizeof(char));
			mxGetString(prhs[imi],argvm[isi],mxGetM(prhs[imi])*mxGetN(prhs[imi])+1);
		}
		isi++;
		imi++;
	}
	argcm=isi;
	mexPrintf("argcm=%d\n",argcm);
	for (i=0; i<argcm; i++)
		mexPrintf("argvm[%d]=\"%s\"\n",i,argvm[i]);

	/* Do the actual computations in a subroutine */

	mexPrintf("Gmapx:\n");
	ierr=Gmapx(&maptabi,
			   argcm,
			   argvm,
			   nvert,
			   nedge2,
			   adjir,
			   adjjc,
			   vli,
			   vwi,
			   ewi,
			   archtyp,
			   napar,
			   api);
	mexPrintf("%s -- Error %d from Gmapx.\n",THISFUNCTION,ierr);

/*  for (i=0; i<nvert; i++)
		mexPrintf("maptabi[%d][0]=%d, maptabi[%d][1]=%d\n",
			 	  i,maptabi[i][0],i,maptabi[i][1]); */

	/* Create matrices for the return arguments */

	if (maptabi) {
		plhs[imo]=mxCreateDoubleMatrix(nvert, 2, mxREAL);
		maptabd = mxGetPr(plhs[imo]);
		k=0;
		for (j=0; j<2; j++)
			for (i=0; i<nvert; i++)
				maptabd[k++]=(double)maptabi[i][j];
		free(maptabi);
	}
	else {
		plhs[imo]=mxCreateDoubleMatrix(0, 2, mxREAL);
	}
	imo++;

	if (api)     free(api);
	if (archtyp) free(archtyp);
	if (ewi)     free(ewi);
	if (vwi)     free(vwi);
	if (vli)     free(vli);
	if (adjjc)   free(adjjc);
	if (adjir)   free(adjir);
	if (argvm) {
		for (i=argcm-1; i>=0; i--)
			free(argvm[i]);
		free(argvm);
	}

	return;
}

void GmapUsage( void )
{

    mexPrintf("\n");
    mexPrintf("Usage: [maptab]=Gmap(adjmat,vertlb,vertwt,edgewt,archtyp,archpar,\n");
    mexPrintf("                         Scotch-specific parameters);\n");
    mexPrintf("\n");

    return;
}

