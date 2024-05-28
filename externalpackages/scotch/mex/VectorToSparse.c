
#define THISFUNCTION "VectorToSparse"

#include <stdio.h>
#include <string.h>    /*  strcasecmp  */
#include <time.h>      /*  clock,time,difftime  */
#include "mex.h"


/* Input Arguments */

#define    IR_IN        prhs[0]
#define    JC_IN        prhs[1]
#define    PR_IN        prhs[2]
#define    M_IN         prhs[3]
#define    N_IN         prhs[4]

/* Output Arguments */

#define    A_OUT        plhs[0]


void VectorToSparseUsage( void );


void mexFunction( int nlhs,
				  mxArray *plhs[],
				  int nrhs,
				  const mxArray *prhs[] )
{
	int     min,nin;
	mwIndex *ir =NULL,*jc =NULL;
	double  *ird=NULL,*jcd=NULL,*pr=NULL,*prd=NULL;
	int     i,mrow;

	/* Check for proper number of arguments */
   
	if      (nrhs == 0 && nlhs == 0) {
		VectorToSparseUsage();
		return;
	}
	else if (nrhs <  2 || nlhs != 1) {
		VectorToSparseUsage();
		mexErrMsgTxt(" ");
	}

	/* Create matrices for the return arguments */

	if (!mxIsNumeric(IR_IN) || !mxIsNumeric(JC_IN)) {
		mexPrintf("%s -- Input matrices IR and JC must be numeric.\n",THISFUNCTION);
		mexErrMsgTxt(" ");
	}

	mrow = 0;
	ird = mxGetPr(IR_IN);
	for (i=0; i<mxGetM(IR_IN)*mxGetN(IR_IN); i++)
		if ((int)ird[i]+1 > mrow)
			mrow=(int)ird[i]+1;

	if (nrhs >= 4 && mxIsNumeric(M_IN) && !mxIsEmpty(M_IN))
		min = mxGetScalar(M_IN);
	else {
		min = mrow;
	}

	if (mrow > min) {
		mexPrintf("%s -- Number of rows specified by M (%d) and IR (%d) is inconsistent.\n",
				  THISFUNCTION,min,mrow);
		mexErrMsgTxt(" ");
	}

	if (nrhs >= 5 && mxIsNumeric(N_IN) && !mxIsEmpty(N_IN))
		nin = mxGetScalar(N_IN);
	else
		nin = mxGetM(JC_IN)*mxGetN(JC_IN)-1;

	if (mxGetM(JC_IN)*mxGetN(JC_IN)-1 != nin) {
		mexPrintf("%s -- Number of columns specified by N (%d) and JC (%d) is inconsistent.\n",
				  THISFUNCTION,nin,mxGetM(JC_IN)*mxGetN(JC_IN)-1);
		mexErrMsgTxt(" ");
	}

	A_OUT = mxCreateSparse(min, nin, mxGetM(IR_IN)*mxGetN(IR_IN), mxREAL);
	if (mxGetM(IR_IN)*mxGetN(IR_IN)) {
		ird = mxGetPr(IR_IN);
		ir  = mxGetIr(A_OUT);
		for (i=0; i<mxGetM(IR_IN)*mxGetN(IR_IN); i++)
			ir[i]=(mwIndex)ird[i];
	}
	if (mxGetM(JC_IN)*mxGetN(JC_IN)) {
		jcd = mxGetPr(JC_IN);
		jc  = mxGetJc(A_OUT);
		for (i=0; i<mxGetM(JC_IN)*mxGetN(JC_IN); i++)
			jc[i]=(mwIndex)jcd[i];
	}

	if (nrhs >= 3 && mxIsNumeric(PR_IN) && !mxIsEmpty(PR_IN)) {
		if (mxGetM(PR_IN)*mxGetN(PR_IN) != mxGetM(IR_IN)*mxGetN(IR_IN)) {
			mexPrintf("%s -- Number of terms specified by IR (%d) and PR (%d) is inconsistent.\n",
					  THISFUNCTION,mxGetM(IR_IN)*mxGetN(IR_IN),mxGetM(PR_IN)*mxGetN(PR_IN));
			mexErrMsgTxt(" ");
		}

		if (mxGetM(PR_IN)*mxGetN(PR_IN)) {
			prd = mxGetPr(PR_IN);
			pr  = mxGetPr(A_OUT);
			for (i=0; i<mxGetM(PR_IN)*mxGetN(PR_IN); i++)
				pr[i]=prd[i];
		}
	}

	else {
		mexPrintf("%s -- Populating sparse matrix terms with ones.\n",THISFUNCTION);
		mexWarnMsgTxt(" ");

		if (mxGetM(IR_IN)*mxGetN(IR_IN)) {
			pr  = mxGetPr(A_OUT);
			for (i=0; i<mxGetM(IR_IN)*mxGetN(IR_IN); i++)
				pr[i]=1.;
		}
	}

	mexPrintf("%s -- Output matrix is of size %d by %d with %d non-zeroes.\n",
			  THISFUNCTION,min,nin,mxGetM(IR_IN)*mxGetN(IR_IN));

	return;
}

void VectorToSparseUsage( void )
{

    mexPrintf("\n");
    mexPrintf("Usage: [a]=VectorToSparse(ir,jc,pr,m,n);\n");
    mexPrintf("\n");

    return;
}

