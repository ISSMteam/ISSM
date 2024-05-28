/*!\file DenseGslSolve.cpp
 * \brief: solve dense matrix system with GSL library
 */

/*Header files: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../shared/shared.h"
#include "../../classes/Params/GenericParam.h"
#include "../../classes/Params/Parameters.h"
#include "../adolc/adolcincludes.h"
#include "../codipack/codipackincludes.h"
#include "./gslincludes.h"

#ifdef _HAVE_GSL_
#include <gsl/gsl_linalg.h>
#endif

/*}}}*/

void DenseGslSolve(IssmPDouble** pX,IssmPDouble* A,IssmPDouble* B, int n){ /*{{{*/

	/*Intermediary: */
	IssmPDouble *X  = xNew<IssmPDouble>(n);
	SolverxSeq(X,A,B,n);

	/*allocate output pointers: */
	*pX=X;
}
/*}}}*/
void DenseGslSolve(IssmPDouble** px,IssmPDouble* Kff,int Kff_M,int Kff_N,IssmPDouble* pf,int pf_M,Parameters* parameters){ /*{{{*/

	/*Intermediary: */

	if(Kff_N!=pf_M)_error_("Right hand side vector of size " << pf_M << ", when matrix is of size " << Kff_M << "-" << Kff_N << " !");
	if(Kff_M!=Kff_N)_error_("Stiffness matrix should be square!");

	IssmPDouble *x  = xNew<IssmPDouble>(Kff_N);
	SolverxSeq(x,Kff,pf,Kff_N);

	/*allocate output pointers: */
	*px=x;
}
/*}}}*/
void SolverxSeq(IssmPDouble *X, IssmPDouble *A, IssmPDouble *B,int n){ /*{{{*/
#ifdef _HAVE_GSL_
	/*GSL Matrices and vectors: */
	int              s;
	gsl_matrix_view  a;
	gsl_vector_view  b,x;
	gsl_permutation *p = NULL;
//	for (int i=0; i<n*n; ++i) std::cout << "SolverxSeq A["<< i << "]=" << A[i] << std::endl;
//	for (int i=0; i<n; ++i) std::cout << "SolverxSeq b["<< i << "]=" << B[i] << std::endl;
	/*A will be modified by LU decomposition. Use copy*/
	double* Acopy = xNew<double>(n*n);
	xMemCpy(Acopy,A,n*n);

	/*Initialize gsl matrices and vectors: */
	a = gsl_matrix_view_array (Acopy,n,n);
	b = gsl_vector_view_array (B,n);
	x = gsl_vector_view_array (X,n);

	/*Run LU and solve: */
	p = gsl_permutation_alloc (n);
	gsl_linalg_LU_decomp (&a.matrix, p, &s);
	gsl_linalg_LU_solve (&a.matrix, p, &b.vector, &x.vector);

	/*Clean up and assign output pointer*/
	xDelete(Acopy);
	gsl_permutation_free(p);
#endif
}
/*}}}*/

#ifdef _HAVE_ADOLC_
int EDF_for_solverx(int n, IssmPDouble *x, int m, IssmPDouble *y){ /*{{{*/
	SolverxSeq(y,x, x+m*m, m); // x is where the matrix starts, x+m*m is where the right-hand side starts
	return 0;
} /*}}}*/
int EDF_fos_forward_for_solverx(int n, IssmPDouble *inVal, IssmPDouble *inDeriv, int m, IssmPDouble *outVal, IssmPDouble *outDeriv) { /*{{{*/
#ifdef _HAVE_GSL_
	//  for (int i=0; i<m*m; ++i) std::cout << "EDF_fos_forward_for_solverx A["<< i << "]=" << inVal[i] << std::endl;
	//  for (int i=0; i<m; ++i) std::cout << "EDF_fos_forward_for_solverx b["<< i << "]=" << inVal[i+m*m] << std::endl;
	// the matrix will be modified by LU decomposition. Use gsl_A copy
	double* Acopy = xNew<double>(m*m);
	xMemCpy(Acopy,inVal,m*m);
	/*Initialize gsl matrices and vectors: */
	gsl_matrix_view gsl_A = gsl_matrix_view_array (Acopy,m,m);
	gsl_vector_view gsl_b = gsl_vector_view_array (inVal+m*m,m); // the right hand side starts at address inVal+m*m
	gsl_permutation *perm_p = gsl_permutation_alloc (m);
	int  signPerm;
	// factorize
	gsl_linalg_LU_decomp (&gsl_A.matrix, perm_p, &signPerm);
	gsl_vector *gsl_x_p = gsl_vector_alloc (m);
	// solve for the value
	gsl_linalg_LU_solve (&gsl_A.matrix, perm_p, &gsl_b.vector, gsl_x_p);
	/*Copy result*/
	xMemCpy(outVal,gsl_vector_ptr(gsl_x_p,0),m);
	gsl_vector_free(gsl_x_p);
	//  for (int i=0; i<m; ++i) std::cout << "EDF_fos_forward_for_solverx x["<< i << "]=" << outVal[i] << std::endl;
	// solve for the derivatives acc. to A * dx = r  with r=db - dA * x
	// compute the RHS
	double* r=xNew<double>(m);
	for (int i=0; i<m; i++) {
		r[i]=inDeriv[m*m+i]; // this is db[i]
		for (int j=0;j<m; j++) {
			r[i]-=inDeriv[i*m+j]*outVal[j]; // this is dA[i][j]*x[j]
		}
	}
	gsl_vector_view gsl_r=gsl_vector_view_array(r,m);
	gsl_vector *gsl_dx_p = gsl_vector_alloc(m);
	gsl_linalg_LU_solve (&gsl_A.matrix, perm_p, &gsl_r.vector, gsl_dx_p);
	xMemCpy(outDeriv,gsl_vector_ptr(gsl_dx_p,0),m);
	gsl_vector_free(gsl_dx_p);
	xDelete(r);
	gsl_permutation_free(perm_p);
	xDelete(Acopy);
#endif
	return 0;
} /*}}}*/
int EDF_fov_forward_for_solverx(int n, IssmPDouble *inVal, int directionCount, IssmPDouble **inDeriv, int m, IssmPDouble *outVal, IssmPDouble **outDeriv) { /*{{{*/
#ifdef _HAVE_GSL_
	// the matrix will be modified by LU decomposition. Use gsl_A copy
	double* Acopy = xNew<double>(m*m);
	xMemCpy(Acopy,inVal,m*m);
	/*Initialize gsl matrices and vectors: */
	gsl_matrix_view gsl_A = gsl_matrix_view_array (Acopy,m,m);
	gsl_vector_view gsl_b = gsl_vector_view_array (inVal+m*m,m); // the right hand side starts at address inVal+m*m
	gsl_permutation *perm_p = gsl_permutation_alloc (m);
	int  signPerm;
	// factorize
	gsl_linalg_LU_decomp (&gsl_A.matrix, perm_p, &signPerm);
	gsl_vector *gsl_x_p = gsl_vector_alloc (m);
	// solve for the value
	gsl_linalg_LU_solve (&gsl_A.matrix, perm_p, &gsl_b.vector, gsl_x_p);
	/*Copy result*/
	xMemCpy(outVal,gsl_vector_ptr(gsl_x_p,0),m);
	gsl_vector_free(gsl_x_p);
	// solve for the derivatives acc. to A * dx = r  with r=db - dA * x
	double* r=xNew<double>(m);
	gsl_vector *gsl_dx_p = gsl_vector_alloc(m);
	for (int dir=0;dir<directionCount;++dir) {
		// compute the RHS
		for (int i=0; i<m; i++) {
			r[i]=inDeriv[m*m+i][dir]; // this is db[i]
			for (int j=0;j<m; j++) {
				r[i]-=inDeriv[i*m+j][dir]*outVal[j]; // this is dA[i][j]*x[j]
			}
		}
		gsl_vector_view gsl_r=gsl_vector_view_array(r,m);
		gsl_linalg_LU_solve (&gsl_A.matrix, perm_p, &gsl_r.vector, gsl_dx_p);
		// reuse r
		xMemCpy(r,gsl_vector_ptr(gsl_dx_p,0),m);
		for (int i=0; i<m; i++) {
			outDeriv[i][dir]=r[i];
		}
	}
	gsl_vector_free(gsl_dx_p);
	xDelete(r);
	gsl_permutation_free(perm_p);
	xDelete(Acopy);
#endif
	return 0;
}
/*}}}*/
int EDF_fos_reverse_for_solverx(int m, double *dp_U, int n, double *dp_Z, double* dp_x, double* dp_y) { /*{{{*/
	// copy to transpose the matrix
	double* transposed=xNew<double>(m*m);
	for (int i=0; i<m; ++i) for (int j=0; j<m; ++j) transposed[j*m+i]=dp_x[i*m+j];
	gsl_matrix_view aTransposed = gsl_matrix_view_array (transposed,m,m);
	// the adjoint of the solution is our right-hand side
	gsl_vector_view x_bar=gsl_vector_view_array(dp_U,m);
	// the last m elements of dp_Z representing the adjoint of the right-hand side we want to compute:
	gsl_vector_view b_bar=gsl_vector_view_array(dp_Z+m*m,m);
	gsl_permutation *perm_p = gsl_permutation_alloc (m);
	int permSign;
	gsl_linalg_LU_decomp (&aTransposed.matrix, perm_p, &permSign);
	gsl_linalg_LU_solve (&aTransposed.matrix, perm_p, &x_bar.vector, &b_bar.vector);
	// now do the adjoint of the matrix
	for (int i=0; i<m; ++i) for (int j=0; j<m; ++j) dp_Z[i*m+j]-=dp_Z[m*m+i]*dp_y[j];
	gsl_permutation_free(perm_p);
	xDelete(transposed);
	return 0;
}
/*}}}*/
int EDF_fov_reverse_for_solverx(int m, int p, double **dpp_U, int n, double **dpp_Z, double* dp_x, double* dp_y) { /*{{{*/
	// copy to transpose the matrix
	double* transposed=xNew<double>(m*m);
	for (int i=0; i<m; ++i) for (int j=0; j<m; ++j) transposed[j*m+i]=dp_x[i*m+j];
	gsl_matrix_view aTransposed = gsl_matrix_view_array (transposed,m,m);
	gsl_permutation *perm_p = gsl_permutation_alloc (m);
	int permSign;
	gsl_linalg_LU_decomp (&aTransposed.matrix, perm_p, &permSign);
	for (int weightsRowIndex=0;weightsRowIndex<p;++weightsRowIndex) {
		// the adjoint of the solution is our right-hand side
		gsl_vector_view x_bar=gsl_vector_view_array(dpp_U[weightsRowIndex],m);
		// the last m elements of dp_Z representing the adjoint of the right-hand side we want to compute:
		gsl_vector_view b_bar=gsl_vector_view_array(dpp_Z[weightsRowIndex]+m*m,m);
		gsl_linalg_LU_solve (&aTransposed.matrix, perm_p, &x_bar.vector, &b_bar.vector);
		// now do the adjoint of the matrix
		for (int i=0; i<m; ++i) for (int j=0; j<m; ++j) dpp_Z[weightsRowIndex][i*m+j]-=dpp_Z[weightsRowIndex][m*m+i]*dp_y[j];
	}
	gsl_permutation_free(perm_p);
	xDelete(transposed);
	return 0;
}
/*}}}*/
void DenseGslSolve(/*output*/ IssmDouble** px,/*stiffness matrix:*/ IssmDouble* Kff, int Kff_M, int Kff_N, /*right hand side load vector: */ IssmDouble* pf, int pf_M, Parameters* parameters){ /*{{{*/

	/*Intermediary: */

	if(Kff_N!=pf_M)_error_("Right hand side vector of size " << pf_M << ", when matrix is of size " << Kff_M << "-" << Kff_N << " !");
	if(Kff_M!=Kff_N)_error_("Stiffness matrix should be square!");

// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
	IssmDouble *x  = xNew<IssmDouble>(Kff_N,"t");
#else
	IssmDouble *x  = xNew<IssmDouble>(Kff_N);
#endif

	SolverxSeq(x,Kff,pf,Kff_N,parameters);

	/*allocate output pointers: */
	*px=x;
}
/*}}}*/
void SolverxSeq(IssmDouble *X,IssmDouble *A,IssmDouble *B,int n, Parameters* parameters){/*{{{*/
	// pack inputs to conform to the EDF-prescribed interface
// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
        IssmDouble*  adoubleEDFin=xNew<IssmDouble>(n*(n+1),"t");
#else
        IssmDouble*  adoubleEDFin=xNew<IssmDouble>(n*(n+1));
#endif  
	// packed inputs, i.e. matrix and right hand side
        for(int i=0; i<n*n;i++)adoubleEDFin[i]    =A[i];      // pack matrix
        for(int i=0; i<n;  i++)adoubleEDFin[i+n*n]=B[i];      // pack the right hand side
	// call the wrapped solver through the registry entry we retrieve from parameters
	call_ext_fct(xDynamicCast<GenericParam<Adolc_edf> * >(parameters->FindParamObject(AdolcParamEnum))->GetParameterValue().myEDF_for_solverx_p,
	             n*(n+1), adoubleEDFin,
	             n, X);
	// for(int i=0; i<n;  i++) {ADOLC_DUMP_MACRO(X[i]);}
	xDelete(adoubleEDFin);
}
/*}}}*/
#endif

#ifdef _HAVE_CODIPACK_
#if _CODIPACK_MAJOR_==2
using Tape = typename IssmDouble::Tape;
using AccessInterface = codi::VectorAccessInterface<typename Tape::Real, typename Tape::Identifier>;
void SolverxSeq_codi_b(Tape* tape,void* data_in, AccessInterface* ra) {/*{{{*/

	/*recast data_in and tape*/
	codi::ExternalFunctionUserData* data = (codi::ExternalFunctionUserData*)data_in;

  IssmDouble::Real* valueATrans;
  IssmDouble::Identifier* indexATrans;
  IssmDouble::Identifier* indexB;
  IssmDouble::Real* valueX;
  IssmDouble::Identifier* indexX;
  int n;

  data->getData(valueATrans);
  data->getData(indexATrans);
  data->getData(indexB);
  data->getData(valueX);
  data->getData(indexX);
  data->getData(n);

  // create the adjoint vector for x and reset the adjoint values on the tape
  IssmDouble::Gradient* adjX = xNew<IssmDouble::Gradient>(n);
  getVectorAdjoint(*tape, indexX, adjX, n);

  IssmDouble::Gradient* sol  = xNew<IssmDouble::Gradient>(n);
  SolverxSeq(sol, valueATrans, adjX, n);

  updateVectorAdjoint(*tape, indexB, sol, n);
  for(int i=0; i<n; ++i) {
    for (int j=0; j<n; ++j) {
      // we access the transposed matrix here because we stored the indices in a transposed way
      updateAdjoint(*tape, indexATrans[i*n+j], -sol[j]*valueX[i]);
    }
  }

  xDelete(sol);
  xDelete(adjX);
}
/*}}}*/
void SolverxSeq_codi_delete(Tape* tape,void* data_in) {/*{{{*/

	/*recast data_in*/
	codi::ExternalFunctionUserData* data = (codi::ExternalFunctionUserData*)data_in;

  IssmDouble::Real* valueATrans;
  IssmDouble::Identifier* indexATrans;
  IssmDouble::Identifier* indexB;
  IssmDouble::Real* valueX;
  IssmDouble::Identifier* indexX;
  int n;

  data->getData(valueATrans);
  data->getData(indexATrans);
  data->getData(indexB);
  data->getData(valueX);
  data->getData(indexX);
  data->getData(n);

  xDelete(valueATrans);
  xDelete(indexATrans);
  xDelete(indexB);
  xDelete(valueX);
  xDelete(indexX);

  delete data;
}
/*}}}*/
void SolverxSeq(IssmDouble *X,IssmDouble *A,IssmDouble *B,int n, Parameters* parameters){/*{{{*/
  IssmDouble::Tape& tape = IssmDouble::getTape();
  codi::ExternalFunctionUserData* dataHandler = NULL;

  if(tape.isActive()) {
    dataHandler = new codi::ExternalFunctionUserData();

    // create the index vector and the double data for A and B
    IssmDouble::Real* valueATrans = xNew<IssmDouble::Real>(n*n);
    IssmDouble::Identifier* indexATrans = xNew<IssmDouble::Identifier>(n*n);

    // read the data for matrix in a transposed fashion
	  for (int i=0; i<n; ++i) {
      for (int j=0; j<n; ++j) {
        getPrimalAndGradData(A[i*n+j], valueATrans[j*n+i], indexATrans[j*n+i]);
      }
    }

    // read the data from B (primal values are not required vor B
    IssmDouble::Identifier* indexB = xNew<IssmDouble::Identifier>(n);
    getVectorGradData(B, indexB, n);

    dataHandler->addData(valueATrans);
    dataHandler->addData(indexATrans);
    dataHandler->addData(indexB);
  }

  // unpack the primal values from the matrix and the vector
  IssmDouble::Real* valueA = xNew<IssmDouble::Real>(n*n);
  IssmDouble::Real* valueB = xNew<IssmDouble::Real>(n);
  // read the data from A and B
  getVectorPrimal(A, valueA, n*n);
  getVectorPrimal(B, valueB, n);

  // create the placeholder for X and solve the system
  IssmDouble::Real* valueX = xNew<IssmDouble::Real>(n);
  SolverxSeq(valueX, valueA, valueB, n);

  // pack the values into x
  setVectorPrimal(X, valueX, n);

  if(tape.isActive()) {
    // create the index vector X and register x as active variables
    IssmDouble::Identifier* indexX = xNew<IssmDouble::Identifier>(n);
    registerVector(X, indexX, n);

    dataHandler->addData(valueX);
    dataHandler->addData(indexX);

    // store other arguments
    dataHandler->addData(n);

	 tape.pushExternalFunction(codi::ExternalFunction<Tape>::create(&SolverxSeq_codi_b,(void*)dataHandler, &SolverxSeq_codi_delete));
  }
  else{
    // if the tape is active valueX is stored in the dataHandler and deleted in the reverse sweep
    xDelete(valueX);
  }

  xDelete(valueB);
  xDelete(valueA);
}
/*}}}*/
#elif _CODIPACK_MAJOR_==1
void SolverxSeq_codi_b(void* tape_in,void* data_in,void* ra) {/*{{{*/

	/*recast data_in and tape*/
	codi::DataStore* data = (codi::DataStore*)data_in;
	//IssmDouble::TapeType& tape = (IssmDouble::TapeType&)tape_in;
	IssmDouble::TapeType& tape = IssmDouble::getGlobalTape();

	IssmDouble::Real* valueATrans;
	IssmDouble::GradientData* indexATrans;
	IssmDouble::GradientData* indexB;
	IssmDouble::Real* valueX;
	IssmDouble::GradientData* indexX;
	int n;

	data->getData(valueATrans);
	data->getData(indexATrans);
	data->getData(indexB);
	data->getData(valueX);
	data->getData(indexX);
	data->getData(n);


	// create the adjoint vector for x and reset the adjoint values on the tape
	IssmDouble::GradientValue* adjX = xNew<IssmDouble::GradientValue>(n);
	getVectorAdjoint(tape, indexX, adjX, n);

	IssmDouble::GradientValue* sol  = xNew<IssmDouble::GradientValue>(n);
	SolverxSeq(sol, valueATrans, adjX, n);

	updateVectorAdjoint(tape, indexB, sol, n);
	for(int i=0; i<n; ++i) {
		for (int j=0; j<n; ++j) {
			// we access the transposed matrix here because we stored the indices in a transposed way
			updateAdjoint(tape, indexATrans[i*n+j], -sol[j]*valueX[i]);
		}
	}

	xDelete(sol);
	xDelete(adjX);
}
/*}}}*/
void SolverxSeq_codi_delete(void* tape_in,void* data_in) {/*{{{*/

	/*recast data_in*/
	codi::DataStore* data = (codi::DataStore*)data_in;

	IssmDouble::Real* valueATrans;
	IssmDouble::GradientData* indexATrans;
	IssmDouble::GradientData* indexB;
	IssmDouble::Real* valueX;
	IssmDouble::GradientData* indexX;
	int n;

	data->getData(valueATrans);
	data->getData(indexATrans);
	data->getData(indexB);
	data->getData(valueX);
	data->getData(indexX);
	data->getData(n);

	xDelete(valueATrans);
	xDelete(indexATrans);
	xDelete(indexB);
	xDelete(valueX);
	xDelete(indexX);

	delete data;
}
/*}}}*/
void SolverxSeq(IssmDouble *X,IssmDouble *A,IssmDouble *B,int n, Parameters* parameters){/*{{{*/
	IssmDouble::TapeType& tape = IssmDouble::getGlobalTape();
	codi::DataStore* dataHandler = NULL;

	if(tape.isActive()) {
		dataHandler = new codi::DataStore();

		// create the index vector and the double data for A and B
		IssmDouble::Real* valueATrans = xNew<IssmDouble::Real>(n*n);
		IssmDouble::GradientData* indexATrans = xNew<IssmDouble::GradientData>(n*n);

		// read the data for matrix in a transposed fashion
		for (int i=0; i<n; ++i) {
			for (int j=0; j<n; ++j) {
				getPrimalAndGradData(A[i*n+j], valueATrans[j*n+i], indexATrans[j*n+i]);
			}
		}

		// read the data from B (primal values are not required vor B
		IssmDouble::GradientData* indexB = xNew<IssmDouble::GradientData>(n);
		getVectorGradData(B, indexB, n);

		dataHandler->addData(valueATrans);
		dataHandler->addData(indexATrans);
		dataHandler->addData(indexB);
	}

	// unpack the primal values from the matrix and the vector
	IssmDouble::Real* valueA = xNew<IssmDouble::Real>(n*n);
	IssmDouble::Real* valueB = xNew<IssmDouble::Real>(n);
	// read the data from A and B
	getVectorPrimal(A, valueA, n*n);
	getVectorPrimal(B, valueB, n);

	// create the placeholder for X and solve the system
	IssmDouble::Real* valueX = xNew<IssmDouble::Real>(n);
	SolverxSeq(valueX, valueA, valueB, n);

	// pack the values into x
	setVectorPrimal(X, valueX, n);

	if(tape.isActive()) {
		// create the index vector X and register x as active variables
		IssmDouble::GradientData* indexX = xNew<IssmDouble::GradientData>(n);
		registerVector(X, indexX, n);

		dataHandler->addData(valueX);
		dataHandler->addData(indexX);

		// store other arguments
		dataHandler->addData(n);

		tape.pushExternalFunctionHandle(&SolverxSeq_codi_b, dataHandler, &SolverxSeq_codi_delete);
	}
	else{
		// if the tape is active valueX is stored in the dataHandler and deleted in the reverse sweep
		xDelete(valueX);
	}

	xDelete(valueB);
	xDelete(valueA);
}
/*}}}*/
#else
#error "_CODIPACK_MAJOR_ not supported"
#endif
void DenseGslSolve(/*output*/ IssmDouble** px,/*stiffness matrix:*/ IssmDouble* Kff, int Kff_M, int Kff_N, /*right hand side load vector: */ IssmDouble* pf, int pf_M, Parameters* parameters){ /*{{{*/

	/*Intermediary: */

	if(Kff_N!=pf_M)_error_("Right hand side vector of size " << pf_M << ", when matrix is of size " << Kff_M << "-" << Kff_N << " !");
	if(Kff_M!=Kff_N)_error_("Stiffness matrix should be square!");

	IssmDouble *x  = xNew<IssmDouble>(Kff_N,"t");

	SolverxSeq(x,Kff,pf,Kff_N,parameters);

	/*allocate output pointers: */
	*px=x;
}
/*}}}*/

#endif
