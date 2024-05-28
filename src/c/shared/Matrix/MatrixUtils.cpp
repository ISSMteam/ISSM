/*!\file: MatrixUtils
 * \brief triple multiply
 */ 

/*Headers*/
/*{{{*/
#include <stdio.h>
#include <sys/types.h>
#include <math.h>
#include <float.h>    /*  DBL_EPSILON  */
#include <cstdarg>
#include <iostream>

#include "./matrix.h"
#include "../Exceptions/exceptions.h"
#include "../MemOps/MemOps.h"
#include "../io/io.h"
/*}}}*/

int TripleMultiply(IssmDouble* a, int nrowa, int ncola, int itrna, IssmDouble* b, int nrowb, int ncolb, int itrnb, IssmDouble* c, int nrowc, int ncolc, int itrnc, IssmDouble* d, int iaddd){/*{{{*/
	/*TripleMultiply    Perform triple matrix product a*b*c+d.*/

	int         idima,idimb,idimc,idimd;
	IssmDouble  dtemp_static[600];
	IssmDouble* dtemp_dynamic = NULL;
	IssmDouble* dtemp         = NULL;
	_assert_(a && b && c && d);

	/*  set up dimensions for triple product  */
	if (!itrna){
		idima=nrowa;
		idimb=ncola;
	}
	else{
		idima=ncola;
		idimb=nrowa;
	}

	if (!itrnb){
		if (nrowb != idimb) _error_("Matrix A and B inner vectors not equal size.");
		idimc=ncolb;
	}
	else{
		if (ncolb != idimb) _error_("Matrix A and B inner vectors not equal size.");
		idimc=nrowb;
	}

	if (!itrnc) {
		if (nrowc != idimc) _error_("Matrix B and C inner vectors not equal size.");
		idimd=ncolc;
	}
	else{
		if (ncolc != idimc) _error_("Matrix B and C inner vectors not equal size.");
		idimd=nrowc;
	}

	/*Depending on the size of incoming matrices, we might need to use a dynamic allocation*/
	if(idima*idimc>600){
		dtemp_dynamic = xNew<IssmDouble>(idima*idimc);
		dtemp         = dtemp_dynamic;
	}
	else{
		dtemp = &dtemp_static[0];
	}

	/*  perform the matrix triple product in the order that minimizes the
		 number of multiplies and the temporary space used, noting that
		 (a*b)*c requires ac(b+d) multiplies and ac IssmDoubles, and a*(b*c)
		 requires bd(a+c) multiplies and bd IssmDoubles (both are the same for
		 a symmetric triple product)  */

	/*  multiply (a*b)*c+d  */
	if (idima*idimc*(idimb+idimd) <= idimb*idimd*(idima+idimc)) {
		MatrixMultiply(a,nrowa,ncola,itrna,b,nrowb,ncolb,itrnb,dtemp,0);
		MatrixMultiply(dtemp,idima,idimc,0,c,nrowc,ncolc,itrnc,d,iaddd);
	}

	/*  multiply a*(b*c)+d  */
	else{
		MatrixMultiply(b,nrowb,ncolb,itrnb,c,nrowc,ncolc,itrnc,dtemp,0);
		MatrixMultiply(a,nrowa,ncola,itrna,dtemp,idimb,idimd,0,d,iaddd);
	}

	/*Cleanup and return*/
	xDelete<IssmDouble>(dtemp_dynamic);
	return 1;
}/*}}}*/
int MatrixMultiply(IssmDouble* a, int nrowa, int ncola, int itrna, IssmDouble* b, int nrowb, int ncolb, int itrnb, IssmDouble* c, int iaddc ){/*{{{*/
	/*MatrixMultiply    Perform matrix multiplication a*b+c.*/
	int noerr=1;
	int i,j,k,ipta,iptb,iptc;
	int nrowc,ncolc,iinca,jinca,iincb,jincb,ntrma,ntrmb,nterm;

	_assert_(a && b && c);

	/*  set up dimensions and increments for matrix a  */
	if (!itrna) {
		nrowc=nrowa;
		ntrma=ncola;
		iinca=ncola;
		jinca=1;
	}
	else {
		nrowc=ncola;
		ntrma=nrowa;
		iinca=1;
		jinca=ncola;
	}

	/*  set up dimensions and increments for matrix b  */
	if (!itrnb) {
		ncolc=ncolb;
		ntrmb=nrowb;
		iincb=ncolb;
		jincb=1;
	}
	else {
		ncolc=nrowb;
		ntrmb=ncolb;
		iincb=1;
		jincb=ncolb;
	}

	if (ntrma != ntrmb) _error_("Matrix A and B inner vectors not equal size");

	nterm=ntrma;

	/*  zero matrix c, if not being added to product  */
	if (!iaddc) for (i=0;i<nrowc*ncolc;i++) c[i]=0.;

	/*  perform the matrix multiplication  */
	iptc=0;
	for (i=0; i<nrowc; i++){
		for (j=0; j<ncolc; j++){
			ipta=i*iinca;
			iptb=j*jincb;

			for (k=0; k<nterm; k++){
				c[iptc]+=a[ipta]*b[iptb];
				ipta+=jinca;
				iptb+=iincb;
			}
			iptc++;
		}
	}

	return noerr;
}/*}}}*/
int MatrixInverse( IssmDouble* a, int ndim, int nrow, IssmDouble* b, int nvec, IssmDouble* pdet ){/*{{{*/
	/* MatrixInverse    Perform matrix inversion and linear equation solution.

		This function uses Gaussian elimination on the original matrix
		augmented by an identity matrix of the same size to calculate
		the inverse (see for example, "Modern Methods of Engineering
		Computation", Sec. 6.4).  By noting how the matrices are
		unpopulated and repopulated, the calculation may be done in place.

		Gaussian elimination is inherently inefficient, and so this is
		intended for small matrices.  */
	int noerr=1;
	int i,j,k,ipt,jpt,irow,icol,ipiv,ncol;
	int *pivrc1,*pivrc2,*pindx;
	IssmDouble pivot,det,dtemp;

	if (!b && nvec) {
		_error_("No right-hand side for nvec=" << nvec << ".");
		noerr=0;
		return noerr;
	}

	/*In debugging mode, check that we are not dealing with simple matrices*/
	_assert_(!(ndim==2 && nrow==2));
	_assert_(!(ndim==3 && nrow==3));

	/*  initialize local variables and arrays  */

	ncol=nrow;
	det=1.;
	pivrc1 =xNew<int>(nrow);
	pivrc2 =xNew<int>(nrow);
	pindx =xNew<int>(nrow);

	/*  loop over the rows/columns of the matrix  */

	for (i=0; i<nrow; i++) {

		/*  search for pivot, finding the term with the greatest magnitude
			 in the rows/columns not yet used  */

		pivot=0.;
		for (j=0; j<nrow; j++)
		 if (!pindx[j])
		  for (k=0; k<ncol; k++)
			if (!pindx[k])
			 if (fabs(a[j*ndim+k]) > fabs(pivot)) {
				 irow=j;
				 icol=k;
				 pivot=a[j*ndim+k];
			 }

		if (fabs(pivot) < DBL_EPSILON) {
			xDelete<int>(pivrc1);
			xDelete<int>(pivrc2);
			xDelete<int>(pindx);
			_error_("Pivot " << pivot << " less than machine epsilon");
			noerr=0;
			return noerr;
		}

		pivrc1[i]=irow;
		pivrc2[i]=icol;

		ipiv=icol;
		pindx[ipiv]++;

		//		_printf_("pivot for i=" << i << ": irow=" << irow << ", icol=" << icol  << ", pindx[" << ipiv << "]=" << pindx[ipiv] << "\n\n");

		/*  switch rows to put pivot element on diagonal, noting that the
			 column stays the same and the determinant changes sign  */

		if (irow != icol) {
			//			_printf_("row switch for i=" << i << ": irow=" << irow << ", icol=" << icol << "\n\n");

			ipt=irow*ndim;
			jpt=icol*ndim;
			for (k=0; k<ncol; k++) {
				dtemp   =a[ipt+k];
				a[ipt+k]=a[jpt+k];
				a[jpt+k]=dtemp;
			}

			ipt=irow*nvec;
			jpt=icol*nvec;
			for (k=0; k<nvec; k++) {
				dtemp   =b[ipt+k];
				b[ipt+k]=b[jpt+k];
				b[jpt+k]=dtemp;
			}

			det=-det;
		}

		/*  divide pivot row by pivot element, noting that the original
			 matrix will have 1 on the diagonal, which will be discarded,
			 and the augmented matrix will start with 1 from the identity
			 matrix and then have 1/pivot, which is part of the inverse.  */

		a[ipiv*ndim+ipiv]=1.;

		ipt=ipiv*ndim;
		for (k=0; k<ncol; k++)
		 a[ipt+k]/=pivot;

		ipt=ipiv*nvec;
		for (k=0; k<nvec; k++)
		 b[ipt+k]/=pivot;

		/*  reduce non-pivot rows such that they will have 0 in the pivot
			 column, which will be discarded, and the augmented matrix will
			 start with 0 from the identity matrix and then have non-zero
			 in the corresponding column, which is part of the inverse.
			 only one column of the augmented matrix is populated at a time,
			 which corresponds to the only column of the original matrix
			 being zeroed, so that the inverse may be done in place.  */

		for (j=0; j<nrow; j++) {
			if (j == ipiv) continue;

			dtemp=a[j*ndim+ipiv];
			a[j*ndim+ipiv]=0.;

			if (fabs(dtemp) > DBL_EPSILON) {
				ipt=j   *ndim;
				jpt=ipiv*ndim;
				for (k=0; k<ncol; k++)
				 a[ipt+k]-=dtemp*a[jpt+k];

				ipt=j   *nvec;
				jpt=ipiv*nvec;
				for (k=0; k<nvec; k++)
				 b[ipt+k]-=dtemp*b[jpt+k];
			}
		}

		/*  for a diagonal matrix, the determinant is the product of the
			 diagonal terms, and so it may be accumulated from the pivots,
			 noting that switching rows changes the sign as above  */

		det*=pivot;
	}

	/*  switch columns back in reverse order, noting that a row switch
		 in the original matrix corresponds to a column switch in the
		 inverse matrix  */

	for (i=0; i<nrow; i++) {
		j=(nrow-1)-i;

		if (pivrc1[j] != pivrc2[j]) {
			irow=pivrc1[j];
			icol=pivrc2[j];

			//			_printf_("column switch back for j=" << j << ": irow=" << irow << ", icol=" << icol << "\n\n");

			ipt=0;
			for (k=0; k<nrow; k++) {
				dtemp      =a[ipt+irow];
				a[ipt+irow]=a[ipt+icol];
				a[ipt+icol]=dtemp;
				ipt+=ndim;
			}
		}
	}

	if (pdet) *pdet=det;
	xDelete<int>(pivrc1);
	xDelete<int>(pivrc2);
	xDelete<int>(pindx);
	return noerr;
}/*}}}*/

void Matrix2x2Determinant(IssmDouble* Adet,IssmDouble* A){/*{{{*/
	/*Compute determinant of a 2x2 matrix*/

	/*det = a*d - c*b*/
	*Adet= A[0]*A[3]-A[2]*A[1];
}
/*}}}*/
void Matrix2x2Invert(IssmDouble* Ainv,IssmDouble* A){/*{{{*/

	/*Intermediaries*/
	IssmDouble det,det_reciprocal;

	/*Compute determinant*/
	Matrix2x2Determinant(&det,A);
	if (fabs(det) < DBL_EPSILON) _error_("Determinant smaller than machine epsilon");

	/*Multiplication is faster than divsion, so we multiply by the reciprocal*/
	det_reciprocal = 1./det;  

	/*Compute invert*/
	Ainv[0]=   A[3]*det_reciprocal; /* =  d/det */
	Ainv[1]= - A[1]*det_reciprocal; /* = -b/det */
	Ainv[2]= - A[2]*det_reciprocal; /* = -c/det */
	Ainv[3]=   A[0]*det_reciprocal; /* =  a/det */

}/*}}}*/
void Matrix2x2Eigen(IssmDouble* plambda1,IssmDouble* plambda2,IssmDouble* pvx, IssmDouble* pvy,IssmDouble a11, IssmDouble a21,IssmDouble a22){/*{{{*/
	/*From symetric matrix (a11,a21;a21,a22), get eigen values lambda1 and lambda2 and one eigen vector v*/

	/*Output*/
	IssmDouble lambda1,lambda2;
	IssmDouble vx,vy;

	/*To get the eigen values, we must solve the following equation:
	 *     | a11 - lambda    a21        |
	 * det |                            | = 0
	 *     | a21             a22-lambda |
	 *
	 * We have to solve the following polynom:
	 *  lamda^2 + ( -a11 -a22)*lambda + (a11*a22-a21*a21) = 0*/

	/*Compute polynom determinant*/
	IssmDouble b=-a11-a22;
	IssmDouble delta=b*b - 4*(a11*a22-a21*a21);

	/*Compute norm of M to avoid round off errors*/
	IssmDouble normM=a11*a11 + a22*a22 + a21*a21;

	/*1: normM too small: eigen values = 0*/
	if(normM<1.e-30){
		lambda1=0.;
		lambda2=0.;
		vx=1.;
		vy=0.;
	}
	/*2: delta is small -> double root*/
	else if (delta < 1.e-5*normM){
		lambda1=-b/2.;
		lambda2=-b/2.;
		vx=1.;
		vy=0.;
	}
	/*3: general case -> two roots*/
	else{
		delta   = sqrt(delta);
		lambda1 = (-b-delta)/2.;
		lambda2 = (-b+delta)/2.;

		/*Now, one must find the eigen vectors. For that we use the following property of the inner product
		 *    <Ax,y> = <x,tAy>
		 * Here, M'(M-lambda*Id) is symmetrical, which gives:
		 *    \forall (x,y)\in R²xR² <M'x,y> = <M'y,x>
		 * And we have the following:
		 *    if y\in Ker(M'), \forall x\in R² <M'x,y> = <x,M'y> = 0
		 * We have shown that
		 *    Im(M') \perp Ker(M')
		 *
		 * To find the eigen vectors of M, we only have to find two vectors
		 * of the image of M' and take their perpendicular as long as they are
		 * not 0.
		 * To do that, we take the images (1,0) and (0,1):
		 *  x1 = (a11 - lambda)      x2 = a21
		 *  y1 = a21                 y2 = (a22-lambda)
		 *
		 * We take the vector that has the larger norm and take its perpendicular.*/

		IssmDouble norm1 = (a11-lambda1)*(a11-lambda1) + a21*a21; 
		IssmDouble norm2 = a21*a21 + (a22-lambda1)*(a22-lambda1);

		if(norm2<norm1){
			norm1=sqrt(norm1);
			vx = - a21/norm1;
			vy = (a11-lambda1)/norm1;
		}
		else{
			norm2=sqrt(norm2);
			vx = - (a22-lambda1)/norm2;
			vy = a21/norm2;
		}
	}

	/*Assign output*/
	*plambda1 = lambda1;
	*plambda2 = lambda2;
	if(pvx) *pvx = vx;
	if(pvy) *pvy = vy;

}/*}}}*/

void Matrix3x3Determinant(IssmDouble* Adet,IssmDouble* A){/*{{{*/
	/*Compute determinant of a 3x3 matrix*/

	/*det = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)*/
	*Adet= A[0]*A[4]*A[8]-A[0]*A[5]*A[7]-A[3]*A[1]*A[8]+A[3]*A[2]*A[7]+A[6]*A[1]*A[5]-A[6]*A[2]*A[4];
}
/*}}}*/
IssmDouble Matrix3x3Determinant(IssmDouble a1,IssmDouble a2,IssmDouble a3, IssmDouble b1,IssmDouble b2,IssmDouble b3, IssmDouble c1,IssmDouble c2,IssmDouble c3){/*{{{*/
	/*Compute determinant of a 3x3 matrix*/

	/*det = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
	 * a b c   a1 a2 a3
	 * d e f   b1 b2 b3
	 * g h i   c1 c2 c3 */
	return a1*b2*c3-a1*b3*c2-b1*a2*c3+b1*a3*c2+c1*a2*b3-c1*a3*b2;
}
/*}}}*/
void Matrix3x3Invert(IssmDouble* Ainv,IssmDouble* A){/*{{{*/

	/*Intermediaries*/
	IssmDouble det,det_reciprocal;

	/*Compute determinant*/
	Matrix3x3Determinant(&det,A);
	if (fabs(det) < DBL_EPSILON) _error_("Determinant smaller than machine epsilon");

	/*Multiplication is faster than divsion, so we multiply by the reciprocal*/
	det_reciprocal = 1./det;  

	/*Compute invert*/
	Ainv[0]=(A[4]*A[8]-A[5]*A[7])*det_reciprocal; /* = (e*i-f*h)/det */
	Ainv[1]=(A[2]*A[7]-A[1]*A[8])*det_reciprocal; /* = (c*h-b*i)/det */
	Ainv[2]=(A[1]*A[5]-A[2]*A[4])*det_reciprocal; /* = (b*f-c*e)/det */
	Ainv[3]=(A[5]*A[6]-A[3]*A[8])*det_reciprocal; /* = (f*g-d*i)/det */
	Ainv[4]=(A[0]*A[8]-A[2]*A[6])*det_reciprocal; /* = (a*i-c*g)/det */
	Ainv[5]=(A[2]*A[3]-A[0]*A[5])*det_reciprocal; /* = (c*d-a*f)/det */
	Ainv[6]=(A[3]*A[7]-A[4]*A[6])*det_reciprocal; /* = (d*h-e*g)/det */
	Ainv[7]=(A[1]*A[6]-A[0]*A[7])*det_reciprocal; /* = (b*g-a*h)/det */
	Ainv[8]=(A[0]*A[4]-A[1]*A[3])*det_reciprocal; /* = (a*e-b*d)/det */
}/*}}}*/
void Matrix3x3Solve(IssmDouble* X,IssmDouble* A,IssmDouble* B){/*{{{*/

	IssmDouble Ainv[3][3];

	Matrix3x3Invert(&Ainv[0][0],A);
	for(int i=0;i<3;i++) X[i]=Ainv[i][0]*B[0] + Ainv[i][1]*B[1] + Ainv[i][2]*B[2];

}/*}}}*/

void Matrix4x4Determinant(IssmDouble* Adet,IssmDouble* A){/*{{{*/
	/*Compute determinant of a 4x4 matrix*/

	IssmDouble a1 = A[0*4+0];
	IssmDouble b1 = A[0*4+1]; 
	IssmDouble c1 = A[0*4+2];
	IssmDouble d1 = A[0*4+3];

	IssmDouble a2 = A[1*4+0];
	IssmDouble b2 = A[1*4+1]; 
	IssmDouble c2 = A[1*4+2];
	IssmDouble d2 = A[1*4+3];

	IssmDouble a3 = A[2*4+0]; 
	IssmDouble b3 = A[2*4+1];
	IssmDouble c3 = A[2*4+2];
	IssmDouble d3 = A[2*4+3];

	IssmDouble a4 = A[3*4+0];
	IssmDouble b4 = A[3*4+1]; 
	IssmDouble c4 = A[3*4+2];
	IssmDouble d4 = A[3*4+3];

	*Adet= a1 * Matrix3x3Determinant(b2, b3, b4, c2, c3, c4, d2, d3, d4)
		  - b1 * Matrix3x3Determinant(a2, a3, a4, c2, c3, c4, d2, d3, d4)
		  + c1 * Matrix3x3Determinant(a2, a3, a4, b2, b3, b4, d2, d3, d4)
		  - d1 * Matrix3x3Determinant(a2, a3, a4, b2, b3, b4, c2, c3, c4);
}
/*}}}*/
void Matrix4x4Adjoint(IssmDouble* Aadj,IssmDouble* A){/*{{{*/

    IssmDouble a1 = A[0*4+0];
    IssmDouble b1 = A[0*4+1]; 
    IssmDouble c1 = A[0*4+2];
    IssmDouble d1 = A[0*4+3];

    IssmDouble a2 = A[1*4+0];
    IssmDouble b2 = A[1*4+1]; 
    IssmDouble c2 = A[1*4+2];
    IssmDouble d2 = A[1*4+3];

    IssmDouble a3 = A[2*4+0];
    IssmDouble b3 = A[2*4+1];
    IssmDouble c3 = A[2*4+2];
    IssmDouble d3 = A[2*4+3];

    IssmDouble a4 = A[3*4+0];
    IssmDouble b4 = A[3*4+1]; 
    IssmDouble c4 = A[3*4+2];
    IssmDouble d4 = A[3*4+3];

    /* Row column labeling reversed since we transpose rows & columns*/
    Aadj[0*4+0]  =   Matrix3x3Determinant(b2, b3, b4, c2, c3, c4, d2, d3, d4);
    Aadj[1*4+0]  = - Matrix3x3Determinant(a2, a3, a4, c2, c3, c4, d2, d3, d4);
    Aadj[2*4+0]  =   Matrix3x3Determinant(a2, a3, a4, b2, b3, b4, d2, d3, d4);
    Aadj[3*4+0]  = - Matrix3x3Determinant(a2, a3, a4, b2, b3, b4, c2, c3, c4);

    Aadj[0*4+1]  = - Matrix3x3Determinant(b1, b3, b4, c1, c3, c4, d1, d3, d4);
    Aadj[1*4+1]  =   Matrix3x3Determinant(a1, a3, a4, c1, c3, c4, d1, d3, d4);
    Aadj[2*4+1]  = - Matrix3x3Determinant(a1, a3, a4, b1, b3, b4, d1, d3, d4);
    Aadj[3*4+1]  =   Matrix3x3Determinant(a1, a3, a4, b1, b3, b4, c1, c3, c4);

    Aadj[0*4+2]  =   Matrix3x3Determinant(b1, b2, b4, c1, c2, c4, d1, d2, d4);
    Aadj[1*4+2]  = - Matrix3x3Determinant(a1, a2, a4, c1, c2, c4, d1, d2, d4);
    Aadj[2*4+2]  =   Matrix3x3Determinant(a1, a2, a4, b1, b2, b4, d1, d2, d4);
    Aadj[3*4+2]  = - Matrix3x3Determinant(a1, a2, a4, b1, b2, b4, c1, c2, c4);

    Aadj[0*4+3]  = - Matrix3x3Determinant(b1, b2, b3, c1, c2, c3, d1, d2, d3);
    Aadj[1*4+3]  =   Matrix3x3Determinant(a1, a2, a3, c1, c2, c3, d1, d2, d3);
    Aadj[2*4+3]  = - Matrix3x3Determinant(a1, a2, a3, b1, b2, b3, d1, d2, d3);
    Aadj[3*4+3]  =   Matrix3x3Determinant(a1, a2, a3, b1, b2, b3, c1, c2, c3);
}/*}}}*/
void Matrix4x4Invert(IssmDouble* Ainv,IssmDouble* A){/*{{{*/

	/*Intermediaries*/
	IssmDouble det,det_reciprocal;

	/*Compute determinant*/
	Matrix4x4Determinant(&det,A);
	if(fabs(det) < DBL_EPSILON) _error_("Determinant smaller than machine epsilon");

	/*Multiplication is faster than division, so we multiply by the reciprocal*/
	det_reciprocal = 1./det;

	/*Compute adjoint matrix*/
	Matrix4x4Adjoint(Ainv,A);

	/*Scalte adjoint matrix to get inverse*/
	for(int i=0;i<4*4;i++) Ainv[i] = Ainv[i]*det_reciprocal;
}/*}}}*/
void Matrix4x4Solve(IssmDouble* X,IssmDouble* A,IssmDouble *B){/*{{{*/
	IssmDouble Ainv[4][4];

	Matrix4x4Invert(&Ainv[0][0],A);
	for(int i=0;i<4;i++) X[i]=Ainv[i][0]*B[0] + Ainv[i][1]*B[1] + Ainv[i][2]*B[2] + Ainv[i][3]*B[3];
}/*}}}*/

void newcell(IssmDouble** pcell, IssmDouble newvalue, bool top, int m){  /*{{{*/

    IssmDouble* cell=NULL;
    IssmDouble* dummy=NULL;

    /*recover pointer: */
    cell=*pcell;

    /*reallocate:*/
    dummy=xNew<IssmDouble>(m+1);

	/*copy data:*/
    if(top){
        dummy[0]=newvalue;
        for(int i=0;i<m;i++)dummy[i+1]=cell[i];
    }
    else{
        dummy[m]=newvalue;
        for(int i=0;i<m;i++)dummy[i]=cell[i];
    }

    /*delete and reassign: */
    xDelete<IssmDouble>(cell); cell=dummy;

    /*assign output pointer:*/
    *pcell=cell;
} /*}}}*/
IssmDouble  cellsum(IssmDouble* cell, int m){ /*{{{*/

	IssmDouble sum=0;

	for(int i=0;i<m;i++)sum+=cell[i];

	return sum;
} /*}}}*/
void celldelete(IssmDouble** pcell, int m, int* indices, int nind){ /*{{{*/

	/*input: */
	IssmDouble* cell=*pcell;

	/*output: */
	IssmDouble* newcell=xNew<IssmDouble>(nind);

	for(int i=0;i<nind;i++){
		newcell[i]=cell[indices[i]];
	}

	/*free allocation:*/
	xDelete<IssmDouble>(cell);

	/*assign output pointers: */
	*pcell=newcell;
} /*}}}*/
void cellsplit(IssmDouble** pcell, int m, int i,IssmDouble scale) { /*{{{*/

	/*input: */
	IssmDouble* cell=*pcell;

	/*output: */
	IssmDouble* newcell=xNew<IssmDouble>(m+1);

	for(int j=0;j<i;j++)newcell[j]=cell[j]; 
	newcell[i]=scale*cell[i];
	newcell[i+1]=scale* cell[i];
	for(int j=i+2;j<m+1;j++)newcell[j]=cell[j-1];

	/*free allocation:*/
	xDelete<IssmDouble>(cell);

	/*assign output pointers: */
	*pcell=newcell;
} /*}}}*/
void cellecho(int numcells, int m, ...) { /*{{{*/

	va_list arguments;                     
	IssmDouble** celllist= NULL;

	/*allocate variable length array: */
	celllist=xNew<IssmDouble*>(numcells); 

	va_start(arguments,m);

	for ( int x = 0; x < numcells; x++ ){
		celllist[x]= va_arg ( arguments, IssmDouble*); 
	}
	va_end ( arguments );                  

	_printf_("Echo of cell: \n");
	for(int i=0;i<m;i++){
		_printf_(i << ": ");
		for (int j=0;j<numcells;j++)_printf_(setprecision(10) << celllist[j][i] << " ");
		_printf_("\n");
	}

	/*deallocate:*/
	xDelete<IssmDouble*>(celllist);

} /*}}}*/
void CholeskyRealPositiveDefinite(IssmDouble* Lchol, IssmDouble* A, int ndim) { /*{{{*/
   /*CholeskyRealPositiveDefinite   computes lower triangular matrix of the Cholesky decomposition of A
   Follows the Cholesky–Banachiewicz algorithm
   Lchol should point to an IssmDouble* of same dimensions as A*/

	/*ensure zero-initialization*/
	for(int i=0;i<(ndim*ndim);i++) Lchol[i]=0;;

	for(int i=0;i<ndim;i++){
		for(int j=0;j<=i;j++){
			IssmDouble sum=0.;
			for(int k=0;k<j;k++) sum += Lchol[i*ndim+k]*Lchol[j*ndim+k];
			
			if(i==j) Lchol[i*ndim+j] = sqrt(A[i*ndim+j]-sum);
			else Lchol[i*ndim+j]     = 1./Lchol[j*ndim+j] * (A[i*ndim+j]-sum);
		}
	}
} /*}}}*/
