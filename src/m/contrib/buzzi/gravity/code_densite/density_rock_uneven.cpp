#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <cassert>
#include <gsl/gsl_multifit.h>
#include "mpi.h"
using namespace std;

class Matrix{/*{{{*/
	private:
		int     M;        /*Number of lines   */
		int     N;        /*Number if Columns */
		double *values;
	public:
		Matrix(int m_in,int n_in){/*{{{*/
			this->M = m_in;
			this->N = n_in;
			this->values = new double[M*N]();
		}/*}}}*/
		~Matrix(){/*{{{*/
			delete [] this->values;
		}/*}}}*/
		void Echo(void){/*{{{*/
			for(int i=0;i<M;i++){
				for(int j=0;j<N;j++){
					cout << " " << this->values[i*N+j];
				}
				cout << endl;
			}
		}/*}}}*/
		void SetValue(int i,int j,double value){/*{{{*/
			this->values[i*N+j] = value;
		}/*}}}*/
		double GetValue(int i,int j){/*{{{*/
			return this->values[i*N+j];
		}/*}}}*/
		void GetSize(int* pM,int* pN){/*{{{*/
			*pM = this->M;
			*pN = this->N;
		}/*}}}*/
		double* GetPointer(void){/*{{{*/
			return this->values;
		}/*}}}*/
		void MatrixSum(Matrix* A,Matrix* B){/*{{{*/
			/*Check that sizes are compatible*/
			int M_B,N_B,M_A,N_A;
			B->GetSize(&M_B,&N_B);
			A->GetSize(&M_A,&N_A);
			assert(this->M==M_B && this->N==N_B);
			assert(this->M==M_A && this->N==N_A);
			for(int i=0;i<this->M;i++){
				for(int j=0;j<this->N;j++){
					this->SetValue(i,j,A->GetValue(i,j) + B->GetValue(i,j));
				}
			}
		}/*}}}*/
		void MatrixDiff(Matrix* A,Matrix* B){/*{{{*/
			/*Check that sizes are compatible*/
			int M_B,N_B,M_A,N_A;
			B->GetSize(&M_B,&N_B);
			A->GetSize(&M_A,&N_A);
			assert(this->M==M_B && this->N==N_B);
			assert(this->M==M_A && this->N==N_A);
			for(int i=0;i<this->M;i++){
				for(int j=0;j<this->N;j++){
					this->SetValue(i,j,A->GetValue(i,j) - B->GetValue(i,j));
				}
			}
		}/*}}}*/
		void MatrixAbs(Matrix* A){/*{{{*/
			for(int i=0;i<this->M;i++){
				for(int j=0;j<this->N;j++){
					this->SetValue(i,j,fabs(A->GetValue(i,j)));
				}
			}
		}/*}}}*/
		void MatrixEqual(Matrix* A){/*{{{*/
			for(int i=0;i<this->M;i++){
				for(int j=0;j<this->N;j++){
					this->SetValue(i,j,A->GetValue(i,j));
				}
			}
		}/*}}}*/
		double MatrixInternSum(){/*{{{*/
			double sum=0;
			for(int i=0;i<this->M;i++){
				for(int j=0;j<this->N;j++){
					sum+=this->GetValue(i,j);
				}
			}
			return sum;
		}/*}}}*/
		void ExtractLine(Matrix* A,int i){/*{{{*/
			/* Check that the size of A is compatible */
			int M_A,N_A;
			A->GetSize(&M_A,&N_A);
			assert(M_A==1 && this->N==N_A);
			for(int j=0;j<this->N;j++){
				A->SetValue(0,j,this->GetValue(i,j));
			}
		}/*}}}*/
		void ExtractColumn(Matrix* A,int j){/*{{{*/
			/* Check that the size of A is compatible */
			int M_A,N_A;
			A->GetSize(&M_A,&N_A);
			assert(N_A==1 && this->M==M_A);
			for(int i=0;i<this->M;i++){
				A->SetValue(i,0,this->GetValue(i,j));
			}
		}/*}}}*/
		void AddNumber(double a){/*{{{*/
			for(int i=0;i<this->M;i++){
				for(int j=0;j<this->N;j++){
					this->SetValue(i,j,this->GetValue(i,j) + a);
				}
			}
		}/*}}}*/
};/*}}}*/

/*Local prototypes{{{*/
void makep(Matrix *Pobs,int nx,int ny, int dx, int dy);
void vec2grid(Matrix *V,Matrix *V1,Matrix *V2,int nx, int ny);
void msplit( Matrix *m, Matrix *m1,Matrix *m2,Matrix *dlevel);
void plouff(Matrix *g,Matrix *Pobs,Matrix *Pp,Matrix * mesh,Matrix *rho,int dx,int dy, int dn,int m,int n,int l,Matrix *evalid,int my_rank,int num_procs);
void vec2gridsimple(Matrix *V,Matrix *V1,int nx, int ny);
void reshape(Matrix* V,Matrix* V1,int nx,int ny);
double misfit(Matrix* m0,Matrix* evalid,Matrix* gobs,Matrix *dlevel,Matrix* Pobs,Matrix* xobs,Matrix* yobs,Matrix* Pp,Matrix* rho1, Matrix* rho2,int dx,int dy,int dn,int nx,int ny, int mx,int my,int my_rank,int num_procs);
void GSLsquarefit(Matrix** pX,Matrix* A,Matrix* B);
/*}}}*/

int main(int argc,char *argv[]){/*{{{*/

	int my_rank,num_procs;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

	/* Define the variables {{{*/

	int    dx     = 1000;   /* prism dimension in x-direction                           */
	int    dy     = 1000;   /* prism dimension in y-direction                           */
	int    mx     = 99;    /* number of prisms in x-direction                          */
	int    my     = 99;    /* number of prisms in y-direction                          */
	int    nx     = 99;    /* number of data points in x-direction                     */
	int    ny     = 99;    /* number of data points in y-direction                     */
	int    dn     = 15000; /* distance for neighbouting prisms for gravity calculation */

	Matrix *Pobs=new Matrix(nx*ny,2); /* data positions */
	makep(Pobs,nx,ny,dx,dy);
	// Pobs->Echo();



	Matrix *Pp=new Matrix(mx*my,2); /* prisms positions */
	makep(Pp,mx,my,dx,dy);
	// Pp->Echo();

	double  rhoi = 917;           /* ice density     */
	double  rhow = 1030;          /* water density   */
	// double  rhos = 2013;		      /* sediment density */

	double rhoc_min=2000.;
	double rhoc_max=3000.;

	Matrix *Rho  = new Matrix(1,2);
	Rho->SetValue(0,0,rhoi);
	Rho->SetValue(0,1,rhow);
	Matrix *rho1  = new Matrix(1,3);
	rho1->SetValue(0,0,rhoi);
	rho1->SetValue(0,1,rhow);
	rho1->SetValue(0,2,rhoc_min);
	Matrix *rho2  = new Matrix(1,2);
	rho2->SetValue(0,0,rhoi-rhoc_min);
	rho2->SetValue(0,1,rhow-rhoc_min);


	Matrix *xobs= new Matrix(ny,nx);
	Matrix *yobs= new Matrix(ny,nx);

	vec2grid(Pobs,xobs,yobs,nx,ny);
	//	xobs->Echo();
	//	yobs->Echo();


	/*}}}*/     
	/* load the data {{{*/


	double inputnumber;

	/* Levels of data acquisition */

	ifstream file0("dataalti.txt");
	Matrix * dlevel= new Matrix(nx*ny,1);
	for(int i=0;i<ny*nx; i++){
		file0 >> inputnumber;
		dlevel->SetValue(i,0,inputnumber);
	}
	file0.close();

	/* Observed gravity anomaly */

	ifstream file1("gravityraw.txt");
	Matrix * gobs= new Matrix(nx*ny,1);
	for(int i=0;i<ny*nx; i++){ 
		file1 >> inputnumber;
		gobs->SetValue(i,0, inputnumber*1e-5);
	}
	file1.close();
	//	gobs->Echo();


	/* id of grid to evaluate misfit */


	ifstream file4("evalid1.txt");
	Matrix * evalid= new Matrix(nx*ny,1);
	for(int s=0;s<nx*ny; s++){ 
		file4 >> inputnumber;
		evalid->SetValue(s,0,inputnumber);
	}
	file4.close();
	//	evalid->Echo();

	/* initial guess of the model */

	ifstream file5("m0_102714contzach.txt");
	Matrix * mesh_ini= new Matrix(mx*my,3);
	for(int s=0;s<mx*my; s++){ 
		for(int j=0;j<3;j++){
			file5 >> inputnumber;
			mesh_ini->SetValue(s,j,inputnumber);
		}
	}
	file5.close();
	//	mesh_ini->Echo();
	/*}}}*/
	/* Test {{{ */


	double rhoc=rhoc_min;
	double rhoc_opti=rhoc_min;
	double E=misfit(mesh_ini,evalid,gobs,dlevel,Pobs,xobs,yobs,Pp,rho1,rho2,dx,dy,dn,nx,ny,mx,my,my_rank,num_procs);
	double E_opti=E;

	for(int i=rhoc_min;i<rhoc_max+1;i++){
		rhoc=i;
		rho1->SetValue(0,2,rhoc);
		rho2->SetValue(0,0,rhoi-rhoc);
		rho2->SetValue(0,1,rhow-rhoc);

		E=misfit(mesh_ini,evalid,gobs,dlevel,Pobs,xobs,yobs,Pp,rho1,rho2,dx,dy,dn,nx,ny,mx,my,my_rank,num_procs);

		if(E<E_opti){
			E_opti=E;
			rhoc_opti=rhoc;
		}
		if(my_rank==0){
			cout<<rhoc<<"  "<<rhoc_opti<<"  "<<E<<"  "<<E_opti<<endl;
		}
	}




	delete Pobs;
	delete Pp;
	delete Rho;
	delete rho1;
	delete rho2;
	delete xobs;
	delete yobs;
	delete mesh_ini;
	delete evalid;
	delete gobs;

	/*}}}*/

	MPI_Finalize();

	return 0;
}/*}}}*/

void GSLsquarefit(Matrix** pX,Matrix* A,Matrix* B){/*{{{*/

	/*GSL Matrices and vectors: */
	int    M,N;
	double chisq;
	/*Get system size*/
	A->GetSize(&M,&N);

	/*Initialize gsl matrices and vectors: */
	gsl_matrix* a = gsl_matrix_alloc(M,N);
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
			gsl_matrix_set (a,i,j,A->GetValue(i,j));
		}
	}
	gsl_vector* b = gsl_vector_alloc(M);
	for(int i=0;i<M;i++){
		gsl_vector_set(b,i,B->GetValue(i,0));
	}

	gsl_vector* x = gsl_vector_alloc(N);
	gsl_matrix* cov = gsl_matrix_alloc(N,N);

	/*Least square fit: */
	gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(M,N);
	gsl_multifit_linear (a, b, x, cov, &chisq, work);
	gsl_multifit_linear_free (work);

	/*Clean up and assign output pointer*/
	Matrix* X = new Matrix(N,1);
	for(int j=0;j<N;j++){
		X->SetValue(j,0,gsl_vector_get(x,j));
	}
	*pX = X;

	gsl_matrix_free(a);
	gsl_vector_free(x);
	gsl_vector_free(b);
	gsl_matrix_free(cov);

}/*}}}*/
void makep(Matrix *Pobs,int nx,int ny, int dx, int dy){/*{{{*/
	for(int i=0;i<ny;i++){
		for(int j=0;j<nx;j++){
			Pobs->SetValue(j+nx*i,0,j*dx);
			Pobs->SetValue(j+nx*i,1,i*dy);
		}
	}
}/*}}}*/
void vec2grid(Matrix *V,Matrix *V1,Matrix *V2,int nx, int ny){/*{{{*/
	for(int i=0;i<ny;i++){
		for (int j=0;j<nx;j++){
			V1->SetValue(i,j, V->GetValue(j+nx*i,0));
			V2->SetValue(i,j, V->GetValue(j+nx*i,1));
		}
	}
}/*}}}*/
void msplit( Matrix *m, Matrix *m1,Matrix *m2,Matrix* dlevel){/*{{{*/
	int sizem1,sizem2;
	m->GetSize(&sizem1,&sizem2);
	for(int i=0;i<sizem1;i++){
		for(int j=0;j<sizem2+1;j++){
			if(j<sizem2){
				m1->SetValue(i,j,1e-10*(sizem2+1-j));
				m2->SetValue(i,j,m->GetValue(i,j));
				if(m->GetValue(i,j)<0){
					m1->SetValue(i,j,m->GetValue(i,j));
					m2->SetValue(i,j,i*1e-10);
				}
				m1->SetValue(i,j,m1->GetValue(i,j)+dlevel->GetValue(i,1));
				m2->SetValue(i,j,m2->GetValue(i,j)+dlevel->GetValue(i,1));
			}
			else{
				m1->SetValue(i,j,1e-10+dlevel->GetValue(i,1));
			}
		}
	}
}/*}}}*/
void plouff(Matrix *g,Matrix *Pobs,Matrix *Pp,Matrix * mesh,Matrix *rho,int dx,int dy, int dn,int m,int n,int l,Matrix *evalid,int my_rank,int num_procs){/*{{{*/
	double gg=6.673e-11;
	int si,sj,id,s;
	double R,Q,P;
	Matrix *xp= new Matrix(1,2);
	Matrix *yp= new Matrix(1,2);
	Matrix *xpp= new Matrix(1,2);
	Matrix *ypp= new Matrix(1,2);
	Matrix *U= new Matrix(l,4);
	Matrix *U1=new Matrix(1,4);
	Matrix *U2=new Matrix(1,4);
	Matrix *gl= new Matrix(1,l-1);
	bool test=true;

	double *glocal=new double[n]();

	for(int c=my_rank;c<n;c+=num_procs){
		glocal[c]=0;
		if(evalid->GetValue(c,0)==1){
			for(int a=0;a<m;a++){
				test=true;
				xp->SetValue(0,0,Pp->GetValue(a,0)-Pobs->GetValue(c,0));
				xp->SetValue(0,1,Pp->GetValue(a,0)-Pobs->GetValue(c,0)+dx);
				if(xp->GetValue(0,0)<0 && xp->GetValue(0,0)<xp->GetValue(0,1) && xp->GetValue(0,0)*xp->GetValue(0,1)>=0){
					xpp->SetValue(0,0,xp->GetValue(0,1));
					xpp->SetValue(0,1,xp->GetValue(0,0));
					xp->MatrixAbs(xpp);
				}
				yp->SetValue(0,0,Pp->GetValue(a,1)-Pobs->GetValue(c,1));
				yp->SetValue(0,1,Pp->GetValue(a,1)-Pobs->GetValue(c,1)+dy);
				if(yp->GetValue(0,0)<0 && yp->GetValue(0,0)<yp->GetValue(0,1) && yp->GetValue(0,0)*yp->GetValue(0,1)>=0){
					ypp->SetValue(0,0,yp->GetValue(0,1));
					ypp->SetValue(0,1,yp->GetValue(0,0));
					yp->MatrixAbs(ypp);
				}
				P=sqrt(xp->GetValue(0,0)*xp->GetValue(0,0)+yp->GetValue(0,0)*yp->GetValue(0,0));
				if(P>dn){
					test=false;
					for(int i=0;i<l-1;i++){
						gl->SetValue(0,i,0);
					}
				}
				if(test==true){
					si=1;
					sj=1;
					id=0;
					for(int i=0;i<2;i++){
						si*=-1;
						for(int j=0;j<2;j++){
							sj*=-1;
							s=si*sj;
							for(int k=0;k<l;k++){
								R=sqrt(xp->GetValue(0,i)*xp->GetValue(0,i)+yp->GetValue(0,j)*yp->GetValue(0,j)+mesh->GetValue(a,k)*mesh->GetValue(a,k));
								Q=atan(xp->GetValue(0,i)*yp->GetValue(0,j)/(mesh->GetValue(a,k)*R));
								U->SetValue(k,id,s*(mesh->GetValue(a,k)*Q-xp->GetValue(0,i)*log(R+yp->GetValue(0,j))-yp->GetValue(0,j)*log(R+xp->GetValue(0,i))));
							}
							id++;
						}
					}
					for(int b=0;b<l-1;b++){
						U->ExtractLine(U1,b);
						U->ExtractLine(U2,b+1);
						gl->SetValue(0,b,rho->GetValue(0,b)*(U1->MatrixInternSum()*(-1)+U2->MatrixInternSum()));
					}
				}
				glocal[c]=glocal[c]+gg*gl->MatrixInternSum();
			}
		}
	}

	MPI_Allreduce(glocal,g->GetPointer(),n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	delete xp;
	delete yp;
	delete xpp;
	delete ypp;
	delete gl;
	delete U;
	delete U1;
	delete U2;
	delete []glocal;
}/*}}}*/
void vec2gridsimple(Matrix *V,Matrix *V1,int nx, int ny){/*{{{*/
	for(int i=0;i<ny;i++){
		for (int j=0;j<nx;j++){
			V1->SetValue(i,j, V->GetValue(j+nx*i,0));
		}
	}
}/*}}}*/
void reshape(Matrix* V,Matrix* V1,int nx,int ny){/*{{{*/
	for (int i=0;i<ny;i++){
		for(int j=0;j<nx;j++){
			V1->SetValue(j+nx*i,0,V->GetValue(i,j));
		}
	}
}/*}}}*/
double misfit(Matrix* m0,Matrix* evalid,Matrix* gobs,Matrix *dlevel,Matrix* Pobs,Matrix* xobs,Matrix* yobs,Matrix* Pp,Matrix* rho1, Matrix* rho2,int dx,int dy,int dn,int nx,int ny, int mx,int my,int my_rank,int num_procs){/*{{{*/
	Matrix* m1=new Matrix(mx*my,4);
	Matrix* m2=new Matrix(mx*my,3);
	Matrix* g1=new Matrix(nx*ny,1);
	Matrix* g2=new Matrix(nx*ny,1);
	Matrix* g=new Matrix(nx*ny,1);
	Matrix* gcalgr=new Matrix(ny,nx);
	Matrix* gcalvec=new Matrix(nx*ny,1);
	Matrix* df=new Matrix(nx*ny,1);
	Matrix* G=new Matrix(nx*ny,3);
	double a=0;
	double b=0;
	double e=0;
	msplit(m0,m1,m2,dlevel);
	plouff(g1,Pobs,Pp,m1,rho1,dx,dy,dn,mx*my,nx*ny,4,evalid, my_rank, num_procs);
	plouff(g2,Pobs,Pp,m2,rho2,dx,dy,dn,mx*my,nx*ny,3,evalid, my_rank, num_procs);
	g->MatrixSum(g1,g2);
	vec2gridsimple(g,gcalgr,nx,ny);
	reshape(gcalgr,gcalvec,nx,ny);
	for (int i=0;i<nx*ny;i++){
		df->SetValue(i,0,evalid->GetValue(i,0)*(gobs->GetValue(i,0)-gcalvec->GetValue(i,0)));
		G->SetValue(i,0,evalid->GetValue(i,0)*Pobs->GetValue(i,0));
		G->SetValue(i,1,evalid->GetValue(i,0)*Pobs->GetValue(i,1));
		G->SetValue(i,2,evalid->GetValue(i,0));
	}
	Matrix* M = NULL;
	GSLsquarefit(&M,G,df);

	for (int i=0;i<ny;i++){
		for(int j=0;j<nx;j++){
			gcalgr->SetValue(i,j,gcalgr->GetValue(i,j)+xobs->GetValue(i,j)*M->GetValue(0,0)+yobs->GetValue(i,j)*M->GetValue(1,0)+M->GetValue(2,0));
		}
	}
	reshape(gcalgr,g,nx,ny);
	for (int i=0;i<nx*ny;i++){
		a=a+fabs(evalid->GetValue(i,0)*(gobs->GetValue(i,0)-g->GetValue(i,0)));
		b=b+fabs(evalid->GetValue(i,0)*(gobs->GetValue(i,0)+g->GetValue(i,0)));
	}
	e=2*a/(a+b);

	delete m1;
	delete m2;
	delete g1;
	delete g2;
	delete g;
	delete gcalgr;
	delete gcalvec;
	delete df;
	delete G;
	delete M;

	return e;
}/*}}}*/
