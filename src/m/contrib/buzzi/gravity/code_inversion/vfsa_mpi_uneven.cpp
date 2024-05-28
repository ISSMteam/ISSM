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
double signe(double a);
void filtergrav(Matrix* A,Matrix* Ain,double ctr,double sd,int mx,int my);
void newmodelgen(Matrix* m0,Matrix* m1,Matrix* bathy,Matrix* icethick,int mx,int my,double T,double ptval,double mmax,double mmax2,double ctr,double sd, Matrix *landmask);
double coolshed(double T0,double k,double c,double D);
/*}}}*/

int main(int argc,char *argv[]){/*{{{*/
	
	int my_rank,num_procs;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

	/* Seed the random number generator {{{*/
		srand (time(NULL));            /*}}}*/
	/* Define the variables {{{*/

	int    dx     = 1000;   /* prism dimension in x-direction                           */
	int    dy     = 1000;   /* prism dimension in y-direction                           */
	int    mx     = 99;    /* number of prisms in x-direction                          */
	int    my     = 99;    /* number of prisms in y-direction                          */
	int    nx     = 99;    /* number of data points in x-direction                     */
	int    ny     = 99;    /* number of data points in y-direction                     */
	int    dn     = 15000; /* distance for neighbouting prisms for gravity calculation */
	double ptval  = 100.;  /* max. amount to perturb model                             */
	double ptval2 = 100.;

	Matrix *Pobs=new Matrix(nx*ny,2); /* data positions */
	makep(Pobs,nx,ny,dx,dy);
	// Pobs->Echo();


	Matrix *Pp=new Matrix(mx*my,2); /* prisms positions */
	makep(Pp,mx,my,dx,dy);
	// Pp->Echo();

	double  rhoi = 917;           /* ice density     */
	double  rhow = 1030;          /* water density   */
	// double  rhos = 2013;		      /* sediment density */
	double  rhoc = 2670;          /* bedrock density */

	Matrix *Rho  = new Matrix(1,2);
	Rho->SetValue(0,0,rhoi);
	Rho->SetValue(0,1,rhow);
	Matrix *rho1  = new Matrix(1,3);
	rho1->SetValue(0,0,rhoi);
	rho1->SetValue(0,1,rhow);
	rho1->SetValue(0,2,rhoc);
	Matrix *rho2  = new Matrix(1,2);
	rho2->SetValue(0,0,rhoi-rhoc);
	rho2->SetValue(0,1,rhow-rhoc);

	double ctr=1;            /* parameter for filtering */
	double sd=0.1;

	Matrix *xobs= new Matrix(ny,nx);
	Matrix *yobs= new Matrix(ny,nx);

	vec2grid(Pobs,xobs,yobs,nx,ny);
	//	xobs->Echo();
	//	yobs->Echo();


	double mmax  = 1000;               /* max value for layer interfaces */
	double mmax2 = 1000;
	double mmax3 = 1000;

	/* control parameter for temperature schedule  */

	double ca=0.9;                    /* for acceptance */
	double cm=0.5;                    /* for model perturbation */

	double T0a          = 0.1;      /* initial temperature for acceptance           */
	double T0m          = 0.9;      /* initial temperature for model perturbation   */
	double D            = 2;        /* dimension of the model                       */
	int    maxconsecrej = 1000;     /* max consecutive rejection                    */
	int    maxsuccess   = 100;      /* max number of success within one temperature */
	double T_min        = 1e-10;    /* stopping temp                                */
	double Tred         = 1;
	double E_min        = -1000000;
	double E_exp        = 0.0291;   /* expected misfit                              */
	int    maxiter      = 10000;
	int    maxtotaliter = 1000000;
	double Tol          = 1e-10;    /* tolerance on misfit                          */
	int    sfreq        = 100;

	/*}}}*/     
	/* load the data {{{*/

	/*landmask */

	ifstream file("landmaskzach.txt");
	Matrix * landmask= new Matrix(nx*ny,1);
	double inputnumber;
	for(int i=0;i<ny*nx; i++){ 
		file >> inputnumber;
		landmask->SetValue(i,0,inputnumber);
	}
	file.close();

	/* Levels of data acquisition */

	ifstream file0("altizach.txt");
	Matrix * dlevel= new Matrix(nx*ny,1);
	for(int i=0;i<ny*nx; i++){
		file0 >> inputnumber;
		dlevel->SetValue(i,0,inputnumber);
	}
	file0.close();

	/* Observed gravity anomaly */

	ifstream file1("gravityzach.txt");
	Matrix * gobs= new Matrix(nx*ny,1);
	for(int i=0;i<ny*nx; i++){ 
		file1 >> inputnumber;
		gobs->SetValue(i,0, inputnumber*1e-5);
	}
	file1.close();
	//	gobs->Echo();

	/* load data about the ice thickness */

	ifstream file2("icethickzach.txt");
	Matrix * icethick= new Matrix(mx*my,1);
	for(int s=0;s<mx*my; s++){ 
		file2 >> inputnumber;
		icethick->SetValue(s,0,inputnumber);
	}
	file2.close();
	//	icethick->Echo();

	/* load the batimethry data */

	ifstream file3("bathymetryzach.txt");
	Matrix * bathy= new Matrix(mx*my,1);
	for(int s=0;s<mx*my; s++){ 
		file3 >> inputnumber;
		bathy->SetValue(s,0,inputnumber);
	}
	file3.close();
	//	bathy->Echo();

	/* id of grid to evaluate misfit */


	ifstream file4("evalidzach.txt");
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
	/* VFSA {{{ */

	/* name of the files to save results */
	std::ofstream savefile1 ("r_zach.txt");
	std::ofstream savefile2("m_zach.txt");

	/* counters initialization */
	int    success   = 0;
	int    finished  = 0;
	int    consec    = 0;
	double Ta        = T0a;
	double Tm        = T0m;
	int    iterT     = 0;   /* iteration within a T      */
	int    total     = 0;   /* total number of iteration */
	int    totaliter = 0;
	int    msave     = 0;
	double E_new;
	double E_final;
	double dE;
	double P;
	double rn;
	Matrix* m_old    = new Matrix(mx *my,3);
	Matrix* m_min    = new Matrix(mx *my,3);
	Matrix* m_new    = new Matrix(mx *my,3);
	m_old->MatrixEqual(mesh_ini);

	/* calculate initial misfit */
	double E_old=misfit(m_old,evalid,gobs,dlevel,Pobs,xobs,yobs,Pp,rho1,rho2,dx,dy,dn,nx,ny,mx,my, my_rank, num_procs);
	/* record initial settings */
	if(!my_rank){
		savefile1 << "P     "<< "Ta    "<< "Tm    "<< "Eold  "<< "totaliter "<< "Tred   "<< endl;
		savefile1 << "nan   "<<  Ta<<"   "<< Tm<<"   "<< E_old<<"     "<< totaliter<<"         "<< Tred <<"  "<< endl;
		savefile2 << totaliter<< endl;
		for(int i=0;i<mx*my;i++){
			savefile2 << m_old->GetValue(i,0)<<"   "<< m_old->GetValue(i,1)<<"   "<< m_old->GetValue(i,2)<<endl;
		}
		savefile2 << "111111111111111111111111111111111111111111111111111111111111111111111111111"<< endl;
	}
	/* beginning of the loop */

	while(finished==0){
		iterT++;
		totaliter++;

		/* stop or reduce T */
		if(iterT>=maxiter || success>maxsuccess){
			if(Ta<T_min || total>maxtotaliter || fabs(E_old)<=Tol){
				finished=1;
				total+=iterT;
				break;
			}
			else{ /* reduce T */
				Ta=coolshed(T0a,Tred,ca,D);
				Tm=coolshed(T0m,Tred,cm,D);
				total+=iterT;
				iterT=0;
				success=1;
				Tred++;
				consec=0;
			}
		}

		/* update model and calculate energy */

		newmodelgen(m_old,m_new,bathy,icethick,mx,my,Tm,ptval,mmax,mmax2,ctr,sd, landmask);  /* new model */
		E_new=misfit(m_new,evalid,gobs,dlevel,Pobs,xobs,yobs,Pp,rho1,rho2,dx,dy,dn,nx,ny,mx,my, my_rank, num_procs); /* new energy */
		dE=E_new-E_old;                                        /* energy difference */

		/* acceptance probability */

		P=exp(-dE/Ta);

		/* stop if energy is lower than specified minimum */
		if (E_new<E_min){
			m_old->MatrixEqual(m_new);
			E_old=E_new;
			break;
		}

		rn=rand()/double (RAND_MAX);

		/* accept new model or not */
		if(dE<=0){
			m_old->MatrixEqual(m_new);
			E_old=E_new;
			E_final=E_old;
			success++;
			consec=0;
			if(!my_rank){
				savefile1 << P<<"   "<<  Ta<<"   "<< Tm<<"   "<< E_old<<"     "<< totaliter<<"         "<< Tred <<"  "<< endl;
			}
			if(Ta<1e-3){
				if(!my_rank){
					savefile2 << totaliter<< endl;
					for(int i=0;i<mx*my;i++){
						savefile2 << m_old->GetValue(i,0)<<"   "<< m_old->GetValue(i,1)<<"   "<< m_old->GetValue(i,2)<<endl;
					}
					savefile2 << "111111111111111111111111111111111111111111111111111111111111111111111111111"<< endl;
				}
			}
		}
		else{
			if(P>rn){
				m_old->MatrixEqual(m_new);
				E_old=E_new;
				success++;
				consec=0;
				if(!my_rank){
					savefile1 << P<<"   "<<  Ta<<"   "<< Tm<<"   "<< E_old<<"     "<< totaliter<<"         "<< Tred <<"  "<< endl;
					if(Ta<1e-3){
						savefile2 << totaliter<< endl;
						for(int i=0;i<mx*my;i++){
							savefile2 << m_old->GetValue(i,0)<<"   "<< m_old->GetValue(i,1)<<"   "<< m_old->GetValue(i,2)<<endl;
						}
						savefile2 << "111111111111111111111111111111111111111111111111111111111111111111111111111"<< endl;
					}
				}
			}
			else{
				consec++;
			}
		}
	}

	m_min->MatrixEqual(m_old);
	if(!my_rank){
		savefile1 << "nan"<<"   "<<  "nan"<<"   "<< "nan"<<"   "<< E_final<<"     "<< "nan"<<"         "<< "nan" <<"  "<< endl;
		savefile2 << " Mesh final"<< endl;
		for(int i=0;i<mx*my;i++){
			savefile2 << m_min->GetValue(i,0)<<"   "<< m_min->GetValue(i,1)<<"   "<< m_min->GetValue(i,2)<<endl;
		}
	}
	savefile1.close();
	savefile2.close();

	delete m_old;
	delete m_min;
	delete m_new;
	delete Pobs;
	delete Pp;
	delete Rho;
	delete rho1;
	delete rho2;
	delete xobs;
	delete yobs;
	delete mesh_ini;
	delete bathy;
	delete icethick;
	delete evalid;
	delete gobs;
	delete landmask;

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
void plouff(Matrix *g,Matrix *Pobs,Matrix *Pp,Matrix * mesh,Matrix *rho,int dx,int dy, int dn,int m,int n,int l,Matrix* evalid,int my_rank,int num_procs){/*{{{*/
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
		if(evalid->GetValue(i,0)==1){
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
void newmodelgen(Matrix* m0,Matrix* m1,Matrix* bathy,Matrix* icethick,int mx,int my,double T,double ptval,double mmax,double mmax2,double ctr,double sd, Matrix *landmask){/*{{{*/
	Matrix* m1gr=new Matrix(my,mx);
	Matrix* m1grsm=new Matrix(my,mx);
	Matrix* m1col=new Matrix(mx*my,1);
	Matrix* m1gr2=new Matrix(my,mx);
	Matrix* m1grsm2=new Matrix(my,mx);
	Matrix* m1col2=new Matrix(mx*my,1);
	Matrix* nptflag= new Matrix(mx*my,1);
	double u=0;
	double y=0;
	m1->MatrixEqual(m0);
	nptflag->MatrixSum(icethick,bathy);
	/* first layer: ice */
	for (int i=0;i<mx*my;i++){
		if(landmask->GetValue(i,0)==2){
			if(nptflag->GetValue(i,0)==0){
				u=double(rand())/double(RAND_MAX);
				y=signe(u-0.5)*T*(pow(1+1/T,fabs(2*u-1))-1);
				m1->SetValue(i,1,m0->GetValue(i,1)+y*ptval);
				m1->SetValue(i,2,m1->GetValue(i,1)+1e-10);
				if(m1->GetValue(i,1)<=m1->GetValue(i,0)){
					m1->SetValue(i,1,m1->GetValue(i,0)+1e-10);
					m1->SetValue(i,2,m1->GetValue(i,1)+1e-10);
				}
				if(m1->GetValue(i,1)>=m1->GetValue(i,0)+mmax){
					m1->SetValue(i,1,m1->GetValue(i,0)+mmax);
					m1->SetValue(i,2,m1->GetValue(i,1)+1e-10);
				}
			}
		}
		else if(landmask->GetValue(i,0)==0){
			if(nptflag->GetValue(i,0)==0){
				u=double(rand())/double(RAND_MAX);
				y=signe(u-0.5)*T*(pow(1+1/T,fabs(2*u-1))-1);
				m1->SetValue(i,2,m0->GetValue(i,2)+y*ptval);
				if(m1->GetValue(i,2)<=m1->GetValue(i,0)){
					m1->SetValue(i,2,m1->GetValue(i,0)+1e-10);
				}
				if(m1->GetValue(i,2)>=m1->GetValue(i,0)+mmax2){
					m1->SetValue(i,2,m1->GetValue(i,0)+mmax2);
				}
			}
		}
		else if(landmask->GetValue(i,0)==3){
			if(nptflag->GetValue(i,0)==0){
				u=double(rand())/double(RAND_MAX);
				y=signe(u-0.5)*T*(pow(1+1/T,fabs(2*u-1))-1);
				m1->SetValue(i,1,m0->GetValue(i,1)+y*ptval);
				if(m1->GetValue(i,1)<=m1->GetValue(i,0)){
					m1->SetValue(i,1,m1->GetValue(i,0)+1e-10);
				}
				if(m1->GetValue(i,1)>=m1->GetValue(i,0)+mmax){
					m1->SetValue(i,1,m1->GetValue(i,0)+mmax);
				}
				u=double(rand())/double(RAND_MAX);
				y=signe(u-0.5)*T*(pow(1+1/T,fabs(2*u-1))-1);
				m1->SetValue(i,2,m0->GetValue(i,2)+y*ptval);
				if(m1->GetValue(i,2)<=m1->GetValue(i,1)){
					m1->SetValue(i,2,m1->GetValue(i,1)+1e-10);
				}
				if(m1->GetValue(i,2)>=m1->GetValue(i,1)+mmax2){
					m1->SetValue(i,2,m1->GetValue(i,1)+mmax2);
				}
			}
		}
	}

	m1->ExtractColumn(m1col,1);
	vec2gridsimple(m1col,m1gr,mx,my);
	filtergrav(m1grsm,m1gr,ctr,sd,mx,my);
	reshape(m1grsm,m1col,mx,my);
	m1->ExtractColumn(m1col2,2);
	vec2gridsimple(m1col2,m1gr2,mx,my);
	filtergrav(m1grsm2,m1gr2,ctr,sd,mx,my);
	reshape(m1grsm2,m1col2,mx,my);

	for (int i=0;i<mx*my;i++){
		if(landmask->GetValue(i,0)==2){
			if(nptflag->GetValue(i,0)==0){
				m1->SetValue(i,1,m1col->GetValue(i,0));
				m1->SetValue(i,2,m1col2->GetValue(i,0));
				if(m1->GetValue(i,1)<=m1->GetValue(i,0)){
					m1->SetValue(i,1,m1->GetValue(i,0)+1e-10);
					m1->SetValue(i,2,m1->GetValue(i,1)+1e-10);
				}
				if(fabs(m1->GetValue(i,2)-m1->GetValue(i,1))>1){
					m1->SetValue(i,2,m1->GetValue(i,1)+1e-10);
				}
			}
			else{
				m1->SetValue(i,1,m0->GetValue(i,1));
				m1->SetValue(i,2,m0->GetValue(i,2));
			}
		}
		else if(landmask->GetValue(i,0)==0){
			if(nptflag->GetValue(i,0)==0){
				m1->SetValue(i,2,m1col2->GetValue(i,0));
				if(m1->GetValue(i,2)<=m1->GetValue(i,0)){
					m1->SetValue(i,2,m1->GetValue(i,0)+1e-10);
				}
				if(fabs(m1->GetValue(i,0)-m1->GetValue(i,1))>1){
					m1->SetValue(i,1,m1->GetValue(i,0)+1e-10);
				}
			}
			else{
				m1->SetValue(i,1,m0->GetValue(i,1));
				m1->SetValue(i,2,m0->GetValue(i,2));
			}
		}
		else if(landmask->GetValue(i,0)==3){
			if(nptflag->GetValue(i,0)==0){
				m1->SetValue(i,1,m1col->GetValue(i,0));
				m1->SetValue(i,2,m1col2->GetValue(i,0));
				if(m1->GetValue(i,1)<=m1->GetValue(i,0)){
					m1->SetValue(i,1,m1->GetValue(i,0)+1e-10);
				}
				if(m1->GetValue(i,2)<=m1->GetValue(i,1)){
					m1->SetValue(i,2,m1->GetValue(i,1)+1e-10);
				}
			}
			else{
				m1->SetValue(i,1,m0->GetValue(i,1));
				m1->SetValue(i,2,m0->GetValue(i,2));
			}
		}
		else {
			if(nptflag->GetValue(i,0)==0){
				if(fabs(m1->GetValue(i,0)-m1->GetValue(i,1))>1){
					m1->SetValue(i,1,m1->GetValue(i,0));
				}
				if(fabs(m1->GetValue(i,0)-m1->GetValue(i,2))>1){
					m1->SetValue(i,2,m1->GetValue(i,0));
				}
			}
			else{
				m1->SetValue(i,1,m0->GetValue(i,1));
				m1->SetValue(i,2,m0->GetValue(i,2));
			}
		}
	}

				/* second layer: water */
//	for (int i=0;i<mx*my;i++){
//		if(bathy->GetValue(i,0)==0){
//			u=double (rand())/ double(RAND_MAX);
//			y=signe(u-0.5)*T*(pow(1+1/T,fabs(2*u-1))-1);
//			m1->SetValue(i,2,m0->GetValue(i,2)+y*ptval);
//			if(m1->GetValue(i,2)<=m1->GetValue(i,1)){
//				m1->SetValue(i,2,m1->GetValue(i,1)+1e-10);
//			}
//			if(m1->GetValue(i,2)>=m1->GetValue(i,1)+mmax2){
//				m1->SetValue(i,2,m1->GetValue(i,1)+mmax2);
//			}
//		}
//	}
//	m1->ExtractColumn(m1col,2);
//	vec2gridsimple(m1col,m1gr,mx,my);
//	filtergrav(m1grsm,m1gr,ctr,sd,mx,my);
//	reshape(m1grsm,m1col,mx,my);
//	for (int i=0;i<mx*my;i++){
//		if(bathy->GetValue(i,0)==0){
//			m1->SetValue(i,2,m1col->GetValue(i,0));
//		}
//		else{
//			m1->SetValue(i,2,m0->GetValue(i,2));
//		}
//		if(m1->GetValue(i,2)<=m1->GetValue(i,1)){
//			m1->SetValue(i,2,m1->GetValue(i,1)+1e-10);
//		}
//	}
	delete m1gr;
	delete m1grsm;
	delete m1col;
	delete m1gr2;
	delete m1grsm2;
	delete m1col2;
	delete nptflag;
}/*}}}*/
double signe(double a){/*{{{*/
	if(a<0){return -1;}
	else{return 1;}
}/*}}}*/
void filtergrav(Matrix* A,Matrix* Ain,double ctr,double sd,int mx,int my){/*{{{*/
	A->MatrixEqual(Ain);
	for (int i=1;i<my-1;i++){
		for(int j=1;j<mx-1;j++){
			A->SetValue(i,j,(ctr*Ain->GetValue(i,j)+sd*(Ain->GetValue(i-1,j)+Ain->GetValue(i+1,j)+Ain->GetValue(i,j-1)+Ain->GetValue(i,j+1)))/(ctr+4*sd));
		}
	}
}/*}}}*/
double coolshed(double T0,double k,double c,double D){/*{{{*/
	double T1=T0*exp(-c*pow(k,1/D));
	return T1;
}/*}}}*/
