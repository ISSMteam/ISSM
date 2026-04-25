/*!\file: love_core.cpp
 * \brief: core of the LOVE numbers solution
 */

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"
#include "petscblaslapack.h"
#ifdef _HAVE_MPLAPACK_
#include <quadmath.h>
#include <iostream>
#include "mpblas__Float128.h"
#include "mplapack__Float128.h"
#include "mplapack_utils__Float128.h"
#endif

template<typename doubletype> IssmDouble DownCastVarToDouble(doubletype var); // pure declaration

template <> IssmDouble DownCastVarToDouble<IssmDouble>(IssmDouble var){
	return var;
}
template <> IssmDouble DownCastVarToDouble<IssmComplex>(IssmComplex var){
	return std::real(var);
}

template<typename doubletype> IssmDouble DownCastImagVarToDouble(doubletype var); // pure declaration

template <> IssmDouble DownCastImagVarToDouble<IssmDouble>(IssmDouble var){
	return var*0;
}
template <> IssmDouble DownCastImagVarToDouble<IssmComplex>(IssmComplex var){
	return std::imag(var);
}

#ifdef _HAVE_MPLAPACK_
template <> IssmDouble DownCastVarToDouble<__float128>(__float128 var){
	return static_cast<IssmDouble>(var);
}
template <> IssmDouble DownCastImagVarToDouble<__float128>(__float128 var){
	return static_cast<IssmDouble>(var*0);
}
__float128 pow(__float128 x, int y){
	return powq(x,y);
}
__float128 pow(__float128 x, double y){
	return powq(x,y);
}
__float128 pow(double x, __float128 y){
	return powq(x,y);
}

ostream& operator<<(ostream& os, __float128 x){
	char buf[128];
	quadmath_snprintf(buf, sizeof(buf), "%.34Qf", x);

    	os << buf;
    	return os;
}
#endif

/*local definitions:*/
template <class doubletype> class LoveVariables{  /*{{{*/

	public: 
		doubletype g0; 
		doubletype r0;
		doubletype* EarthMass; 
		int nyi, ifreq, nfreq; 
		int starting_layer;
		int* deg_layer_delete;
		int* nstep;
		doubletype* mu;
		doubletype* la;

		LoveVariables(){  /*{{{*/
			g0=0;
			r0=0;
			EarthMass=NULL;
			mu=NULL;
			la=NULL;
			nyi=0;
			nfreq=0;
			ifreq=0;
			starting_layer=0;
			nstep=NULL;
			deg_layer_delete=NULL;
		} /*}}}*/
		LoveVariables(doubletype* EarthMassin,doubletype g0in,doubletype r0in,int nyiin,int starting_layerin, int* deg_layer_deletein, int* nstepin){  /*{{{*/
			EarthMass=EarthMassin;
			g0=g0in;
			r0=r0in;
			nyi=nyiin;
			starting_layer=starting_layerin;
			deg_layer_delete=deg_layer_deletein;
			nstep=nstepin;
			mu=NULL;
			la=NULL;
			nfreq=0;
			ifreq=0;
		} /*}}}*/
		~LoveVariables(){
			xDelete<int>(deg_layer_delete);
			xDelete<int>(nstep);
			if (EarthMass) xDelete<doubletype>(EarthMass);
			if(mu)	xDelete<doubletype>(mu);
			if(la)	xDelete<doubletype>(la);
		};
}; /*}}}*/

template <class doubletype> class LoveNumbers{  /*{{{*/

	public:	
		doubletype* H;
		doubletype* K;
		doubletype* L;
		doubletype* Kernels;
		int sh_nmin, sh_nmax, nfreq, nkernels, lower_row, nfreqtotal; 

		LoveNumbers(){  /*{{{*/
			H=NULL;
			K=NULL;
			L=NULL;
			Kernels=NULL;
			sh_nmin=0;
			sh_nmax=0;
			nfreq=0;
			nkernels=0;
		} /*}}}*/
		LoveNumbers(int sh_nminin, int sh_nmaxin, int nfreqin, int lower_rowin,int nfreqtotalin,Matlitho* matlitho){  /*{{{*/
			sh_nmax=sh_nmaxin;
			sh_nmin=sh_nminin;
			nfreq=nfreqin;
			lower_row=lower_rowin;
			nfreqtotal=nfreqtotalin;
			nkernels=(sh_nmax+1)*(matlitho->numlayers+1)*6;
			H=xNewZeroInit<doubletype>(nfreq*(sh_nmax+1));
			K=xNewZeroInit<doubletype>(nfreq*(sh_nmax+1));
			L=xNewZeroInit<doubletype>(nfreq*(sh_nmax+1));
			Kernels=xNewZeroInit<doubletype>(nfreq*nkernels);
		} /*}}}*/
		void DownCastToDouble(LoveNumbers<IssmDouble>* LoveDouble){
			for(int i=0;i<(sh_nmax+1)*nfreq;i++){
				LoveDouble->H[i]=DownCastVarToDouble<doubletype>(H[i]);
				LoveDouble->K[i]=DownCastVarToDouble<doubletype>(K[i]);
				LoveDouble->L[i]=DownCastVarToDouble<doubletype>(L[i]);
			}
			for(int i=0;i<nkernels*nfreq;i++){
				LoveDouble->Kernels[i]=DownCastVarToDouble<doubletype>(Kernels[i]);
			}

		}
		void DownCastImagToDouble(LoveNumbers<IssmDouble>* LoveImag){
			for(int i=0;i<(sh_nmax+1)*nfreq;i++){
				LoveImag->H[i]=DownCastImagVarToDouble<doubletype>(H[i]);
				LoveImag->K[i]=DownCastImagVarToDouble<doubletype>(K[i]);
				LoveImag->L[i]=DownCastImagVarToDouble<doubletype>(L[i]);
			}
			for(int i=0;i<nkernels*nfreq;i++){
				LoveImag->Kernels[i]=DownCastImagVarToDouble<doubletype>(Kernels[i]);
			}

		}
		void LoveMPI_Gather(LoveNumbers<doubletype>* Love_local, int lower_row){
			int* recvcounts=xNew<int>(IssmComm::GetSize());
			int* displs=xNew<int>(IssmComm::GetSize());
			int  rc;
			int  offset;
			int nf_local = Love_local->nfreq;

			/*Deal H, K, L first, as they share the same size*/
			rc=(sh_nmax+1)*nf_local;
			offset=(sh_nmax+1)*lower_row;
			ISSM_MPI_Allgather(&rc,1,ISSM_MPI_INT,recvcounts,1,ISSM_MPI_INT,IssmComm::GetComm());
			ISSM_MPI_Allgather(&offset,1,ISSM_MPI_INT,displs,1,ISSM_MPI_INT,IssmComm::GetComm());
			ISSM_MPI_Allgatherv(Love_local->H, rc, ISSM_MPI_DOUBLE, H, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());
			ISSM_MPI_Allgatherv(Love_local->K, rc, ISSM_MPI_DOUBLE, K, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());
			ISSM_MPI_Allgatherv(Love_local->L, rc, ISSM_MPI_DOUBLE, L, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());

			/* deal with love kernels now */
			rc=nf_local*nkernels;
			offset=lower_row*nkernels;
			ISSM_MPI_Allgather(&rc,1,ISSM_MPI_INT,recvcounts,1,ISSM_MPI_INT,IssmComm::GetComm());
			ISSM_MPI_Allgather(&offset,1,ISSM_MPI_INT,displs,1,ISSM_MPI_INT,IssmComm::GetComm());
			ISSM_MPI_Allgatherv(Love_local->Kernels, rc, ISSM_MPI_DOUBLE, Kernels, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());

			xDelete<int>(recvcounts);
			xDelete<int>(displs);
		}
		void Broadcast(void){
			//Intended only for IssmDouble type
			Vector<IssmDouble>* vH;
			Vector<IssmDouble>* vK;
			Vector<IssmDouble>* vL;

			vH= new Vector<IssmDouble>((sh_nmax+1)*nfreqtotal);
			for(int i=0;i<(sh_nmax+1)*nfreq;i++) vH->SetValue(lower_row*(sh_nmax+1)+i,DownCastVarToDouble<doubletype>(H[i]),INS_VAL);
			xDelete(H);
			vH->Assemble();
			H=vH->ToMPISerial();
			delete vH;

			vK= new Vector<IssmDouble>((sh_nmax+1)*nfreqtotal);
			for(int i=0;i<(sh_nmax+1)*nfreq;i++) vK->SetValue(lower_row*(sh_nmax+1)+i,DownCastVarToDouble<doubletype>(K[i]),INS_VAL);
			xDelete(K);
			vK->Assemble();
			K=vK->ToMPISerial();
			delete vK;

			vL= new Vector<IssmDouble>((sh_nmax+1)*nfreqtotal);
			for(int i=0;i<(sh_nmax+1)*nfreq;i++) vL->SetValue(lower_row*(sh_nmax+1)+i,DownCastVarToDouble<doubletype>(L[i]),INS_VAL);
			xDelete(L);
			vL->Assemble();
			L=vL->ToMPISerial();
			delete vL;
		}
		void KernelBroadcast(void){
			Vector<IssmDouble>* vKernels;
			vKernels= new Vector<IssmDouble>(nkernels*nfreqtotal);
			for(int i=0;i<nkernels*nfreq;i++){
				vKernels->SetValue(lower_row*nkernels+i,DownCastVarToDouble<doubletype>(Kernels[i]),INS_VAL);
			}
			xDelete(Kernels);
			vKernels->Assemble();
			Kernels=vKernels->ToMPISerial();
			delete vKernels;
		}
		void Copy(LoveNumbers<doubletype>* LoveDup){
			for(int i=0;i<(sh_nmax+1)*nfreq;i++){
				H[i]=LoveDup->H[i];
				K[i]=LoveDup->K[i];
				L[i]=LoveDup->L[i];
			}
			for(int i=0;i<nkernels*nfreq;i++){
				Kernels[i]=LoveDup->Kernels[i];
			}
		}
		~LoveNumbers(){
			xDelete<doubletype>(H);
			xDelete<doubletype>(K);
			xDelete<doubletype>(L);
			xDelete<doubletype>(Kernels);
		};
};

 /*}}}*/

/*self contained support routines used by cores below:*/
template<typename doubletype> doubletype                 angular_frequency(IssmDouble frequency); //pure declaration
template <> IssmDouble                     angular_frequency<IssmDouble>(IssmDouble frequency){ /*{{{*/
	return 2.0*PI*frequency;
} /*}}}*/
template <> IssmComplex                    angular_frequency<IssmComplex>(IssmDouble frequency){ /*{{{*/
	IssmComplex value=reCast<IssmComplex>(complex<double>(0,1))*2.0*PI*reCast<IssmComplex>(frequency);
	return value;
} /*}}}*/
#ifdef _HAVE_MPLAPACK_
template <> __float128                     angular_frequency<__float128>(IssmDouble frequency){ /*{{{*/
	return 2.0*PI*frequency;
} /*}}}*/
#endif
template<typename doubletype> void                       allgesv(int* pnyi, int* pnrhs, doubletype* yilocal, int* plda, int* ipiv, doubletype* rhslocal, int* pldb, int* pinfo);
template <> void                           allgesv<IssmDouble>(int* pnyi, int* pnrhs, IssmDouble* yilocal, int* plda, int* ipiv, IssmDouble* rhslocal, int* pldb, int* pinfo){ /*{{{*/
	dgesv_(pnyi, pnrhs, yilocal, plda, ipiv, rhslocal, pldb, pinfo);
} /*}}}*/
template <> void                           allgesv<IssmComplex>(int* pnyi, int* pnrhs, IssmComplex* yilocal, int* plda, int* ipiv, IssmComplex* rhslocal, int* pldb, int* pinfo){ /*{{{*/
	//_error_("zgesv_ not linked correctly yet! ");
	Zgesvx(pnyi, pnrhs, yilocal, plda, ipiv, rhslocal, pldb, pinfo);
} /*}}}*/
#ifdef _HAVE_MPLAPACK_
template <> void                           allgesv<__float128>(int* pnyi, int* pnrhs, __float128* yilocal, int* plda, int* ipiv, __float128* rhslocal, int* pldb, int* pinfo){ /*{{{*/
	mplapackint nyi=*pnyi;
	mplapackint nrhs=*pnrhs;
	mplapackint lda=*plda;
	mplapackint* qipiv=NULL;
	mplapackint ldb=*pldb;
	mplapackint info=0;
	qipiv=xNewZeroInit<mplapackint>(*pnyi);
	
	Rgesv(nyi, nrhs, yilocal, lda, qipiv, rhslocal, ldb, info);

	for (int i;i=0;i<*pnyi) ipiv[i]=qipiv[i];
	*pinfo=info;
	xDelete<mplapackint>(qipiv);
} /*}}}*/
#endif

template<typename doubletype> doubletype   factorial(int n){ /*{{{*/
	doubletype prod=1;
	for (int i=2;i<n+1;i++) prod*=i;
	return prod;
}/*}}}*/
template<typename doubletype> doubletype   n_C_r(int n, int r){ /*{{{*/ 
	//n choose r
	int primes[169] = 
	{2,    3,    5,    7,   11,   13,   17,   19,   23,   29,  
		31,   37,   41,   43,   47,   53,   59,   61,   67,   71,  
		73,   79,   83,   89,   97,  101,  103,  107,  109,  113,  
		127,  131,  137,  139,  149,  151,  157,  163,  167,  173,  
		179,  181,  191,  193,  197,  199,  211,  223,  227,  229,  
		233,  239,  241,  251,  257,  263,  269,  271,  277,  281,  
		283,  293,  307,  311,  313,  317,  331,  337,  347,  349,  
		353,  359,  367,  373,  379,  383,  389,  397,  401,  409,  
		419,  421,  431,  433,  439,  443,  449,  457,  461,  463,  
		467,  479,  487,  491,  499,  503,  509,  521,  523,  541,  
		547,  557,  563,  569,  571,  577,  587,  593,  599,  601,  
		607,  613,  617,  619,  631,  641,  643,  647,  653,  659,  
		661,  673,  677,  683,  691,  701,  709,  719,  727,  733,  
		739,  743,  751,  757,  761,  769,  773,  787,  797,  809,  
		811,  821,  823,  827,  829,  839,  853,  857,  859,  863,  
		877,  881,  883,  887,  907,  911,  919,  929,  937,  941,  
		947,  953,  967,  971,  977,  983,  991,  997, 1009};
	int num, den;
	num = 1;
	den = 1;

	for (int i=0;i<r;i++){
		num = num*(n-i);
		den = den*(i+1);
		if (i>0) {
			// Divide out common prime factors
			for (int k=0;k<169;k++){ //169 is the length of the prime vector here
				if ( i % primes[k] == 0) { // modulo
					num = num/primes[k];
					den = den/primes[k];
				}
			}
		}
	}

	doubletype res;        
	return res = num/den;
}/*}}}*/
template<typename doubletype> doubletype*  postwidder_coef(int NTit){ /*{{{*/
	//Coefficients of the Post-Widder method through Saltzer summation for inverse Laplace transform:
	//The Mth iteration estimate will be: f(t)_M = sum_{k=1:2*M}(xi_[M,k] * f(s_k))
	//The method is based on equations (2), (6), (7) in: 
	//Valko PP, Abate J. Comparison of sequence accelerators for the Gaver method of numerical Laplace transform inversion. Computational Mathematics and Applications. (2004)
	//Note that the coefficients xi lack the factor s=k*log(2)/t. 
	//That is because we are computing the heaviside response of the system rather than its impulse response, 
	//and Laplace_Transform(Heaviside(t).*f(t)) = f(s)/s. So s cancels out in the sum for f(t)_M.

	doubletype* xi=xNewZeroInit<doubletype>(2*NTit*NTit);
	int indxi;
	for (int M=1;M<NTit+1;M++){
		for (int k=1;k<2*M+1;k++){
			indxi=(M-1)*(2*NTit)+k-1;
			for (int j=floor((k+1)/2);j<min(k,M)+1;j++){
				xi[indxi]+=pow(j,M+1.0)/factorial<doubletype>(M)*n_C_r<doubletype>(M,j)*n_C_r<doubletype>(2*j,j)*n_C_r<doubletype>(j,k-j);
			}
			xi[indxi]*=pow(-1.0,k+M)/k;
		}
	}
	return xi;
}/*}}}*/
template<typename doubletype> void         postwidder_transform(doubletype* Lovet, doubletype* Lovef,int d, int t, int sh_nmax,int NTit, doubletype* xi, FemModel* femmodel){ /*{{{*/
	//Computes Lovet for time step t and degree d from the PW coefficients xi and the corresponding 2*NTit frequency samples in Lovef

	int indxi, indf;
	doubletype PW_test;
	IssmDouble PW_threshold;
	femmodel->parameters->FindParam(&PW_threshold,LovePostWidderThresholdEnum);

	indf=(t*2*NTit)*(sh_nmax+1)+d;
	doubletype* LoveM = NULL;
	doubletype zero = 0;


	// test variation across frequencies tested, something with little frequency dependence is not worth going through PW tranform
	// that transform would also be numerically unstable
	PW_test = abs((Lovef[indf+(2*NTit-1)*(sh_nmax+1)]-Lovef[indf])/Lovef[indf]); 

	if (PW_test==zero){ //elastic or fluid response: Love(t) = Love(s), we can do an early return
		Lovet[t*(sh_nmax+1)+d]=Lovef[indf];
		return;
	}

	LoveM=xNewZeroInit<doubletype>(NTit);

	for (int M=1;M<NTit+1;M++){
		LoveM[M-1]=0.0;
		for (int k=1;k<2*M+1;k++){
			indxi=(M-1)*(2*NTit)+k-1;
			LoveM[M-1]+=xi[indxi]*Lovef[indf+(k-1)*(sh_nmax+1)];
		}

		//Make sure we are not getting into numerical instability
		//Diverging once: ok, we'll give that the benefit of the doubt, it can be an inflexion point in the convergence series
		//Diverging twice in a row: we are definitely propagating numerical error: revert to the last stable value and exit
		if (M>3){ 
			if ( abs(LoveM[M-1]-LoveM[M-2]) > abs(LoveM[M-2]-LoveM[M-3]) &&
					abs(LoveM[M-2]-LoveM[M-3]) > abs(LoveM[M-3]-LoveM[M-4]) ){
				Lovet[t*(sh_nmax+1)+d]=LoveM[M-3];
				xDelete<doubletype>(LoveM);
				return;
			}
		}
	}
	Lovet[t*(sh_nmax+1)+d]=LoveM[NTit-1];
	xDelete<doubletype>(LoveM);
}/*}}}*/

template <typename doubletype> doubletype HypergeomTableLookup(doubletype z1, doubletype alpha, IssmDouble* h1, IssmDouble* z, int nz, int nalpha){/*{{{*/
	int iz1, iz2;	
	doubletype hf,h00,h10, h01, h11, za, zd, ha, hb,hc,hd, m0,m1,t;
	doubletype dalpha=1.0/(nalpha-1); // alpha table resolution given 0 <= alpha <= 1
	int        ialpha  = static_cast<int>(DownCastVarToDouble(alpha/dalpha));
	doubletype lincoef = alpha/dalpha - std::floor(DownCastVarToDouble(alpha/dalpha));//linear fraction in [0;1] for alpha interpolation
	iz1=nz;
	for (int i=0;i<nz;i++){
		if (abs(z[i])>abs(z1)) {
			iz1=i-1;
			break;
		}
	}

	if (iz1<0){
		//1-hf for very small abs(z) tends to 0, and is very log-linear with respect to log(z), so we can simply extrapolate the value of hf via the loglog slope
		hf=(1.0-lincoef)*h1[ialpha*nz+0]+lincoef*h1[(ialpha+1)*nz+0];
		hf=1.0- (1.0-hf)*pow(10.0,(log10(abs(z1))-log10(abs(z[0]))));
		//hf=1.0;
	}
	else if (iz1==nz){
		//hf for very large abs(z) tends to 0, and is very log-linear with respect to log(z), so we can simply extrapolate the value of hf via the loglog slope
		hf=(1.0-lincoef)*h1[ialpha*nz+nz-1]+lincoef*h1[(ialpha+1)*nz+nz-1];
		hf=hf *pow(10.0,-(log10(abs(z1))-log10(abs(z[nz-1]))));

		//hf=0;
	}
	else{ //cubic spline interpolation
		//edge cases: extrapolate 1 point
		if (iz1==0){
			za=2.0*z[0]-z[1];
			ha=(1.0-lincoef)*h1[ialpha*nz+0]+lincoef*h1[(ialpha+1)*nz+0];
			ha=1.0- (1.0-ha) *pow(10.0,log10(abs(za))-log10(abs(z[0])));
		} 
		else {
			za=z[iz1-1];
			ha=(1.0-lincoef)*h1[ialpha*nz+iz1-1] + lincoef*h1[(ialpha+1)*nz+iz1-1];
		}

		if (iz1==nz-2){
			zd=2.0*z[nz-1]-z[nz-2];
			hd=(1.0-lincoef)*h1[ialpha*nz+nz-1]+lincoef*h1[(ialpha+1)*nz+nz-1];
			hd=hd *pow(10.0,-(log10(abs(zd))-log10(abs(z[nz-1]))));
		} 
		else {
			zd=z[iz1+2];
			hd=(1.0-lincoef)*h1[ialpha*nz+iz1+2]+lincoef*h1[(ialpha+1)*nz+iz1+2];
		}

		hb=(1.0-lincoef)*h1[ialpha*nz+iz1] +lincoef*h1[(ialpha+1)*nz+iz1];
		hc=(1.0-lincoef)*h1[ialpha*nz+iz1+1] +lincoef*h1[(ialpha+1)*nz+iz1+1];

		//left derivative
		m0= 0.5*(z[iz1+1]-z[iz1])*((hc-hb)/(z[iz1+1]-z[iz1]) + (hb-ha)/(z[iz1]-za));
		//right derivative
		m1= 0.5*(z[iz1+1]-z[iz1])*((hd-hc)/(zd-z[iz1+1]) + (hc-hb)/(z[iz1+1]-z[iz1]));

		//interpolation abscissa
		t=(z1-z[iz1])/(z[iz1+1]-z[iz1]);
		
		//cubic spline functions
		h00=2.*pow(t,3)-3.*pow(t,2)+1.;
		h10=pow(t,3)-2.*pow(t,2)+t;
		h01=-2.*pow(t,3)+3.*pow(t,2);
		h11=pow(t,3)-pow(t,2);

		hf=h00*hb + h10*m0 + h01*hc + h11*m1;
	}
	return hf;

}/*}}}*/

template <typename doubletype> doubletype muEBM(int layer_index, doubletype omega, Matlitho* matlitho, FemModel* femmodel); //pure declaration
template <> IssmComplex muEBM<IssmComplex>(int layer_index, IssmComplex omega, Matlitho* matlitho, FemModel* femmodel){/*{{{*/
	// Initialization
	int nz, nalpha, dummy1, dummy2;
	IssmComplex mu;
	IssmDouble* z=NULL;
	IssmDouble* h1=NULL;
	IssmDouble* h2=NULL;
	IssmComplex hf11, hf12, hf21, hf22;
	IssmDouble  factor=0;
	IssmComplex z1, z2;
	IssmComplex U1, U2;
	IssmComplex j=reCast<IssmComplex>(complex<double>(0,1));
	//Matlitho parameters
	IssmDouble alpha=matlitho->ebm_alpha[layer_index];
	IssmDouble delta=matlitho->ebm_delta[layer_index];
	IssmDouble taul=matlitho->ebm_taul[layer_index];
	IssmDouble tauh=matlitho->ebm_tauh[layer_index];
	IssmDouble vi=matlitho->viscosity[layer_index];
	IssmDouble mu0=matlitho->lame_mu[layer_index];
	//fetch hypergeometric function tables and parameters
	femmodel->parameters->FindParam(&nz,LoveHypergeomNZEnum);
	femmodel->parameters->FindParam(&nalpha,LoveHypergeomNAlphaEnum);
	femmodel->parameters->FindParam(&z,&dummy1,LoveHypergeomZEnum);
	femmodel->parameters->FindParam(&h1,&dummy1,&dummy2,LoveHypergeomTable1Enum);
	femmodel->parameters->FindParam(&h2,&dummy1,&dummy2,LoveHypergeomTable2Enum);
	omega=omega/j;

	z1= -pow(omega*tauh,2.0);
	z2= -pow(omega*taul,2.0);
	//Table1 h1 should be 2F1([1 1+alpha], [2+alpha/2], z)
	//Table2 h2 should be 2F1([1 0.5+alpha], [1.5+alpha/2], z)
	hf11=HypergeomTableLookup<IssmComplex>(z1, alpha, h1, z, nz, nalpha);
	hf21=HypergeomTableLookup<IssmComplex>(z1, alpha, h2, z, nz, nalpha);
	hf12=HypergeomTableLookup<IssmComplex>(z2, alpha, h1, z, nz, nalpha);
	hf22=HypergeomTableLookup<IssmComplex>(z2, alpha, h2, z, nz, nalpha);

	//Ivins et al (2020) p11
	U1=(pow(tauh,alpha)-pow(taul,alpha))/alpha-pow(omega,2.0)/(2.0+alpha)*(pow(tauh,2.0+alpha)*hf11-pow(taul,2.0+alpha)*hf12);
	U2=(pow(tauh,1.0+alpha)*hf21-pow(taul,1.0+alpha)*hf22)/(1.0+alpha);

	factor= alpha*delta/(pow(tauh,alpha)-pow(taul,alpha));
	U1=(1.0+factor*U1);
	U2=(factor*omega*U2 +mu0/vi/omega);
	mu=mu0*(U1+j*U2)/(pow(U1,2.0)+pow(U2,2.0));
	omega=omega*j;

	xDelete<IssmDouble>(z);
	xDelete<IssmDouble>(h1);
	xDelete<IssmDouble>(h2);
	return mu;
}/*}}}*/

template <> IssmDouble muEBM<IssmDouble>(int layer_index, IssmDouble omega, Matlitho* matlitho, FemModel* femmodel){/*{{{*/
	// Initialization
	int nz, nalpha, dummy1, dummy2;
	IssmDouble mu;
	IssmDouble* z=NULL;
	IssmDouble* h1=NULL;
	IssmDouble hf11, hf12;
	IssmDouble  factor, B, D, z1, z2;
	//Matlitho parameters
	IssmDouble alpha=matlitho->ebm_alpha[layer_index];
	IssmDouble delta=matlitho->ebm_delta[layer_index];
	IssmDouble taul=matlitho->ebm_taul[layer_index];
	IssmDouble tauh=matlitho->ebm_tauh[layer_index];
	IssmDouble vi=matlitho->viscosity[layer_index];
	IssmDouble mu0=matlitho->lame_mu[layer_index];
	//fetch hypergeometric function tables and parameters
	femmodel->parameters->FindParam(&nz,LoveHypergeomNZEnum);
	femmodel->parameters->FindParam(&nalpha,LoveHypergeomNAlphaEnum);
	femmodel->parameters->FindParam(&z,&dummy1,LoveHypergeomZEnum);
	femmodel->parameters->FindParam(&h1,&dummy1,&dummy2,LoveHypergeomTable1Enum);

	z1=-omega*tauh;
	z2=-omega*taul;
	//Table1 h1 should be 2F1([1 1+alpha], [2+alpha], z)
	hf11=HypergeomTableLookup<IssmDouble>(z1, alpha, h1, z, nz, nalpha);
	hf12=HypergeomTableLookup<IssmDouble>(z2, alpha, h1, z, nz, nalpha);

	//Ivins et al. (2022) p1979
	factor= alpha*delta/(pow(tauh,alpha)-pow(taul,alpha));
	B= factor/(1.0+alpha) *mu0/vi * (pow(tauh,1.0+alpha)*hf11 - pow(taul,1.0+alpha)*hf12);
	D= omega*vi/mu0* 1.0/(1.0+omega*vi/mu0*(1.0+delta) -pow(omega*vi/mu0,2.0)*B);

	xDelete<IssmDouble>(z);
	xDelete<IssmDouble>(h1);
	return mu=mu0*D;
}/*}}}*/
#ifdef _HAVE_MPLAPACK_
template <> __float128 muEBM<__float128>(int layer_index, __float128 omega, Matlitho* matlitho, FemModel* femmodel){/*{{{*/
	// Initialization
	int nz, nalpha, dummy1, dummy2;
	IssmDouble* z=NULL;
	IssmDouble* h1=NULL;
	__float128 mu;
	__float128 hf11, hf12;
	__float128  factor, B, D, z1, z2;
	//Matlitho parameters
	__float128 alpha=matlitho->ebm_alpha[layer_index];
	__float128 delta=matlitho->ebm_delta[layer_index];
	__float128 taul=matlitho->ebm_taul[layer_index];
	__float128 tauh=matlitho->ebm_tauh[layer_index];
	__float128 vi=matlitho->viscosity[layer_index];
	__float128 mu0=matlitho->lame_mu[layer_index];
	//fetch hypergeometric function tables and parameters
	femmodel->parameters->FindParam(&nz,LoveHypergeomNZEnum);
	femmodel->parameters->FindParam(&nalpha,LoveHypergeomNAlphaEnum);
	femmodel->parameters->FindParam(&z,&dummy1,LoveHypergeomZEnum);
	femmodel->parameters->FindParam(&h1,&dummy1,&dummy2,LoveHypergeomTable1Enum);

	z1=-(omega*tauh);
	z2=-(omega*taul);
	//Table1 h1 should be 2F1([1 1+alpha], [2+alpha], z)
	hf11=HypergeomTableLookup<__float128>(z1, alpha, h1, z, nz, nalpha);
	hf12=HypergeomTableLookup<__float128>(z2, alpha, h1, z, nz, nalpha);

	//Ivins et al. (2022) p1979
	//Note: therein, mu(s') = s'*mu~(s'); s'=omega*tauM=omega*vi/mu0
	factor= alpha*delta/(pow(tauh,alpha)-pow(taul,alpha));
	B= factor/(1.0q+alpha) *mu0/vi * (pow(tauh,1.0q+alpha)*hf11 - pow(taul,1.0q+alpha)*hf12);
	D= omega*vi/mu0* 1.0q/(1.0q+omega*vi/mu0*(1.0q+delta) -pow(omega*vi/mu0,2.0q)*B);

	xDelete<IssmDouble>(z);
	xDelete<IssmDouble>(h1);
	return mu=mu0*D;
}/*}}}*/
#endif
template <typename doubletype> void        GetEarthRheology(doubletype* pla, doubletype* pmu, int layer_index, doubletype omega,  Matlitho* matlitho, FemModel* femmodel){ /*{{{*/

	//returns lame parameters (material rigity) lambda and mu for the right frequency and layer
	doubletype mu,la;

	doubletype vi=matlitho->viscosity[layer_index];
	doubletype mu0=matlitho->lame_mu[layer_index];
	doubletype la0=matlitho->lame_lambda[layer_index];
	int rheo=matlitho->rheologymodel[layer_index];
	doubletype zero = 0;

	if (vi!=zero && omega!=zero){ //take into account viscosity in the rigidity if the material isn't a perfect fluid
		doubletype ka=la0 + 2.0/3.0*mu0; //Bulk modulus
		if (rheo==2){//EBM
			mu=muEBM<doubletype>(layer_index, omega, matlitho, femmodel);
			la=ka-2.0/3.0*mu;
		} 
		else if (rheo==1){//Burgers
			doubletype vi2=matlitho->burgers_viscosity[layer_index];
			doubletype mu2=matlitho->burgers_mu[layer_index];

			mu=mu0*omega*(omega+mu2/vi2)/((omega+mu2/vi2)*(omega+mu0/vi)+mu0/vi2*omega);
			la=ka-2.0/3.0*mu;
		}
		else{//Maxwell
			la = (la0 + mu0*ka/vi/omega)/(1.0 + mu0/vi/omega);
			mu = mu0/(1.0+mu0/vi/omega);
		}
	}
	else{//Otherwise return the elastic value
	la=la0;
	mu=mu0;
	}

	*pla=la;
	*pmu=mu;

} /*}}}*/

template <typename doubletype> void        EarthRheology(LoveVariables<doubletype>* vars, IssmDouble* frequencies, int nfreq,  Matlitho* matlitho, FemModel* femmodel){/*{{{*/
	doubletype omega;
	//reset pointers to NULL if this function was previously called
	if(vars->mu)	xDelete<doubletype>(vars->mu);
	if(vars->la)	xDelete<doubletype>(vars->la);
	//precompute rheology at the requested frequencies
	vars->mu=xNewZeroInit<doubletype>(nfreq*matlitho->numlayers);
	vars->la=xNewZeroInit<doubletype>(nfreq*matlitho->numlayers);
	vars->nfreq=nfreq;
	for (int i=0;i<matlitho->numlayers;i++){
		for (int fr=0;fr<nfreq;fr++){
			omega=angular_frequency<doubletype>(frequencies[fr]);
			GetEarthRheology<doubletype>(&vars->la[i*nfreq+fr], &vars->mu[i*nfreq+fr], i,omega,matlitho, femmodel);
		}
	}
}/*}}}*/

template <typename doubletype> doubletype	GetGravity(doubletype r2, int layer_index, FemModel* femmodel, Matlitho* matlitho,LoveVariables<doubletype>* vars){ /*{{{*/
	//computes gravity at radius r2
	doubletype* EarthMass;
	doubletype g, GG;
	IssmDouble GGp;

	EarthMass=vars->EarthMass;
	femmodel->parameters->FindParam(&GGp,LoveGravitationalConstantEnum);
	GG=GGp;
	doubletype ro=matlitho->density[layer_index];
	doubletype M=0;
	doubletype r1=0;
	if (layer_index==0){
		M=4.0/3.0*PI*ro*pow(r2,3.0);
	}
	else{ 
		r1=matlitho->radius[layer_index];
		M=EarthMass[layer_index-1]+4.0/3.0*PI*ro*(pow(r2,3.0)-pow(r1,3.0));
	}
	return	g= GG*M/pow(r2,2.0);
}/*}}}*/
template <typename doubletype> void        fill_yi_prefactor(doubletype* yi_prefactor, int* pdeg, doubletype* pomega, FemModel* femmodel, Matlitho* matlitho, LoveVariables<doubletype>* vars){ /*{{{*/
	//precalculates partial derivative factors for function yi_derivatives
	doubletype ra=matlitho->radius[matlitho->numlayers];
	doubletype  g0,r0,mu0;
	IssmDouble mu0p, GG;
	int nstep,nsteps,nindex, starting_layer;

	femmodel->parameters->FindParam(&mu0p,LoveMu0Enum);
	femmodel->parameters->FindParam(&GG,LoveGravitationalConstantEnum);

	g0=vars->g0;
	r0=vars->r0;
	mu0=mu0p;
	starting_layer=vars->starting_layer;

	doubletype frh,frhg0,fgr0,fgr,fn,rm0,rlm,flm;
	doubletype xmin,xmax,x,dr;
	doubletype g,ro;
	bool       issolid;

	if (pomega) { //frequency and degree dependent terms /*{{{*/
		doubletype la,mu;
		doubletype f[12];
		int deg=*pdeg;
		doubletype omega=*pomega;	
		fn=deg*(deg+1.0);

		for (int layer_index=starting_layer;layer_index<matlitho->numlayers;layer_index++){
			nstep=vars->nstep[layer_index];
			nsteps=0;
			for (int i=0;i<layer_index;i++)	nsteps+=vars->nstep[i];

			ro=matlitho->density[layer_index];
			issolid=matlitho->issolid[layer_index];
			if(issolid==1){
				//GetEarthRheology<doubletype>(&la, &mu,layer_index,omega,matlitho);   
				mu=vars->mu[layer_index*vars->nfreq+vars->ifreq];
				la=vars->la[layer_index*vars->nfreq+vars->ifreq];

				/*_______Expressions*/
				flm=(la+2.0*mu);
				rlm=(3.0*la+2.0*mu)/(la+2.0*mu);
				rm0=mu/mu0;
				frh=ro*ra/mu0;

				f[0]=(-2.0*la/flm);
				f[1]=mu0/flm;
				f[2]=(la*fn/flm);
				f[3]=rm0*rlm;
				//f[4]=-ro*pow(omega,2.0)*ra*ra/mu0;
				f[4]=0;
				f[5]=(-4.0*mu/flm);
				f[6]=fn*frh;
				f[7]=-(2.0*rm0*rlm)*fn;
				f[8]=1.0/rm0;
				f[9]=-2.0*rm0*rlm;
				f[10]=-la/flm;
				f[11]=2.0*rm0*(la*(2.0*fn-1.0)+2.0*mu*(fn-1.0))/flm;

				xmin=matlitho->radius[layer_index]/ra;
				xmax=(matlitho->radius[layer_index+1])/ra;
				dr = (xmax -xmin)/reCast<doubletype>(nstep);
				x=xmin;

				for (int n=0;n<nstep;n++){

					g=GetGravity<doubletype>(x*ra,layer_index,femmodel,matlitho,vars);
					nindex=nsteps*36+n*36;
					yi_prefactor[nindex+ 0*6+0]= f[0]/x;                      // in dy[0*6+0]
					yi_prefactor[nindex+ 0*6+1]= f[1];                        // in dy[0*6+1]
					yi_prefactor[nindex+ 0*6+2]= f[2]/x;                      // in dy[0*6+2]
					yi_prefactor[nindex+ 1*6+0]= 4.0*(-frh*g+f[3]/x)/x + f[4];// in dy[1*6+0]
					yi_prefactor[nindex+ 1*6+1]= f[5]/x;                      // in dy[1*6+1]
					yi_prefactor[nindex+ 1*6+2]= (f[6]*g+f[7]/x)/x;           // in dy[1*6+2]
					yi_prefactor[nindex+ 2*6+3]= f[8];                        // in dy[2*6+3]
					yi_prefactor[nindex+ 3*6+0]= (frh*g+f[9]/x)/x;            // in dy[3*6+0]
					yi_prefactor[nindex+ 3*6+1]= f[10]/x;                     // in dy[3*6+1]
					yi_prefactor[nindex+ 3*6+2]= f[11]/(x*x) + f[4];          // in dy[3*6+2]
					x=x+dr;
				}
			}
		}
		/*}}}*/
	} else if (pdeg) { // degree dependent terms /*{{{*/
		int deg=*pdeg;
		fn=(deg*(deg+1.0));

		for (int layer_index=starting_layer;layer_index<matlitho->numlayers;layer_index++){
			nstep=vars->nstep[layer_index];
			nsteps=0;
			for (int i=0;i<layer_index;i++)	nsteps+=vars->nstep[i];

			ro=matlitho->density[layer_index];
			issolid=matlitho->issolid[layer_index];

			/*_______Expressions*/
			fgr=4.0*PI*GG*ro*ra;

			xmin=matlitho->radius[layer_index]/ra;
			xmax=(matlitho->radius[layer_index+1])/ra;
			dr = (xmax -xmin)/reCast<doubletype>(nstep);
			x=xmin;

			for (int n=0;n<nstep;n++){
				nindex=nsteps*36+n*36;
				g=GetGravity<doubletype>(x*ra,layer_index,femmodel,matlitho,vars);

				if(issolid){
					yi_prefactor[nindex+ 1*6+3]= fn/x;                  // in dy[1*6+3]
					yi_prefactor[nindex+ 5*6+2]= -(fgr/g0*fn)/x;        // in dy[5*6+2]
					yi_prefactor[nindex+ 5*6+4]= fn/(x*x);		     // in dy[5*6+4]
				} else {
					yi_prefactor[nindex+ 1*6+0]= (-4.0*(fgr/g)+fn/x)/x; // in dy[1*6+0] liquid layer
				}
				x=x+dr;
			}
		}
		/*}}}*/
	} else { // static terms /*{{{*/
		for (int layer_index=starting_layer;layer_index<matlitho->numlayers;layer_index++){
			nstep=vars->nstep[layer_index];
			nsteps=0;
			for (int i=0;i<layer_index;i++)	nsteps+=vars->nstep[i];

			ro=matlitho->density[layer_index];
			issolid=matlitho->issolid[layer_index];

			/*_______Expressions*/
			frhg0=ro*g0*ra/mu0;
			fgr=4.0*PI*GG*ro*ra;

			xmin=matlitho->radius[layer_index]/ra;
			xmax=(matlitho->radius[layer_index+1])/ra;
			dr = (xmax -xmin)/reCast<doubletype>(nstep);
			x=xmin;

			for (int n=0;n<nstep;n++){
				g=GetGravity<doubletype>(x*ra,layer_index,femmodel,matlitho,vars);
				nindex=nsteps*36+n*36;
				if(issolid){
					yi_prefactor[nindex+ 1*6+5]= -frhg0;       // in dy[1*6+5]
					yi_prefactor[nindex+ 2*6+0]= -1.0/x;       // in dy[2*6+0]
					yi_prefactor[nindex+ 2*6+2]= 1.0/x;        // in dy[2*6+2]
					yi_prefactor[nindex+ 3*6+3]= -3.0/x;       // in dy[3*6+3]
					yi_prefactor[nindex+ 3*6+4]= -frhg0/x;     // in dy[3*6+4]
					yi_prefactor[nindex+ 4*6+0]= fgr/g0;       // in dy[4*6+0]
					yi_prefactor[nindex+ 4*6+5]= 1.0;          // in dy[4*6+5]
					yi_prefactor[nindex+ 5*6+5]= -2.0/x;       // in dy[5*6+5]
				} else {
					yi_prefactor[nindex+ 0*6+0]= fgr/g;        // in dy[0*6+0] liquid layer
					yi_prefactor[nindex+ 0*6+1]= 1.0;          // in dy[0*6+1] liquid layer
					yi_prefactor[nindex+ 1*6+1]= -2.0/x-fgr/g; // in dy[1*6+1] liquid layer
				}
				x=x+dr;
			}
		}
		/*}}}*/
	}
}/*}}}*/
template <typename doubletype> void        yi_derivatives(doubletype* dydx, doubletype* y, int layer_index, int n, doubletype* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<doubletype>* vars){ /*{{{*/
	//computes yi derivatives at r=radius[layer_index]+ n/nstep*(radius[layer_index+1]-radius[layer_index])

	bool issolid=matlitho->issolid[layer_index];
	int iy,id,ny, nindex, nstep, nsteps;
	//femmodel->parameters->FindParam(&nstep,LoveIntStepsPerLayerEnum);

	nstep=vars->nstep[layer_index];
	nsteps=0;
	for (int i=0;i<layer_index;i++) nsteps+=vars->nstep[i];

	/*{{{*/ /* For reference:
			   flm=(la+2.0*mu);
			   rlm=(3.0*la+2.0*mu)/(la+2.0*mu);
			   rm0=mu/mu0;
			   rg0=g/g0;
			   frh=ro*g*ra/mu0;
			   fgr=4.0*PI*GG*ro*ra/g0;
			   fn=(deg*(deg+1.0));

			   if(issolid==1){
			   ny = 6;

			   dy[0*6+0]= (-2.0*la/flm)/x;
			   dy[0*6+1]= mu0/flm;
			   dy[0*6+2]= (la*fn/flm)/x;
			   dy[0*6+3]= 0.0;
			   dy[0*6+4]= 0.0;
			   dy[0*6+5]= 0.0;

			   dy[1*6+0]=  4.0*(-frh+rm0*rlm/x)/x + ro*pow(omega,2.0)*ra/mu0;
			   dy[1*6+1]=(-4.0*mu/flm)/x;
			   dy[1*6+2]= fn*(frh-2.0*rm0*rlm/x)/x;
			   dy[1*6+3]= fn/x;
			   dy[1*6+4]= 0.0;
			   dy[1*6+5]= -frh/rg0;

			   dy[2*6+0]= -1.0/x;
			   dy[2*6+1]= 0.0;
			   dy[2*6+2]= 1.0/x;
			   dy[2*6+3]= 1/rm0;
			   dy[2*6+4]= 0.0;
			   dy[2*6+5]= 0.0;

			   dy[3*6+0]= (frh-2.0*rm0*rlm/x)/x;
			   dy[3*6+1]= ( -la/flm)/x;
			   dy[3*6+2]= (2.0*rm0*(la*(2.0*fn-1.0)+2.0*mu*(fn-1.0))/flm)/(x*x) + ro*pow(omega,2.0)*ra/mu0;
			   dy[3*6+3]= -3.0/x;
			   dy[3*6+4]= -(frh/rg0)/x;
			   dy[3*6+5]= 0.0;

			   dy[4*6+0]= fgr;
			   dy[4*6+1]= 0.0;
			   dy[4*6+2]= 0.0;
			   dy[4*6+3]= 0.0;
			   dy[4*6+4]= 0.0;
			   dy[4*6+5]= 1.0;

			   dy[5*6+0]= 0.0;
			   dy[5*6+1]= 0.0;
			   dy[5*6+2]= -(fgr*fn)/x;
			   dy[5*6+3]= 0.0;
			   dy[5*6+4]= fn/(x*x);
			   dy[5*6+5]= -2.0/x;

			   } else {
			   ny = 2;

			   dy[0*6+0]= fgr/rg0;
			   dy[0*6+1]= 1.0;
			   dy[1*6+0]= (-4.0*(fgr/rg0)+fn/x)/x;
			   dy[1*6+1]= -2.0/x-fgr/rg0;

			   }
	*/ /*}}}*/
	nindex=nsteps*36+n*36;

	if(issolid==1){
		ny = 6;
	} else {
		ny = 2;
	}

	for (id=0;id<ny;id++){
		dydx[id]=0.0;
		for (iy=0;iy<ny;iy++){
			dydx[id]+=yi_prefactor[nindex+id*6+iy]*y[iy];
		}
	}
	return;
}/*}}}*/
template <typename doubletype> void        propagate_yi_euler(doubletype* y, doubletype xmin, doubletype xmax, int layer_index, doubletype* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<doubletype>* vars){ /*{{{*/
	//computes this: if we have y[j]=1.0 and y[!j]=0.0 at the bottom of the layer, what is y at the top of the layer?
	//euler method
	int nstep;
	nstep=vars->nstep[layer_index];

	doubletype* dydx=xNewZeroInit<doubletype>(6);
	doubletype dr = (xmax -xmin)/reCast<doubletype>(nstep);
	doubletype x=xmin;
	for(int i = 0;i<nstep;i++){
		yi_derivatives<doubletype>(dydx,y,layer_index, i,yi_prefactor,femmodel,matlitho,vars);
		for (int j=0;j<6;j++){
			y[j]+=dydx[j]*dr;
		}
		x = x + dr;
	}
	xDelete<doubletype>(dydx);
}/*}}}*/
template <typename doubletype> void        propagate_yi_RK2(doubletype* y, doubletype xmin, doubletype xmax, int layer_index, doubletype* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<doubletype>* vars){ /*{{{*/
	//computes this: if we have y[j]=1.0 and y[!j]=0.0 at the bottom of the layer, what is y at the top of the layer?
	//Implements Runge-Kutta 2nd order (midpoint method)
	int nstep;
	nstep=vars->nstep[layer_index];

	doubletype k1[6]={0};
	doubletype k2[6]={0};
	doubletype k3[6]={0};
	doubletype k4[6]={0};
	doubletype y1[6]={0};
	doubletype y2[6]={0};
	doubletype y3[6]={0};

	doubletype dr = (xmax -xmin)/reCast<doubletype>(nstep);
	doubletype x=xmin;

	for(int i = 0;i<nstep/2;i++){
		yi_derivatives<doubletype>(k1,y,layer_index, 2*i,yi_prefactor,femmodel,matlitho,vars);
		for (int j=0;j<6;j++) {y1[j]=y[j]+k1[j]*dr;}
		yi_derivatives<doubletype>(k2,y1,layer_index, 2*i+1,yi_prefactor,femmodel,matlitho,vars);		

		for (int j=0;j<6;j++){
			y[j]+=k2[j]*2.0*dr;
		}
		x = x + 2.0*dr;
	}
}/*}}}*/
	template <typename doubletype> void        propagate_yi_RK4(doubletype* y, doubletype xmin, doubletype xmax, int layer_index, doubletype* yi_prefactor, FemModel* femmodel, Matlitho* matlitho,LoveVariables<doubletype>* vars){ /*{{{*/
	//computes this: if we have y[j]=1.0 and y[!j]=0.0 at the bottom of the layer, what is y at the top of the layer?
	//Implements Runge-Kutta 4th order
	int nstep;
	nstep=vars->nstep[layer_index];

	doubletype k1[6]={0};
	doubletype k2[6]={0};
	doubletype k3[6]={0};
	doubletype k4[6]={0};
	doubletype y1[6]={0};
	doubletype y2[6]={0};
	doubletype y3[6]={0};

	doubletype dr = (xmax -xmin)/reCast<doubletype>(nstep);
	doubletype x=xmin;

	for(int i = 0;i<nstep/2-1;i++){
		yi_derivatives<doubletype>(k1,y,layer_index, 2*i,yi_prefactor,femmodel,matlitho,vars);
		for (int j=0;j<6;j++) {y1[j]=y[j]+k1[j]*dr;}
		yi_derivatives<doubletype>(k2,y1,layer_index, 2*i+1,yi_prefactor,femmodel,matlitho,vars);
		for (int j=0;j<6;j++) {y2[j]=y[j]+k2[j]*dr;}
		yi_derivatives<doubletype>(k3,y2,layer_index, 2*i+1,yi_prefactor,femmodel,matlitho,vars);
		for (int j=0;j<6;j++) {y3[j]=y[j]+k3[j]*2.0*dr;}
		yi_derivatives<doubletype>(k4,y3,layer_index, 2*i+2,yi_prefactor,femmodel,matlitho,vars);		

		for (int j=0;j<6;j++){
			y[j]+=(k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j])/3.0*dr;		
		}
	}

	//Last step: we don't know the derivative at xmax, so we have to use the midpoint method for the last step
	int i=nstep/2-1;
	yi_derivatives<doubletype>(k1,y,layer_index, 2*i,yi_prefactor,femmodel,matlitho,vars);
	for (int j=0;j<6;j++) {y1[j]=y[j]+k1[j]*dr;}
	yi_derivatives<doubletype>(k2,y1,layer_index, 2*i+1,yi_prefactor,femmodel,matlitho,vars);	
				

	for (int j=0;j<6;j++){
		y[j]+=k2[j]*2.0*dr;		
	}

	x = x + 2.0*dr;

}/*}}}*/
template <typename doubletype> void        Innersphere_boundaryconditions(doubletype* yi, int layer_index, int degree, doubletype omega, FemModel* femmodel, Matlitho* matlitho, LoveVariables<doubletype>* vars){ /*{{{*/
	//fills the boundary conditions at the bottom of layer[layer_index] in yi[0:2][0:5]

	int nyi;
	doubletype r = matlitho->radius[layer_index];
	doubletype ra=matlitho->radius[matlitho->numlayers];
	doubletype  g0,r0,mu0, GG;
	IssmDouble mu0p, GGp;
	doubletype deg = reCast<doubletype>(degree);

	femmodel->parameters->FindParam(&mu0p,LoveMu0Enum);
	femmodel->parameters->FindParam(&GGp,LoveGravitationalConstantEnum);

	g0=vars->g0;
	r0=vars->r0;
	mu0=mu0p;
	GG=GGp;
	nyi=vars->nyi;


	doubletype g=GetGravity<doubletype>(r,layer_index,femmodel,matlitho,vars);
	doubletype la,mu,ro;
	
	int i=layer_index-1;
	if (layer_index==0) i=layer_index;

	ro=matlitho->density[i];

	//elastic values
	la=matlitho->lame_lambda[i];
	mu=matlitho->lame_mu[i];
	doubletype Kappa=(la+2.0/3.0*mu);

	//update to viscoelastic values
	mu=vars->mu[i*vars->nfreq+vars->ifreq];
	la = Kappa-2.0/3.0*mu; 

	doubletype cst = 4.0*PI*GG*ro;
	doubletype r2=pow(r,2.0);

	//Greff-Lefftz and Legros 1997, p701, analytical solution for incompressible elastic layer for y3, y4, y5 ensuring they =0 at r=0
	//These equations are then divided by r^n for numerical stability at higher degrees

	//all terms in y1 y2 y6 are 0 in that layer because they are of the type r^l with l<0 and would diverge at the origin, that's why we only need 3 equations

	yi[0+nyi*0]=1.0*r/ra;
	yi[0+nyi*1]=1.0/(r*ra);
	yi[0+nyi*2]=0.0;

	yi[1+nyi*0]=(2.0*mu*(deg-1.0-3.0/deg) + cst/3.0*ro*r2)/mu0;
	yi[1+nyi*1]=(2.0*mu*(deg-1.0)/r2 + cst/3.0*ro)/mu0;
	yi[1+nyi*2]=-ro/mu0;

	yi[2+nyi*0]=(deg+3.0)/(deg*(deg+1.0))*r/ra;
	yi[2+nyi*1]=1.0/(deg*r*ra);
	yi[2+nyi*2]=0.0;

	yi[3+nyi*0]=2.0*mu*(deg+2.0)/((deg+1.0)*mu0);
	yi[3+nyi*1]=2.0*mu*(deg-1.0)/(deg*r2*mu0);
	yi[3+nyi*2]=0.0;

	yi[4+nyi*0]=0.0;
	yi[4+nyi*1]=0.0;
	yi[4+nyi*2]=1.0/(g0*ra);

	yi[5+nyi*0]=-cst*r/g0;
	yi[5+nyi*1]=-cst/(r*g0);
	yi[5+nyi*2]=deg/(r*g0);


	/*doubletype vp2 = (la + 2.0*mu)/ro;
	yi[0+nyi*0]=1.0*r/ra;
	yi[0+nyi*1]=1.0/(r*ra);
	yi[0+nyi*2]=0.0;

	yi[1+nyi*0]=(2.0*mu*(deg-1.0-3.0/deg) + cst/3.0*ro*r2)/mu0;
	yi[1+nyi*1]=(2.0*mu*(deg-1.0)/r2 + cst/3.0*ro)/mu0;
	yi[1+nyi*2]=-ro/mu0;

	yi[2+nyi*0]=(deg+3.0)/(deg*(deg+1.0))*r/ra;
	yi[2+nyi*1]=1.0/(deg*r*ra);
	yi[2+nyi*2]=0.0;

	yi[3+nyi*0]=2.0*mu*(deg+2.0)/((deg+1.0)*mu0);
	yi[3+nyi*1]=2.0*mu*(deg-1.0)/(deg*r2*mu0);
	yi[3+nyi*2]=0.0;

	yi[4+nyi*0]=0.0;
	yi[4+nyi*1]=0.0;
	yi[4+nyi*2]=1.0/(g0*ra);

	yi[5+nyi*0]=-cst*r/g0;
	yi[5+nyi*1]=-cst/(r*g0);
	yi[5+nyi*2]=deg/(r*g0);*/



}/*}}}*/
template <typename doubletype> void        build_yi_system(doubletype* yi, int deg, doubletype omega, doubletype* yi_prefactor, FemModel* femmodel, Matlitho* matlitho,LoveVariables<doubletype>* vars) { /*{{{*/

	doubletype	g0,r0,mu0,x,ro1, GG;
	int		nyi,starting_layer, nstep;
	doubletype 	xmin,xmax,one,ro,g, ra;
	IssmDouble 	mu0p, GGp;
	bool 		debug;
	int ny,is,ii,jj;
	int scheme;
	int ici = 0;   // Index of current interface 
	int cmb=0;

	femmodel->parameters->FindParam(&cmb,LoveCoreMantleBoundaryEnum);
	femmodel->parameters->FindParam(&mu0p,LoveMu0Enum);
	femmodel->parameters->FindParam(&GGp,LoveGravitationalConstantEnum);
	femmodel->parameters->FindParam(&debug,LoveDebugEnum);
	femmodel->parameters->FindParam(&scheme,LoveIntegrationSchemeEnum);

	g0=vars->g0;
	r0=vars->r0;
	nyi=vars->nyi;
	starting_layer=vars->starting_layer;
	mu0=mu0p;
	GG=GGp;
	ra=matlitho->radius[matlitho->numlayers];

	for (int i=0;i<6*(matlitho->numlayers+1);i++){
		for(int j=0;j<6*(matlitho->numlayers+1);j++){
			yi[i+6*(matlitho->numlayers+1)*j]=0.0;
		}
	}

	doubletype ystart[6];
	for (int k=0;k<6;k++) ystart[k]=0.0;		



	for (int i = starting_layer; i<matlitho->numlayers;i++){ 
		ici=i-starting_layer;
		xmin=matlitho->radius[i]/ra;
		xmax=(matlitho->radius[i+1])/ra;

		if (matlitho->issolid[i]){
			ny = 6;
			is = 0;
			one= 1.0;
		} else {	
			ny = 2;
			is = 4;
			one= -1.0;
		}

		for (int j = 0;j<ny;j++){
			for (int k=0;k<6;k++){ystart[k]=0.0;}
			ystart[j]= 1.0;

			// Numerical Integration 
			if (debug) propagate_yi_euler<doubletype>(&ystart[0], xmin, xmax, i, yi_prefactor,femmodel, matlitho, vars);
			else {
				if (scheme==0) propagate_yi_euler<doubletype>(&ystart[0], xmin, xmax, i, yi_prefactor,femmodel, matlitho, vars);
				else if (scheme==1) propagate_yi_RK2<doubletype>(&ystart[0], xmin, xmax, i, yi_prefactor,femmodel, matlitho, vars);
				else if (scheme==2) propagate_yi_RK4<doubletype>(&ystart[0], xmin, xmax, i, yi_prefactor,femmodel, matlitho, vars);
				else _error_("Love core error: integration scheme not found");
			}
			// Boundary Condition matrix - propagation part 
			ii = 6*(ici+1)+is;
			jj = 6*(ici+1)+j+is-3;
			for (int kk=0;kk<ny;kk++){
				yi[(ii+kk)+nyi*jj] = ystart[kk]*one;
			}
		}

		// Boundary Condition matrix - solid regions
		if (matlitho->issolid[i]){
			one = -1.0;
			if (i>0) if (!matlitho->issolid[i-1]) one = 1.0;
			for (int j=0;j<6;j++){
				yi[(j+6*ici)+ nyi*(j+6*ici+3)] = one;
			}
		} else { // Boundary Condition matrix - liquid regions
			ro1=matlitho->density[i];
			g=GetGravity<doubletype>(matlitho->radius[i], i, femmodel,matlitho,vars);
			ii = 6*ici;
			jj = 6*ici+3;
			yi[ii+nyi*(jj)] = -1.0;
			yi[ii+nyi*(jj+4)] = -g0/g;
			yi[(ii+1)+nyi*(jj)]=-ro1*g*ra/mu0;
			yi[(ii+2)+nyi*(jj+1)]=-1.0;
			yi[(ii+5)+nyi*(jj)]= 4.0*PI*GG*ro1*ra/g0;
			yi[(ii+4)+nyi*(jj+4)]=-1.0;
			yi[(ii+5)+nyi*(jj+5)]=-1.0;
			g=GetGravity<doubletype>(matlitho->radius[i+1], i,femmodel,matlitho,vars);
			ii = 6*(ici+1);

			yi[ii+nyi*(jj+2)]=-1.0;
			yi[ii+nyi*(jj+4)]=yi[(ii+4)+nyi*(jj+4)]*g0/g; // yi(17,14) solution integration 1 of z5 CMB
			yi[ii+nyi*(jj+5)]=yi[(ii+4)+nyi*(jj+5)]*g0/g; // yi(17,15) solution integration 2 of z5 CMB
			// yi(13,..) y1 CMB
			yi[(ii+1)+nyi*(jj+2)]=-ro1*g*ra/mu0;
			yi[(ii+2)+nyi*(jj+3)]=-1.0;
			yi[(ii+5)+nyi*(jj+2)]= 4.0*PI*GG*ro1*ra/g0;
		}	
		ici = ici+1;
	}

	//-- Internal sphere: integration starts here rather than r=0 for numerical reasons

	Innersphere_boundaryconditions<doubletype>(yi, starting_layer, deg, omega, femmodel, matlitho,vars);

	//-- Surface conditions
	yi[(nyi-6)+nyi*(nyi-3)]=-1.0;
	yi[(nyi-4)+nyi*(nyi-2)]=-1.0;
	yi[(nyi-2)+nyi*(nyi-1)]=-1.0;
	yi[(nyi-1)+nyi*(nyi-1)]=deg+1.0;
//y4,y3,y2,y1,y6,y5
//y4,y6,y2,y1,y3,y5
	//-- Degree 1 special case
	if(deg==1){
		for (int i=0;i<nyi;i++){
			yi[(nyi-1)+nyi*i]=0.0;
		}
		yi[(nyi-1)+nyi*(nyi-1)]=1.0;
	}

}/*}}}*/
template <typename doubletype> void        yi_boundary_conditions(doubletype* yi_righthandside, int degree, FemModel* femmodel, Matlitho* matlitho,LoveVariables<doubletype>* vars, int forcing_type){ /*{{{*/

	doubletype  g0,r0,mu0,ra,rb,rc;
	int nyi,icb,cmb,starting_layer;
	doubletype* EarthMass;
	IssmDouble mu0p;
	doubletype deg = reCast<doubletype>(degree);

	g0=vars->g0;
	r0=vars->r0;
	nyi=vars->nyi;
	starting_layer=vars->starting_layer;
	EarthMass=vars->EarthMass;

	femmodel->parameters->FindParam(&mu0p,LoveMu0Enum);
	femmodel->parameters->FindParam(&icb,LoveInnerCoreBoundaryEnum);
	femmodel->parameters->FindParam(&cmb,LoveCoreMantleBoundaryEnum);

	mu0=mu0p;
	// In Case of a Inner Core - Outer Core - Mantle planet and Boundary conditions on these 3 interfaces
	ra=matlitho->radius[matlitho->numlayers];	
	rb=0;
	rc=0;
	if (forcing_type<=4){
		rc=matlitho->radius[icb];
	} 
	else if (forcing_type<=8){
		rb=matlitho->radius[cmb];
	}

	doubletype ro_mean=EarthMass[matlitho->numlayers-1]/(4.0/3.0*PI*pow(ra,3.0));

	for (int i=0;i<(matlitho->numlayers+1)*6;i++) yi_righthandside[i]=0.0;

	switch (forcing_type) {

		//-- forcings at the Inner Core Boundary
		case 1:	//'ICB --Volumetric Potential'
			yi_righthandside[6*icb+5]=(deg)/(rc*g0);
			yi_righthandside[6*icb+4]=1.0/(ra*g0);
			break;
		case 2: //'ICB --Pressure'
			yi_righthandside[6*icb+1]=-ro_mean/mu0;
			break;
		case 3://'ICB --Loading'
			yi_righthandside[6*icb+1]=-ro_mean*(2.0*deg+1.0)/(3.0*mu0)*ra/rc;
			yi_righthandside[6*icb+5]= (2.0*deg+1.0)/(rc*g0);
			break;
		case 4://'ICB --Tangential Traction'
			yi_righthandside[6*icb+3]= ro_mean/mu0;
			break;

			//--forcings at the Core Mantle Boundary
		case 5://'CMB --Volumetric Potential'
			yi_righthandside[6*cmb+1]=-ro_mean/mu0*ra/rb;
			yi_righthandside[6*cmb+5]= (2.0*deg+1.0)/(rb*g0);
			break;
		case 6://'CMB --Pressure'
			yi_righthandside[6*cmb+1]=-ro_mean/mu0;
			break;
		case 7://'CMB --Loading'
			yi_righthandside[6*cmb+1]=-ro_mean*(2.0*deg+1.0)/(3.0*mu0)*ra/rb;
			yi_righthandside[6*cmb+5]= (2.0*deg+1.0)/(rb*g0);
			break;
		case 8://'CMB --Tangential Traction'
			yi_righthandside[6*cmb+3]=-ro_mean/mu0;
			break;

			//--forcings at the surface
		case 9://'SURF--Volumetric Potential'
			if (degree>1) yi_righthandside[nyi-1]=(2.0*deg+1.0)/(ra*g0);
			break;
		case 10://'SURF--Pressure'
			yi_righthandside[nyi-5]=-ro_mean/mu0;
			break;
		case 11://'SURF--Loading'
			yi_righthandside[nyi-5]=-ro_mean*(2.0*deg+1.0)/(3.0*mu0);
			if (degree>1) yi_righthandside[nyi-1]= (2.0*deg+1.0)/(ra*g0);
			break;
		case 12://'SURF--Tangential Traction'
			yi_righthandside[nyi-3]= ro_mean/mu0;
			break;
		default:
			_error_("love core error: forcing_type not supported yet");
	}
}/*}}}*/
template <typename doubletype> void        solve_yi_system(doubletype* loveh, doubletype* lovel, doubletype* lovek, int deg, doubletype omega, IssmDouble* frequencies, doubletype* yi, doubletype* rhs, FemModel* femmodel, Matlitho* matlitho, LoveVariables<doubletype>* vars, bool verbosecpu){ /*{{{*/

	doubletype  g0,r0,mu0;
	//IssmDouble* frequencies;
	int nyi,starting_layer, dummy,cmb;
	bool allow_layer_deletion, debug;
	doubletype* EarthMass=NULL;
	IssmDouble mu0p,loveratio,underflow_tol;

	g0=vars->g0;
	r0=vars->r0;
	nyi=vars->nyi;
	starting_layer=vars->starting_layer;
	EarthMass=vars->EarthMass;

	femmodel->parameters->FindParam(&mu0p,LoveMu0Enum);
	femmodel->parameters->FindParam(&allow_layer_deletion,LoveAllowLayerDeletionEnum);
	femmodel->parameters->FindParam(&underflow_tol,LoveUnderflowTolEnum);
	femmodel->parameters->FindParam(&debug,LoveDebugEnum);
	femmodel->parameters->FindParam(&cmb,LoveCoreMantleBoundaryEnum);
	//femmodel->parameters->FindParam(&frequencies,&dummy,LoveFrequenciesEnum);
	mu0=mu0p;
	doubletype ra=matlitho->radius[matlitho->numlayers];
	bool exit=false;
	int lda,ldb;


	for(;!exit;){ //cycles of: attempt to solve the yi system, then delete a layer if necessary
		lda=nyi;
		ldb=nyi;
		doubletype*  yilocal=xNew<doubletype>(nyi*nyi); // we will need to redeclare these inside here as nyi changes
		doubletype*  rhslocal=xNew<doubletype>(nyi);

		//we need to do a local copy of yi,rhs to send them to LAPACK with the appropriate size and to keep the original matrices in case layers are deleted and we need to try again
		for (int i=0;i<nyi;i++){ 
			rhslocal[i]=rhs[i];
			for (int j=0;j<nyi;j++){
				yilocal[i+j*nyi]=yi[i+j*nyi];
			}
		}
		
		if (debug){
			IssmDouble*  yidebug=xNew<IssmDouble>(nyi*nyi);
			IssmDouble*  rhsdebug=xNew<IssmDouble>(nyi);
			for (int i=0;i<nyi;i++){ 
				rhsdebug[i]=DownCastVarToDouble(rhs[i]);
				for (int j=0;j<nyi;j++){
					yidebug[i+j*nyi]=DownCastVarToDouble(yi[i+j*nyi]);
				}
			}
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveYiEnum,yidebug,nyi,nyi,0,0));
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveRhsEnum,rhsdebug,nyi,1,0,0));
			xDelete<IssmDouble>(yidebug);
			xDelete<IssmDouble>(rhsdebug);
		}

		//-- Resolution
		int* ipiv=xNewZeroInit<int>(nyi); //pivot index vector
		int  info = 0;// error checker
		int  nrhs=1; // number of right hand size columns

		allgesv<doubletype>(&nyi, &nrhs, yilocal, &lda, ipiv, rhslocal, &ldb, &info);

		xDelete<int>(ipiv);

			/*_printf_("i j yi[i+nyi*j] rhs[i]");
			for (int i=0;i<nyi;i++){
					_printf_(i<<" "<<rhs[i]<<"\n");
			}

			for (int i=0;i<nyi;i++){
				for (int j=0;j<nyi;j++){
					_printf_(i<<" "<<j<<" "<<yi[i+nyi*j]<<" "<<rhs[i]<<"\n");
				}
			}
			_error_("love core warning in DGESV : LAPACK linear equation solver couldn't resolve the system");*/

		if(VerboseSolution() && verbosecpu && info!=0){ 
			_printf_("i j yi[i+nyi*j] rhs[i]\n");
			for (int i=0;i<nyi;i++){
				for (int j=0;j<nyi;j++){
					_printf_(i<<" "<<j<<" "<<yi[i+nyi*j]<<" "<<rhs[i]<<"\n");
				}
			}
			_error_("deg " << deg << " love core error in DGESV : LAPACK linear equation solver couldn't resolve the system");
		}


		*loveh = rhslocal[nyi-3]*ra*g0;
		*lovel = rhslocal[nyi-2]*ra*g0;
		*lovek = rhslocal[nyi-1]*ra*g0;

		doubletype loveh1 = rhslocal[3];
		doubletype lovel1 = rhslocal[5];
		doubletype lovek1 = rhslocal[7] - pow(matlitho->radius[starting_layer]/ra,deg)/(g0*ra);

		doubletype loveh1s = rhslocal[nyi-3];
		doubletype lovel1s = rhslocal[nyi-2];
		doubletype lovek1s = rhslocal[nyi-1] - 1.0/(g0*ra);
		doubletype zero    = 0;

		loveratio = abs(loveh1/loveh1s); //ratio of center to surface love numbers, determines if we should remove layers
		if (abs(lovel1/lovel1s) < loveratio) loveratio = abs(lovel1/lovel1s); 
		if (abs(lovek1/lovek1s) < loveratio) loveratio = abs(lovek1/lovek1s);

		if (debug) goto save_results;

		if (!allow_layer_deletion || nyi<=12 || omega!=angular_frequency<doubletype>(frequencies[0]) || deg==0){ 
			goto save_results;
			/*We are not allowed to delete layers, or there is only one layer left. We also don't want to delete 
			  layers in the middle of a loop on frequencies, as that can lead to a jump that would compromise the 
			  inverse laplace transform.*/
		}

		if (omega==zero){ // if running elastic love_numbers, record at which degree we must delete layers, this way we synch layer deletion between cpus next time we calculate love numbers
			//We need to delete a layer and try again if the ratio between deepest love number to surface love number is too low (risk of underflow) or garbage
			if (loveratio<=underflow_tol || xIsNan(loveratio) || xIsInf(loveratio)) {
				vars->deg_layer_delete[starting_layer]=deg;
				if(VerboseSolution() && verbosecpu){
					_printf_("\n   Degree: " << deg <<", surface/Depth Love number ratio small: "<<loveratio<<"\n");
					_printf_("    Changing the interface where integration starts\n");
					_printf_("    New start interface: r="<< matlitho->radius[starting_layer+1] <<"m\n\n");
				}
			}
		}

		if (deg==vars->deg_layer_delete[starting_layer]){ // if we are at the degree where we should delete the current layer, proceed to delete the bottom layer
			//if (omega!=0 && VerboseSolution()  && verbosecpu) _printf_(", deleting layer " << starting_layer << "\n");
			nyi-=6;
			starting_layer+=1;
			vars->nyi=nyi;
			vars->starting_layer=starting_layer;

			for (int i=0;i<nyi;i++){//shift everything down by 1 layer
				rhs[i]=rhs[i+6];
				for (int j=0;j<nyi;j++){
					yi[j+i*nyi]=yi[j+6+(i+6)*(nyi+6)];
				}
			}

	Innersphere_boundaryconditions<doubletype>(yi, starting_layer, deg, omega, femmodel, matlitho,vars);
		} else { //we are ready to save the outputs and break the main loop

save_results:
			for (int i=0;i<nyi;i++){
				rhs[i]=rhslocal[i];
				for (int j=0;j<nyi;j++){
					yi[j+i*nyi]=yilocal[j+i*nyi];
				}
			}

			//make sure we can't output numbers from deleted layers
			for (int i=nyi;i<(matlitho->numlayers+1)*6;i++){ 
				rhs[i]=0.0;
				for (int j=0;j<(matlitho->numlayers+1)*6;j++){
					yi[j+i*(matlitho->numlayers+1)*6]=0.0;
				}
			}
			for (int i=0;i<nyi;i++){
				for (int j=nyi;j<(matlitho->numlayers+1)*6;j++){
					yi[j+i*(matlitho->numlayers+1)*6]=0.0;
				}
			}

			exit = true;
		}
		xDelete<doubletype>(yilocal);
		xDelete<doubletype>(rhslocal);
	}
	//xDelete<IssmDouble>(frequencies);	

}/*}}}*/
template <typename doubletype> void        love_freq_to_temporal(LoveNumbers<doubletype>* Lovet, LoveNumbers<doubletype>* Tidalt, doubletype* pmtf_colineart, doubletype* pmtf_orthot, LoveNumbers<doubletype>* Lovef,LoveNumbers<doubletype>* Tidalf, IssmDouble* frequencies, FemModel* femmodel, bool verbosecpu){ /*{{{*/
	//Transforms all frequency-dependent love numbers into time-dependent love numbers
	int nfreq,sh_nmax,sh_nmin,indxi,indf, NTit, forcing_type, nt;
	IssmDouble kf,Omega,moi_e,moi_p,alpha;
	doubletype* pmtf_colinearf=NULL;
	doubletype* pmtf_orthof=NULL;
	bool chandler_wobble=false;

	nfreq=Lovef->nfreq;
	sh_nmin=Lovef->sh_nmin;
	sh_nmax=Lovef->sh_nmax;

	femmodel->parameters->FindParam(&NTit,LoveNTemporalIterationsEnum);

	//Parameters for the rotationnal feedback
	femmodel->parameters->FindParam(&chandler_wobble,LoveChandlerWobbleEnum);
	femmodel->parameters->FindParam(&forcing_type,LoveForcingTypeEnum);
	femmodel->parameters->FindParam(&kf,TidalLoveK2SecularEnum);
	femmodel->parameters->FindParam(&Omega,RotationalAngularVelocityEnum);
	femmodel->parameters->FindParam(&moi_e,RotationalEquatorialMoiEnum);
	femmodel->parameters->FindParam(&moi_p,RotationalPolarMoiEnum);

	nt=Lovet->nfreq;
	doubletype* xi=postwidder_coef<doubletype>(NTit);

	if(VerboseSolution()  && verbosecpu) _printf_("   Inverse Laplace Transform... ");

	for (int d=sh_nmin;d<sh_nmax+1;d++){
		for (int t=0;t<nt;t++){
			postwidder_transform<doubletype>(Lovet->H,Lovef->H,d,t,sh_nmax,NTit,xi,femmodel);
			postwidder_transform<doubletype>(Lovet->K,Lovef->K,d,t,sh_nmax,NTit,xi,femmodel);
			postwidder_transform<doubletype>(Lovet->L,Lovef->L,d,t,sh_nmax,NTit,xi,femmodel);
		}
	}

	if(VerboseSolution()  && verbosecpu) _printf_("done!\n");

	if (forcing_type==11){ //Let's retrieve the functions necessary for the rotational_feedback
		if(VerboseSolution()  && verbosecpu) _printf_("     Transforming PMTF and tidal love numbers... ");
		pmtf_colinearf = xNewZeroInit<doubletype>(3*nfreq);
		pmtf_orthof = xNewZeroInit<doubletype>(3*nfreq);
		int d=2;
		doubletype s,R1,R2;
		alpha=(moi_p-moi_e)/moi_e; //Earth flattening

		if (chandler_wobble){ //Chandler Wobble is untested yet
			for (int fr=0;fr<nfreq;fr++){		
				s=angular_frequency<doubletype>(frequencies[fr]);
				R1=alpha*Omega*(1.0-Tidalf->K[fr*3+d]/kf);
				R2=1.0+alpha*Tidalf->K[fr*3+d]/kf;
				pmtf_colinearf[fr*3+d]=alpha*(1.0+Lovef->K[fr*(sh_nmax+1)+d])*(Omega*R1-pow(s,2)*R2)/(pow(R1,2.0)+pow(s*R2,2.0));
				pmtf_orthof[fr*3+d]=alpha*(1.0+Lovef->K[fr*(sh_nmax+1)+d])*s*Omega*(1.0+alpha)/(pow(R1,2.0)+pow(s*R2,2.0));
			}
		}
		else {
			for (int fr=0;fr<nfreq;fr++){		
				pmtf_colinearf[fr*3+d]=(1.0+Lovef->K[fr*(sh_nmax+1)+d])/(1.0-Tidalf->K[fr*3+d]/kf);
				pmtf_orthof[fr*3+d]=0.0;
			}
		}
		for (int t=0;t<nt;t++){
			postwidder_transform<doubletype>(Tidalt->H,Tidalf->H,2,t,2,NTit,xi,femmodel);
			postwidder_transform<doubletype>(Tidalt->K,Tidalf->K,2,t,2,NTit,xi,femmodel);
			postwidder_transform<doubletype>(Tidalt->L,Tidalf->L,2,t,2,NTit,xi,femmodel);
			postwidder_transform<doubletype>(pmtf_colineart,pmtf_colinearf,2,t,2,NTit,xi,femmodel);
			postwidder_transform<doubletype>(pmtf_orthot,pmtf_orthof,2,t,2,NTit,xi,femmodel);
		}
		xDelete<doubletype>(pmtf_colinearf);
		xDelete<doubletype>(pmtf_orthof);
		if(VerboseSolution() && verbosecpu) _printf_("done!\n");
	}

	xDelete<doubletype>(xi);
}/*}}}*/

template <typename doubletype> void        compute_love_numbers(LoveNumbers<doubletype>* Lovef, LoveNumbers<doubletype>* Elastic, int forcing_type, int sh_cutoff, IssmDouble* frequencies, FemModel* femmodel, Matlitho* matlitho, LoveVariables<doubletype>* vars, bool verbosecpu){

	int nsteps, kernel_index,kernel_indexe,deleted_layer_offset, deg, sh_nmin, sh_nmax, nfreq;
	doubletype  lovek, loveh, lovel, loveratio;
	doubletype  omega;
	doubletype* yi_prefactor=NULL;
	doubletype* yi_righthandside=NULL;
	doubletype* yi=NULL;
	doubletype  underflow_tol;
	IssmDouble dr;
	bool freq_skip, istemporal;
	int cmb=0;
	int nyi_init=0;

	//femmodel->parameters->FindParam(&nstep,LoveIntStepsPerLayerEnum);
	femmodel->parameters->FindParam(&istemporal,LoveIsTemporalEnum);
	femmodel->parameters->FindParam(&cmb,LoveCoreMantleBoundaryEnum);

	nfreq=Lovef->nfreq;
	sh_nmin=Lovef->sh_nmin;
	sh_nmax=Lovef->sh_nmax;
	if (Elastic==NULL) sh_cutoff=sh_nmax;

	// reset deleted layers in case we have called this function before;
	vars->starting_layer=0;
	vars->nyi=6*(matlitho->numlayers-vars->starting_layer+1);
	nyi_init=6*(matlitho->numlayers+1);
	nsteps=0;
	for (int i=0;i<matlitho->numlayers;i++)	nsteps+=vars->nstep[i];

	//yi_prefactor=xNewZeroInit<doubletype>(6*6*nstep*matlitho->numlayers);
	yi_prefactor=xNewZeroInit<doubletype>(6*6*nsteps);
	yi_righthandside=xNewZeroInit<doubletype>(nyi_init);
	yi=xNewZeroInit<doubletype>(nyi_init*nyi_init);

	//precompute yi coefficients that do not depend on degree or frequency
	fill_yi_prefactor<doubletype>(yi_prefactor, NULL, NULL,femmodel, matlitho,vars);

	if (VerboseSolution() && Elastic  && verbosecpu) _printf_("\n");

	for(int deg=0;deg<2;deg++){ // calculation is in the center of mass reference frame, neutralize degree 0 and 1 mass changes, i.e 1+k=0
		for (int fr=0;fr<nfreq;fr++){
			Lovef->K[fr*(sh_nmax+1)+deg]=-1.0;
		}
	}

	for(int deg=sh_nmin;deg<sh_cutoff+1;deg++){
		if (VerboseSolution() && Elastic && verbosecpu) {
			_printf_("\r   Degree: " << deg << "/" << sh_nmax << "    ");
		}

		//precompute yi coefficients that depend on degree but not frequency
		fill_yi_prefactor<doubletype>(yi_prefactor, &deg, NULL,femmodel, matlitho,vars); 

		for (int fr=0;fr<nfreq;fr++){
			omega=angular_frequency<doubletype>(frequencies[fr]);
			vars->ifreq=fr;
			//precompute yi coefficients that depend on degree and frequency
			fill_yi_prefactor<doubletype>(yi_prefactor, &deg,&omega,femmodel, matlitho,vars);

			//solve the system
			yi_boundary_conditions<doubletype>(yi_righthandside,deg,femmodel,matlitho,vars,forcing_type);
			build_yi_system<doubletype>(yi,deg,omega,yi_prefactor,femmodel,matlitho,vars);
			solve_yi_system<doubletype>(&loveh,&lovel,&lovek, deg, omega, frequencies, yi, yi_righthandside,femmodel, matlitho,vars,verbosecpu && !Elastic);

			Lovef->H[fr*(sh_nmax+1)+deg]=loveh;
			Lovef->K[fr*(sh_nmax+1)+deg]=lovek-1.0;
			Lovef->L[fr*(sh_nmax+1)+deg]=lovel;
			deleted_layer_offset=(matlitho->numlayers+1)*6-vars->nyi;// =6 per deleted layer
			kernel_index=fr*(sh_nmax+1)*(matlitho->numlayers+1)*6 + deg*(matlitho->numlayers+1)*6 + deleted_layer_offset;
			for (int i=0;i<vars->nyi;i++){
				Lovef->Kernels[kernel_index+i]=yi_righthandside[i];
			}

		}
	}

	if (Elastic) { // if elastic values were provided, we copy elastic love numbers above the cutoff degree instead of computing them
		for(int deg=sh_cutoff+1;deg<sh_nmax+1;deg++){
			if (VerboseSolution() && Elastic  && verbosecpu) {
				if (deg==sh_nmax || deg%100==0)	_printf_("\r   Degree: " << deg << "/" << Lovef->sh_nmax << "    ");
			}
			for (int fr=0;fr<nfreq;fr++){
				// just copy the elastic values
				Lovef->H[fr*(sh_nmax+1)+deg]=Elastic->H[deg];
				Lovef->K[fr*(sh_nmax+1)+deg]=Elastic->K[deg];
				Lovef->L[fr*(sh_nmax+1)+deg]=Elastic->L[deg];
				deleted_layer_offset=(matlitho->numlayers+1)*6-vars->nyi;// =6 per deleted layer
				kernel_index=fr*(sh_nmax+1)*(matlitho->numlayers+1)*6 + deg*(matlitho->numlayers+1)*6 + deleted_layer_offset;
				kernel_indexe=deg*(matlitho->numlayers+1)*6 + deleted_layer_offset;
				for (int i=0;i<vars->nyi;i++){
					Lovef->Kernels[kernel_index+i]=Elastic->Kernels[kernel_indexe+i];
				}
			}
		}
	}

	if (VerboseSolution() && Elastic  && verbosecpu) _printf_("\n");
	xDelete<doubletype>(yi);
	xDelete<doubletype>(yi_righthandside);
	xDelete<doubletype>(yi_prefactor);
}/*}}}*/

/*templated cores:*/
template <typename doubletype> LoveVariables<doubletype>*	love_init(FemModel* femmodel, Matlitho* matlitho, bool verbosecpu){/*{{{*/

	/*initialize Planet_Mass(r) for efficient computation of gravity, value of surface gravity and inital size of the yi equation system*/

	bool        verbosemod = (int)VerboseModule();
	int         numlayers  = matlitho->numlayers;
	int 	    minsteps;
	doubletype* r=NULL;
	doubletype  r1,r2,ro, GG;
	IssmDouble GGp;
	IssmDouble dr;

	/*outputs:*/
	doubletype* EarthMass=NULL;
	doubletype  g0,r0;
	int         nyi,starting_layer,cmb;
	int*	    deg_layer_delete;
	int*	    nstep;

	
	femmodel->parameters->FindParam(&GGp,LoveGravitationalConstantEnum);
	femmodel->parameters->FindParam(&minsteps, LoveMinIntegrationStepsEnum);
	femmodel->parameters->FindParam(&dr, LoveMaxIntegrationdrEnum);
	femmodel->parameters->FindParam(&cmb,LoveCoreMantleBoundaryEnum);
	GG=GGp;
	EarthMass=xNewZeroInit<doubletype>(numlayers+1);
	deg_layer_delete=xNewZeroInit<int>(numlayers);

	r=xNewZeroInit<doubletype>(numlayers+1);
	nstep=xNewZeroInit<int>(numlayers);
	for (int i=0;i<numlayers+1;i++){
		r[i] = matlitho->radius[i];
		if (i<numlayers) {
			// nstep[i] is the largest even integer such that (radius[i+1]-radius[i])/nstep[i]<dr
			nstep[i]=ceil((matlitho->radius[i+1]-matlitho->radius[i])/dr/2)*2;
			if (nstep[i]<minsteps) nstep[i]=minsteps;
		}
	}

	for (int i=0;i<numlayers;i++){
		r2 = r[i+1];
		ro = matlitho->density[i];
		if (i==0){
			EarthMass[i] = ro*pow(r2,3.0)*4.0*PI/3.0;
		}else{
			r1=r[i];
			EarthMass[i] = EarthMass[i-1] + ro*(pow(r2,3.0)-pow(r1,3.0))*4.0*PI/3.0;;
		}
	}
	g0=EarthMass[numlayers-1]*GG/pow(r[numlayers],2.0);
	r0=r[numlayers];
	starting_layer=0;
	nyi=6*(numlayers-starting_layer+1);


	if(VerboseSolution() && verbosecpu){
		_printf_("     Surface gravity: " << g0 << " m.s^-2\n");
		_printf_("     Mean density: " << EarthMass[numlayers-1]/(4.0/3.0*PI*pow(r0,3.0)) << " kg.m^-3\n");
	}

	xDelete<doubletype>(r);
	return new LoveVariables<doubletype>(EarthMass,g0,r0,nyi,starting_layer,deg_layer_delete,nstep);

} /*}}}*/
template <typename doubletype> void        love_core_template(FemModel* femmodel){ /*{{{*/

	Matlitho*   matlitho=NULL;
	int         nfreq, NTit,nt, forcing_type,dummy, sh_cutoff;
	int         sh_nmin,sh_nmax,kernel_index,deleted_layer_offset;
	bool        allow_layer_deletion,love_kernels, istemporal, freq_skip;
	bool        verbosemod = (int)VerboseModule();
	IssmDouble *frequencies = NULL;
	IssmDouble *frequencies_local=NULL;
	bool        save_results;
	bool        complex_computation;
	bool	    quad_precision;
	bool	    verbosecpu=false;

	doubletype  omega;
	doubletype  lovek, loveh, lovel, loveratio;
	IssmDouble pw_threshold, pw_test_h, pw_test_l,pw_test_k;

	/* parallel computing */
	LoveNumbers<doubletype>* Lovef_local=NULL;
	LoveNumbers<doubletype>* Tidalf_local=NULL;

	/*elastic & fluid love numbers*/
	LoveNumbers<doubletype>* Elastic=NULL;
	IssmDouble* frequencies_elastic=NULL;
	LoveNumbers<doubletype>* Fluid=NULL;
	IssmDouble* frequencies_fluid=NULL;

	LoveVariables<doubletype>* vars=NULL;

	/*recover materials parameters: there is only one Matlitho, chase it down the hard way:*/
	for (Object* & object: femmodel->materials->objects){
		Material* material=xDynamicCast<Material*>(object);
		if(material->ObjectEnum()==MatlithoEnum){
			matlitho=xDynamicCast<Matlitho*>(material);
			break;
		}
	}
	_assert_(matlitho);

	/*recover parameters: */
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&nfreq,LoveNfreqEnum);
	femmodel->parameters->FindParam(&frequencies,&dummy,LoveFrequenciesEnum); _assert_(nfreq==dummy);
	femmodel->parameters->FindParam(&sh_nmax,LoveShNmaxEnum);
	femmodel->parameters->FindParam(&sh_nmin,LoveShNminEnum);
	femmodel->parameters->FindParam(&allow_layer_deletion,LoveAllowLayerDeletionEnum);
	femmodel->parameters->FindParam(&love_kernels,LoveKernelsEnum);
	femmodel->parameters->FindParam(&forcing_type,LoveForcingTypeEnum);
	femmodel->parameters->FindParam(&istemporal,LoveIsTemporalEnum);
	femmodel->parameters->FindParam(&complex_computation,LoveComplexComputationEnum);
	femmodel->parameters->FindParam(&quad_precision,LoveQuadPrecisionEnum);
	femmodel->parameters->FindParam(&pw_threshold,LovePostWidderThresholdEnum);
	if (istemporal)	femmodel->parameters->FindParam(&NTit,LoveNTemporalIterationsEnum);

	Elastic= new LoveNumbers<doubletype>(sh_nmin,sh_nmax,1,1,1,matlitho);
	Fluid= new LoveNumbers<doubletype>(sh_nmin,sh_nmax,1,1,1,matlitho);
	//distribute frequencies for parallel computation /*{{{*/
	int nt_local, nf_local, lower_row, upper_row;
	if (istemporal){ 
		//temporal love numbers are obtained via blocks of 2*NTit frequencies samples
		//here we are making sure no block is split between different cpus, which would make the inverse laplace transform impossible
		nt=nfreq/2/NTit;
		nt_local=DetermineLocalSize(nt,IssmComm::GetComm());
		nf_local=nt_local*2*NTit; // number of local frequencies
		GetOwnershipBoundariesFromRange(&lower_row,&upper_row,nt_local,IssmComm::GetComm());
		lower_row*=2*NTit;
		upper_row*=2*NTit;
	}
	else{
		//in this case frequency samples are completely independent so we can split them evenly across cpus
		nf_local=DetermineLocalSize(nfreq,IssmComm::GetComm());
		GetOwnershipBoundariesFromRange(&lower_row,&upper_row,nf_local,IssmComm::GetComm());
	}

	if (lower_row==0) verbosecpu=true; //let only cpu1 be verbose
	if(VerboseSolution() && verbosecpu) _printf0_("   computing LOVE numbers\n");
	vars=love_init<doubletype>(femmodel,matlitho,verbosecpu);

	frequencies_local=xNewZeroInit<IssmDouble>(nf_local);
	for (int fr=0;fr<nf_local;fr++){
		frequencies_local[fr]=frequencies[lower_row+fr];
	}
	Lovef_local= new LoveNumbers<doubletype>(sh_nmin,sh_nmax,nf_local,lower_row,nfreq, matlitho);
	Tidalf_local= new LoveNumbers<doubletype>(2,2,nf_local,lower_row,nfreq, matlitho);

	/*}}}*/
	frequencies_elastic=xNewZeroInit<IssmDouble>(1);
	frequencies_fluid=xNewZeroInit<IssmDouble>(1);
	for (int fr=0;fr<nfreq;fr++){ // find the lowest non-zero frequency requested
		if (frequencies_fluid[0]==0) frequencies_fluid[0]=frequencies[fr];
		else if(frequencies[fr]!=0 && frequencies_fluid[0]>frequencies[fr]) frequencies_fluid[0]=frequencies[fr];
	}

	// run elastic and fluid love numbers

	EarthRheology<doubletype>(vars,frequencies_elastic,1,matlitho,femmodel);

	if(VerboseSolution() && verbosecpu) _printf_("     elastic\n");

	compute_love_numbers<doubletype>(Elastic, NULL, forcing_type, sh_nmax,frequencies_elastic, femmodel, matlitho, vars,verbosecpu);

	if (nfreq>1){
		EarthRheology<doubletype>(vars,frequencies_fluid,1,matlitho,femmodel);

		compute_love_numbers<doubletype>(Fluid, NULL, forcing_type, sh_nmax,frequencies_fluid, femmodel, matlitho, vars,verbosecpu);
		sh_cutoff=sh_nmax;
		for (int deg=100;deg<sh_nmax+1;deg++){
			pw_test_h=abs((Fluid->H[deg]-Elastic->H[deg])/Elastic->H[deg]);
			pw_test_k=abs((Fluid->K[deg]-Elastic->K[deg])/Elastic->K[deg]);
			pw_test_l=abs((Fluid->L[deg]-Elastic->L[deg])/Elastic->L[deg]);
			if (pw_test_h<pw_threshold && pw_test_k<pw_threshold && pw_test_l<pw_threshold){
				sh_cutoff=deg;
				if(VerboseSolution() && verbosecpu){
					_printf_("   Degree: " << deg << "/" << sh_nmax << "    ");
					_printf_("      found negligible variation across frequencies, will copy elastic values after this degree\n");
					_printf_("      Delta_h/h=" << pw_test_h << "; Delta_k/k="<< pw_test_k << "; Delta_l/l=" << pw_test_l << "; threshold set to " << pw_threshold << "\n");
				}
				break;
			}
		}
	} 
	else sh_cutoff=sh_nmax; 

	delete Fluid;

	//Requested forcing_type
	if (nfreq>1){ // if we are not running just elastic love numbers
		if(VerboseSolution() && verbosecpu){
			if (forcing_type==11) _printf_("     loading\n");
			else if(forcing_type==9) _printf_("     tidal\n");
			else _printf_("     love\n");
		}
		EarthRheology<doubletype>(vars,frequencies_local,nf_local,matlitho,femmodel);
		compute_love_numbers<doubletype>(Lovef_local, Elastic, forcing_type, sh_cutoff, frequencies_local, femmodel, matlitho, vars,verbosecpu);
	}
	else{
		Lovef_local->Copy(Elastic);
	}
	/*}}}*/

	//Take care of rotationnal feedback love numbers, if relevant /*{{{*/
	if (forcing_type==11 && sh_nmin<=2 && sh_nmax>=2){ // if forcing is surface loading and we have degree 2
		if(VerboseSolution() && verbosecpu) _printf_("     tidal\n");
		int tidal_forcing_type=9;
		//no need to call EarthRheology, we already have the right one
		compute_love_numbers<doubletype>(Tidalf_local, NULL,tidal_forcing_type=9, 2,frequencies_local, femmodel, matlitho, vars,verbosecpu);
	}
	/*}}}*/

	//Temporal love numbers
	if (istemporal && !complex_computation){
		/*Initialize*/
		/*Downcast arrays to be exported in parameters*/
		IssmDouble*  pmtf_colineartDouble=NULL;
		IssmDouble*  pmtf_orthotDouble=NULL;

		/* parallel computing */
		LoveNumbers<IssmDouble>* LovefDouble_local=NULL;
		LoveNumbers<IssmDouble>* LovetDouble_local=NULL;
		LoveNumbers<IssmDouble>* TidaltDouble_local=NULL;
		LoveNumbers<IssmDouble>* TidaltDouble=NULL;
		IssmDouble*  pmtf_colineartDouble_local=NULL;
		IssmDouble*  pmtf_orthotDouble_local=NULL;

		doubletype*  pmtf_colineart_local=NULL;
		doubletype*  pmtf_orthot_local=NULL;	
		LoveNumbers<doubletype>* Lovet_local=NULL;
		LoveNumbers<doubletype>* Tidalt_local=NULL;	

		Lovet_local= new LoveNumbers<doubletype>(sh_nmin,sh_nmax,nt_local,lower_row/2/NTit,nt,matlitho);
		Tidalt_local= new LoveNumbers<doubletype>(2,2,nt_local,lower_row/2/NTit,nt,matlitho);	
		pmtf_colineart_local=xNewZeroInit<doubletype>(3*nt_local);
		pmtf_orthot_local=xNewZeroInit<doubletype>(3*nt_local);

		love_freq_to_temporal<doubletype>(Lovet_local,Tidalt_local,pmtf_colineart_local,pmtf_orthot_local,Lovef_local,Tidalf_local,frequencies_local,femmodel,verbosecpu);

		if(VerboseSolution() && verbosecpu) _printf_("   Assembling parallel vectors...");
		//delete Lovef_local;
		delete Tidalf_local;
		//Lovet
		LovetDouble_local= new LoveNumbers<IssmDouble>(sh_nmin,sh_nmax,nt_local,lower_row/2/NTit,nt,matlitho);
		Lovet_local->DownCastToDouble(LovetDouble_local);
		delete Lovet_local;
		LovetDouble_local->Broadcast();

		//Lovef
		LovefDouble_local= new LoveNumbers<IssmDouble>(sh_nmin,sh_nmax,nf_local,lower_row,nfreq,matlitho);
		Lovef_local->DownCastToDouble(LovefDouble_local);
		delete Lovef_local;
		LovefDouble_local->Broadcast();	

		if (forcing_type==11 && sh_nmin<=2 && sh_nmax>=2){			
			TidaltDouble_local= new LoveNumbers<IssmDouble>(2,2,nt_local,lower_row/2/NTit,nt,matlitho);
			TidaltDouble= new LoveNumbers<IssmDouble>(2,2,nt,lower_row/2/NTit,nt,matlitho);
			Tidalt_local->DownCastToDouble(TidaltDouble_local);
			delete Tidalt_local;
			TidaltDouble->LoveMPI_Gather(TidaltDouble_local, lower_row/2/NTit);
			delete TidaltDouble_local;
			//TidaltDouble_local->Broadcast();
		}

		//pmtf:
		pmtf_colineartDouble_local=xNew<IssmDouble>(nt_local);
		pmtf_orthotDouble_local=xNew<IssmDouble>(nt_local);
		/*Downcast*/ /*{{{*/
		for(int i=0;i<nt_local;i++){
			pmtf_colineartDouble_local[i]=DownCastVarToDouble<doubletype>(pmtf_colineart_local[i*3+2]);
			pmtf_orthotDouble_local[i]=DownCastVarToDouble<doubletype>(pmtf_orthot_local[i*3+2]);
		}
		/*}}}*/	
		xDelete<doubletype>(pmtf_colineart_local);
		xDelete<doubletype>(pmtf_orthot_local);
		pmtf_colineartDouble=xNew<IssmDouble>(nt);
		pmtf_orthotDouble=xNew<IssmDouble>(nt);

		int* recvcounts=xNew<int>(IssmComm::GetSize());
		int* displs=xNew<int>(IssmComm::GetSize());
		int  rc;
		int  offset;
		rc=nt_local;
		offset=lower_row/2/NTit;
		ISSM_MPI_Allgather(&rc,1,ISSM_MPI_INT,recvcounts,1,ISSM_MPI_INT,IssmComm::GetComm());
		ISSM_MPI_Allgather(&offset,1,ISSM_MPI_INT,displs,1,ISSM_MPI_INT,IssmComm::GetComm());
		ISSM_MPI_Allgatherv(pmtf_colineartDouble_local, rc, ISSM_MPI_DOUBLE, pmtf_colineartDouble, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());
		ISSM_MPI_Allgatherv(pmtf_orthotDouble_local, rc, ISSM_MPI_DOUBLE, pmtf_orthotDouble, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());
		xDelete<int>(recvcounts);
		xDelete<int>(displs);

		xDelete<IssmDouble>(pmtf_colineartDouble_local);
		xDelete<IssmDouble>(pmtf_orthotDouble_local);
		/*}}}*/	

		if(VerboseSolution() && verbosecpu) _printf_("done\n");

		if(VerboseSolution() && verbosecpu) _printf_("   saving results\n");

		/* Add to parameters */ /*{{{*/
		if(forcing_type==9){ //tidal loading
			femmodel->parameters->AddObject(new DoubleMatParam(TidalLoveHEnum,LovetDouble_local->H,(sh_nmax+1)*nt,1));
			femmodel->parameters->AddObject(new DoubleMatParam(TidalLoveKEnum,LovetDouble_local->K,(sh_nmax+1)*nt,1));
			femmodel->parameters->AddObject(new DoubleMatParam(TidalLoveLEnum,LovetDouble_local->L,(sh_nmax+1)*nt,1));
				}
		else if(forcing_type==11){ //surface loading
			femmodel->parameters->AddObject(new DoubleMatParam(LoadLoveHEnum,LovetDouble_local->H,(sh_nmax+1)*nt,1));
			femmodel->parameters->AddObject(new DoubleMatParam(LoadLoveKEnum,LovetDouble_local->K,(sh_nmax+1)*nt,1));
			femmodel->parameters->AddObject(new DoubleMatParam(LoadLoveLEnum,LovetDouble_local->L,(sh_nmax+1)*nt,1));
			femmodel->parameters->AddObject(new DoubleMatParam(TidalLoveHEnum,TidaltDouble->H,3*nt,1));
			femmodel->parameters->AddObject(new DoubleMatParam(TidalLoveKEnum,TidaltDouble->K,3*nt,1));
			femmodel->parameters->AddObject(new DoubleMatParam(TidalLoveLEnum,TidaltDouble->L,3*nt,1));
			femmodel->parameters->AddObject(new DoubleMatParam(LovePolarMotionTransferFunctionColinearEnum,pmtf_colineartDouble,nt,1));
			femmodel->parameters->AddObject(new DoubleMatParam(LovePolarMotionTransferFunctionOrthogonalEnum,pmtf_orthotDouble,nt,1));
		}
		/*}}}*/	
	
		/*Add into external results*/
		femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveKtEnum,LovetDouble_local->K,nt,sh_nmax+1,0,0));
		femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveHtEnum,LovetDouble_local->H,nt,sh_nmax+1,0,0));
		femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveLtEnum,LovetDouble_local->L,nt,sh_nmax+1,0,0));
		femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveKfEnum,LovefDouble_local->K,nfreq,sh_nmax+1,0,0));
		femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveHfEnum,LovefDouble_local->H,nfreq,sh_nmax+1,0,0));
		femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveLfEnum,LovefDouble_local->L,nfreq,sh_nmax+1,0,0));
		if(forcing_type==11){ 
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveTidalKtEnum,TidaltDouble->K,nt,3,0,0));
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveTidalHtEnum,TidaltDouble->H,nt,3,0,0));
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveTidalLtEnum,TidaltDouble->L,nt,3,0,0));
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LovePMTF1tEnum,pmtf_colineartDouble,nt,1,0,0));
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LovePMTF2tEnum,pmtf_orthotDouble,nt,1,0,0));
		}
		/*Only when love_kernels is on*/
		if (love_kernels==1) {
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveKernelsEnum,LovefDouble_local->Kernels,nfreq,(sh_nmax+1)*(matlitho->numlayers+1)*6,0,0));
		}

		xDelete<IssmDouble>(pmtf_colineartDouble);
		xDelete<IssmDouble>(pmtf_orthotDouble);
		delete LovetDouble_local;
		delete LovefDouble_local;
		delete TidaltDouble;

	}
	else{
		LoveNumbers<IssmDouble>* LovefDouble=NULL;
		LoveNumbers<IssmDouble>* LovefDouble_local=NULL;
		LovefDouble= new LoveNumbers<IssmDouble>(sh_nmin,sh_nmax,nfreq,lower_row,nfreq,matlitho);
		LovefDouble_local= new LoveNumbers<IssmDouble>(sh_nmin,sh_nmax,nf_local,lower_row,nfreq,matlitho);

		LoveNumbers<IssmDouble>* TidalfDouble=NULL;
		LoveNumbers<IssmDouble>* TidalfDouble_local=NULL;
		TidalfDouble= new LoveNumbers<IssmDouble>(2,2,nfreq,lower_row,nfreq,matlitho);
		TidalfDouble_local= new LoveNumbers<IssmDouble>(2,2,nf_local,lower_row,nfreq,matlitho);

		LoveNumbers<IssmDouble>* LovefImag=NULL;
		LoveNumbers<IssmDouble>* LovefImag_local=NULL;
		LoveNumbers<IssmDouble>* TidalfImag=NULL;
		LoveNumbers<IssmDouble>* TidalfImag_local=NULL;

		Lovef_local->DownCastToDouble(LovefDouble_local);
		Tidalf_local->DownCastToDouble(TidalfDouble_local);

		//LovefDouble_local->Broadcast();
		//TidalfDouble_local->Broadcast();

		/*MPI_Gather*/
		if (nfreq>1){
			LovefDouble->LoveMPI_Gather(LovefDouble_local, lower_row);
			if (forcing_type==11 && sh_nmin<=2 && sh_nmax>=2){
				TidalfDouble->LoveMPI_Gather(TidalfDouble_local, lower_row);		
			}
		}
		else{
			Elastic->DownCastToDouble(LovefDouble);
			Tidalf_local->DownCastToDouble(TidalfDouble);
		}

		if (complex_computation){
			LovefImag= new LoveNumbers<IssmDouble>(sh_nmin,sh_nmax,nfreq,lower_row,nfreq,matlitho);
			LovefImag_local= new LoveNumbers<IssmDouble>(sh_nmin,sh_nmax,nf_local,lower_row,nfreq,matlitho);
			TidalfImag= new LoveNumbers<IssmDouble>(2,2,nfreq,lower_row,nfreq,matlitho);
			TidalfImag_local= new LoveNumbers<IssmDouble>(2,2,nf_local,lower_row,nfreq,matlitho);

			Lovef_local->DownCastImagToDouble(LovefImag_local);
			//LovefImag_local->Broadcast();
			LovefImag->LoveMPI_Gather(LovefImag_local, lower_row);

			if (forcing_type==11 && sh_nmin<=2 && sh_nmax>=2){
				Tidalf_local->DownCastImagToDouble(TidalfImag_local);
				//TidalfImag_local->Broadcast();
				TidalfImag->LoveMPI_Gather(TidalfImag_local, lower_row);		
			}
		}

		/*Add into parameters:*/
		if(forcing_type==9){ //tidal loading
			femmodel->parameters->AddObject(new DoubleMatParam(TidalLoveHEnum,LovefDouble->H,(sh_nmax+1)*nfreq,1));
			femmodel->parameters->AddObject(new DoubleMatParam(TidalLoveKEnum,LovefDouble->K,(sh_nmax+1)*nfreq,1));
			femmodel->parameters->AddObject(new DoubleMatParam(TidalLoveLEnum,LovefDouble->L,(sh_nmax+1)*nfreq,1));
		}
		else if(forcing_type==11){ //surface loading
			femmodel->parameters->AddObject(new DoubleMatParam(LoadLoveHEnum,LovefDouble->H,(sh_nmax+1)*nfreq,1));
			femmodel->parameters->AddObject(new DoubleMatParam(LoadLoveKEnum,LovefDouble->K,(sh_nmax+1)*nfreq,1));
			femmodel->parameters->AddObject(new DoubleMatParam(LoadLoveLEnum,LovefDouble->L,(sh_nmax+1)*nfreq,1));
		}

		/*Add into external results:*/
		if (complex_computation){
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveKfEnum,LovefDouble->K,nfreq,sh_nmax+1,0,0));
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveHfEnum,LovefDouble->H,nfreq,sh_nmax+1,0,0));
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveLfEnum,LovefDouble->L,nfreq,sh_nmax+1,0,0));
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveKfiEnum,LovefImag->K,nfreq,sh_nmax+1,0,0));
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveHfiEnum,LovefImag->H,nfreq,sh_nmax+1,0,0));
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveLfiEnum,LovefImag->L,nfreq,sh_nmax+1,0,0));
			///*Only when love_kernels is on*/
			//if (love_kernels==1) {
			//	femmodel->results->AddObject(new GenericExternalResult<IssmComplex*>(femmodel->results->Size()+1,LoveKernelsEnum,Lovef->Kernels,nfreq,(sh_nmax+1)*(matlitho->numlayers+1)*6,0,0));
			//}
		}
		else{
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveKfEnum,LovefDouble->K,nfreq,sh_nmax+1,0,0));
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveHfEnum,LovefDouble->H,nfreq,sh_nmax+1,0,0));
			femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveLfEnum,LovefDouble->L,nfreq,sh_nmax+1,0,0));
			/*Only when love_kernels is on*/
			if (love_kernels==1) {
				femmodel->results->AddObject(new GenericExternalResult<IssmDouble*>(femmodel->results->Size()+1,LoveKernelsEnum,LovefDouble->Kernels,nfreq,(sh_nmax+1)*(matlitho->numlayers+1)*6,0,0));
			}
		}

		delete Lovef_local;
		delete LovefDouble;
		delete LovefDouble_local;
		delete Tidalf_local;
		delete TidalfDouble;
		delete TidalfDouble_local;
		if (complex_computation){
			delete LovefImag_local;
			delete TidalfImag_local;
			delete LovefImag;
			delete TidalfImag;
		}
	}
	/*Free resources:*/
	xDelete<IssmDouble>(frequencies);
	xDelete<IssmDouble>(frequencies_local);
	xDelete<IssmDouble>(frequencies_elastic);
	xDelete<IssmDouble>(frequencies_fluid);

	delete Elastic;
	delete vars;

} /*}}}*/

/*cores and template instantiations:*/
/*template instantiations :{{{*/
// IssmDouble
template void love_core_template<IssmDouble>(FemModel* femmodel);
template LoveVariables<IssmDouble>*	love_init<IssmDouble>(FemModel* femmodel, Matlitho* matlitho,bool verbosecpu);
template void        fill_yi_prefactor<IssmDouble>(IssmDouble* yi_prefactor, int* pdeg, IssmDouble* pomega, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmDouble>* vars);
template void        GetEarthRheology<IssmDouble>(IssmDouble* pla, IssmDouble* pmu, int layer_index, IssmDouble omega,  Matlitho* matlitho, FemModel* femmodel);
template IssmDouble	GetGravity<IssmDouble>(IssmDouble r2, int layer_index, FemModel* femmodel, Matlitho* matlitho,LoveVariables<IssmDouble>* vars);
template void        yi_boundary_conditions<IssmDouble>(IssmDouble* yi_righthandside, int deg, FemModel* femmodel, Matlitho* matlitho,LoveVariables<IssmDouble>* vars, int forcing_type);
template void        yi_derivatives<IssmDouble>(IssmDouble* dydx, IssmDouble* y, int layer_index, int n, IssmDouble* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmDouble>* vars);
template void        propagate_yi_RK2<IssmDouble>(IssmDouble* y, IssmDouble xmin, IssmDouble xmax, int layer_index, IssmDouble* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmDouble>* vars);
template void        propagate_yi_RK4<IssmDouble>(IssmDouble* y, IssmDouble xmin, IssmDouble xmax, int layer_index, IssmDouble* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmDouble>* vars);
template void        propagate_yi_euler<IssmDouble>(IssmDouble* y, IssmDouble xmin, IssmDouble xmax, int layer_index, IssmDouble* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmDouble>* vars);
template void        Innersphere_boundaryconditions<IssmDouble>(IssmDouble* yi, int layer_index, int deg, IssmDouble omega, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmDouble>* vars);
template void        build_yi_system<IssmDouble>(IssmDouble* yi, int deg, IssmDouble omega, IssmDouble* yi_prefactor, FemModel* femmodel, Matlitho* matlitho,LoveVariables<IssmDouble>* vars);
template void        solve_yi_system<IssmDouble>(IssmDouble* loveh, IssmDouble* lovel, IssmDouble* lovek, int deg, IssmDouble omega, IssmDouble* frequencies, IssmDouble* yi, IssmDouble* rhs, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmDouble>* vars,bool verbosecpu);
template void	     compute_love_numbers<IssmDouble>(LoveNumbers<IssmDouble>* Lovef, LoveNumbers<IssmDouble>* Elastic, int forcing_type, int sh_cutoff,IssmDouble* frequencies, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmDouble>* vars, bool verbosecpu);
template IssmDouble  factorial<IssmDouble>(int n);
template IssmDouble* postwidder_coef<IssmDouble>(int NTit);
template IssmDouble  n_C_r<IssmDouble>(int n, int r);
template void         postwidder_transform<IssmDouble>(IssmDouble* Lovet, IssmDouble* Lovef,int d, int t, int sh_nmax,int NTit, IssmDouble* xi, FemModel* femmodel);
template void        EarthRheology<IssmDouble>(LoveVariables<IssmDouble>* vars, IssmDouble* frequencies, int nfreq,  Matlitho* matlitho, FemModel* femmodel);
template IssmDouble HypergeomTableLookup(IssmDouble z1, IssmDouble alpha, IssmDouble* h1, IssmDouble* z, int nz, int nalpha);

//IssmComplex
template void love_core_template<IssmComplex>(FemModel* femmodel);
template LoveVariables<IssmComplex>*	love_init<IssmComplex>(FemModel* femmodel, Matlitho* matlitho,bool verboscpu);
template void        fill_yi_prefactor<IssmComplex>(IssmComplex* yi_prefactor, int* pdeg, IssmComplex* pomega, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmComplex>* vars);
template void        GetEarthRheology<IssmComplex>(IssmComplex* pla, IssmComplex* pmu, int layer_index, IssmComplex omega,  Matlitho* matlitho, FemModel* femmodel);
template IssmComplex	GetGravity<IssmComplex>(IssmComplex r2, int layer_index, FemModel* femmodel, Matlitho* matlitho,LoveVariables<IssmComplex>* vars);
template void        yi_boundary_conditions<IssmComplex>(IssmComplex* yi_righthandside, int deg, FemModel* femmodel, Matlitho* matlitho,LoveVariables<IssmComplex>* vars, int forcing_type);
template void        yi_derivatives<IssmComplex>(IssmComplex* dydx, IssmComplex* y, int layer_index, int n, IssmComplex* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmComplex>* vars);
template void        propagate_yi_RK2<IssmComplex>(IssmComplex* y, IssmComplex xmin, IssmComplex xmax, int layer_index, IssmComplex* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmComplex>* vars);
template void        propagate_yi_RK4<IssmComplex>(IssmComplex* y, IssmComplex xmin, IssmComplex xmax, int layer_index, IssmComplex* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmComplex>* vars);
template void        propagate_yi_euler<IssmComplex>(IssmComplex* y, IssmComplex xmin, IssmComplex xmax, int layer_index, IssmComplex* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmComplex>* vars);
template void        Innersphere_boundaryconditions<IssmComplex>(IssmComplex* yi, int layer_index, int deg, IssmComplex omega, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmComplex>* vars);
template void        build_yi_system<IssmComplex>(IssmComplex* yi, int deg, IssmComplex omega, IssmComplex* yi_prefactor, FemModel* femmodel, Matlitho* matlitho,LoveVariables<IssmComplex>* vars);
template void        solve_yi_system<IssmComplex>(IssmComplex* loveh, IssmComplex* lovel, IssmComplex* lovek, int deg, IssmComplex omega, IssmDouble* frequencies, IssmComplex* yi, IssmComplex* rhs, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmComplex>* vars,bool verbosecpu);
template void	     compute_love_numbers<IssmComplex>(LoveNumbers<IssmComplex>* Lovef, LoveNumbers<IssmComplex>* Elastic, int forcing_type, int sh_cutoff, IssmDouble* frequencies, FemModel* femmodel, Matlitho* matlitho, LoveVariables<IssmComplex>* vars, bool verbosecpu);
template void        EarthRheology<IssmComplex>(LoveVariables<IssmComplex>* vars, IssmDouble* frequencies, int nfreq,  Matlitho* matlitho, FemModel* femmodel);
template IssmComplex HypergeomTableLookup(IssmComplex z1, IssmComplex alpha, IssmDouble* h1, IssmDouble* z, int nz, int nalpha);

//__float128
#ifdef _HAVE_MPLAPACK_
template void love_core_template<__float128>(FemModel* femmodel);
template LoveVariables<__float128>*	love_init<__float128>(FemModel* femmodel, Matlitho* matlitho, bool verbosecpu);
template void        fill_yi_prefactor<__float128>(__float128* yi_prefactor, int* pdeg, __float128* pomega, FemModel* femmodel, Matlitho* matlitho, LoveVariables<__float128>* vars);
template void        GetEarthRheology<__float128>(__float128* pla, __float128* pmu, int layer_index, __float128 omega,  Matlitho* matlitho, FemModel* femmodel);
template __float128	GetGravity<__float128>(__float128 r2, int layer_index, FemModel* femmodel, Matlitho* matlitho,LoveVariables<__float128>* vars);
template void        yi_boundary_conditions<__float128>(__float128* yi_righthandside, int deg, FemModel* femmodel, Matlitho* matlitho,LoveVariables<__float128>* vars, int forcing_type);
template void        yi_derivatives<__float128>(__float128* dydx, __float128* y, int layer_index, int n, __float128* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<__float128>* vars);
template void        propagate_yi_RK2<__float128>(__float128* y, __float128 xmin, __float128 xmax, int layer_index, __float128* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<__float128>* vars);
template void        propagate_yi_RK4<__float128>(__float128* y, __float128 xmin, __float128 xmax, int layer_index, __float128* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<__float128>* vars);
template void        propagate_yi_euler<__float128>(__float128* y, __float128 xmin, __float128 xmax, int layer_index, __float128* yi_prefactor, FemModel* femmodel, Matlitho* matlitho, LoveVariables<__float128>* vars);
template void        Innersphere_boundaryconditions<__float128>(__float128* yi, int layer_index, int deg, __float128 omega, FemModel* femmodel, Matlitho* matlitho, LoveVariables<__float128>* vars);
template void        build_yi_system<__float128>(__float128* yi, int deg, __float128 omega, __float128* yi_prefactor, FemModel* femmodel, Matlitho* matlitho,LoveVariables<__float128>* vars);
template void        solve_yi_system<__float128>(__float128* loveh, __float128* lovel, __float128* lovek, int deg, __float128 omega, IssmDouble* frequencies, __float128* yi, __float128* rhs, FemModel* femmodel, Matlitho* matlitho, LoveVariables<__float128>* vars,bool verbosecpu);
template void	     compute_love_numbers<__float128>(LoveNumbers<__float128>* Lovef, LoveNumbers<__float128>* Elastic, int forcing_type, int sh_cutoff, IssmDouble* frequencies, FemModel* femmodel, Matlitho* matlitho, LoveVariables<__float128>* vars, bool verbosecpu);
template __float128  factorial<__float128>(int n);
template __float128* postwidder_coef<__float128>(int NTit);
template __float128  n_C_r<__float128>(int n, int r);
template void        postwidder_transform<__float128>(__float128* Lovet, __float128* Lovef,int d, int t, int sh_nmax,int NTit, __float128* xi, FemModel* femmodel);
template void        EarthRheology<__float128>(LoveVariables<__float128>* vars, IssmDouble* frequencies, int nfreq,  Matlitho* matlitho, FemModel* femmodel);
template __float128 HypergeomTableLookup(__float128 z1, __float128 alpha, IssmDouble* h1, IssmDouble* z, int nz, int nalpha);
#endif

/*}}}*/
void           love_core(FemModel* femmodel){ /*{{{*/
	bool        complex_computation;
	bool        quad_precision;

	femmodel->parameters->FindParam(&complex_computation,LoveComplexComputationEnum);
	femmodel->parameters->FindParam(&quad_precision,LoveQuadPrecisionEnum);

	if(complex_computation) love_core_template<IssmComplex>(femmodel);
#ifdef _HAVE_MPLAPACK_
	else if(quad_precision) love_core_template<__float128>(femmodel);
#endif
	else                    love_core_template<IssmDouble>(femmodel);

} /*}}}*/
