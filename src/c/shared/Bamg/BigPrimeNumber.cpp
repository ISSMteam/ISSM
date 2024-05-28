#include "./Abs.h"
#include "./extrema.h"

namespace bamg {

	long BigPrimeNumber(long n){
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/AGoodNumberPrimeWith)*/

		/*list of big prime numbers*/
		const long BigPrimeNumber[] ={ 567890359L,
			567890431L,  567890437L,  567890461L,  567890471L,
			567890483L,  567890489L,  567890497L,  567890507L,
			567890591L,  567890599L,  567890621L,  567890629L , 0};

		/*initialize o and pi*/
		long o  = 0;
		long pi = BigPrimeNumber[1];

		/*loop until BigPrimeNumber[i]==0 (end of BigPrimeNumber)*/
		for (int i=0; BigPrimeNumber[i]; i++){

			/*compute r, remainder of the division of BigPrimeNumber[i] by n*/
			long r = BigPrimeNumber[i] % n;

			/*compute oo = min ( r , n-r , |n - 2r|, |n-3r|)*/
			long oo = Min(Min(r,n-r),Min(Abs(n-2*r),Abs(n-3*r)));
			if(o < oo){
				o  = oo;
				pi = BigPrimeNumber[i];
			}
		}
		return pi; 
	}
}
