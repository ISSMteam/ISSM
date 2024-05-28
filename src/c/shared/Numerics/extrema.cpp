/*!\file:  extrema.cpp
 * \brief min and max functions
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./types.h"

#ifndef _HAVE_CODIPACK_// already defined in codipack headers
IssmDouble min(IssmDouble a,IssmDouble b){
	if (a<b)return a;
	else return b;
}
#endif
int min(int a,int b){
	if (a<b)return a;
	else return b;
}
#ifndef _HAVE_CODIPACK_// already defined in codipack headers
IssmDouble max(IssmDouble a,IssmDouble b){
	if (a>b)return a;
	else return b;
}
#endif
int max(int a,int b){
	if (a>b)return a;
	else return b;
}

#ifdef _HAVE_AD_
IssmPDouble  min(IssmPDouble a,IssmPDouble b){
	if (a<b)return a;
	else return b;
}
IssmPDouble  max(IssmPDouble a,IssmPDouble b){
	if (a>b)return a;
	else return b;
}
#endif
