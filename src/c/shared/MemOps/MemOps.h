/* \file MemOps.h
 * \brief: header file for memory allocations as well as templated new/delete checking for non-null pointers
 */

#ifndef _MEM_OPS_H_
#define _MEM_OPS_H_

#include <cassert>
#include <cstring> /*for memcpy*/

/* memory management of types T with non-trivial constructors require C++ style memory management*/
#define USE_CXX_MEMORY_MANAGMENT_FOR_NON_POD_TYPES
/* but for speed one may alternatively use C memory management but can do so safely only for T that are at most 
 * plain old data structures (POD)*/
#ifndef USE_CXX_MEMORY_MANAGMENT_FOR_NON_POD_TYPES
#include <cstdlib>
#endif 

static char const DEFCONTIG = 'f';

/* AD (mostly ADOLC) is sensitive to calls to ensurecontiguous. These changes limit its use.*/
template <class T> T* xNew(unsigned int size, const char* const contig = &DEFCONTIG){
#ifdef USE_CXX_MEMORY_MANAGMENT_FOR_NON_POD_TYPES
  T* aT_p=new T[size];
  assert(aT_p);
  return aT_p;
#else
  T* aT_p=(T*)malloc(size*sizeof(T));
  assert(aT_p);
  return aT_p;
#endif  
}
template <class T> T** xNew(unsigned int dim1, unsigned int dim2) { /*{{{*/
#ifdef USE_CXX_MEMORY_MANAGMENT_FOR_NON_POD_TYPES
  T* buf=xNew<T>(dim1*dim2);
  T** aT_pp =new T*[dim1];
  assert(aT_pp );
  for (unsigned int i=0;i<dim1;++i) {
    aT_pp [i]=buf;
    buf+=dim2;
  }
  return aT_pp ;
#else
  T* buf=(T*)malloc(dim1*dim2*sizeof(T));
  assert(buf );
  T** aT_pp =(T**)malloc(dim1*sizeof(T*));
  assert(aT_pp );
  for (unsigned int i=0;i<dim1;++i) {
    aT_pp [i]=buf;
    buf+=dim2;
  }
  return aT_pp ;
#endif
}/*}}}*/
// AD (mostly ADOLC) is sensitive to calls to ensurecontiguous. These changes limit its use.
template <class T> T* xNewZeroInit(unsigned int size,const char* const contig = &DEFCONTIG){
#ifdef USE_CXX_MEMORY_MANAGMENT_FOR_NON_POD_TYPES
#ifdef _HAVE_AD_
  T* aT_p=xNew<T>(size,contig);
#else
  T* aT_p=xNew<T>(size);
#endif
  for(unsigned int i=0; i<size;++i) aT_p[i]=(T)0;
  return aT_p;
#else
  T* aT_p=(T*)calloc(size,sizeof(T));
  assert(aT_p);
  return aT_p;
#endif
}
template <class T> T** xNewZeroInit(unsigned int dim1, unsigned int dim2) {/*{{{*/
#ifdef USE_CXX_MEMORY_MANAGMENT_FOR_NON_POD_TYPES
  T** aT_pp=xNew<T>(dim1,dim2);
  for (unsigned int i=0; i<dim1*dim2;++i) (*aT_pp)[i]=(T)0;
  return aT_pp;
#else
  T* buf=(T*)calloc(dim1*dim2*sizeof(T));
  assert(buf );
  T** aT_pp =(T**)malloc(dim1*sizeof(T*));
  assert(aT_pp );
  for (unsigned int i=0;i<dim1;++i) {
    aT_pp [i]=buf;
    buf+=dim2;
  }
  return aT_pp ;
#endif
}/*}}}*/
template <class T> void xDelete(T**& aT_pp) {/*{{{*/
  if (aT_pp) {
#ifdef USE_CXX_MEMORY_MANAGMENT_FOR_NON_POD_TYPES
    delete [](*aT_pp);
    delete [](aT_pp);
#else
    free((void*)*aT_pp)
    free((void**)aT_pp);
#endif
  }
  aT_pp=0;
}/*}}}*/
template <class T> void xDelete(T*& aT_p) {/*{{{*/
  if (aT_p) 
#ifdef USE_CXX_MEMORY_MANAGMENT_FOR_NON_POD_TYPES
    delete []aT_p;
#else
    free((void*)aT_p);
#endif
  aT_p=0;
}/*}}}*/
template <class T> T* xReNew(T* old, unsigned int old_size, unsigned int size) {/*{{{*/
#ifdef USE_CXX_MEMORY_MANAGMENT_FOR_NON_POD_TYPES
	T* aT_p=0;
	if (!old) { // no old memory
		if (size)  
		 aT_p=xNew<T>(size); // according to realloc behavior in manual page 
	}
	else { // have old memory
		if (!size)  // but 0 size
		 xDelete<T>(old); // according to realloc behavior in manual page
		else { // non-zero size
			assert(old_size); // have old memory - need to have old_size set or this call is bad
			// allocate new, delete old; ; even for the case when size is 
			// less than old_size we can't just keep the memory unchanged 
			// because otherwise classes that have ctors/dtors with side-effects 
			// may misbehave, for example classes with static instance/operations counters. 
			aT_p=xNew<T>(size);
			unsigned int iMax=(old_size<size)?old_size:size;
			for (unsigned int i=0; i<iMax;++i) { 
				// we need to copy the items by explicit assignments
				aT_p[i]=old[i];
			}
			xDelete<T>(old);
		}
	}
	return aT_p;
#else
	T* aT_p=0;
	aT_p=(T*)realloc((void*)old,size*sizeof(T));
	if (size) 
	 assert(aT_p); // according to realloc behavior in manual page
	return aT_p;
#endif 
}/*}}}*/
template <class T>  T* xMemCpy(T* dest, const T* src, unsigned int size) {/*{{{*/
  assert(dest); assert(src);
  #if defined(_HAVE_ADOLC_) || defined(_HAVE_CODIPACK_)
  for (int i=0; i<size;++i) dest[i]=src[i];
  #else
  memcpy(dest,src,size*sizeof(T));
  #endif
  return dest;
};
/*}}}*/

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)
#include "../Numerics/types.h"
template <> adouble*  xNew(unsigned int size, const char* const contig);
#endif

#endif
