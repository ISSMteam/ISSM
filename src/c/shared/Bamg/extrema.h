#ifndef _EXTREMA_H_
#define _EXTREMA_H_

namespace bamg {

	template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
	template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}
	template<class T> inline T Max3 (const T &a,const T & b,const T & c){return Max(Max(a,b),c);}

}

#endif
