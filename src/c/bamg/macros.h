#ifndef _BAMGMACROS_H
#define _BAMGMACROS_H

#include "./typedefs.h"

namespace bamg {

	const double Pi =3.141592653589793238462643383279502884197169399375105820974944592308;
	const float  fPi=3.141592653589793238462643383279502884197169399375105820974944592308;
	const  int   IsVertexOnGeom = 8;
	const  int   IsVertexOnVertex = 16;
	const  int   IsVertexOnEdge = 32;
	static const short VerticesOfTriangularEdge[3][2] = {{1,2},{2,0},{0,1}};
	static const short EdgesVertexTriangle[3][2] = {{1,2},{2,0},{0,1}};
	static const short OppositeVertex[3] = {0,1,2};
	static const short OppositeEdge[3] =  {0,1,2};
	static const short NextEdge[3] = {1,2,0};
	static const short PreviousEdge[3] = {2,0,1};
	static const short NextVertex[3] = {1,2,0};
	static const short PreviousVertex[3] = {2,0,1};
}

#endif
