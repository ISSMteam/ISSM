/*!\file VertexCoordinatesx
 * \brief: compute a vector x,y and z of vertex coordinates by
 * marching through all our vertices.
 */

#include "./VertexCoordinatesx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void VertexCoordinatesx( IssmDouble** px, IssmDouble** py, IssmDouble** pz,Vertices* vertices,bool spherical) {

	/*figure out how many vertices we have: */
	int numberofvertices=vertices->NumberOfVertices();

	Vector<IssmDouble>* vx=new Vector<IssmDouble>(numberofvertices);
	Vector<IssmDouble>* vy=new Vector<IssmDouble>(numberofvertices);
	Vector<IssmDouble>* vz=new Vector<IssmDouble>(numberofvertices);

	/*march through our vertices: */
	for(Object* & object : vertices->objects){
		Vertex* vertex=(Vertex*)object;
		vertex->VertexCoordinates(vx,vy,vz,spherical);
	}

	/*Assemble*/
	vx->Assemble();
	vy->Assemble();
	vz->Assemble();

	/*serialize: */
	IssmDouble* x=vx->ToMPISerial();
	IssmDouble* y=vy->ToMPISerial();
	IssmDouble* z=vz->ToMPISerial();

	/*Free resources: */
	delete vx;
	delete vy;
	delete vz;

	/*output: */
	if(px) *px=x;
	else xDelete<IssmDouble>(x);
	if(py) *py=y;
	else xDelete<IssmDouble>(y);
	if(pz) *pz=z;
	else xDelete<IssmDouble>(z);
}
