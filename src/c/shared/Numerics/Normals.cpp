/*!\file:  Normals.cpp
 * \brief Normal vectors for sections.
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./types.h"
#include "../shared.h"

void LineSectionNormal(IssmDouble* normal, IssmDouble* xyz_list) {
	IssmDouble vector[2];
	vector[0]=xyz_list[1*3+0] - xyz_list[0*3+0];
	vector[1]=xyz_list[1*3+1] - xyz_list[0*3+1];

    normal[0] = vector[1];
    normal[1] = -vector[0];

    //normalize
	IssmDouble norm=sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
    norm = max(norm, 1e-10); // avoid NaN
    normal[0] /= norm;
    normal[1] /= norm;
}

void TriangleFacetNormal(IssmDouble* normal, IssmDouble* xyz_list) {
    IssmDouble AB[3];
	IssmDouble AC[3];

	AB[0]=xyz_list[1*3+0] - xyz_list[0*3+0];
	AB[1]=xyz_list[1*3+1] - xyz_list[0*3+1];
	AB[2]=xyz_list[1*3+2] - xyz_list[0*3+2];
	AC[0]=xyz_list[2*3+0] - xyz_list[0*3+0];
	AC[1]=xyz_list[2*3+1] - xyz_list[0*3+1];
	AC[2]=xyz_list[2*3+2] - xyz_list[0*3+2];

	cross(normal,AB,AC);
	
    //normalize
    IssmDouble norm=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    norm = max(norm, 1e-10); // avoid NaN
	for(int i=0;i<3;i++) normal[i] /= norm;
}
