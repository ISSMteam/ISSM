/*!\file PentaRef.cpp
 * \brief: implementation of the PentaRef object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
/*}}}*/

/*Element macros*/
#define NUMNODESP0    1
#define NUMNODESP1    6
#define NUMNODESP1_2d 3
#define NUMNODESP1b   7
#define NUMNODESP1xP2 9
#define NUMNODESP1xP3 12
#define NUMNODESP1xP4 15
#define NUMNODESP2xP1 12
#define NUMNODESP2    18
#define NUMNODESP2b   19
#define NUMNODESP2xP4 30
#define NUMNODESMAX   30

/*Object constructors and destructor*/
PentaRef::PentaRef(){/*{{{*/
}
/*}}}*/
PentaRef::~PentaRef(){/*{{{*/
}
/*}}}*/

/*Reference Element numerics*/
void PentaRef::BasalNodeIndices(int* pnumindices,int** pindices,int finiteelement){/*{{{*/

	/*Output*/
	int  numindices;
	int* indices = NULL;

	switch(finiteelement){
		case P1Enum: case P1DGEnum:
			numindices = 3;
			indices    = xNew<int>(numindices);
			indices[0] = 0;
			indices[1] = 1;
			indices[2] = 2;
			break;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			numindices = 3;
			indices    = xNew<int>(numindices);
			indices[0] = 0;
			indices[1] = 1;
			indices[2] = 2;
			break;
		case P2xP1Enum:
			numindices = 6;
			indices    = xNew<int>(numindices);
			indices[0] = 0;
			indices[1] = 1;
			indices[2] = 2;
			indices[3] = 6;
			indices[4] = 7;
			indices[5] = 8;
			break;
		case P1xP2Enum:
			numindices = 3;
			indices    = xNew<int>(numindices);
			indices[0] = 0;
			indices[1] = 1;
			indices[2] = 2;
			break;
		case P1xP3Enum:
			numindices = 3;
			indices    = xNew<int>(numindices);
			indices[0] = 0;
			indices[1] = 1;
			indices[2] = 2;
			break;
		case P1xP4Enum:
			numindices = 3;
			indices    = xNew<int>(numindices);
			indices[0] = 0;
			indices[1] = 1;
			indices[2] = 2;
			break;
		case P2Enum:
			numindices = 6;
			indices    = xNew<int>(numindices);
			indices[0] = 0;
			indices[1] = 1;
			indices[2] = 2;
			indices[3] = 9;
			indices[4] = 10;
			indices[5] = 11;
			break;
		case P2bubbleEnum:
			numindices = 6;
			indices    = xNew<int>(numindices);
			indices[0] = 0;
			indices[1] = 1;
			indices[2] = 2;
			indices[3] = 9;
			indices[4] = 10;
			indices[5] = 11;
			break;
		case P2xP4Enum:
			numindices = 6;
			indices    = xNew<int>(numindices);
			indices[0] = 0;
			indices[1] = 1;
			indices[2] = 2;
			indices[3] = 9;
			indices[4] = 10;
			indices[5] = 11;
			break;
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

	/*Assign output pointer*/
	*pnumindices = numindices;
	*pindices    = indices;
}
/*}}}*/
void PentaRef::GetInputDerivativeValue(IssmDouble* p, IssmDouble* plist,IssmDouble* xyz_list, Gauss* gauss,int finiteelement){/*{{{*/
	/*From node values of parameter p (p_list[0], p_list[1], p_list[2],
	 * p_list[3], p_list[4] and p_list[4]), return parameter derivative value at
	 * gaussian point specified by gauss_coord:
	 *   dp/dx=p_list[0]*dh1/dx+p_list[1]*dh2/dx+p_list[2]*dh3/dx+p_list[3]*dh4/dx+p_list[4]*dh5/dx+p_list[5]*dh6/dx;
	 *   dp/dy=p_list[0]*dh1/dy+p_list[1]*dh2/dy+p_list[2]*dh3/dy+p_list[3]*dh4/dy+p_list[4]*dh5/dy+p_list[5]*dh6/dy;
	 *   dp/dz=p_list[0]*dh1/dz+p_list[1]*dh2/dz+p_list[2]*dh3/dz+p_list[3]*dh4/dz+p_list[4]*dh5/dz+p_list[5]*dh6/dz;
	 *
	 *   p is a vector of size 3x1 already allocated.
	 *
	 * WARNING: For a significant gain in performance, it is better to use
	 * static memory allocation instead of dynamic.
	 */

	/*Allocate derivatives of basis functions*/
	IssmDouble  dbasis[3*NUMNODESMAX];

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->NumberofNodes(finiteelement);
	_assert_(numnodes<=NUMNODESMAX);

	/*Get basis functions derivatives at this point*/
	GetNodalFunctionsDerivatives(&dbasis[0],xyz_list,gauss,finiteelement);

	/*Calculate parameter for this Gauss point*/
	IssmDouble dpx=0.;
	IssmDouble dpy=0.;
	IssmDouble dpz=0.;
	for(int i=0;i<numnodes;i++) dpx += dbasis[0*numnodes+i]*plist[i];
	for(int i=0;i<numnodes;i++) dpy += dbasis[1*numnodes+i]*plist[i];
	for(int i=0;i<numnodes;i++) dpz += dbasis[2*numnodes+i]*plist[i];

	/*Assign values*/
	p[0]=dpx;
	p[1]=dpy;
	p[2]=dpz;
}
/*}}}*/
void PentaRef::GetInputValue(IssmDouble* pvalue,IssmDouble* plist,Gauss* gauss,int finiteelement){/*{{{*/
	/* WARNING: For a significant gain in performance, it is better to use
	 * static memory allocation instead of dynamic.*/

	/*Allocate basis functions*/
	IssmDouble  basis[NUMNODESMAX];

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->NumberofNodes(finiteelement);
	_assert_(numnodes<=NUMNODESMAX);

	/*Get basis functions at this point*/
	GetNodalFunctions(&basis[0],gauss,finiteelement);

	/*Calculate parameter for this Gauss point*/
	IssmDouble value =0.;
	for(int i=0;i<numnodes;i++) value += basis[i]*plist[i];

	/*Assign output pointer*/
	*pvalue = value;
}
/*}}}*/
void PentaRef::GetJacobian(IssmDouble* J, IssmDouble* xyz_list,Gauss* gauss_in){/*{{{*/
	/*The Jacobian is constant over the element, discard the gaussian points. 
	 * J is assumed to have been allocated of size 2x2.*/

	IssmDouble A1,A2,A3;  // area coordinates
	IssmDouble xi,eta,zi; // parametric coordinates
	IssmDouble x1,x2,x3,x4,x5,x6;
	IssmDouble y1,y2,y3,y4,y5,y6;
	IssmDouble z1,z2,z3,z4,z5,z6;
	IssmDouble j_const_reciprocal; // SQRT3/12.0

	/*Cast gauss to GaussPenta*/
	_assert_(gauss_in->Enum()==GaussPentaEnum);
	GaussPenta* gauss = xDynamicCast<GaussPenta*>(gauss_in);

	/*Figure out xi,eta and zi (parametric coordinates), for this gaussian point: */
	A1  = gauss->coord1;
	A2  = gauss->coord2;
	A3  = gauss->coord3;
	xi  = A2-A1;
	eta = SQRT3*A3;
	zi  = gauss->coord4;

	x1=xyz_list[3*0+0];
	x2=xyz_list[3*1+0];
	x3=xyz_list[3*2+0];
	x4=xyz_list[3*3+0];
	x5=xyz_list[3*4+0];
	x6=xyz_list[3*5+0];

	y1=xyz_list[3*0+1];
	y2=xyz_list[3*1+1];
	y3=xyz_list[3*2+1];
	y4=xyz_list[3*3+1];
	y5=xyz_list[3*4+1];
	y6=xyz_list[3*5+1];

	z1=xyz_list[3*0+2];
	z2=xyz_list[3*1+2];
	z3=xyz_list[3*2+2];
	z4=xyz_list[3*3+2];
	z5=xyz_list[3*4+2];
	z6=xyz_list[3*5+2];

	j_const_reciprocal=SQRT3/12;

	J[3*0+0] = 0.25*(x1-x2-x4+x5)*zi+0.25*(-x1+x2-x4+x5);
	J[3*1+0] = j_const_reciprocal*(x1+x2-2*x3-x4-x5+2*x6)*zi+j_const_reciprocal*(-x1-x2+2*x3-x4-x5+2*x6);
	J[3*2+0] = j_const_reciprocal*(x1+x2-2*x3-x4-x5+2*x6)*eta+0.25*(x1-x2-x4+x5)*xi +0.25*(-x1+x5-x2+x4);

	J[3*0+1] = 0.25*(y1-y2-y4+y5)*zi+0.25*(-y1+y2-y4+y5);
	J[3*1+1] = j_const_reciprocal*(y1+y2-2*y3-y4-y5+2*y6)*zi+j_const_reciprocal*(-y1-y2+2*y3-y4-y5+2*y6);
	J[3*2+1] = j_const_reciprocal*(y1+y2-2*y3-y4-y5+2*y6)*eta+0.25*(y1-y2-y4+y5)*xi+0.25*(y4-y1+y5-y2);

	J[3*0+2] = 0.25*(z1-z2-z4+z5)*zi+0.25*(-z1+z2-z4+z5);
	J[3*1+2] = j_const_reciprocal*(z1+z2-2*z3-z4-z5+2*z6)*zi+j_const_reciprocal*(-z1-z2+2*z3-z4-z5+2*z6);
	J[3*2+2] = j_const_reciprocal*(z1+z2-2*z3-z4-z5+2*z6)*eta+0.25*(z1-z2-z4+z5)*xi+0.25*(-z1+z5-z2+z4);
}
/*}}}*/
void PentaRef::GetJacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*On a penta, Jacobian varies according to coordinates. We need to get the Jacobian, and take 
	 * the determinant of it: */
	IssmDouble J[3][3];

	/*Get Jacobian*/
	GetJacobian(&J[0][0],xyz_list,gauss);

	/*Get Determinant*/
	Matrix3x3Determinant(Jdet,&J[0][0]);
	if(*Jdet<0) _error_("negative jacobian determinant!");

}
/*}}}*/
void PentaRef::GetJacobianInvert(IssmDouble* Jinv, IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	/*Jacobian*/
	IssmDouble J[3][3];

	/*Call Jacobian routine to get the jacobian:*/
	GetJacobian(&J[0][0], xyz_list, gauss);

	/*Invert Jacobian matrix: */
	Matrix3x3Invert(Jinv,&J[0][0]);
}
/*}}}*/
void PentaRef::GetNodalFunctions(IssmDouble* basis,Gauss* gauss_in,int finiteelement){/*{{{*/
	/*This routine returns the values of the nodal functions  at the gaussian point.*/

	_assert_(basis);

	/*Cast gauss to GaussPenta*/
	_assert_(gauss_in->Enum()==GaussPentaEnum);
	GaussPenta* gauss = xDynamicCast<GaussPenta*>(gauss_in);

	/*Get current coordinates in reference element*/
	IssmDouble zeta=gauss->coord4;

	switch(finiteelement){
		case P0Enum: 
			basis[0]=1.;
			return;
		case P1Enum: case P1DGEnum:
			basis[0]=gauss->coord1*(1.-zeta)/2.;
			basis[1]=gauss->coord2*(1.-zeta)/2.;
			basis[2]=gauss->coord3*(1.-zeta)/2.;
			basis[3]=gauss->coord1*(1.+zeta)/2.;
			basis[4]=gauss->coord2*(1.+zeta)/2.;
			basis[5]=gauss->coord3*(1.+zeta)/2.;
			return;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			basis[0]=gauss->coord1*(1.-zeta)/2.;
			basis[1]=gauss->coord2*(1.-zeta)/2.;
			basis[2]=gauss->coord3*(1.-zeta)/2.;
			basis[3]=gauss->coord1*(1.+zeta)/2.;
			basis[4]=gauss->coord2*(1.+zeta)/2.;
			basis[5]=gauss->coord3*(1.+zeta)/2.;
			basis[6]=27.*gauss->coord1*gauss->coord2*gauss->coord3*(1.+zeta)*(1.-zeta);
			return;
		#ifndef _HAVE_AD_ /*speed up CoDiPack Compilation by hiding higher order elements*/
		case P2xP1Enum:
			/*Corner nodes*/
			basis[ 0]=gauss->coord1*(2.*gauss->coord1-1.)*(1.-zeta)/2.;
			basis[ 1]=gauss->coord2*(2.*gauss->coord2-1.)*(1.-zeta)/2.;
			basis[ 2]=gauss->coord3*(2.*gauss->coord3-1.)*(1.-zeta)/2.;
			basis[ 3]=gauss->coord1*(2.*gauss->coord1-1.)*(1.+zeta)/2.;
			basis[ 4]=gauss->coord2*(2.*gauss->coord2-1.)*(1.+zeta)/2.;
			basis[ 5]=gauss->coord3*(2.*gauss->coord3-1.)*(1.+zeta)/2.;
			/*mid-sides of triangles*/
			basis[ 6]=4.*gauss->coord3*gauss->coord2*(1.-zeta)/2.;
			basis[ 7]=4.*gauss->coord3*gauss->coord1*(1.-zeta)/2.;
			basis[ 8]=4.*gauss->coord1*gauss->coord2*(1.-zeta)/2.;
			basis[ 9]=4.*gauss->coord3*gauss->coord2*(1.+zeta)/2.;
			basis[10]=4.*gauss->coord3*gauss->coord1*(1.+zeta)/2.;
			basis[11]=4.*gauss->coord1*gauss->coord2*(1.+zeta)/2.;
			return;
		case P1xP2Enum:
			/*Corner nodes*/
			basis[ 0]=gauss->coord1*zeta*(zeta-1.)/2.;
			basis[ 1]=gauss->coord2*zeta*(zeta-1.)/2.;
			basis[ 2]=gauss->coord3*zeta*(zeta-1.)/2.;
			basis[ 3]=gauss->coord1*zeta*(zeta+1.)/2.;
			basis[ 4]=gauss->coord2*zeta*(zeta+1.)/2.;
			basis[ 5]=gauss->coord3*zeta*(zeta+1.)/2.;
			/*mid-sides of quads*/
			basis[ 6]=gauss->coord1*(1.-zeta*zeta);
			basis[ 7]=gauss->coord2*(1.-zeta*zeta);
			basis[ 8]=gauss->coord3*(1.-zeta*zeta);
			return;
		case P2Enum:
			/*Corner nodes*/
			basis[ 0]=gauss->coord1*(2.*gauss->coord1-1.)*zeta*(zeta-1.)/2.;
			basis[ 1]=gauss->coord2*(2.*gauss->coord2-1.)*zeta*(zeta-1.)/2.;
			basis[ 2]=gauss->coord3*(2.*gauss->coord3-1.)*zeta*(zeta-1.)/2.;
			basis[ 3]=gauss->coord1*(2.*gauss->coord1-1.)*zeta*(zeta+1.)/2.;
			basis[ 4]=gauss->coord2*(2.*gauss->coord2-1.)*zeta*(zeta+1.)/2.;
			basis[ 5]=gauss->coord3*(2.*gauss->coord3-1.)*zeta*(zeta+1.)/2.;
			/*mid-sides of quads*/
			basis[ 6]=gauss->coord1*(2.*gauss->coord1-1.)*(1.-zeta*zeta);
			basis[ 7]=gauss->coord2*(2.*gauss->coord2-1.)*(1.-zeta*zeta);
			basis[ 8]=gauss->coord3*(2.*gauss->coord3-1.)*(1.-zeta*zeta);
			/*mid-sides of triangles*/
			basis[ 9]=4.*gauss->coord3*gauss->coord2*zeta*(zeta-1.)/2.;
			basis[10]=4.*gauss->coord3*gauss->coord1*zeta*(zeta-1.)/2.;
			basis[11]=4.*gauss->coord1*gauss->coord2*zeta*(zeta-1.)/2.;
			basis[12]=4.*gauss->coord3*gauss->coord2*zeta*(zeta+1.)/2.;
			basis[13]=4.*gauss->coord3*gauss->coord1*zeta*(zeta+1.)/2.;
			basis[14]=4.*gauss->coord1*gauss->coord2*zeta*(zeta+1.)/2.;
			/*quad faces*/
			basis[15]=4.*gauss->coord3*gauss->coord2*(1.-zeta*zeta);
			basis[16]=4.*gauss->coord3*gauss->coord1*(1.-zeta*zeta);
			basis[17]=4.*gauss->coord1*gauss->coord2*(1.-zeta*zeta);
			return;
		case P2bubbleEnum: case P2bubblecondensedEnum:
			/*Corner nodes*/
			basis[ 0]=gauss->coord1*(2.*gauss->coord1-1.)*zeta*(zeta-1.)/2.;
			basis[ 1]=gauss->coord2*(2.*gauss->coord2-1.)*zeta*(zeta-1.)/2.;
			basis[ 2]=gauss->coord3*(2.*gauss->coord3-1.)*zeta*(zeta-1.)/2.;
			basis[ 3]=gauss->coord1*(2.*gauss->coord1-1.)*zeta*(zeta+1.)/2.;
			basis[ 4]=gauss->coord2*(2.*gauss->coord2-1.)*zeta*(zeta+1.)/2.;
			basis[ 5]=gauss->coord3*(2.*gauss->coord3-1.)*zeta*(zeta+1.)/2.;
			/*mid-sides of quads*/
			basis[ 6]=gauss->coord1*(2.*gauss->coord1-1.)*(1.-zeta*zeta);
			basis[ 7]=gauss->coord2*(2.*gauss->coord2-1.)*(1.-zeta*zeta);
			basis[ 8]=gauss->coord3*(2.*gauss->coord3-1.)*(1.-zeta*zeta);
			/*mid-sides of triangles*/
			basis[ 9]=4.*gauss->coord3*gauss->coord2*zeta*(zeta-1.)/2.;
			basis[10]=4.*gauss->coord3*gauss->coord1*zeta*(zeta-1.)/2.;
			basis[11]=4.*gauss->coord1*gauss->coord2*zeta*(zeta-1.)/2.;
			basis[12]=4.*gauss->coord3*gauss->coord2*zeta*(zeta+1.)/2.;
			basis[13]=4.*gauss->coord3*gauss->coord1*zeta*(zeta+1.)/2.;
			basis[14]=4.*gauss->coord1*gauss->coord2*zeta*(zeta+1.)/2.;
			/*quad faces*/
			basis[15]=4.*gauss->coord3*gauss->coord2*(1.-zeta*zeta);
			basis[16]=4.*gauss->coord3*gauss->coord1*(1.-zeta*zeta);
			basis[17]=4.*gauss->coord1*gauss->coord2*(1.-zeta*zeta);
			/*bubble*/
			basis[18]=27.*gauss->coord1*gauss->coord2*gauss->coord3*(1.+zeta)*(1.-zeta);
			return;
		case P2xP4Enum :
			/*Corner nodes*/
			basis[ 0]=gauss->coord1*(2.*gauss->coord1-1.)*(2./3.)*(zeta-1.)*(zeta-0.5 )*(zeta)*(zeta+0.5);
			basis[ 1]=gauss->coord2*(2.*gauss->coord2-1.)*(2./3.)*(zeta-1.)*(zeta-0.5 )*(zeta)*(zeta+0.5);
			basis[ 2]=gauss->coord3*(2.*gauss->coord3-1.)*(2./3.)*(zeta-1.)*(zeta-0.5 )*(zeta)*(zeta+0.5);
			basis[ 3]=gauss->coord1*(2.*gauss->coord1-1.)*(2./3.)*(zeta-0.5)*(zeta)*(zeta+0.5)*(zeta +1.);
			basis[ 4]=gauss->coord2*(2.*gauss->coord2-1.)*(2./3.)*(zeta-0.5)*(zeta)*(zeta+0.5)*(zeta +1.);
			basis[ 5]=gauss->coord3*(2.*gauss->coord3-1.)*(2./3.)*(zeta-0.5)*(zeta)*(zeta+0.5)*(zeta +1.);
			/*mid-sides of quads*/
			basis[ 6]=gauss->coord1*(2.*gauss->coord1-1)*4.*(zeta-1.)*(zeta-0.5)*(zeta+0.5)*(zeta+1.);
			basis[ 7]=gauss->coord2*(2.*gauss->coord2-1)*4.*(zeta-1.)*(zeta-0.5)*(zeta+0.5)*(zeta+1.);
			basis[ 8]=gauss->coord3*(2.*gauss->coord3-1)*4.*(zeta-1.)*(zeta-0.5)*(zeta+0.5)*(zeta+1.);
			/*mid-sides of triangles*/
			basis[ 9]=4.*gauss->coord2*gauss->coord3*(2./3.)*(zeta-1.)*(zeta-0.5)*zeta*(zeta+0.5);
			basis[10]=4.*gauss->coord1*gauss->coord3*(2./3.)*(zeta-1.)*(zeta-0.5)*zeta*(zeta+0.5);
			basis[11]=4.*gauss->coord1*gauss->coord2*(2./3.)*(zeta-1.)*(zeta-0.5)*zeta*(zeta+0.5);
			basis[12]=4.*gauss->coord2*gauss->coord3*(2./3.)*(zeta-0.5)*zeta*(zeta+0.5)*(zeta+1.);
			basis[13]=4.*gauss->coord1*gauss->coord3*(2./3.)*(zeta-0.5)*zeta*(zeta+0.5)*(zeta+1.);
			basis[14]=4.*gauss->coord1*gauss->coord2*(2./3.)*(zeta-0.5)*zeta*(zeta+0.5)*(zeta+1.);
			/*quarter-sides of quads*/
			basis[15]=gauss->coord1*(2.*gauss->coord1-1.)*(-8./3.)*(zeta-1.0)*(zeta-0.5)*zeta*(zeta+1.);
			basis[16]=gauss->coord2*(2.*gauss->coord2-1.)*(-8./3.)*(zeta-1.0)*(zeta-0.5)*zeta*(zeta+1.);
			basis[17]=gauss->coord3*(2.*gauss->coord3-1.)*(-8./3.)*(zeta-1.0)*(zeta-0.5)*zeta*(zeta+1.);
			basis[18]=gauss->coord1*(2.*gauss->coord1-1.)*(-8./3.)*(zeta-1.0)*zeta*(zeta+0.5)*(zeta+1.);
			basis[19]=gauss->coord2*(2.*gauss->coord2-1.)*(-8./3.)*(zeta-1.0)*zeta*(zeta+0.5)*(zeta+1.);
			basis[20]=gauss->coord3*(2.*gauss->coord3-1.)*(-8./3.)*(zeta-1.0)*zeta*(zeta+0.5)*(zeta+1.);
			/* mid-sides of interior triangles*/
			basis[21]=4.*gauss->coord2*gauss->coord3*(-8./3.)*(zeta-1.0)*(zeta-0.5)*zeta*(zeta+1.);
			basis[22]=4.*gauss->coord1*gauss->coord3*(-8./3.)*(zeta-1.0)*(zeta-0.5)*zeta*(zeta+1.);
			basis[23]=4.*gauss->coord1*gauss->coord2*(-8./3.)*(zeta-1.0)*(zeta-0.5)*zeta*(zeta+1.);
			basis[24]=4.*gauss->coord2*gauss->coord3*4.*(zeta-1.0)*(zeta-0.5)*(zeta+0.5)*(zeta+1.);
			basis[25]=4.*gauss->coord1*gauss->coord3*4.*(zeta-1.0)*(zeta-0.5)*(zeta+0.5)*(zeta+1.);
			basis[26]=4.*gauss->coord1*gauss->coord2*4.*(zeta-1.0)*(zeta-0.5)*(zeta+0.5)*(zeta+1.);
			basis[27]=4.*gauss->coord2*gauss->coord3*(-8./3.)*(zeta-1.0)*zeta*(zeta+0.5)*(zeta+1.);
			basis[28]=4.*gauss->coord1*gauss->coord3*(-8./3.)*(zeta-1.0)*zeta*(zeta+0.5)*(zeta+1.);
			basis[29]=4.*gauss->coord1*gauss->coord2*(-8./3.)*(zeta-1.0)*zeta*(zeta+0.5)*(zeta+1.);
			return;
		case P1xP3Enum :
			/*Corner nodes*/
			basis[ 0]=-(9.)/(16.)*gauss->coord1*(zeta-1)*(zeta-1./3.)*(zeta+1./3.);	
			basis[ 1]=-(9.)/(16.)*gauss->coord2*(zeta-1)*(zeta-1./3.)*(zeta+1./3.);	
			basis[ 2]=-(9.)/(16.)*gauss->coord3*(zeta-1)*(zeta-1./3.)*(zeta+1./3.);	
			basis[ 3]=(9.)/(16.)*gauss->coord1*(zeta-1./3.)*(zeta+1./3.)*(zeta+1.);
			basis[ 4]=(9.)/(16.)*gauss->coord2*(zeta-1./3.)*(zeta+1./3.)*(zeta+1.);
			basis[ 5]=(9.)/(16.)*gauss->coord3*(zeta-1./3.)*(zeta+1./3.)*(zeta+1.);
			/*third-sides of quads*/
			basis[ 6]=(27.)/(16.)*gauss->coord1*(zeta-1)*(zeta-1./3.)*(zeta+1.);
			basis[ 7]=(27.)/(16.)*gauss->coord2*(zeta-1)*(zeta-1./3.)*(zeta+1.);
			basis[ 8]=(27.)/(16.)*gauss->coord3*(zeta-1)*(zeta-1./3.)*(zeta+1.);
			basis[ 9]=-(27.)/(16.)*gauss->coord1*(zeta-1)*(zeta+1./3.)*(zeta+1.);
			basis[10]=-(27.)/(16.)*gauss->coord2*(zeta-1)*(zeta+1./3.)*(zeta+1.);
			basis[11]=-(27.)/(16.)*gauss->coord3*(zeta-1)*(zeta+1./3.)*(zeta+1.);
			return;
		case P1xP4Enum :
			/*Corner nodes*/
			basis[ 0]=gauss->coord1*(2./3.)*(zeta-1.)*(zeta-0.5 )*(zeta)*(zeta+0.5);
			basis[ 1]=gauss->coord2*(2./3.)*(zeta-1.)*(zeta-0.5 )*(zeta)*(zeta+0.5);
			basis[ 2]=gauss->coord3*(2./3.)*(zeta-1.)*(zeta-0.5 )*(zeta)*(zeta+0.5);
			basis[ 3]=gauss->coord1*(2./3.)*(zeta-0.5)*(zeta)*(zeta+0.5)*(zeta +1.);
			basis[ 4]=gauss->coord2*(2./3.)*(zeta-0.5)*(zeta)*(zeta+0.5)*(zeta +1.);
			basis[ 5]=gauss->coord3*(2./3.)*(zeta-0.5)*(zeta)*(zeta+0.5)*(zeta +1.);
			/*mid-sides of quads (center of vertical edges)*/
			basis[ 6]=gauss->coord1*4.*(zeta-1.)*(zeta-0.5)*(zeta+0.5)*(zeta+1.);
			basis[ 7]=gauss->coord2*4.*(zeta-1.)*(zeta-0.5)*(zeta+0.5)*(zeta+1.);
			basis[ 8]=gauss->coord3*4.*(zeta-1.)*(zeta-0.5)*(zeta+0.5)*(zeta+1.);
			/*quarter-sides of quads (-0.5 and +0.5 of vertical edges)*/
			basis[ 9]=gauss->coord1*(-8./3.)*(zeta-1.0)*(zeta-0.5)*zeta*(zeta+1.);
			basis[10]=gauss->coord2*(-8./3.)*(zeta-1.0)*(zeta-0.5)*zeta*(zeta+1.);
			basis[11]=gauss->coord3*(-8./3.)*(zeta-1.0)*(zeta-0.5)*zeta*(zeta+1.);
			basis[12]=gauss->coord1*(-8./3.)*(zeta-1.0)*zeta*(zeta+0.5)*(zeta+1.);
			basis[13]=gauss->coord2*(-8./3.)*(zeta-1.0)*zeta*(zeta+0.5)*(zeta+1.);
			basis[14]=gauss->coord3*(-8./3.)*(zeta-1.0)*zeta*(zeta+0.5)*(zeta+1.);
			return;
		#endif
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}
}
/*}}}*/
void PentaRef::GetNodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list, Gauss* gauss,int finiteelement){/*{{{*/

	/*This routine returns the values of the nodal functions derivatives  (with respect to the 
	 * actual coordinate system): */
	IssmDouble    Jinv[3][3];

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->NumberofNodes(finiteelement);

	/*Get nodal functions derivatives in reference triangle*/
	IssmDouble dbasis_ref[3*NUMNODESMAX];
	GetNodalFunctionsDerivativesReference(dbasis_ref,gauss,finiteelement);

	/*Get Jacobian invert: */
	GetJacobianInvert(&Jinv[0][0], xyz_list, gauss);

	/*Build dbasis: 
	 *
	 * [dhi/dx]= Jinv'*[dhi/dr]
	 * [dhi/dy]        [dhi/ds]
	 * [dhi/dz]        [dhi/dzeta]
	 */

	for(int i=0;i<numnodes;i++){
		dbasis[numnodes*0+i]=Jinv[0][0]*dbasis_ref[0*numnodes+i]+Jinv[0][1]*dbasis_ref[1*numnodes+i]+Jinv[0][2]*dbasis_ref[2*numnodes+i];
		dbasis[numnodes*1+i]=Jinv[1][0]*dbasis_ref[0*numnodes+i]+Jinv[1][1]*dbasis_ref[1*numnodes+i]+Jinv[1][2]*dbasis_ref[2*numnodes+i];
		dbasis[numnodes*2+i]=Jinv[2][0]*dbasis_ref[0*numnodes+i]+Jinv[2][1]*dbasis_ref[1*numnodes+i]+Jinv[2][2]*dbasis_ref[2*numnodes+i];
	}
}
/*}}}*/
void PentaRef::GetNodalFunctionsDerivativesReference(IssmDouble* dbasis,Gauss* gauss_in,int finiteelement){/*{{{*/

	/*This routine returns the values of the nodal functions derivatives  (with respect to the 
	 * natural coordinate system) at the gaussian point. */

	_assert_(dbasis && gauss_in);

	/*Cast gauss to GaussPenta*/
	_assert_(gauss_in->Enum()==GaussPentaEnum);
	GaussPenta* gauss = xDynamicCast<GaussPenta*>(gauss_in);

	/*Get current coordinates in reference element*/
	IssmDouble zeta=gauss->coord4;

	switch(finiteelement){
		case P0Enum: 
			/*Zero derivative*/
			dbasis[NUMNODESP0*0+0]   = 0.;
			dbasis[NUMNODESP0*1+0]   = 0.;
			dbasis[NUMNODESP0*2+0]   = 0.;
			return;
		case P1Enum: case P1DGEnum:
			/*Nodal function 1*/
			dbasis[NUMNODESP1*0+0]   = (zeta-1.)/4.;
			dbasis[NUMNODESP1*1+0]   = SQRT3/12.*(zeta-1.);
			dbasis[NUMNODESP1*2+0]   = -.5*gauss->coord1;
			/*Nodal function 2*/
			dbasis[NUMNODESP1*0+1]   = (1.-zeta)/4.;
			dbasis[NUMNODESP1*1+1]   = SQRT3/12.*(zeta-1);
			dbasis[NUMNODESP1*2+1]   = -.5*gauss->coord2;
			/*Nodal function 3*/
			dbasis[NUMNODESP1*0+2]   = 0.;
			dbasis[NUMNODESP1*1+2]   = SQRT3/6.*(1.-zeta);
			dbasis[NUMNODESP1*2+2]   = -.5*gauss->coord3;
			/*Nodal function 4*/
			dbasis[NUMNODESP1*0+3]   = -(1.+zeta)/4.;
			dbasis[NUMNODESP1*1+3]   = -SQRT3/12.*(1.+zeta);
			dbasis[NUMNODESP1*2+3]   = .5*gauss->coord1;
			/*Nodal function 5*/
			dbasis[NUMNODESP1*0+4]   = (1.+zeta)/4.;
			dbasis[NUMNODESP1*1+4]   = -SQRT3/12.*(1.+zeta);
			dbasis[NUMNODESP1*2+4]   = .5*gauss->coord2;
			/*Nodal function 6*/
			dbasis[NUMNODESP1*0+5]   = 0.;
			dbasis[NUMNODESP1*1+5]   = SQRT3/6.*(1.+zeta);
			dbasis[NUMNODESP1*2+5]   = .5*gauss->coord3;
			return;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			/*Nodal function 1*/
			dbasis[NUMNODESP1b*0+0]   = (zeta-1.)/4.;
			dbasis[NUMNODESP1b*1+0]   = SQRT3/12.*(zeta-1.);
			dbasis[NUMNODESP1b*2+0]   = -.5*gauss->coord1;
			/*Nodal function 2*/
			dbasis[NUMNODESP1b*0+1]   = (1.-zeta)/4.;
			dbasis[NUMNODESP1b*1+1]   = SQRT3/12.*(zeta-1);
			dbasis[NUMNODESP1b*2+1]   = -.5*gauss->coord2;
			/*Nodal function 3*/
			dbasis[NUMNODESP1b*0+2]   = 0.;
			dbasis[NUMNODESP1b*1+2]   = SQRT3/6.*(1.-zeta);
			dbasis[NUMNODESP1b*2+2]   = -.5*gauss->coord3;
			/*Nodal function 4*/
			dbasis[NUMNODESP1b*0+3]   = -(1.+zeta)/4.;
			dbasis[NUMNODESP1b*1+3]   = -SQRT3/12.*(1.+zeta);
			dbasis[NUMNODESP1b*2+3]   = .5*gauss->coord1;
			/*Nodal function 5*/
			dbasis[NUMNODESP1b*0+4]   = (1.+zeta)/4.;
			dbasis[NUMNODESP1b*1+4]   = -SQRT3/12.*(1.+zeta);
			dbasis[NUMNODESP1b*2+4]   = .5*gauss->coord2;
			/*Nodal function 6*/
			dbasis[NUMNODESP1b*0+5]   = 0.;
			dbasis[NUMNODESP1b*1+5]   = SQRT3/6.*(1.+zeta);
			dbasis[NUMNODESP1b*2+5]   = .5*gauss->coord3;
			/*Nodal function 7*/
			dbasis[NUMNODESP1b*0+6] = 27.*(1.+zeta)*(1.-zeta)*(-.5*gauss->coord2*gauss->coord3 + .5*gauss->coord1*gauss->coord3);
			dbasis[NUMNODESP1b*1+6] = 27.*(1.+zeta)*(1.-zeta)*SQRT3*(-1./6.*gauss->coord2*gauss->coord3 - 1./6.*gauss->coord1*gauss->coord3 +1./3.*gauss->coord1*gauss->coord2);
			dbasis[NUMNODESP1b*2+6] = -54*gauss->coord1*gauss->coord2*gauss->coord3*zeta;
			return;
		#ifndef _HAVE_AD_ /*speed up CoDiPack Compilation by hiding higher order elements*/
		case P2xP1Enum:
			/*Nodal function 1*/
			dbasis[NUMNODESP2xP1*0+0 ] = .5*(1.-zeta)*(-2.*gauss->coord1 + 0.5);
			dbasis[NUMNODESP2xP1*1+0 ] = .5*(1.-zeta)*(-2.*SQRT3/3.*gauss->coord1 + SQRT3/6.);
			dbasis[NUMNODESP2xP1*2+0 ] = -.5*gauss->coord1*(2.*gauss->coord1-1.);
			/*Nodal function 2*/
			dbasis[NUMNODESP2xP1*0+1 ] = .5*(1.-zeta)*(+2.*gauss->coord2 - 0.5);
			dbasis[NUMNODESP2xP1*1+1 ] = .5*(1.-zeta)*(-2.*SQRT3/3.*gauss->coord2 + SQRT3/6.);
			dbasis[NUMNODESP2xP1*2+1 ] = -.5*gauss->coord2*(2.*gauss->coord2-1.);
			/*Nodal function 3*/
			dbasis[NUMNODESP2xP1*0+2 ] = 0.;
			dbasis[NUMNODESP2xP1*1+2 ] = .5*(1.-zeta)*(4.*SQRT3/3.*gauss->coord3 - SQRT3/3.);
			dbasis[NUMNODESP2xP1*2+2 ] = -.5*gauss->coord3*(2.*gauss->coord3-1.);
			/*Nodal function 4*/
			dbasis[NUMNODESP2xP1*0+3 ] = .5*(1.+zeta)*(-2.*gauss->coord1 + 0.5);
			dbasis[NUMNODESP2xP1*1+3 ] = .5*(1.+zeta)*(-2.*SQRT3/3.*gauss->coord1 + SQRT3/6.);
			dbasis[NUMNODESP2xP1*2+3 ] = .5*gauss->coord1*(2.*gauss->coord1-1.);
			/*Nodal function 5*/
			dbasis[NUMNODESP2xP1*0+4 ] = .5*(1.+zeta)*(+2.*gauss->coord2 - 0.5);
			dbasis[NUMNODESP2xP1*1+4 ] = .5*(1.+zeta)*(-2.*SQRT3/3.*gauss->coord2 + SQRT3/6.);
			dbasis[NUMNODESP2xP1*2+4 ] = .5*gauss->coord2*(2.*gauss->coord2-1.);
			/*Nodal function 6*/
			dbasis[NUMNODESP2xP1*0+5 ] = 0.;
			dbasis[NUMNODESP2xP1*1+5 ] = .5*(1.+zeta)*(4.*SQRT3/3.*gauss->coord3 - SQRT3/3.);
			dbasis[NUMNODESP2xP1*2+5 ] = .5*gauss->coord3*(2.*gauss->coord3-1.);

			/*Nodal function 7*/
			dbasis[NUMNODESP2xP1*0+6 ] = (1.-zeta)*gauss->coord3;
			dbasis[NUMNODESP2xP1*1+6 ] = .5*(1.-zeta)*(+4.*SQRT3/3.*gauss->coord2 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2xP1*2+6 ] = -2.*gauss->coord3*gauss->coord2;
			/*Nodal function 8*/
			dbasis[NUMNODESP2xP1*0+7 ] = -(1.-zeta)*gauss->coord3;
			dbasis[NUMNODESP2xP1*1+7 ] = .5*(1.-zeta)*(+4.*SQRT3/3.*gauss->coord1 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2xP1*2+7 ] = -2.*gauss->coord3*gauss->coord1;
			/*Nodal function 9*/
			dbasis[NUMNODESP2xP1*0+8 ] = (1.-zeta)*(gauss->coord1-gauss->coord2);
			dbasis[NUMNODESP2xP1*1+8 ] = .5*(1.-zeta)*(-2.*SQRT3/3.*(gauss->coord1+gauss->coord2));
			dbasis[NUMNODESP2xP1*2+8 ] = -2.*gauss->coord1*gauss->coord2;
			/*Nodal function 10*/
			dbasis[NUMNODESP2xP1*0+9 ] = (1.+zeta)*gauss->coord3;
			dbasis[NUMNODESP2xP1*1+9 ] = .5*(1.+zeta)*(+4.*SQRT3/3.*gauss->coord2 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2xP1*2+9 ] = 2.*gauss->coord3*gauss->coord2;
			/*Nodal function 11*/
			dbasis[NUMNODESP2xP1*0+10] = -(1.+zeta)*gauss->coord3;
			dbasis[NUMNODESP2xP1*1+10] = .5*(1.+zeta)*(+4.*SQRT3/3.*gauss->coord1 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2xP1*2+10] = 2.*gauss->coord3*gauss->coord1;
			/*Nodal function 12*/
			dbasis[NUMNODESP2xP1*0+11] = (1.+zeta)*(gauss->coord1-gauss->coord2);
			dbasis[NUMNODESP2xP1*1+11] = .5*(1.+zeta)*(-2.*SQRT3/3.*(gauss->coord1+gauss->coord2));
			dbasis[NUMNODESP2xP1*2+11] = 2.*gauss->coord1*gauss->coord2;
			return;
		#endif
		#ifndef _HAVE_AD_ /*speed up CoDiPack Compilation by hiding higher order elements*/
		case P1xP2Enum:
			/*Nodal function 1*/
			dbasis[NUMNODESP1xP2*0+0]   = -zeta*(zeta-1.)/4.;
			dbasis[NUMNODESP1xP2*1+0]   = -SQRT3/12.*zeta*(zeta-1.);
			dbasis[NUMNODESP1xP2*2+0]   = .5*(2.*zeta-1.)*gauss->coord1;
			/*Nodal function 2*/
			dbasis[NUMNODESP1xP2*0+1]   = zeta*(zeta-1.)/4.;
			dbasis[NUMNODESP1xP2*1+1]   = -SQRT3/12.*zeta*(zeta-1);
			dbasis[NUMNODESP1xP2*2+1]   = .5*(2.*zeta-1.)*gauss->coord2;
			/*Nodal function 3*/
			dbasis[NUMNODESP1xP2*0+2]   = 0.;
			dbasis[NUMNODESP1xP2*1+2]   = SQRT3/6.*zeta*(zeta-1.);
			dbasis[NUMNODESP1xP2*2+2]   = .5*(2.*zeta-1.)*gauss->coord3;
			/*Nodal function 4*/
			dbasis[NUMNODESP1xP2*0+3]   = -zeta*(zeta+1)/4.;
			dbasis[NUMNODESP1xP2*1+3]   = -SQRT3/12.*zeta*(zeta+1.);
			dbasis[NUMNODESP1xP2*2+3]   = .5*(2.*zeta+1.)*gauss->coord1;
			/*Nodal function 5*/
			dbasis[NUMNODESP1xP2*0+4]   = zeta*(zeta+1.)/4.;
			dbasis[NUMNODESP1xP2*1+4]   = -SQRT3/12.*zeta*(zeta+1.);
			dbasis[NUMNODESP1xP2*2+4]   = .5*(2.*zeta+1.)*gauss->coord2;
			/*Nodal function 6*/
			dbasis[NUMNODESP1xP2*0+5]   = 0.;
			dbasis[NUMNODESP1xP2*1+5]   = SQRT3/6.*zeta*(zeta+1.);
			dbasis[NUMNODESP1xP2*2+5]   = .5*(2.*zeta+1.)*gauss->coord3;

			/*Nodal function 7*/
			dbasis[NUMNODESP1xP2*0+6 ] = -0.5*(1.-zeta*zeta);
			dbasis[NUMNODESP1xP2*1+6 ] = -SQRT3/6.*(1.-zeta*zeta);
			dbasis[NUMNODESP1xP2*2+6 ] = -2.*zeta*gauss->coord1;
			/*Nodal function 8*/
			dbasis[NUMNODESP1xP2*0+7 ] = 0.5*(1.-zeta*zeta);
			dbasis[NUMNODESP1xP2*1+7 ] = -SQRT3/6.*(1.-zeta*zeta);
			dbasis[NUMNODESP1xP2*2+7 ] = -2.*zeta*gauss->coord2;
			/*Nodal function 9*/
			dbasis[NUMNODESP1xP2*0+8 ] = 0.;
			dbasis[NUMNODESP1xP2*1+8 ] = SQRT3/3.*(1.-zeta*zeta);
			dbasis[NUMNODESP1xP2*2+8 ] = -2.*zeta*gauss->coord3;
			return;
		#endif
		#ifndef _HAVE_AD_ /*speed up CoDiPack Compilation by hiding higher order elements*/
		case P2Enum:
			/*Nodal function 1*/
			dbasis[NUMNODESP2*0+0 ] = .5*zeta*(zeta-1.)*(-2.*gauss->coord1 + 0.5);
			dbasis[NUMNODESP2*1+0 ] = .5*zeta*(zeta-1.)*(-2.*SQRT3/3.*gauss->coord1 + SQRT3/6.);
			dbasis[NUMNODESP2*2+0 ] = .5*(2.*zeta-1.)*gauss->coord1*(2.*gauss->coord1-1.);
			/*Nodal function 2*/
			dbasis[NUMNODESP2*0+1 ] = .5*zeta*(zeta-1.)*(+2.*gauss->coord2 - 0.5);
			dbasis[NUMNODESP2*1+1 ] = .5*zeta*(zeta-1.)*(-2.*SQRT3/3.*gauss->coord2 + SQRT3/6.);
			dbasis[NUMNODESP2*2+1 ] = .5*(2.*zeta-1.)*gauss->coord2*(2.*gauss->coord2-1.);
			/*Nodal function 3*/
			dbasis[NUMNODESP2*0+2 ] = 0.;
			dbasis[NUMNODESP2*1+2 ] = .5*zeta*(zeta-1.)*(4.*SQRT3/3.*gauss->coord3 - SQRT3/3.);
			dbasis[NUMNODESP2*2+2 ] = .5*(2.*zeta-1.)*gauss->coord3*(2.*gauss->coord3-1.);
			/*Nodal function 4*/
			dbasis[NUMNODESP2*0+3 ] = .5*zeta*(zeta+1.)*(-2.*gauss->coord1 + 0.5);
			dbasis[NUMNODESP2*1+3 ] = .5*zeta*(zeta+1.)*(-2.*SQRT3/3.*gauss->coord1 + SQRT3/6.);
			dbasis[NUMNODESP2*2+3 ] = .5*(2.*zeta+1.)*gauss->coord1*(2.*gauss->coord1-1.);
			/*Nodal function 5*/
			dbasis[NUMNODESP2*0+4 ] = .5*zeta*(zeta+1.)*(+2.*gauss->coord2 - 0.5);
			dbasis[NUMNODESP2*1+4 ] = .5*zeta*(zeta+1.)*(-2.*SQRT3/3.*gauss->coord2 + SQRT3/6.);
			dbasis[NUMNODESP2*2+4 ] = .5*(2.*zeta+1.)*gauss->coord2*(2.*gauss->coord2-1.);
			/*Nodal function 6*/
			dbasis[NUMNODESP2*0+5 ] = 0.;
			dbasis[NUMNODESP2*1+5 ] = .5*zeta*(zeta+1.)*(4.*SQRT3/3.*gauss->coord3 - SQRT3/3.);
			dbasis[NUMNODESP2*2+5 ] = .5*(2.*zeta+1.)*gauss->coord3*(2.*gauss->coord3-1.);

			/*Nodal function 7*/
			dbasis[NUMNODESP2*0+6 ] = (-2.*gauss->coord1 + 0.5)*(1.-zeta*zeta);
			dbasis[NUMNODESP2*1+6 ] = (-2.*SQRT3/3.*gauss->coord1 + SQRT3/6.)*(1.-zeta*zeta);
			dbasis[NUMNODESP2*2+6 ] = -2.*zeta*gauss->coord1*(2.*gauss->coord1-1.);
			/*Nodal function 8*/
			dbasis[NUMNODESP2*0+7 ] = (+2.*gauss->coord2 - 0.5)*(1.-zeta*zeta);
			dbasis[NUMNODESP2*1+7 ] = (-2.*SQRT3/3.*gauss->coord2 + SQRT3/6.)*(1.-zeta*zeta);
			dbasis[NUMNODESP2*2+7 ] = -2.*zeta*gauss->coord2*(2.*gauss->coord2-1.);
			/*Nodal function 9*/
			dbasis[NUMNODESP2*0+8 ] = 0.;
			dbasis[NUMNODESP2*1+8 ] = (+4.*SQRT3/3.*gauss->coord3 - SQRT3/3.)*(1.-zeta*zeta);
			dbasis[NUMNODESP2*2+8 ] = -2.*zeta*gauss->coord3*(2.*gauss->coord3-1.);

			/*Nodal function 10*/
			dbasis[NUMNODESP2*0+9 ] = zeta*(zeta-1.)*gauss->coord3;
			dbasis[NUMNODESP2*1+9 ] = .5*zeta*(zeta-1.)*(+4.*SQRT3/3.*gauss->coord2 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2*2+9 ] = 2.*gauss->coord3*gauss->coord2*(2.*zeta-1.);
			/*Nodal function 11*/
			dbasis[NUMNODESP2*0+10] = -zeta*(zeta-1.)*gauss->coord3;
			dbasis[NUMNODESP2*1+10] = .5*zeta*(zeta-1.)*(+4.*SQRT3/3.*gauss->coord1 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2*2+10] = 2.*gauss->coord3*gauss->coord1*(2.*zeta-1.);
			/*Nodal function 12*/
			dbasis[NUMNODESP2*0+11] = zeta*(zeta-1.)*(gauss->coord1-gauss->coord2);
			dbasis[NUMNODESP2*1+11] = .5*zeta*(zeta-1.)*(-2.*SQRT3/3.*(gauss->coord1+gauss->coord2));
			dbasis[NUMNODESP2*2+11] = 2.*gauss->coord1*gauss->coord2*(2.*zeta-1.);
			/*Nodal function 13*/
			dbasis[NUMNODESP2*0+12] = zeta*(zeta+1.)*gauss->coord3;
			dbasis[NUMNODESP2*1+12] = .5*zeta*(zeta+1.)*(+4.*SQRT3/3.*gauss->coord2 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2*2+12] = 2.*gauss->coord3*gauss->coord2*(2.*zeta+1.);
			/*Nodal function 14*/
			dbasis[NUMNODESP2*0+13] = -zeta*(zeta+1.)*gauss->coord3;
			dbasis[NUMNODESP2*1+13] = .5*zeta*(zeta+1.)*(+4.*SQRT3/3.*gauss->coord1 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2*2+13] = 2.*gauss->coord3*gauss->coord1*(2.*zeta+1.);
			/*Nodal function 15*/
			dbasis[NUMNODESP2*0+14] = zeta*(zeta+1.)*(gauss->coord1-gauss->coord2);
			dbasis[NUMNODESP2*1+14] = .5*zeta*(zeta+1.)*(-2.*SQRT3/3.*(gauss->coord1+gauss->coord2));
			dbasis[NUMNODESP2*2+14] = 2.*gauss->coord1*gauss->coord2*(2.*zeta+1.);

			/*Nodal function 16*/
			dbasis[NUMNODESP2*0+15] = 2.*gauss->coord3*(1.-zeta*zeta);
			dbasis[NUMNODESP2*1+15] = (1.-zeta*zeta)*(+4.*SQRT3/3.*gauss->coord2 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2*2+15] = -2.*zeta*4.*gauss->coord3*gauss->coord2;
			/*Nodal function 17*/
			dbasis[NUMNODESP2*0+16] = -2.*gauss->coord3*(1.-zeta*zeta);
			dbasis[NUMNODESP2*1+16] = (1.-zeta*zeta)*(+4.*SQRT3/3.*gauss->coord1 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2*2+16] = -2.*zeta*4.*gauss->coord3*gauss->coord1;
			/*Nodal function 18*/
			dbasis[NUMNODESP2*0+17] = 2.*(gauss->coord1-gauss->coord2)*(1.-zeta*zeta);
			dbasis[NUMNODESP2*1+17] = (1.-zeta*zeta)*(-2.*SQRT3/3.*(gauss->coord1+gauss->coord2));
			dbasis[NUMNODESP2*2+17] = -2.*zeta*4.*gauss->coord1*gauss->coord2;
			return;
		case P2bubbleEnum: case P2bubblecondensedEnum:
			/*Nodal function 1*/
			dbasis[NUMNODESP2b*0+0 ] = .5*zeta*(zeta-1.)*(-2.*gauss->coord1 + 0.5);
			dbasis[NUMNODESP2b*1+0 ] = .5*zeta*(zeta-1.)*(-2.*SQRT3/3.*gauss->coord1 + SQRT3/6.);
			dbasis[NUMNODESP2b*2+0 ] = .5*(2.*zeta-1.)*gauss->coord1*(2.*gauss->coord1-1.);
			/*Nodal function 2*/
			dbasis[NUMNODESP2b*0+1 ] = .5*zeta*(zeta-1.)*(+2.*gauss->coord2 - 0.5);
			dbasis[NUMNODESP2b*1+1 ] = .5*zeta*(zeta-1.)*(-2.*SQRT3/3.*gauss->coord2 + SQRT3/6.);
			dbasis[NUMNODESP2b*2+1 ] = .5*(2.*zeta-1.)*gauss->coord2*(2.*gauss->coord2-1.);
			/*Nodal function 3*/
			dbasis[NUMNODESP2b*0+2 ] = 0.;
			dbasis[NUMNODESP2b*1+2 ] = .5*zeta*(zeta-1.)*(4.*SQRT3/3.*gauss->coord3 - SQRT3/3.);
			dbasis[NUMNODESP2b*2+2 ] = .5*(2.*zeta-1.)*gauss->coord3*(2.*gauss->coord3-1.);
			/*Nodal function 4*/
			dbasis[NUMNODESP2b*0+3 ] = .5*zeta*(zeta+1.)*(-2.*gauss->coord1 + 0.5);
			dbasis[NUMNODESP2b*1+3 ] = .5*zeta*(zeta+1.)*(-2.*SQRT3/3.*gauss->coord1 + SQRT3/6.);
			dbasis[NUMNODESP2b*2+3 ] = .5*(2.*zeta+1.)*gauss->coord1*(2.*gauss->coord1-1.);
			/*Nodal function 5*/
			dbasis[NUMNODESP2b*0+4 ] = .5*zeta*(zeta+1.)*(+2.*gauss->coord2 - 0.5);
			dbasis[NUMNODESP2b*1+4 ] = .5*zeta*(zeta+1.)*(-2.*SQRT3/3.*gauss->coord2 + SQRT3/6.);
			dbasis[NUMNODESP2b*2+4 ] = .5*(2.*zeta+1.)*gauss->coord2*(2.*gauss->coord2-1.);
			/*Nodal function 6*/
			dbasis[NUMNODESP2b*0+5 ] = 0.;
			dbasis[NUMNODESP2b*1+5 ] = .5*zeta*(zeta+1.)*(4.*SQRT3/3.*gauss->coord3 - SQRT3/3.);
			dbasis[NUMNODESP2b*2+5 ] = .5*(2.*zeta+1.)*gauss->coord3*(2.*gauss->coord3-1.);

			/*Nodal function 7*/
			dbasis[NUMNODESP2b*0+6 ] = (-2.*gauss->coord1 + 0.5)*(1.-zeta*zeta);
			dbasis[NUMNODESP2b*1+6 ] = (-2.*SQRT3/3.*gauss->coord1 + SQRT3/6.)*(1.-zeta*zeta);
			dbasis[NUMNODESP2b*2+6 ] = -2.*zeta*gauss->coord1*(2.*gauss->coord1-1.);
			/*Nodal function 8*/
			dbasis[NUMNODESP2b*0+7 ] = (+2.*gauss->coord2 - 0.5)*(1.-zeta*zeta);
			dbasis[NUMNODESP2b*1+7 ] = (-2.*SQRT3/3.*gauss->coord2 + SQRT3/6.)*(1.-zeta*zeta);
			dbasis[NUMNODESP2b*2+7 ] = -2.*zeta*gauss->coord2*(2.*gauss->coord2-1.);
			/*Nodal function 9*/
			dbasis[NUMNODESP2b*0+8 ] = 0.;
			dbasis[NUMNODESP2b*1+8 ] = (+4.*SQRT3/3.*gauss->coord3 - SQRT3/3.)*(1.-zeta*zeta);
			dbasis[NUMNODESP2b*2+8 ] = -2.*zeta*gauss->coord3*(2.*gauss->coord3-1.);

			/*Nodal function 10*/
			dbasis[NUMNODESP2b*0+9 ] = zeta*(zeta-1.)*gauss->coord3;
			dbasis[NUMNODESP2b*1+9 ] = .5*zeta*(zeta-1.)*(+4.*SQRT3/3.*gauss->coord2 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2b*2+9 ] = 2.*gauss->coord3*gauss->coord2*(2.*zeta-1.);
			/*Nodal function 11*/
			dbasis[NUMNODESP2b*0+10] = -zeta*(zeta-1.)*gauss->coord3;
			dbasis[NUMNODESP2b*1+10] = .5*zeta*(zeta-1.)*(+4.*SQRT3/3.*gauss->coord1 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2b*2+10] = 2.*gauss->coord3*gauss->coord1*(2.*zeta-1.);
			/*Nodal function 12*/
			dbasis[NUMNODESP2b*0+11] = zeta*(zeta-1.)*(gauss->coord1-gauss->coord2);
			dbasis[NUMNODESP2b*1+11] = .5*zeta*(zeta-1.)*(-2.*SQRT3/3.*(gauss->coord1+gauss->coord2));
			dbasis[NUMNODESP2b*2+11] = 2.*gauss->coord1*gauss->coord2*(2.*zeta-1.);
			/*Nodal function 13*/
			dbasis[NUMNODESP2b*0+12] = zeta*(zeta+1.)*gauss->coord3;
			dbasis[NUMNODESP2b*1+12] = .5*zeta*(zeta+1.)*(+4.*SQRT3/3.*gauss->coord2 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2b*2+12] = 2.*gauss->coord3*gauss->coord2*(2.*zeta+1.);
			/*Nodal function 14*/
			dbasis[NUMNODESP2b*0+13] = -zeta*(zeta+1.)*gauss->coord3;
			dbasis[NUMNODESP2b*1+13] = .5*zeta*(zeta+1.)*(+4.*SQRT3/3.*gauss->coord1 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2b*2+13] = 2.*gauss->coord3*gauss->coord1*(2.*zeta+1.);
			/*Nodal function 15*/
			dbasis[NUMNODESP2b*0+14] = zeta*(zeta+1.)*(gauss->coord1-gauss->coord2);
			dbasis[NUMNODESP2b*1+14] = .5*zeta*(zeta+1.)*(-2.*SQRT3/3.*(gauss->coord1+gauss->coord2));
			dbasis[NUMNODESP2b*2+14] = 2.*gauss->coord1*gauss->coord2*(2.*zeta+1.);

			/*Nodal function 16*/
			dbasis[NUMNODESP2b*0+15] = 2.*gauss->coord3*(1.-zeta*zeta);
			dbasis[NUMNODESP2b*1+15] = (1.-zeta*zeta)*(+4.*SQRT3/3.*gauss->coord2 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2b*2+15] = -2.*zeta*4.*gauss->coord3*gauss->coord2;
			/*Nodal function 17*/
			dbasis[NUMNODESP2b*0+16] = -2.*gauss->coord3*(1.-zeta*zeta);
			dbasis[NUMNODESP2b*1+16] = (1.-zeta*zeta)*(+4.*SQRT3/3.*gauss->coord1 - 2.*SQRT3/3.*gauss->coord3);
			dbasis[NUMNODESP2b*2+16] = -2.*zeta*4.*gauss->coord3*gauss->coord1;
			/*Nodal function 18*/
			dbasis[NUMNODESP2b*0+17] = 2.*(gauss->coord1-gauss->coord2)*(1.-zeta*zeta);
			dbasis[NUMNODESP2b*1+17] = (1.-zeta*zeta)*(-2.*SQRT3/3.*(gauss->coord1+gauss->coord2));
			dbasis[NUMNODESP2b*2+17] = -2.*zeta*4.*gauss->coord1*gauss->coord2;

			/*Nodal function 19*/
			dbasis[NUMNODESP2b*0+18] = 27.*(1.+zeta)*(1.-zeta)*(-.5*gauss->coord2*gauss->coord3 + .5*gauss->coord1*gauss->coord3);
			dbasis[NUMNODESP2b*1+18] = 27.*(1.+zeta)*(1.-zeta)*SQRT3*(-1./6.*gauss->coord2*gauss->coord3 - 1./6.*gauss->coord1*gauss->coord3 +1./3.*gauss->coord1*gauss->coord2);
			dbasis[NUMNODESP2b*2+18] = -54*gauss->coord1*gauss->coord2*gauss->coord3*zeta;
			return;
		case P2xP4Enum :
			/*Nodal function 1*/
			dbasis[NUMNODESP2xP4*0+0 ] = (-2* gauss->coord1 + 0.5 ) *(2./3.) *(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5) ; 
			dbasis[NUMNODESP2xP4*1+0 ] = (-((2.*SQRT3)/(3.))*gauss->coord1 + (SQRT3/6.) )*(2./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5) ; 
			dbasis[NUMNODESP2xP4*2+0 ] =  gauss->coord1 *(2.* gauss->coord1 -1)* 2./3.*( (2.*zeta-1)*(zeta -0.5)*(zeta +0.5) + 2.* zeta *zeta *(zeta -1.)); 
			/*Nodal function 2*/
			dbasis[NUMNODESP2xP4*0+1 ] = (2.*gauss->coord2 - 0.5 ) *(2./3.) *(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5) ; 
			dbasis[NUMNODESP2xP4*1+1 ] = (-((2.*SQRT3)/(3.))*gauss->coord2 + (SQRT3/6.) )*(2./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5) ; 
			dbasis[NUMNODESP2xP4*2+1 ] = gauss->coord2*(2.*gauss->coord2 -1.)* 2./3.* ((2.*zeta-1.)*(zeta -0.5)*(zeta +0.5) + 2. * zeta *zeta*(zeta -1.)); 
			/*Nodal function 3*/
			dbasis[NUMNODESP2xP4*0+2 ] = 0. ; 
			dbasis[NUMNODESP2xP4*1+2 ] = (((4.*SQRT3)/(3.))*gauss->coord3 - (SQRT3)/(3.))*(2./3.)*(zeta -1.)*(zeta-0.5)*(zeta)*(zeta+0.5); 
			dbasis[NUMNODESP2xP4*2+2 ] = gauss->coord3*(2.* gauss->coord3 -1.)* 2./3.*( (2.*zeta-1.)*(zeta -0.5)*(zeta +0.5) + 2.*zeta *zeta*(zeta -1.)); 
			/*Nodal function 4*/
			dbasis[NUMNODESP2xP4*0+3 ] = (-2.* gauss->coord1 + 0.5 ) *(2./3.) *(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP2xP4*1+3 ] = (-((2.*SQRT3)/(3.)) *gauss->coord1 + (SQRT3)/(6.) ) *(2./3.) *(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP2xP4*2+3 ] = gauss->coord1*(2.*gauss->coord1 -1.)* 2./3.*( (2.*zeta+1.)*(zeta -0.5)*(zeta +0.5) + 2.*zeta *zeta*(zeta +1.)); 
			/*Nodal function 5*/
			dbasis[NUMNODESP2xP4*0+4 ] = (2*gauss->coord2 - 0.5 ) *(2./3.) *(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP2xP4*1+4 ] = (-((2.*SQRT3)/(3.)) *gauss->coord2 + (SQRT3/6.))*(2./3.)*(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. );
			dbasis[NUMNODESP2xP4*2+4 ] = gauss->coord2 *(2.*gauss->coord2 -1.)* 2./3.*( (2.*zeta+1.)*(zeta -0.5)*(zeta +0.5) + 2.*zeta *zeta*(zeta +1.)); 
			/*Nodal function 6*/
			dbasis[NUMNODESP2xP4*0+5 ] = 0. ; 
			dbasis[NUMNODESP2xP4*1+5 ] = (((4.*SQRT3)/(3.))*gauss->coord3 - SQRT3/3. ) *(2./3.) *(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP2xP4*2+5 ] = gauss->coord3 *(2.*gauss->coord3 -1.)* 2./3.*( (2.*zeta+1.)*(zeta -0.5)*(zeta +0.5) + 2.*zeta *zeta*(zeta +1)); 
			/*Nodal function 7*/
			dbasis[NUMNODESP2xP4*0+6 ] =  (-2.* gauss->coord1 + 0.5 ) * 4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*1+6 ] =  (-((2.*SQRT3)/(3.)) *gauss->coord1 + (SQRT3)/(6.) )* 4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1.)  ; 
			dbasis[NUMNODESP2xP4*2+6 ] =  gauss->coord1*(2.*gauss->coord1-1)* 4.*( 4.*zeta *zeta*zeta - (5./2.)*zeta ); 
			/*Nodal function 8*/
			dbasis[NUMNODESP2xP4*0+7 ] =  (2*gauss->coord2 - 0.5 )* 4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*1+7 ] =  (-((2.*SQRT3)/(3.)) *gauss->coord2 + (SQRT3)/(6.)) * 4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. )  ; 
			dbasis[NUMNODESP2xP4*2+7 ] =  gauss->coord2*(2.*gauss->coord2-1)* 4.*( 4.*zeta *zeta*zeta - (5./2.)*zeta ); 
			/*Nodal function 9*/
			dbasis[NUMNODESP2xP4*0+8 ] = 0. ; 
			dbasis[NUMNODESP2xP4*1+8 ] = (((4.*SQRT3)/(3.))*gauss->coord3 - SQRT3/3. ) * 4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. )  ; 
			dbasis[NUMNODESP2xP4*2+8 ] = gauss->coord3*(2.*gauss->coord3-1)* 4.*( 4.*zeta *zeta*zeta - (5./2.)*zeta ); 
			/*Nodal function 10*/
			dbasis[NUMNODESP2xP4*0+9 ] = 2.*gauss->coord3 * 2./3.*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5) ; 
			dbasis[NUMNODESP2xP4*1+9 ] = (4.* SQRT3/3.* gauss->coord2- 2.*SQRT3/3. *gauss->coord3) *(2./3.) *(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5);
			dbasis[NUMNODESP2xP4*2+9 ] = 4.* gauss->coord2 * gauss->coord3 *(2./3.)*((2.*zeta-1.)*(zeta -0.5)*(zeta +0.5) + 2.* zeta*zeta*(zeta -1.)); 
			/*Nodal function 11*/
			dbasis[NUMNODESP2xP4*0+10] = -2.* gauss->coord3* 2./3.*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5) ; 
			dbasis[NUMNODESP2xP4*1+10] = (4.*SQRT3/3.*gauss->coord1 - 2.*SQRT3/3.*gauss->coord3) * (2./3.) *(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5);
			dbasis[NUMNODESP2xP4*2+10] = 4.* gauss->coord3*gauss->coord1 *(2./3.)*( (2*zeta-1.)*(zeta -0.5)*(zeta +0.5) + 2* zeta * zeta*(zeta -1));
			/*Nodal function 12*/
			dbasis[NUMNODESP2xP4*0+11] = 2.* (gauss->coord1- gauss->coord2)* (2./3.) *(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5) ; 
			dbasis[NUMNODESP2xP4*1+11] = -2.* SQRT3/3.*(gauss->coord2 +gauss->coord1) *(2./3.) *(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5);
			dbasis[NUMNODESP2xP4*2+11] = 4.* gauss->coord1*gauss->coord2 *(2./3.) *( (2.*zeta-1)*(zeta -0.5)*(zeta +0.5) + 2* zeta* zeta*(zeta -1));
			/*Nodal function 13*/
			dbasis[NUMNODESP2xP4*0+12] = 2.* gauss->coord3 * 2./3.*(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*1+12] = (4.*SQRT3/3.* gauss->coord2 - 2.*SQRT3/3. *gauss->coord3) *(2./3.) *(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. );
			dbasis[NUMNODESP2xP4*2+12] = 4.*gauss->coord2*gauss->coord3 *(2./3.)*( (2.*zeta+1.)*(zeta -0.5)*(zeta +0.5) + 2*zeta*zeta*(zeta +1.)); 
			/*Nodal function 14*/
			dbasis[NUMNODESP2xP4*0+13] = - 2.*gauss->coord3* 2./3.*(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*1+13] = (4.*SQRT3/3.*gauss->coord1- 2.*SQRT3/3.*gauss->coord3) * (2./3.) *(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. );
			dbasis[NUMNODESP2xP4*2+13] = 4.*gauss->coord3*gauss->coord1 *(2./3.)*( (2.*zeta+1.)*(zeta -0.5)*(zeta +0.5) + 2.*zeta *zeta*(zeta +1.)); 
			/*Nodal function 15*/
			dbasis[NUMNODESP2xP4*0+14] = 2.* (gauss->coord1- gauss->coord2)* (2./3.) *(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*1+14] = -2.* SQRT3/3.*(gauss->coord2 +gauss->coord1) *(2./3.) *(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. );
			dbasis[NUMNODESP2xP4*2+14] = 4.*gauss->coord1*gauss->coord2 *(2./3.)*( (2.*zeta+1.)*(zeta -0.5)*(zeta +0.5) + 2* zeta *zeta*(zeta +1.)); 
			/*Nodal function 16*/
			dbasis[NUMNODESP2xP4*0+15] = (-2.* gauss->coord1 + 0.5 )* (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*1+15] = (-2.*SQRT3/3. *gauss->coord1 + SQRT3/6.) * (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+15] = gauss->coord1*(2.*gauss->coord1-1) * (-8./3.)*((2.*zeta -1.)*(zeta-0.5)*(zeta +1.) +zeta*(zeta -1.)*( 2.*zeta + 0.5)); 
			/*Nodal function 17*/
			dbasis[NUMNODESP2xP4*0+16] = (2*gauss->coord2 - 0.5 )* (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*1+16] = (-2.*SQRT3/3. *gauss->coord2 + SQRT3/6.)* (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+16] = gauss->coord2*(2.*gauss->coord2-1) * (-8./3.)*((2.*zeta -1.)*(zeta-0.5)*(zeta +1.) +zeta *(zeta -1.)*( 2.*zeta + 0.5)); 
			/*Nodal function 18*/
			dbasis[NUMNODESP2xP4*0+17] = 0. ; 
			dbasis[NUMNODESP2xP4*1+17] = (4.*SQRT3/3.*gauss->coord3 - SQRT3/3. )* (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+17] = gauss->coord3*(2*gauss->coord3-1) * (-8./3.)*((2.*zeta-1.)*(zeta-0.5)*(zeta +1.) +zeta *(zeta -1.)*( 2.*zeta + 0.5));
			/*Nodal function 19*/
			dbasis[NUMNODESP2xP4*0+18] = (-2.* gauss->coord1 + 0.5 ) * (-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*1+18] =  (-2.*SQRT3/3. *gauss->coord1 + SQRT3/6. )* (-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+18] = gauss->coord1*(2.*gauss->coord1-1) * (-8./3.)*((2.*zeta -1. ) *(zeta+0.5)* (zeta +1.) +  zeta* (zeta -1.)*( 2.*zeta + 3./2.));
			/*Nodal function 20*/
			dbasis[NUMNODESP2xP4*0+19] = (2*gauss->coord2 - 0.5 )* (-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*1+19] = (-2.*SQRT3/3.*gauss->coord2 + SQRT3/6.) * (-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+19] = gauss->coord2*(2.*gauss->coord2-1)* (-8./3.)*((2.*zeta -1. )*(zeta+0.5)*(zeta +1.) +  zeta*(zeta -1.)*( 2.*zeta + 3./2.)); 
			/*Nodal function 21*/
			dbasis[NUMNODESP2xP4*0+20] = 0 ; 
			dbasis[NUMNODESP2xP4*1+20] = (4.*SQRT3/3.*gauss->coord3 - SQRT3/3.)* (-8./3.)*(zeta - 1.)*(zeta + 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+20] = gauss->coord3*(2*gauss->coord3-1) * (-8./3.)*((2. *zeta -1. )*(zeta+0.5)*(zeta +1.) +  zeta*(zeta -1.)*( 2.*zeta + 3./2.));
			/*Nodal function 22*/
			dbasis[NUMNODESP2xP4*0+21] = 2. *gauss->coord3 * (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ); 
			dbasis[NUMNODESP2xP4*1+21] = (4.* SQRT3/3.*gauss->coord2- 2.*SQRT3/3.*gauss->coord3) * (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+21] = 4.*gauss->coord2 *gauss->coord3* (-8./3.)*((2.*zeta -1. )*(zeta-0.5)*(zeta +1.) +  zeta*(zeta -1.)*( 2.*zeta + 0.5)); 
			/*Nodal function 23*/
			dbasis[NUMNODESP2xP4*0+22] = -2. *gauss->coord3 *( -8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ); 
			dbasis[NUMNODESP2xP4*1+22] = (4.* SQRT3/3.*gauss->coord1- 2.*SQRT3/3.*gauss->coord3 )*(-8./3.)*(zeta -1.)*(zeta - 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+22] = 4.*gauss->coord1*gauss->coord3* (-8./3.)*((2.*zeta -1. )*(zeta-0.5)*(zeta +1.) +  zeta*(zeta -1.)*( 2.*zeta + 0.5)); 
			/*Nodal function 24*/
			dbasis[NUMNODESP2xP4*0+23] = 2.*(gauss->coord1- gauss->coord2) * (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ); 
			dbasis[NUMNODESP2xP4*1+23] = -2.*SQRT3/3.*(gauss->coord2+ gauss->coord1) * (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+23] = 4.*gauss->coord1*gauss->coord2* (-8./3.)*((2.*zeta -1. )* (zeta-0.5) *(zeta +1.) +  zeta* (zeta -1.)*( 2.*zeta + 0.5));
			/*Nodal function 25*/
			dbasis[NUMNODESP2xP4*0+24] = 2. *gauss->coord3 *4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP2xP4*1+24] = (4.*SQRT3/3.*gauss->coord2 - 2.* SQRT3/3. *gauss->coord3) *4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+24] = 4.* gauss->coord2 * gauss->coord3* 4.*( 4.* zeta*zeta*zeta - 5./2. *zeta ); 
			/*Nodal function 26*/
			dbasis[NUMNODESP2xP4*0+25] = -2. *gauss->coord3*4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP2xP4*1+25] = (4.*SQRT3/3.*gauss->coord1- 2.*SQRT3/3.*gauss->coord3 )*4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+25] = 4. * gauss->coord1 * gauss->coord3 *4.*( 4. *zeta*zeta*zeta - 5./2.* zeta );
			/*Nodal function 27*/
			dbasis[NUMNODESP2xP4*0+26] = 2.*( gauss->coord1-gauss->coord2) * 4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP2xP4*1+26] = -2.* SQRT3/3.*(gauss->coord1+ gauss->coord2 )*4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+26] = 4. *gauss->coord1 *gauss->coord2 *4.*( 4.* zeta*zeta*zeta - 5./2.*zeta ); 
			/*Nodal function 28*/
			dbasis[NUMNODESP2xP4*0+27] = 2.* gauss->coord3* (-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP2xP4*1+27] = (4.*SQRT3/3.*gauss->coord2- 2.*SQRT3/3.*gauss->coord3) * (-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+27] = 4.* gauss->coord2*gauss->coord3* (-8./3.)*((2.*zeta -1. )*(zeta+0.5)*(zeta +1.) +zeta*(zeta -1.)*( 2.*zeta + 3./2.)); 
			/*Nodal function 29*/
			dbasis[NUMNODESP2xP4*0+28] = -2.*gauss->coord3 *(-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP2xP4*1+28] = (4.*SQRT3/3.*gauss->coord1- 2.*SQRT3/3.*gauss->coord3) * (-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+28] = 4.* gauss->coord1*gauss->coord3* (-8./3.)*((2.*zeta -1. )*(zeta+0.5)*(zeta +1.) +zeta*(zeta -1.)*( 2.*zeta + 3./2.)); 
			/*Nodal function 30*/
			dbasis[NUMNODESP2xP4*0+29] = 2.*(gauss->coord1-gauss->coord2)* (-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP2xP4*1+29] = -2.*SQRT3/3.*(gauss->coord1+gauss->coord2) * (-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP2xP4*2+29] = 4.* gauss->coord1*gauss->coord2 * (-8./3.)*((2.*zeta -1. )*(zeta+0.5)*(zeta +1) +zeta*(zeta -1.)*( 2.*zeta + 3./2.)); 
			return;
		case P1xP3Enum :
			/*Nodal function 1*/
			dbasis[NUMNODESP1xP3*0+0 ] =  (9./32.)*(zeta-1)*(zeta-1./3.)*(zeta+1./3.);
			dbasis[NUMNODESP1xP3*1+0 ] = ((3.*SQRT3)/32.)*(zeta-1)*(zeta-1./3.)*(zeta+1./3.);
			dbasis[NUMNODESP1xP3*2+0 ] =- (9./16.)* gauss->coord1 *( 2. *zeta *( zeta -1. ) + ( zeta - (1./3.) )*( zeta + (1./3.) ));
			/*Nodal function 2*/
			dbasis[NUMNODESP1xP3*0+1 ] = - (9./32.)*(zeta-1)*(zeta-1./3.)*(zeta+1./3.);
			dbasis[NUMNODESP1xP3*1+1 ] = ((3.*SQRT3)/32.) *(zeta-1)*(zeta-1./3.)*(zeta+1./3.);
			dbasis[NUMNODESP1xP3*2+1 ] =- (9./16.)*gauss->coord2 *( 2.* zeta* ( zeta -1. ) + ( zeta - (1./3.) )*( zeta + (1./3.) ));
			/*Nodal function 3*/
			dbasis[NUMNODESP1xP3*0+2 ] =  0.;
			dbasis[NUMNODESP1xP3*1+2 ] = - ((3.*SQRT3)/16.)*(zeta-1)*(zeta-1./3.)*(zeta+1./3.);
			dbasis[NUMNODESP1xP3*2+2 ] = - (9./16.)* gauss->coord3* ( 2. *zeta *( zeta -1. ) + ( zeta - (1./3.) )*( zeta + (1./3.) ));
			/*Nodal function 4*/	 
			dbasis[NUMNODESP1xP3*0+3 ] = -  (9./32.)*(zeta-1./3.)*(zeta+1./3.)*(zeta+1.);
			dbasis[NUMNODESP1xP3*1+3 ] =  -((3.*SQRT3)/32.) *(zeta-1./3.)*(zeta+1./3.)*(zeta+1.);
			dbasis[NUMNODESP1xP3*2+3 ] = (9./16.)* gauss->coord1*( 2.* zeta *( zeta +1. ) + ( zeta - (1./3.) )*( zeta + (1./3.) ));
			/*Nodal function 5*/	
			dbasis[NUMNODESP1xP3*0+4 ] =   (9./32.)* (zeta-1./3.)*(zeta+1./3.)*(zeta+1.);
			dbasis[NUMNODESP1xP3*1+4 ] = - ((3.*SQRT3)/32.) *(zeta-1./3.)*(zeta+1./3.)*(zeta+1.);
			dbasis[NUMNODESP1xP3*2+4 ] = (9./16.)* gauss->coord2* ( 2.* zeta *( zeta +1. ) + ( zeta - (1./3.) )*( zeta + (1./3.) ));
			/*Nodal function 6*/	
			dbasis[NUMNODESP1xP3*0+5 ] =  0.;
			dbasis[NUMNODESP1xP3*1+5 ] =  ((3.*SQRT3)/16.)  *(zeta-1./3.)*(zeta+1./3.)*(zeta+1.);
			dbasis[NUMNODESP1xP3*2+5 ] =  (9./16.)* gauss->coord3 *( 2.* zeta * ( zeta  + 1. ) + ( zeta - (1./3.) )*( zeta + (1./3.) ));
			/*Nodal function 7*/	
			dbasis[NUMNODESP1xP3*0+6 ] = -  (27./32.) *(zeta-1)*(zeta-1./3.)*(zeta+1.);
			dbasis[NUMNODESP1xP3*1+6 ] = -  (9.*SQRT3/32.) *(zeta-1)*(zeta-1./3.)*(zeta+1.);
			dbasis[NUMNODESP1xP3*2+6 ] =  gauss->coord1*(27./16.)*( 2.* zeta *( zeta - (1./3.)) + ( zeta - 1. )*( zeta + 1. ));
			/*Nodal function 8*/	
			dbasis[NUMNODESP1xP3*0+7 ] =  (27./32.) *(zeta-1)*(zeta-1./3.)*(zeta+1.);
			dbasis[NUMNODESP1xP3*1+7 ] = -((9.*SQRT3)/32.) *(zeta-1)*(zeta-1./3.)*(zeta+1.);
			dbasis[NUMNODESP1xP3*2+7 ] =  gauss->coord2*(27./16.)*( 2.* zeta *( zeta - (1./3.)) + ( zeta - 1. )*( zeta + 1. ));
			/*Nodal function 9*/	
			dbasis[NUMNODESP1xP3*0+8 ] = 0.;
			dbasis[NUMNODESP1xP3*1+8 ] =  ((9.*SQRT3)/16.) *(zeta-1.)*(zeta-1./3.)*(zeta+1.);
			dbasis[NUMNODESP1xP3*2+8 ] =  gauss->coord3*(27./16.)*( 2. *zeta *( zeta - (1./3.)) + ( zeta - 1. )*( zeta + 1. ));
			/*Nodal function 10*/	
			dbasis[NUMNODESP1xP3*0+9 ] = (27./32.) *(zeta-1.)*(zeta+1./3.)*(zeta+1.);
			dbasis[NUMNODESP1xP3*1+9 ] = ((9.*SQRT3)/32.) *(zeta-1.)*(zeta+1./3.)*(zeta+1.);
			dbasis[NUMNODESP1xP3*2+9 ] =  -gauss->coord1 *(27./16.)*( 2* zeta *( zeta + (1./3.)) + ( zeta - 1. )*( zeta + 1. ));
			/*Nodal function 11*/	
			dbasis[NUMNODESP1xP3*0+10] = - (27./32.) *(zeta-1)*(zeta+1./3.)*(zeta+1);
			dbasis[NUMNODESP1xP3*1+10] = ((9.*SQRT3)/32.)  *(zeta-1.)*(zeta+1./3.)*(zeta+1);
			dbasis[NUMNODESP1xP3*2+10] = -gauss->coord2 *(27./16.) *( 2.* zeta *( zeta + (1./3.)) + ( zeta - 1. )*( zeta + 1. ));
			/*Nodal function 12*/	
			dbasis[NUMNODESP1xP3*0+11] = 0.;
			dbasis[NUMNODESP1xP3*1+11] = -((9.*SQRT3)/16.) *(zeta-1.)*(zeta+1./3.)*(zeta+1);
			dbasis[NUMNODESP1xP3*2+11] = -gauss->coord3 *(27./16.)*( 2.* zeta *( zeta + (1./3.)) + ( zeta - 1. )*( zeta + 1. ));
			return;
		case P1xP4Enum :
			/*Nodal function 1*/
			dbasis[NUMNODESP1xP4*0+0 ] = -0.5*(2./3.)     *(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5);
			dbasis[NUMNODESP1xP4*1+0 ] = -SQRT3/6.*(2./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5); 
			dbasis[NUMNODESP1xP4*2+0 ] =  gauss->coord1 * 2./3.*( (2.*zeta-1)*(zeta -0.5)*(zeta +0.5) + 2.* zeta *zeta *(zeta -1.)); 
			/*Nodal function 2*/
			dbasis[NUMNODESP1xP4*0+1 ] = +0.5*(2./3.)     *(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5);
			dbasis[NUMNODESP1xP4*1+1 ] = -SQRT3/6.*(2./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta+0.5);
			dbasis[NUMNODESP1xP4*2+1 ] = gauss->coord2* 2./3.* ((2.*zeta-1.)*(zeta -0.5)*(zeta +0.5) + 2. * zeta *zeta*(zeta -1.)); 
			/*Nodal function 3*/
			dbasis[NUMNODESP1xP4*0+2 ] = 0. ; 
			dbasis[NUMNODESP1xP4*1+2 ] = SQRT3/3.*(2./3.)*(zeta -1.)*(zeta-0.5)*(zeta)*(zeta+0.5); 
			dbasis[NUMNODESP1xP4*2+2 ] = gauss->coord3* 2./3.*( (2.*zeta-1.)*(zeta -0.5)*(zeta +0.5) + 2.*zeta *zeta*(zeta -1.)); 
			/*Nodal function 4*/
			dbasis[NUMNODESP1xP4*0+3 ] = -0.5 *(2./3.) *(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP1xP4*1+3 ] = -SQRT3/6.*(2./3.) *(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP1xP4*2+3 ] = gauss->coord1* 2./3.*( (2.*zeta+1.)*(zeta -0.5)*(zeta +0.5) + 2.*zeta *zeta*(zeta +1.)); 
			/*Nodal function 5*/
			dbasis[NUMNODESP1xP4*0+4 ] = +0.5 *    (2./3.) *(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP1xP4*1+4 ] = -SQRT3/6.*(2./3.)*(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. );
			dbasis[NUMNODESP1xP4*2+4 ] = gauss->coord2 * 2./3.*( (2.*zeta+1.)*(zeta -0.5)*(zeta +0.5) + 2.*zeta *zeta*(zeta +1.)); 
			/*Nodal function 6*/
			dbasis[NUMNODESP1xP4*0+5 ] = 0. ; 
			dbasis[NUMNODESP1xP4*1+5 ] = SQRT3/3.*(2./3.) *(zeta - 0.5)*(zeta)*(zeta+0.5)*(zeta +1. ); 
			dbasis[NUMNODESP1xP4*2+5 ] = gauss->coord3 * 2./3.*( (2.*zeta+1.)*(zeta -0.5)*(zeta +0.5) + 2.*zeta *zeta*(zeta +1)); 

			/*Nodal function 7*/
			dbasis[NUMNODESP1xP4*0+6 ] = -0.5 * 4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP1xP4*1+6 ] = -SQRT3/6.* 4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1.)  ; 
			dbasis[NUMNODESP1xP4*2+6 ] = gauss->coord1* 4.*( 4.*zeta *zeta*zeta - (5./2.)*zeta ); 
			/*Nodal function 8*/
			dbasis[NUMNODESP1xP4*0+7 ] = +0.5* 4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP1xP4*1+7 ] = -SQRT3/6.* 4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. )  ; 
			dbasis[NUMNODESP1xP4*2+7 ] = gauss->coord2* 4.*( 4.*zeta *zeta*zeta - (5./2.)*zeta ); 
			/*Nodal function 9*/
			dbasis[NUMNODESP1xP4*0+8 ] = 0. ; 
			dbasis[NUMNODESP1xP4*1+8 ] = SQRT3/3. * 4.*(zeta - 1.)*(zeta - 0.5)*(zeta+0.5)*(zeta +1. )  ; 
			dbasis[NUMNODESP1xP4*2+8 ] = gauss->coord3* 4.*( 4.*zeta *zeta*zeta - (5./2.)*zeta ); 

			/*Nodal function 10*/
			dbasis[NUMNODESP1xP4*0+9 ] = -0.5* (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP1xP4*1+9 ] = -SQRT3/6.* (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP1xP4*2+9 ] = gauss->coord1* (-8./3.)*((2.*zeta -1.)*(zeta-0.5)*(zeta +1.) +zeta*(zeta -1.)*( 2.*zeta + 0.5)); 
			/*Nodal function 11*/
			dbasis[NUMNODESP1xP4*0+10] = +0.5* (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP1xP4*1+10] = -SQRT3/6.* (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP1xP4*2+10] = gauss->coord2* (-8./3.)*((2.*zeta -1.)*(zeta-0.5)*(zeta +1.) +zeta *(zeta -1.)*( 2.*zeta + 0.5)); 
			/*Nodal function 12*/
			dbasis[NUMNODESP1xP4*0+11] = 0. ; 
			dbasis[NUMNODESP1xP4*1+11] = SQRT3/3.* (-8./3.)*(zeta - 1.)*(zeta - 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP1xP4*2+11] = gauss->coord3* (-8./3.)*((2.*zeta-1.)*(zeta-0.5)*(zeta +1.) +zeta *(zeta -1.)*( 2.*zeta + 0.5));
			/*Nodal function 13*/
			dbasis[NUMNODESP1xP4*0+12] = -0.5* (-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP1xP4*1+12] = -SQRT3/6.* (-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP1xP4*2+12] = gauss->coord1* (-8./3.)*((2.*zeta -1. ) *(zeta+0.5)* (zeta +1.) +  zeta* (zeta -1.)*( 2.*zeta + 3./2.));
			/*Nodal function 14*/
			dbasis[NUMNODESP1xP4*0+13] = +0.5* (-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP1xP4*1+13] = -SQRT3/6.* (-8./3.)*(zeta - 1.)*(zeta)*(zeta+0.5)*(zeta +1. ) ; 
			dbasis[NUMNODESP1xP4*2+13] = gauss->coord2* (-8./3.)*((2.*zeta -1. )*(zeta+0.5)*(zeta +1.) +  zeta*(zeta -1.)*( 2.*zeta + 3./2.)); 
			/*Nodal function 15*/
			dbasis[NUMNODESP1xP4*0+14] = 0 ; 
			dbasis[NUMNODESP1xP4*1+14] = SQRT3/3.* (-8./3.)*(zeta - 1.)*(zeta + 0.5)*(zeta)*(zeta +1. ) ; 
			dbasis[NUMNODESP1xP4*2+14] = gauss->coord3* (-8./3.)*((2. *zeta -1. )*(zeta+0.5)*(zeta +1.) +  zeta*(zeta -1.)*( 2.*zeta + 3./2.));

			return;
		#endif
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

}
/*}}}*/
void PentaRef::GetQuadJacobianDeterminant(IssmDouble* Jdet,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*This routine returns the values of the nodal functions  at the gaussian point.*/

	IssmDouble x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;

	x1=xyz_list[0*3+0];
	y1=xyz_list[0*3+1];
	z1=xyz_list[0*3+2];
	x2=xyz_list[1*3+0];
	y2=xyz_list[1*3+1];
	z2=xyz_list[1*3+2];
	x3=xyz_list[2*3+0];
	y3=xyz_list[2*3+1];
	z3=xyz_list[2*3+2];
	x4=xyz_list[3*3+0];
	y4=xyz_list[3*3+1];
	z4=xyz_list[3*3+2];

	/*Jdet = (Area of the trapezoid)/(Area trapezoid ref) with AreaRef = 4*/
	/*Area of a trabezoid = altitude * (base1 + base2)/2 */
	*Jdet= sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)) * (z4-z1 + z3-z2)/8.;
	if(*Jdet<0.) _error_("negative jacobian determinant!");

}
/*}}}*/
void PentaRef::GetSegmentJacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*The Jacobian determinant is constant over the element, discard the gaussian points. 
	 * J is assumed to have been allocated of size 2x2.*/

	IssmDouble x1=xyz_list[3*0+0];
	IssmDouble y1=xyz_list[3*0+1];
	IssmDouble z1=xyz_list[3*0+2];
	IssmDouble x2=xyz_list[3*1+0];
	IssmDouble y2=xyz_list[3*1+1];
	IssmDouble z2=xyz_list[3*1+2];

	*Jdet=.5*sqrt(pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
	if(*Jdet<0) _error_("negative jacobian determinant!");

}
/*}}}*/
void PentaRef::GetTriaJacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*The Jacobian determinant is constant over the element, discard the gaussian points. 
	 * J is assumed to have been allocated of size 2x2.*/

	IssmDouble x1=xyz_list[3*0+0];
	IssmDouble y1=xyz_list[3*0+1];
	IssmDouble z1=xyz_list[3*0+2];
	IssmDouble x2=xyz_list[3*1+0];
	IssmDouble y2=xyz_list[3*1+1];
	IssmDouble z2=xyz_list[3*1+2];
	IssmDouble x3=xyz_list[3*2+0];
	IssmDouble y3=xyz_list[3*2+1];
	IssmDouble z3=xyz_list[3*2+2];

	/*Jdet = norm( AB ^ AC ) / (2 * area of the reference triangle), with areaRef=sqrt(3) */
	*Jdet=SQRT3/6.*pow(pow(((y2-y1)*(z3-z1)-(z2-z1)*(y3-y1)),2)+pow(((z2-z1)*(x3-x1)-(x2-x1)*(z3-z1)),2)+pow(((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1)),2),0.5);
	if(*Jdet<0) _error_("negative jacobian determinant!");
}
/*}}}*/
void PentaRef::VerticalSegmentIndicesBase(int** pindices,int* pnumseg,int finiteelement){/*{{{*/

	/*Output*/
	int  numseg;
	int* indices = NULL;

	switch(finiteelement){
		case P1Enum: case P1DGEnum:
			numseg = 3;
			indices = xNew<int>(numseg*2);
			indices[0*2 + 0] = 0; indices[0*2 + 1] = 3;
			indices[1*2 + 0] = 1; indices[1*2 + 1] = 4;
			indices[2*2 + 0] = 2; indices[2*2 + 1] = 5;
			break;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			numseg = 3;
			indices = xNew<int>(numseg*2);
			indices[0*2 + 0] = 0; indices[0*2 + 1] = 3;
			indices[1*2 + 0] = 1; indices[1*2 + 1] = 4;
			indices[2*2 + 0] = 2; indices[2*2 + 1] = 5;
			break;
		case P2xP1Enum:
			numseg = 6;
			indices = xNew<int>(numseg*2);
			indices[0*2 + 0] = 0; indices[0*2 + 1] = 3;
			indices[1*2 + 0] = 1; indices[1*2 + 1] = 4;
			indices[2*2 + 0] = 2; indices[2*2 + 1] = 5;
			indices[3*2 + 0] = 6; indices[3*2 + 1] = 9;
			indices[4*2 + 0] = 7; indices[4*2 + 1] = 10;
			indices[5*2 + 0] = 8; indices[5*2 + 1] = 11;
			break;
		case P1xP2Enum:
			numseg = 3;
			indices = xNew<int>(numseg*2);
			indices[0*2 + 0] = 0; indices[0*2 + 1] = 6;
			indices[1*2 + 0] = 1; indices[1*2 + 1] = 7;
			indices[2*2 + 0] = 2; indices[2*2 + 1] = 8;
			break;
		case P1xP3Enum:
			numseg = 3;
			indices = xNew<int>(numseg*2);
			indices[0*2 + 0] = 0; indices[0*2 + 1] = 6;
			indices[1*2 + 0] = 1; indices[1*2 + 1] = 7;
			indices[2*2 + 0] = 2; indices[2*2 + 1] = 8;
			break;
		case P1xP4Enum:
			numseg = 3;
			indices = xNew<int>(numseg*2);
			indices[0*2 + 0] = 0; indices[0*2 + 1] = 6;
			indices[1*2 + 0] = 1; indices[1*2 + 1] = 7;
			indices[2*2 + 0] = 2; indices[2*2 + 1] = 8;
			break;
		case P2Enum:
			numseg = 6;
			indices = xNew<int>(numseg*2);
			indices[0*2 + 0] = 0;  indices[0*2 + 1] = 6;
			indices[1*2 + 0] = 1;  indices[1*2 + 1] = 7;
			indices[2*2 + 0] = 2;  indices[2*2 + 1] = 8;
			indices[3*2 + 0] = 9;  indices[3*2 + 1] = 15;
			indices[4*2 + 0] = 10; indices[4*2 + 1] = 16;
			indices[5*2 + 0] = 11; indices[5*2 + 1] = 17;
			break;
		case P2bubbleEnum:
			numseg = 6;
			indices = xNew<int>(numseg*2);
			indices[0*2 + 0] = 0;  indices[0*2 + 1] = 6;
			indices[1*2 + 0] = 1;  indices[1*2 + 1] = 7;
			indices[2*2 + 0] = 2;  indices[2*2 + 1] = 8;
			indices[3*2 + 0] = 9;  indices[3*2 + 1] = 15;
			indices[4*2 + 0] = 10; indices[4*2 + 1] = 16;
			indices[5*2 + 0] = 11; indices[5*2 + 1] = 17;
			break;
		case P2xP4Enum:
			numseg = 6;
			indices = xNew<int>(numseg*2);
			indices[0*2 + 0] = 0;  indices[0*2 + 1] = 6;
			indices[1*2 + 0] = 1;  indices[1*2 + 1] = 7;
			indices[2*2 + 0] = 2;  indices[2*2 + 1] = 8;
			indices[3*2 + 0] = 9;  indices[3*2 + 1] = 15;
			indices[4*2 + 0] = 10; indices[4*2 + 1] = 16;
			indices[5*2 + 0] = 11; indices[5*2 + 1] = 17;
			break;
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

	/*Assign output pointer*/
	*pnumseg   = numseg;
	*pindices  = indices;
}
/*}}}*/
int  PentaRef::NumberofNodes(int finiteelement){/*{{{*/

	switch(finiteelement){
		case NoneEnum:              return 0;
		case P0Enum:                return NUMNODESP0;
		case P1Enum:                return NUMNODESP1;
		case P1DGEnum:              return NUMNODESP1;
		case P1bubbleEnum:          return NUMNODESP1b;
		case P1bubblecondensedEnum: return NUMNODESP1b;
		case P2Enum:                return NUMNODESP2;
		case P2bubbleEnum:          return NUMNODESP2b;
		case P2bubblecondensedEnum: return NUMNODESP2b;
		case P2xP1Enum:             return NUMNODESP2xP1;
		case P1xP2Enum:             return NUMNODESP1xP2;
		case P2xP4Enum:             return NUMNODESP2xP4;
		case P1xP3Enum:             return NUMNODESP1xP3;
		case P1xP4Enum:             return NUMNODESP1xP4;
		case P1P1Enum:              return NUMNODESP1*2;
		case P1P1GLSEnum:           return NUMNODESP1*2;
		case MINIcondensedEnum:     return NUMNODESP1b+NUMNODESP1;
		case MINIEnum:              return NUMNODESP1b+NUMNODESP1;
		case TaylorHoodEnum:        return NUMNODESP2+NUMNODESP1;
		case LATaylorHoodEnum:      return NUMNODESP2;
		case OneLayerP4zEnum:       return NUMNODESP2xP4+NUMNODESP1;
		case CrouzeixRaviartEnum:   return NUMNODESP2b+NUMNODESP1;
		case LACrouzeixRaviartEnum: return NUMNODESP2b;
		default:       _error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

	return -1;
}
/*}}}*/
int  PentaRef::PressureInterpolation(int fe_stokes){/*{{{*/

	switch(fe_stokes){
		case P1P1Enum:              return P1Enum;
		case P1P1GLSEnum:           return P1Enum;
		case MINIcondensedEnum:     return P1Enum;
		case MINIEnum:              return P1Enum;
		case TaylorHoodEnum:        return P1Enum;
		case LATaylorHoodEnum:      return NoneEnum;
		case OneLayerP4zEnum:       return P1Enum;
		case CrouzeixRaviartEnum:   return P1DGEnum;
		case LACrouzeixRaviartEnum: return NoneEnum;
		default:       _error_("Element type "<<EnumToStringx(fe_stokes)<<" not supported yet");
	}

	return -1;
}
/*}}}*/
void PentaRef::SurfaceNodeIndices(int* pnumindices,int** pindices,int finiteelement){/*{{{*/

	/*Output*/
	int  numindices;
	int* indices = NULL;

	switch(finiteelement){
		case P1Enum: case P1DGEnum:
			numindices = 3;
			indices    = xNew<int>(numindices);
			indices[0] = 3;
			indices[1] = 4;
			indices[2] = 5;
			break;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			numindices = 3;
			indices    = xNew<int>(numindices);
			indices[0] = 3;
			indices[1] = 4;
			indices[2] = 5;
			break;
		case P2xP1Enum:
			numindices = 6;
			indices    = xNew<int>(numindices);
			indices[0] = 3;
			indices[1] = 4;
			indices[2] = 5;
			indices[3] = 9;
			indices[4] = 10;
			indices[5] = 11;
			break;
		case P1xP2Enum:
		case P1xP3Enum:
		case P1xP4Enum:
			numindices = 3;
			indices    = xNew<int>(numindices);
			indices[0] = 3;
			indices[1] = 4;
			indices[2] = 5;
			break;
		case P2Enum:
			numindices = 6;
			indices    = xNew<int>(numindices);
			indices[0] = 3;
			indices[1] = 4;
			indices[2] = 5;
			indices[3] = 12;
			indices[4] = 13;
			indices[5] = 14;
			break;
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

	/*Assign output pointer*/
	*pnumindices = numindices;
	*pindices    = indices;
}
/*}}}*/
int  PentaRef::TensorInterpolation(int fe_stokes){/*{{{*/

	switch(fe_stokes){
		case XTaylorHoodEnum:    return P1DGEnum;
		default: _error_("Element type "<<EnumToStringx(fe_stokes)<<" not supported yet");
	}

	return -1;
}
/*}}}*/
int  PentaRef::VelocityInterpolation(int fe_stokes){/*{{{*/

	switch(fe_stokes){
		case P1P1Enum:              return P1Enum;
		case P1P1GLSEnum:           return P1Enum;
		case MINIcondensedEnum:     return P1bubbleEnum;
		case MINIEnum:              return P1bubbleEnum;
		case TaylorHoodEnum:        return P2Enum;
		case LATaylorHoodEnum:      return P2Enum;
		case OneLayerP4zEnum:       return P2xP4Enum;
		case CrouzeixRaviartEnum:   return P2bubbleEnum;
		case LACrouzeixRaviartEnum: return P2bubbleEnum;
		default:       _error_("Element type "<<EnumToStringx(fe_stokes)<<" not supported yet");
	}

	return -1;
}
/*}}}*/
