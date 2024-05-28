/*!\file SegRef.c
 * \brief: implementation of the SegRef object
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
#define NUMNODESP0  1
#define NUMNODESP1  2
#define NUMNODESMAX 2

/*Object constructors and destructor*/
SegRef::SegRef(){/*{{{*/
}
/*}}}*/
SegRef::~SegRef(){/*{{{*/
}
/*}}}*/

/*Reference Element numerics*/
void SegRef::GetInputDerivativeValue(IssmDouble* p, IssmDouble* plist,IssmDouble* xyz_list, GaussSeg* gauss,int finiteelement){/*{{{*/

	/*From node values of parameter p (plist[0],plist[1]), return parameter derivative value at gaussian 
	 * point specified by gauss_basis:
	 *   dp/dx=plist[0]*dh1/dx+plist[1]*dh2/dx
	 *
	 * p is a vector already allocated.
	 *
	 * WARNING: For a significant gain in performance, it is better to use
	 * static memory allocation instead of dynamic.
	 */

	/*Allocate derivatives of basis functions*/
	IssmDouble  dbasis[1*NUMNODESMAX];

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->NumberofNodes(finiteelement);
	_assert_(numnodes<=NUMNODESMAX);

	/*Get basis functions derivatives at this point*/
	GetNodalFunctionsDerivatives(&dbasis[0],xyz_list,gauss,finiteelement);

	/*Calculate parameter for this Gauss point*/
	IssmDouble dpx=0.;
	for(int i=0;i<numnodes;i++) dpx += dbasis[0*numnodes+i]*plist[i];

	/*Assign values*/
	p[0]=dpx;
}
/*}}}*/
void SegRef::GetInputValue(IssmDouble* p, IssmDouble* plist, GaussSeg* gauss,int finiteelement){/*{{{*/
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
	*p = value;
}
/*}}}*/
void SegRef::GetJacobian(IssmDouble* J, IssmDouble* xyz_list,GaussSeg* gauss){/*{{{*/
	/*The Jacobian is constant over the element, discard the gaussian points. 
	 * J is assumed to have been allocated of size 1*/

	IssmDouble x1=xyz_list[3*0+0];
	IssmDouble x2=xyz_list[3*1+0];

	*J=.5*fabs(x2-x1);
}
/*}}}*/
void SegRef::GetJacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,GaussSeg* gauss){/*{{{*/
	/*The Jacobian determinant is constant over the element, discard the gaussian points. 
	 * J is assumed to have been allocated of size 1.*/

	/*Call Jacobian routine to get the jacobian:*/
	GetJacobian(Jdet, xyz_list, gauss);
	if(*Jdet<0) _error_("negative jacobian determinant!");

}
/*}}}*/
void SegRef::GetJacobianInvert(IssmDouble* Jinv, IssmDouble* xyz_list,GaussSeg* gauss){/*{{{*/

	/*Jacobian*/
	IssmDouble J;

	/*Call Jacobian routine to get the jacobian:*/
	GetJacobian(&J, xyz_list, gauss);

	/*Invert Jacobian matrix: */
	*Jinv = 1./J;
}
/*}}}*/
void SegRef::GetNodalFunctions(IssmDouble* basis,GaussSeg* gauss,int finiteelement){/*{{{*/
	/*This routine returns the values of the nodal functions  at the gaussian point.*/

	_assert_(basis);

	switch(finiteelement){
		case P0Enum:
			basis[0]=1.;
			return;
		case P1Enum: case P1DGEnum:
			basis[0]=(1.-gauss->coord1)/2.;
			basis[1]=(1.+gauss->coord1)/2.;
			return;
		case P2Enum:
			basis[0]=(gauss->coord1-1.)*gauss->coord1/2.;
			basis[1]=gauss->coord1*(1.+gauss->coord1)/2.;
			basis[2]=(1.-gauss->coord1)*(1.+gauss->coord1);
			return;
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}
}
/*}}}*/
void SegRef::GetNodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list, GaussSeg* gauss,int finiteelement){/*{{{*/

	/*This routine returns the values of the nodal functions derivatives  (with respect to the 
	 * actual coordinate system): */
	IssmDouble    Jinv;

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->NumberofNodes(finiteelement);

	/*Get nodal functions derivatives in reference triangle*/
	IssmDouble dbasis_ref[1*NUMNODESMAX];
	GetNodalFunctionsDerivativesReference(dbasis_ref,gauss,finiteelement); 

	/*Get Jacobian invert: */
	GetJacobianInvert(&Jinv, xyz_list, gauss);

	/*Build dbasis: 
	 * [dhi/dx]= Jinv*[dhi/dr]
	 */
	for(int i=0;i<numnodes;i++){
		dbasis[i] = Jinv*dbasis_ref[i];
	}
}
/*}}}*/
void SegRef::GetNodalFunctionsDerivativesReference(IssmDouble* dbasis,GaussSeg* gauss,int finiteelement){/*{{{*/
	/*This routine returns the values of the nodal functions derivatives  (with respect to the 
	 * natural coordinate system) at the gaussian point. */

	_assert_(dbasis && gauss);

	switch(finiteelement){
		case P0Enum:
			/*Nodal function 1*/
			dbasis[0] = 0.;
			break;
		case P1Enum: case P1DGEnum:
			/*Nodal function 1*/
			dbasis[0] = -0.5;
			/*Nodal function 2*/
			dbasis[1] = 0.5;
			return;
		case P2Enum:
			/*Nodal function 1*/
			dbasis[0] = (gauss->coord1-1.)/2. + gauss->coord1/2.;
			dbasis[1] = (1.+gauss->coord1)/2. + gauss->coord1/2.;
			dbasis[2] = -2.*gauss->coord1;
			return;
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

}
/*}}}*/
int  SegRef::NumberofNodes(int finiteelement){/*{{{*/

	switch(finiteelement){
		case P0Enum:                return NUMNODESP0;
		case P1Enum:                return NUMNODESP1;
		case P1DGEnum:              return NUMNODESP1;
		default: _error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

	return -1;
}
/*}}}*/
