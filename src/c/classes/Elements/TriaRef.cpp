/*!\file TriaRef.c
 * \brief: implementation of the TriaRef object
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
#define NUMNODESP1  3
#define NUMNODESP1b 4
#define NUMNODESP2  6
#define NUMNODESP2b 7
#define NUMNODESMAX 7

/*Object constructors and destructor*/
TriaRef::TriaRef(){/*{{{*/
}
/*}}}*/
TriaRef::~TriaRef(){/*{{{*/
}
/*}}}*/

/*Reference Element numerics*/
void TriaRef::GetInputDerivativeValue(IssmDouble* p, IssmDouble* plist,IssmDouble* xyz_list, Gauss* gauss,int finiteelement){/*{{{*/
	/*From node values of parameter p (plist[0],plist[1],plist[2]), return parameter derivative value at gaussian
	 * point specified by gauss_basis:
	 *   dp/dx=plist[0]*dh1/dx+plist[1]*dh2/dx+plist[2]*dh3/dx
	 *   dp/dx=plist[0]*dh1/dx+plist[1]*dh2/dx+plist[2]*dh3/dx
	 *
	 * p is a vector already allocated.
	 *
	 * WARNING: For a significant gain in performance, it is better to use
	 * static memory allocation instead of dynamic.
	 */

	/*Allocate derivatives of basis functions*/
	IssmDouble  dbasis[2*NUMNODESMAX];

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->NumberofNodes(finiteelement);
	_assert_(numnodes<=NUMNODESMAX);

	/*Get basis functions derivatives at this point*/
	GetNodalFunctionsDerivatives(&dbasis[0],xyz_list,gauss,finiteelement);

	/*Calculate parameter for this Gauss point*/
	IssmDouble dpx=0.;
	IssmDouble dpy=0.;
	for(int i=0;i<numnodes;i++) dpx += dbasis[0*numnodes+i]*plist[i];
	for(int i=0;i<numnodes;i++) dpy += dbasis[1*numnodes+i]*plist[i];

	/*Assign values*/
	p[0]=dpx;
	p[1]=dpy;

}
/*}}}*/
void TriaRef::GetInputValue(IssmDouble* p, IssmDouble* plist, Gauss* gauss,int finiteelement){/*{{{*/
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
void TriaRef::GetJacobian(IssmDouble* J, IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*The Jacobian is constant over the element, discard the gaussian points.
	 * J is assumed to have been allocated of size 2x2.*/

	IssmDouble x1 = xyz_list[3*0+0];
	IssmDouble y1 = xyz_list[3*0+1];
	IssmDouble x2 = xyz_list[3*1+0];
	IssmDouble y2 = xyz_list[3*1+1];
	IssmDouble x3 = xyz_list[3*2+0];
	IssmDouble y3 = xyz_list[3*2+1];

	J[2*0+0] = 0.5*(x2-x1);
	J[2*1+0] = SQRT3/6.0*(2*x3-x1-x2);
	J[2*0+1] = 0.5*(y2-y1);
	J[2*1+1] = SQRT3/6.0*(2*y3-y1-y2);
}
/*}}}*/
void TriaRef::GetJacobianDeterminant(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*The Jacobian determinant is constant over the element, discard the gaussian points.
	 * J is assumed to have been allocated of size 2x2.*/
	IssmDouble J[2][2];

	/*Get Jacobian*/
	GetJacobian(&J[0][0],xyz_list,gauss);

	/*Get Determinant*/
	Matrix2x2Determinant(Jdet,&J[0][0]);
	if(*Jdet<0) _error_("negative jacobian determinant!");

}
/*}}}*/
void TriaRef::GetJacobianDeterminant3D(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*The Jacobian determinant is constant over the element, discard the gaussian points.
	 * J is assumed to have been allocated of size 2x2.*/
	IssmDouble J[2][2];

	/*Get Jacobian*/
	GetJacobian(&J[0][0],xyz_list,gauss);

	/*Get Determinant*/
	Matrix2x2Determinant(Jdet,&J[0][0]);

}
/*}}}*/
void TriaRef::GetJacobianInvert(IssmDouble*  Jinv, IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	/*Jacobian*/
	IssmDouble J[2][2];

	/*Call Jacobian routine to get the jacobian:*/
	GetJacobian(&J[0][0], xyz_list, gauss);

	/*Invert Jacobian matrix: */
	Matrix2x2Invert(Jinv,&J[0][0]);

}
/*}}}*/
void TriaRef::GetNodalFunctions(IssmDouble* basis,Gauss* gauss_in,int finiteelement){/*{{{*/
	/*This routine returns the values of the nodal functions  at the gaussian point.*/

	_assert_(basis);

	/*Cast gauss to GaussTria*/
	_assert_(gauss_in->Enum()==GaussTriaEnum);
	GaussTria* gauss = xDynamicCast<GaussTria*>(gauss_in);

	switch(finiteelement){
		case NoneEnum:
			return;
		case P0Enum: case P0DGEnum:
			basis[0]=1.;
			return;
		case P1Enum: case P1DGEnum:
			basis[0]=gauss->coord1;
			basis[1]=gauss->coord2;
			basis[2]=gauss->coord3;
			return;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			/*Corner nodes*/
			basis[0]=gauss->coord1;
			basis[1]=gauss->coord2;
			basis[2]=gauss->coord3;
			/*bubble*/
			basis[3]=27.*gauss->coord1*gauss->coord2*gauss->coord3;
			return;
		case P2Enum:
			/*Corner nodes*/
			basis[0]=gauss->coord1*(2.*gauss->coord1-1.);
			basis[1]=gauss->coord2*(2.*gauss->coord2-1.);
			basis[2]=gauss->coord3*(2.*gauss->coord3-1.);
			/*Mid-sides*/
			basis[3]=4.*gauss->coord3*gauss->coord2;
			basis[4]=4.*gauss->coord3*gauss->coord1;
			basis[5]=4.*gauss->coord1*gauss->coord2;
			return;
		case P2bubbleEnum: case P2bubblecondensedEnum:
			/*Corner nodes*/
			basis[0]=gauss->coord1*(2.*gauss->coord1-1.);
			basis[1]=gauss->coord2*(2.*gauss->coord2-1.);
			basis[2]=gauss->coord3*(2.*gauss->coord3-1.);
			/*Mid-sides*/
			basis[3]=4.*gauss->coord3*gauss->coord2;
			basis[4]=4.*gauss->coord3*gauss->coord1;
			basis[5]=4.*gauss->coord1*gauss->coord2;
			/*bubble*/
			basis[6]=27.*gauss->coord1*gauss->coord2*gauss->coord3;
			return;
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}
}
/*}}}*/
void TriaRef::GetNodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list, Gauss* gauss,int finiteelement){/*{{{*/

	/*This routine returns the values of the nodal functions derivatives  (with respect to the
	 * actual coordinate system): */
	IssmDouble    Jinv[2][2];

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->NumberofNodes(finiteelement);

	/*Get nodal functions derivatives in reference triangle*/
	IssmDouble dbasis_ref[2*NUMNODESMAX];
	GetNodalFunctionsDerivativesReference(dbasis_ref,gauss,finiteelement);

	/*Get Jacobian invert: */
	GetJacobianInvert(&Jinv[0][0], xyz_list, gauss);

	/*Build dbasis:
	 * [dhi/dx]= Jinv*[dhi/dr]
	 * [dhi/dy]       [dhi/ds]
	 */
	for(int i=0;i<numnodes;i++){
		dbasis[numnodes*0+i] = Jinv[0][0]*dbasis_ref[0*numnodes+i]+Jinv[0][1]*dbasis_ref[1*numnodes+i];
		dbasis[numnodes*1+i] = Jinv[1][0]*dbasis_ref[0*numnodes+i]+Jinv[1][1]*dbasis_ref[1*numnodes+i];
	}

}
/*}}}*/
void TriaRef::GetSegmentJacobianDeterminant(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	/*The Jacobian determinant is constant over the element, discard the gaussian points.
	 * J is assumed to have been allocated*/

	IssmDouble x1 = xyz_list[3*0+0];
	IssmDouble y1 = xyz_list[3*0+1];
	IssmDouble x2 = xyz_list[3*1+0];
	IssmDouble y2 = xyz_list[3*1+1];

	*Jdet = .5*sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
	if(*Jdet<0) _error_("negative jacobian determinant!");

}
/*}}}*/
void TriaRef::GetSegmentNodalFunctions(IssmDouble* basis,Gauss* gauss,int index1,int index2,int finiteelement){/*{{{*/
	/*This routine returns the values of the nodal functions  at the gaussian point.*/

	_assert_(index1>=0 && index1<3);
	_assert_(index2>=0 && index2<3);

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->NumberofNodes(finiteelement);

	/*Get nodal functions*/
	IssmDouble* triabasis=xNew<IssmDouble>(numnodes);
	GetNodalFunctions(triabasis,gauss,finiteelement);

	switch(finiteelement){
		case P0Enum: case P0DGEnum:
			basis[0]=triabasis[0];
			xDelete<IssmDouble>(triabasis);
			return;
		case P1Enum: case P1DGEnum:
			basis[0]=triabasis[index1];
			basis[1]=triabasis[index2];
			xDelete<IssmDouble>(triabasis);
			return;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			basis[0]=triabasis[index1];
			basis[1]=triabasis[index2];
			xDelete<IssmDouble>(triabasis);
			return;
		case P2Enum:
			_assert_(index2<index1);
			basis[0]=triabasis[index1];
			basis[1]=triabasis[index2];
			basis[2]=triabasis[3+index2-1];
			xDelete<IssmDouble>(triabasis);
			return;
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

	/*Clean up*/
	xDelete<IssmDouble>(triabasis);
}
/*}}}*/
void TriaRef::GetSegmentNodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list_tria,Gauss* gauss,int index1,int index2,int finiteelement){/*{{{*/
	/*This routine returns the values of the nodal functions  at the gaussian point.*/

	_assert_(index1>=0 && index1<3);
	_assert_(index2>=0 && index2<3);

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->NumberofNodes(finiteelement);

	/*Get nodal functions*/
	IssmDouble* dtriabasis=xNew<IssmDouble>(2*numnodes);
	GetNodalFunctionsDerivatives(dtriabasis,xyz_list_tria,gauss,finiteelement);

	switch(finiteelement){
		case P1Enum: case P1DGEnum:
			dbasis[2*0+0] = dtriabasis[numnodes*0+index1];
			dbasis[2*0+1] = dtriabasis[numnodes*1+index1];
			dbasis[2*1+0] = dtriabasis[numnodes*0+index2];
			dbasis[2*1+1] = dtriabasis[numnodes*1+index2];
			break;
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

	/*Clean up*/
	xDelete<IssmDouble>(dtriabasis);
}
/*}}}*/
void TriaRef::GetNodalFunctionsDerivativesReference(IssmDouble* dbasis,Gauss* gauss_in,int finiteelement){/*{{{*/
	/*This routine returns the values of the nodal functions derivatives  (with respect to the
	 * natural coordinate system) at the gaussian point. */

	_assert_(dbasis && gauss_in);

	/*Cast gauss to GaussTria*/
	_assert_(gauss_in->Enum()==GaussTriaEnum);
	GaussTria* gauss = xDynamicCast<GaussTria*>(gauss_in);

	switch(finiteelement){
		case P0Enum: case P0DGEnum:
			/*Nodal function 1*/
			dbasis[NUMNODESP0*0+0] = 0.;
			dbasis[NUMNODESP0*1+0] = 0.;
			return;
		case P1Enum: case P1DGEnum:
			/*Nodal function 1*/
			dbasis[NUMNODESP1*0+0] = -0.5;
			dbasis[NUMNODESP1*1+0] = -SQRT3/6.;
			/*Nodal function 2*/
			dbasis[NUMNODESP1*0+1] = 0.5;
			dbasis[NUMNODESP1*1+1] = -SQRT3/6.;
			/*Nodal function 3*/
			dbasis[NUMNODESP1*0+2] = 0;
			dbasis[NUMNODESP1*1+2] = SQRT3/3.;
			return;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			/*Nodal function 1*/
			dbasis[NUMNODESP1b*0+0] = -0.5;
			dbasis[NUMNODESP1b*1+0] = -SQRT3/6.;
			/*Nodal function 2*/
			dbasis[NUMNODESP1b*0+1] = 0.5;
			dbasis[NUMNODESP1b*1+1] = -SQRT3/6.;
			/*Nodal function 3*/
			dbasis[NUMNODESP1b*0+2] = 0;
			dbasis[NUMNODESP1b*1+2] = SQRT3/3.;
			/*Nodal function 4*/
			dbasis[NUMNODESP1b*0+3] = 27.*(-.5*gauss->coord2*gauss->coord3 + .5*gauss->coord1*gauss->coord3);
			dbasis[NUMNODESP1b*1+3] = 27.*SQRT3*(-1./6.*gauss->coord2*gauss->coord3 - 1./6.*gauss->coord1*gauss->coord3 +1./3.*gauss->coord1*gauss->coord2);
			return;
		case P2Enum:
			/*Nodal function 1*/
			dbasis[NUMNODESP2*0+0] = -2.*gauss->coord1 + 0.5;
			dbasis[NUMNODESP2*1+0] = -2.*SQRT3/3.*gauss->coord1 + SQRT3/6.;
			/*Nodal function 2*/
			dbasis[NUMNODESP2*0+1] = +2.*gauss->coord2 - 0.5;
			dbasis[NUMNODESP2*1+1] = -2.*SQRT3/3.*gauss->coord2 + SQRT3/6.;
			/*Nodal function 3*/
			dbasis[NUMNODESP2*0+2] = 0.;
			dbasis[NUMNODESP2*1+2] = +4.*SQRT3/3.*gauss->coord3 - SQRT3/3.;
			/*Nodal function 4*/
			dbasis[NUMNODESP2*0+3] = +2.*gauss->coord3;
			dbasis[NUMNODESP2*1+3] = +4.*SQRT3/3.*gauss->coord2 - 2.*SQRT3/3.*gauss->coord3;
			/*Nodal function 5*/
			dbasis[NUMNODESP2*0+4] = -2.*gauss->coord3;
			dbasis[NUMNODESP2*1+4] = +4.*SQRT3/3.*gauss->coord1 - 2.*SQRT3/3.*gauss->coord3;
			/*Nodal function 6*/
			dbasis[NUMNODESP2*0+5] = 2.*(gauss->coord1-gauss->coord2);
			dbasis[NUMNODESP2*1+5] = -2.*SQRT3/3.*(gauss->coord1+gauss->coord2);
			return;
		case P2bubbleEnum: case P2bubblecondensedEnum:
			/*Nodal function 1*/
			dbasis[NUMNODESP2b*0+0] = -2.*gauss->coord1 + 0.5;
			dbasis[NUMNODESP2b*1+0] = -2.*SQRT3/3.*gauss->coord1 + SQRT3/6.;
			/*Nodal function 2*/
			dbasis[NUMNODESP2b*0+1] = +2.*gauss->coord2 - 0.5;
			dbasis[NUMNODESP2b*1+1] = -2.*SQRT3/3.*gauss->coord2 + SQRT3/6.;
			/*Nodal function 3*/
			dbasis[NUMNODESP2b*0+2] = 0.;
			dbasis[NUMNODESP2b*1+2] = +4.*SQRT3/3.*gauss->coord3 - SQRT3/3.;
			/*Nodal function 4*/
			dbasis[NUMNODESP2b*0+3] = +2.*gauss->coord3;
			dbasis[NUMNODESP2b*1+3] = +4.*SQRT3/3.*gauss->coord2 - 2.*SQRT3/3.*gauss->coord3;
			/*Nodal function 5*/
			dbasis[NUMNODESP2b*0+4] = -2.*gauss->coord3;
			dbasis[NUMNODESP2b*1+4] = +4.*SQRT3/3.*gauss->coord1 - 2.*SQRT3/3.*gauss->coord3;
			/*Nodal function 6*/
			dbasis[NUMNODESP2b*0+5] = 2.*(gauss->coord1-gauss->coord2);
			dbasis[NUMNODESP2b*1+5] = -2.*SQRT3/3.*(gauss->coord1+gauss->coord2);
			/*Nodal function 7*/
			dbasis[NUMNODESP2b*0+6] = 27.*(-.5*gauss->coord2*gauss->coord3 + .5*gauss->coord1*gauss->coord3);
			dbasis[NUMNODESP2b*1+6] = 27.*SQRT3*(-1./6.*gauss->coord2*gauss->coord3 - 1./6.*gauss->coord1*gauss->coord3 +1./3.*gauss->coord1*gauss->coord2);
			return;
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

}
/*}}}*/
void TriaRef::NodeOnEdgeIndices(int* pnumindices,int** pindices,int index,int finiteelement){/*{{{*/

	/*Output*/
	int  numindices;
	int* indices = NULL;

	switch(finiteelement){
		case P1Enum: case P1DGEnum: case P1bubbleEnum: case P1bubblecondensedEnum:
			numindices = 2;
			indices    = xNew<int>(numindices);
			switch(index){
				case 0:
					indices[0] = 1;
					indices[1] = 2;
					break;
				case 1:
					indices[0] = 2;
					indices[1] = 0;
					break;
				case 2:
					indices[0] = 0;
					indices[1] = 1;
					break;
				default:
					_error_("Edge index provided ("<<index<<") is not between 0 and 2");
			}
			break;
		case P2Enum:
			numindices = 3;
			indices    = xNew<int>(numindices);
			switch(index){
				case 0:
					indices[0] = 1;
					indices[1] = 2;
					indices[2] = 3;
					break;
				case 1:
					indices[0] = 2;
					indices[1] = 0;
					indices[2] = 4;
					break;
				case 2:
					indices[0] = 0;
					indices[1] = 1;
					indices[2] = 5;
					break;
				default:
					_error_("Edge index provided ("<<index<<") is not between 0 and 2");
			}
			break;
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

	/*Assign output pointer*/
	*pnumindices = numindices;
	*pindices    = indices;
}
/*}}}*/
int  TriaRef::NumberofNodes(int finiteelement){/*{{{*/

	switch(finiteelement){
		case NoneEnum:                return 0;
		case P0Enum:                  return NUMNODESP0;
		case P0DGEnum:                return NUMNODESP0;
		case P1Enum:                  return NUMNODESP1;
		case P1DGEnum:                return NUMNODESP1;
		case P1bubbleEnum:            return NUMNODESP1b;
		case P1bubblecondensedEnum:   return NUMNODESP1b;
		case P2Enum:                  return NUMNODESP2;
		case P2bubbleEnum:            return NUMNODESP2b;
		case P2bubblecondensedEnum:   return NUMNODESP2b;
		case P1P1Enum:                return NUMNODESP1*2;
		case P1P1GLSEnum:             return NUMNODESP1*2;
		case MINIcondensedEnum:       return NUMNODESP1b+NUMNODESP1;
		case MINIEnum:                return NUMNODESP1b+NUMNODESP1;
		case TaylorHoodEnum:          return NUMNODESP2+NUMNODESP1;
		case LATaylorHoodEnum:        return NUMNODESP2;
		case XTaylorHoodEnum:         return NUMNODESP2+NUMNODESP1;
		case CrouzeixRaviartEnum:     return NUMNODESP2b+NUMNODESP1;
		case LACrouzeixRaviartEnum:   return NUMNODESP2b;
		default: _error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

	return -1;
}
/*}}}*/
int  TriaRef::PressureInterpolation(int fe_stokes){/*{{{*/

	switch(fe_stokes){
		case P1P1Enum:              return P1Enum;
		case P1P1GLSEnum:           return P1Enum;
		case MINIcondensedEnum:     return P1Enum;
		case MINIEnum:              return P1Enum;
		case TaylorHoodEnum:        return P1Enum;
		case LATaylorHoodEnum:      return NoneEnum;
		case XTaylorHoodEnum:       return P1Enum;
		case CrouzeixRaviartEnum:   return P1DGEnum;
		case LACrouzeixRaviartEnum: return NoneEnum;
		default:       _error_("Element type "<<EnumToStringx(fe_stokes)<<" not supported yet");
	}

	return -1;
}
/*}}}*/
int  TriaRef::TensorInterpolation(int fe_stokes){/*{{{*/
	/*This routine returns the values of the nodal functions  at the gaussian point.*/

	switch(fe_stokes){
		case XTaylorHoodEnum: return P1DGEnum;
		default: _error_("Element type "<<EnumToStringx(fe_stokes)<<" not supported yet");
	}
}
/*}}}*/
int  TriaRef::VelocityInterpolation(int fe_stokes){/*{{{*/

	switch(fe_stokes){
		case P1P1Enum:              return P1Enum;
		case P1P1GLSEnum:           return P1Enum;
		case MINIcondensedEnum:     return P1bubbleEnum;
		case MINIEnum:              return P1bubbleEnum;
		case TaylorHoodEnum:        return P2Enum;
		case LATaylorHoodEnum:      return P2Enum;
		case XTaylorHoodEnum:       return P2Enum;
		case CrouzeixRaviartEnum:   return P2bubbleEnum;
		case LACrouzeixRaviartEnum: return P2bubbleEnum;
		default:       _error_("Element type "<<EnumToStringx(fe_stokes)<<" not supported yet");
	}

	return -1;
}
/*}}}*/
