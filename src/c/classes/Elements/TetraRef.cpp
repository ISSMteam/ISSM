/*!\file TetraRef.c
 * \brief: implementation of the TetraRef object
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
#define NUMNODESP1  4
#define NUMNODESP1b 5
#define NUMNODESP2  10
#define NUMNODESMAX 10

/*Object constructors and destructor*/
TetraRef::TetraRef(){/*{{{*/
}
/*}}}*/
TetraRef::~TetraRef(){/*{{{*/
}
/*}}}*/

/*Reference Element numerics*/
void TetraRef::GetInputDerivativeValue(IssmDouble* p, IssmDouble* plist,IssmDouble* xyz_list, GaussTetra* gauss,int finiteelement){/*{{{*/
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
void TetraRef::GetInputValue(IssmDouble* p, IssmDouble* plist, Gauss* gauss,int finiteelement){/*{{{*/
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
void TetraRef::GetJacobian(IssmDouble* J, IssmDouble* xyz_list,GaussTetra* gauss){/*{{{*/
	/*The Jacobian is constant over the element, discard the gaussian points. 
	 * J is assumed to have been allocated of size 1*/

	IssmDouble x1=xyz_list[3*0+0];
	IssmDouble x2=xyz_list[3*1+0];
	IssmDouble x3=xyz_list[3*2+0];
	IssmDouble x4=xyz_list[3*3+0];

	IssmDouble y1=xyz_list[3*0+1];
	IssmDouble y2=xyz_list[3*1+1];
	IssmDouble y3=xyz_list[3*2+1];
	IssmDouble y4=xyz_list[3*3+1];

	IssmDouble z1=xyz_list[3*0+2];
	IssmDouble z2=xyz_list[3*1+2];
	IssmDouble z3=xyz_list[3*2+2];
	IssmDouble z4=xyz_list[3*3+2];

	J[3*0+0] = x2-x1;
	J[3*0+1] = y2-y1;
	J[3*0+2] = z2-z1;

	J[3*1+0] = x3-x1;
	J[3*1+1] = y3-y1;
	J[3*1+2] = z3-z1;

	J[3*2+0] = x4-x1;
	J[3*2+1] = y4-y1;
	J[3*2+2] = z4-z1;
}
/*}}}*/
void TetraRef::GetJacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,GaussTetra* gauss){/*{{{*/
	/*The Jacobian determinant is constant over the element, discard the gaussian points. 
	 * J is assumed to have been allocated of size 2x2.*/
	IssmDouble J[3][3];

	/*Call Jacobian routine to get the jacobian:*/
	GetJacobian(&J[0][0],xyz_list, gauss);

	/*Get Determinant*/
	Matrix3x3Determinant(Jdet,&J[0][0]);
	if(*Jdet<0) _error_("negative jacobian determinant!");

}
/*}}}*/
void TetraRef::GetJacobianDeterminantFace(IssmDouble*  Jdet, IssmDouble* xyz_list,GaussTetra* gauss){/*{{{*/
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
void TetraRef::GetJacobianInvert(IssmDouble* Jinv, IssmDouble* xyz_list,GaussTetra* gauss){/*{{{*/

	/*Jacobian*/
	IssmDouble J[3][3];

	/*Call Jacobian routine to get the jacobian:*/
	GetJacobian(&J[0][0], xyz_list, gauss);

	/*Invert Jacobian matrix: */
	Matrix3x3Invert(Jinv,&J[0][0]);
}
/*}}}*/
void TetraRef::GetNodalFunctions(IssmDouble* basis,Gauss* gauss_in,int finiteelement){/*{{{*/
	/*This routine returns the values of the nodal functions  at the gaussian point.*/

	_assert_(basis);

	/*Cast gauss to GaussTetra*/
	_assert_(gauss_in->Enum()==GaussTetraEnum);
	GaussTetra* gauss = xDynamicCast<GaussTetra*>(gauss_in);

	switch(finiteelement){
		case P0Enum:
			basis[0]=1.;
			return;
		case P1Enum: case P1DGEnum:
			basis[0]=gauss->coord1;
			basis[1]=gauss->coord2;
			basis[2]=gauss->coord3;
			basis[3]=gauss->coord4;
			return;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			/*Corner nodes*/
			basis[0]=gauss->coord1;
			basis[1]=gauss->coord2;
			basis[2]=gauss->coord3;
			basis[3]=gauss->coord4;
			/*bubble*/
			basis[4]=256.*gauss->coord1*gauss->coord2*gauss->coord3*gauss->coord4;
			return;
		case P2Enum:
			/*Vertices*/
			basis[0]=gauss->coord1*(2.*gauss->coord1-1.);
			basis[1]=gauss->coord2*(2.*gauss->coord2-1.);
			basis[2]=gauss->coord3*(2.*gauss->coord3-1.);
			basis[3]=gauss->coord4*(2.*gauss->coord4-1.);
			/*Edges*/
			basis[4]=4.*gauss->coord2*gauss->coord3;
			basis[5]=4.*gauss->coord1*gauss->coord3;
			basis[6]=4.*gauss->coord1*gauss->coord2;
			basis[7]=4.*gauss->coord2*gauss->coord4;
			basis[8]=4.*gauss->coord3*gauss->coord4;
			basis[9]=4.*gauss->coord1*gauss->coord4;
			return;
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}
}
/*}}}*/
void TetraRef::GetNodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list, GaussTetra* gauss,int finiteelement){/*{{{*/

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
void TetraRef::GetNodalFunctionsDerivativesReference(IssmDouble* dbasis,GaussTetra* gauss,int finiteelement){/*{{{*/
	/*This routine returns the values of the nodal functions derivatives  (with respect to the 
	 * natural coordinate system) at the gaussian point. */

	_assert_(dbasis && gauss);

	switch(finiteelement){
		case P0Enum:
			/*Nodal function 1*/
			dbasis[NUMNODESP0*0+0] = 0.;
			dbasis[NUMNODESP0*1+0] = 0.;
			dbasis[NUMNODESP0*2+0] = 0.;
			return;
		case P1Enum: case P1DGEnum:
			dbasis[NUMNODESP1*0+0] = -1.;
			dbasis[NUMNODESP1*1+0] = -1.;
			dbasis[NUMNODESP1*2+0] = -1.;

			dbasis[NUMNODESP1*0+1] = 1.;
			dbasis[NUMNODESP1*1+1] = 0.;
			dbasis[NUMNODESP1*2+1] = 0.;

			dbasis[NUMNODESP1*0+2] = 0.;
			dbasis[NUMNODESP1*1+2] = 1.;
			dbasis[NUMNODESP1*2+2] = 0.;

			dbasis[NUMNODESP1*0+3] = 0.;
			dbasis[NUMNODESP1*1+3] = 0.;
			dbasis[NUMNODESP1*2+3] = 1.;
			return;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			dbasis[NUMNODESP1b*0+0] = -1.;
			dbasis[NUMNODESP1b*1+0] = -1.;
			dbasis[NUMNODESP1b*2+0] = -1.;

			dbasis[NUMNODESP1b*0+1] = 1.;
			dbasis[NUMNODESP1b*1+1] = 0.;
			dbasis[NUMNODESP1b*2+1] = 0.;

			dbasis[NUMNODESP1b*0+2] = 0.;
			dbasis[NUMNODESP1b*1+2] = 1.;
			dbasis[NUMNODESP1b*2+2] = 0.;

			dbasis[NUMNODESP1b*0+3] = 0.;
			dbasis[NUMNODESP1b*1+3] = 0.;
			dbasis[NUMNODESP1b*2+3] = 1.;

			dbasis[NUMNODESP1b*0+4] = 256.*(-gauss->coord2*gauss->coord3*gauss->coord4+gauss->coord1*gauss->coord3*gauss->coord4);
			dbasis[NUMNODESP1b*1+4] = 256.*(-gauss->coord2*gauss->coord3*gauss->coord4+gauss->coord1*gauss->coord2*gauss->coord4);
			dbasis[NUMNODESP1b*2+4] = 256.*(-gauss->coord2*gauss->coord3*gauss->coord4+gauss->coord1*gauss->coord2*gauss->coord3);
			return;
		case P2Enum:
			dbasis[NUMNODESP2*0+0] = -4.*gauss->coord1+1.;
			dbasis[NUMNODESP2*1+0] = -4.*gauss->coord1+1.;
			dbasis[NUMNODESP2*2+0] = -4.*gauss->coord1+1.;

			dbasis[NUMNODESP2*0+1] = 4.*gauss->coord2-1.;
			dbasis[NUMNODESP2*1+1] = 0.;
			dbasis[NUMNODESP2*2+1] = 0.;

			dbasis[NUMNODESP2*0+2] = 0.;
			dbasis[NUMNODESP2*1+2] = 4.*gauss->coord3-1.;
			dbasis[NUMNODESP2*2+2] = 0.;

			dbasis[NUMNODESP2*0+3] = 0.;
			dbasis[NUMNODESP2*1+3] = 0.;
			dbasis[NUMNODESP2*2+3] = 4.*gauss->coord4-1.;

			dbasis[NUMNODESP2*0+4] = 4.*gauss->coord3;
			dbasis[NUMNODESP2*1+4] = 4.*gauss->coord2;
			dbasis[NUMNODESP2*2+4] = 0.;

			dbasis[NUMNODESP2*0+5] = -4.*gauss->coord3;
			dbasis[NUMNODESP2*1+5] = 4.*(gauss->coord1-gauss->coord3);
			dbasis[NUMNODESP2*2+5] = -4.*gauss->coord3;

			dbasis[NUMNODESP2*0+6] = 4.*(gauss->coord1-gauss->coord2);
			dbasis[NUMNODESP2*1+6] = -4.*gauss->coord2;
			dbasis[NUMNODESP2*2+6] = -4.*gauss->coord2;

			dbasis[NUMNODESP2*0+7] = 4.*gauss->coord4;
			dbasis[NUMNODESP2*1+7] = 0.;
			dbasis[NUMNODESP2*2+7] = 4.*gauss->coord2;

			dbasis[NUMNODESP2*0+8] = 0.;
			dbasis[NUMNODESP2*1+8] = 4.*gauss->coord4;
			dbasis[NUMNODESP2*2+8] = 4.*gauss->coord3;

			dbasis[NUMNODESP2*0+9] = -4.*gauss->coord4;
			dbasis[NUMNODESP2*1+9] = -4.*gauss->coord4;
			dbasis[NUMNODESP2*2+9] = 4.*(gauss->coord1-gauss->coord4);
			return;
		default:
			_error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

}
/*}}}*/
int  TetraRef::NumberofNodes(int finiteelement){/*{{{*/

	switch(finiteelement){
		case P0Enum:                return NUMNODESP0;
		case P1Enum:                return NUMNODESP1;
		case P1DGEnum:              return NUMNODESP1;
		case P1bubbleEnum:          return NUMNODESP1b;
		case P1bubblecondensedEnum: return NUMNODESP1b;
		case P2Enum:                return NUMNODESP2;
		case P1P1Enum:              return NUMNODESP1*2;
		case P1P1GLSEnum:           return NUMNODESP1*2;
		case MINIcondensedEnum:     return NUMNODESP1b+NUMNODESP1;
		case MINIEnum:              return NUMNODESP1b+NUMNODESP1;
		case TaylorHoodEnum:        return NUMNODESP2+NUMNODESP1;
		case LATaylorHoodEnum:      return NUMNODESP2;
		case XTaylorHoodEnum:       return NUMNODESP2+NUMNODESP1;
		default: _error_("Element type "<<EnumToStringx(finiteelement)<<" not supported yet");
	}

	return -1;
}
/*}}}*/
int  TetraRef::PressureInterpolation(int fe_stokes){/*{{{*/

	switch(fe_stokes){
		case P1P1Enum:          return P1Enum;
		case P1P1GLSEnum:       return P1Enum;
		case MINIcondensedEnum: return P1Enum;
		case MINIEnum:          return P1Enum;
		case TaylorHoodEnum:    return P1Enum;
		case LATaylorHoodEnum:  return NoneEnum;
		case XTaylorHoodEnum:   return P1Enum;
		default:       _error_("Element type "<<EnumToStringx(fe_stokes)<<" not supported yet");
	}

	return -1;
}/*}}}*/
int  TetraRef::TensorInterpolation(int fe_stokes){/*{{{*/
	/*This routine returns the values of the nodal functions  at the gaussian point.*/

	switch(fe_stokes){
		case XTaylorHoodEnum: return P1DGEnum;
		default: _error_("Element type "<<EnumToStringx(fe_stokes)<<" not supported yet");
	}
}
/*}}}*/
int  TetraRef::VelocityInterpolation(int fe_stokes){/*{{{*/

	switch(fe_stokes){
		case P1P1Enum:          return P1Enum;
		case P1P1GLSEnum:       return P1Enum;
		case MINIcondensedEnum: return P1bubbleEnum;
		case MINIEnum:          return P1bubbleEnum;
		case TaylorHoodEnum:    return P2Enum;
		case LATaylorHoodEnum:  return P2Enum;
		case XTaylorHoodEnum:   return P2Enum;
		default:       _error_("Element type "<<EnumToStringx(fe_stokes)<<" not supported yet");
	}

	return -1;
}
/*}}}*/
