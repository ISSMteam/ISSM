
/*!\file:  TetraRef.h
 * \brief abstract class for handling Tetra oriented routines, like nodal functions, 
 * strain rate generation, etc ...
 */ 

#ifndef _TETRAREF_H_
#define _TETRAREF_H_

class GaussTetra;

class TetraRef{

	public: 
		TetraRef();
		~TetraRef();

		void GetInputDerivativeValue(IssmDouble* p, IssmDouble* plist,IssmDouble* xyz_list, GaussTetra* gauss,int finiteelement);
		void GetInputValue(IssmDouble* p, IssmDouble* plist, Gauss* gauss,int finiteelement);
		void GetJacobian(IssmDouble* J, IssmDouble* xyz_list,GaussTetra* gauss);
		void GetJacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,GaussTetra* gauss);
		void GetJacobianDeterminantFace(IssmDouble*  Jdet, IssmDouble* xyz_list,GaussTetra* gauss);
		void GetJacobianInvert(IssmDouble* Jinv, IssmDouble* xyz_list,GaussTetra* gauss);
		void GetNodalFunctions(IssmDouble* basis,Gauss* gauss_in,int finiteelement);
		void GetNodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list, GaussTetra* gauss,int finiteelement);
		void GetNodalFunctionsDerivativesReference(IssmDouble* dbasis,GaussTetra* gauss,int finiteelement);
		void Marshall(MarshallHandle* marshallhandle){ /*do nothing */};
		int  NumberofNodes(int finiteelement);
		int  PressureInterpolation(int fe_stokes);
		int  TensorInterpolation(int fe_stokes);
		int  VelocityInterpolation(int fe_stokes);
};
#endif
