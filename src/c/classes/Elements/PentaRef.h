/*!\file:  PentaRef.h
 * \brief abstract class for handling Penta oriented routines, like nodal functions, 
 * strain rate generation, etc ...
 */ 

#ifndef _PENTAREF_H_
#define _PENTAREF_H_

class Gauss;
class PentaRef{

	public: 
		PentaRef();
		~PentaRef();

		/*Numerics*/
		void BasalNodeIndices(int* pnumindices,int** pindices,int finiteelement);
		void GetInputDerivativeValue(IssmDouble* pvalues, IssmDouble* plist,IssmDouble* xyz_list, Gauss* gauss,int finiteelement);
		void GetInputValue(IssmDouble* pvalue,IssmDouble* plist, Gauss* gauss,int finiteelement);
		void GetJacobian(IssmDouble* J, IssmDouble* xyz_list,Gauss* gauss);
		void GetJacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,Gauss* gauss);
		void GetJacobianInvert(IssmDouble*  Jinv, IssmDouble* xyz_list,Gauss* gauss);
		void GetLprimeFSSSA(IssmDouble* LprimeFSSSA, IssmDouble* xyz_list, Gauss* gauss);
		void GetNodalFunctions(IssmDouble* basis, Gauss* gauss,int finiteelement);
		void GetNodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss,int finiteelement);
		void GetNodalFunctionsDerivativesReference(IssmDouble* dbasis,Gauss* gauss,int finiteelement);
		void GetQuadJacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,Gauss* gauss);
		void GetSegmentJacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,Gauss* gauss);
		void GetTriaJacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,Gauss* gauss);
		void VerticalSegmentIndicesBase(int** pindices,int* pnumseg,int finiteelement);
		void Marshall(MarshallHandle* marshallhandle){ /*do nothing */};
		int  NumberofNodes(int finiteelement);
		int  PressureInterpolation(int fe_stokes);
		void SurfaceNodeIndices(int* pnumindices,int** pindices,int finiteelement);
		int  TensorInterpolation(int fe_stokes);
		int  VelocityInterpolation(int fe_stokes);
};
#endif
