/*!\file:  TriaRef.h
 * \brief abstract class for handling Tria oriented routines, like nodal functions, 
 * strain rate generation, etc ...
 */ 

#ifndef _TRIAREF_H_
#define _TRIAREF_H_

class Gauss;

class TriaRef{

	public: 
		TriaRef();
		~TriaRef();

		/*Numerics*/
		void GetInputDerivativeValue(IssmDouble* pp, IssmDouble* plist,IssmDouble* xyz_list, Gauss* gauss,int finiteelement);
		void GetInputValue(IssmDouble* pp, IssmDouble* plist, Gauss* gauss,int finiteelement);
		void GetJacobian(IssmDouble* J, IssmDouble* xyz_list,Gauss* gauss);
		void GetJacobianDeterminant(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss);
		void GetJacobianDeterminant3D(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss);
		void GetJacobianInvert(IssmDouble*  Jinv, IssmDouble* xyz_list,Gauss* gauss);
		void GetNodalFunctions(IssmDouble* basis,Gauss* gauss,int finiteelement);
		void GetNodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list, Gauss* gauss,int finiteelement);
		void GetNodalFunctionsDerivativesReference(IssmDouble* dbasis,Gauss* gauss,int finiteelement);
		void GetSegmentJacobianDeterminant(IssmDouble* Jdet, IssmDouble* xyz_list,Gauss* gauss);
		void GetSegmentNodalFunctions(IssmDouble* basis,Gauss* gauss, int index1,int index2,int finiteelement);
		void GetSegmentNodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list_tria,Gauss* gauss, int index1,int index2,int finiteelement);
		void Marshall(MarshallHandle* marshallhandle){ /*do nothing */};
		void NodeOnEdgeIndices(int* pnumindices,int** pindices,int index,int finiteelement);
		int  NumberofNodes(int finiteelement);
		int  PressureInterpolation(int fe_stokes);
		int  TensorInterpolation(int fe_stokes);
		int  VelocityInterpolation(int fe_stokes);

};
#endif
