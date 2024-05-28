
/*!\file:  SegRef.h
 * \brief abstract class for handling Seg oriented routines, like nodal functions, 
 * strain rate generation, etc ...
 */ 

#ifndef _SEGREF_H_
#define _SEGREF_H_

class GaussSeg;

class SegRef{

	public: 
		SegRef();
		~SegRef();

		void GetInputDerivativeValue(IssmDouble* p, IssmDouble* plist,IssmDouble* xyz_list, GaussSeg* gauss,int finiteelement);
		void GetInputValue(IssmDouble* p, IssmDouble* plist, GaussSeg* gauss,int finiteelement);
		void GetJacobian(IssmDouble* J, IssmDouble* xyz_list,GaussSeg* gauss);
		void GetJacobianDeterminant(IssmDouble*  Jdet, IssmDouble* xyz_list,GaussSeg* gauss);
		void GetJacobianInvert(IssmDouble* Jinv, IssmDouble* xyz_list,GaussSeg* gauss);
		void GetNodalFunctions(IssmDouble* basis,GaussSeg* gauss,int finiteelement);
		void GetNodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list, GaussSeg* gauss,int finiteelement);
		void GetNodalFunctionsDerivativesReference(IssmDouble* dbasis,GaussSeg* gauss,int finiteelement);
		void Marshall(MarshallHandle* marshallhandle){ /*do nothing */};
		int  NumberofNodes(int finiteelement);
};
#endif
