/*!\file:  Load.h
 * \brief abstract class for Load object
 * This class is a place holder for the Icefront  and the Penpair loads.
 * It is derived from Load, so DataSets can contain them.
 */ 

#ifndef _LOAD_H_
#define _LOAD_H_

/*Headers:*/
class Node;
template <class doublematrix> class Matrix;
template <class doubletype> class Vector;
class Elements;
class Loads;
class Nodes;
class Vertices;
class Materials;
class Parameters;
#include "../../datastructures/datastructures.h"

class Load: public Object{

	public: 
		virtual       ~Load(){};
		virtual void  Configure(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters)=0;
		virtual void  CreateJacobianMatrix(Matrix<IssmDouble>* Jff)=0;
		virtual void  CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs)=0;
		virtual void  CreatePVector(Vector<IssmDouble>* pf)=0;
		virtual void  GetNodesLidList(int* lidlist)=0;
		virtual void  GetNodesSidList(int* sidlist)=0;
		virtual int   GetNumberOfNodes(void)=0;
		virtual bool  IsPenalty(void)=0;
		virtual void  PenaltyCreateJacobianMatrix(Matrix<IssmDouble>* Jff,IssmDouble kmax)=0;
		virtual void  PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs, IssmDouble kmax)=0;
		virtual void  PenaltyCreatePVector(Vector<IssmDouble>* pf, IssmDouble kmax)=0;
		virtual void  ResetHooks()=0;
		virtual void  SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters)=0;
		virtual void  SetwiseNodeConnectivity(int* d_nz,int* o_nz,Node* node,bool* flags,int* flagsindices,int* flagsindices_counter,int set1_enum,int set2_enum)=0;
};
#endif
