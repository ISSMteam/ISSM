/*!\file Penpair.h
 * \brief: header file for penpair object */

#ifndef _PENPAIR_H_
#define _PENPAIR_H_

/*Headers:*/
/*{{{*/
#include "./Load.h"
#include "../Node.h"
#include "../Elements/Element.h"

class Element;
/*}}}*/

class Penpair: public Load{

	private: 
		int          id;
		Hook        *hnodes;          //hook to 2 nodes
		Node       **nodes;
		Parameters  *parameters;      //pointer to solution parameters

	public:

		/*Penpair constructors, destructors: {{{*/
		Penpair();
		Penpair(int penpair_id,int* penpair_node_ids);
		~Penpair();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object*  copy();
		void     DeepEcho();
		void     Echo();
		int      Id(); 
		void     Marshall(MarshallHandle* marshallhandle);
		int      ObjectEnum();
		/*}}}*/
		/*Update virtual functions resolution: {{{*/
		void  InputUpdateFromConstant(IssmDouble constant, int name);
		void  InputUpdateFromConstant(int constant, int name);
		void  InputUpdateFromConstant(bool constant, int name);
		void  InputUpdateFromIoModel(int index, IoModel* iomodel){_error_("not implemented yet");};
		void  InputUpdateFromMatrixDakota(IssmDouble* matrix, int nrow, int ncols,int name, int type){_error_("Not implemented yet!");}
		void  InputUpdateFromVector(IssmDouble* vector, int name, int type);
		void  InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type){_error_("Not implemented yet!");}
		/*}}}*/
			/*Load virtual functions definitions: {{{*/
		void  Configure(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);
		void  CreateJacobianMatrix(Matrix<IssmDouble>* Jff);
		void  CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs);
		void  CreatePVector(Vector<IssmDouble>* pf);
		void  GetNodesLidList(int* lidlist);
		void  GetNodesSidList(int* sidlist);
		int   GetNumberOfNodes(void);
		bool  IsPenalty(void);
		void  PenaltyCreateJacobianMatrix(Matrix<IssmDouble>* Jff,IssmDouble kmax);
		void  PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff,Matrix<IssmDouble>* Kfs,IssmDouble kmax);
		void  PenaltyCreatePVector(Vector<IssmDouble>* pf, IssmDouble kmax);
		void  ResetHooks();
		void  SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);
		void  SetwiseNodeConnectivity(int* d_nz,int* o_nz,Node* node,bool* flags,int* flagsindices,int* flagsindices_counter,int set1_enum,int set2_enum);
		/*}}}*/
			/*Penpair management: {{{*/
		ElementMatrix* PenaltyCreateKMatrixMasstransport(IssmDouble kmax);
		ElementMatrix* PenaltyCreateKMatrixStressbalanceFS(IssmDouble kmax);
		ElementMatrix* PenaltyCreateKMatrixStressbalanceHoriz(IssmDouble kmax);
		ElementMatrix* PenaltyCreateKMatrixStressbalanceSSAHO(IssmDouble kmax);
		/*}}}*/
};

#endif  /* _PENPAIR_H_ */
