/*!\file Moulin.h
 * \brief: header file for pengrid object */

#ifndef _MOULIN_H_
#define _MOULIN_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include "./Load.h"
class Hook;
class Inputs;
class Parameters;
class IoModel;
/*}}}*/

class Moulin: public Load{

	private: 

		int id;

		/*Hooks*/
		Hook* hnode;     //hook to 1 node
		Hook* hvertex;   //hook to 1 vertex
		Hook* helement;  //hook to 1 element

		/*Corresponding fields*/
		Node    *node;
		Vertex  *vertex;
		Element *element;

		Parameters* parameters; //pointer to solution parameters

	public:

		/*Moulin constructors, destructors {{{*/
		Moulin();
		Moulin(int index, int id, IoModel* iomodel);
		~Moulin();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void  Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum();
		/*}}}*/
		/*Update virtual functions resolution: {{{*/
		void  InputUpdateFromConstant(IssmDouble constant, int name);
		void  InputUpdateFromConstant(int constant, int name);
		void  InputUpdateFromConstant(bool constant, int name);
		void  InputUpdateFromIoModel(int index, IoModel* iomodel){_error_("not implemented yet");};
		void  InputUpdateFromMatrixDakota(IssmDouble* matrix ,int nrows, int ncols, int name, int type);
		void  InputUpdateFromVector(IssmDouble* vector, int name, int type);
		void  InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type);
		/*}}}*/
		/*Load virtual functions definitions: {{{*/
		void  Configure(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);
		void  CreateJacobianMatrix(Matrix<IssmDouble>* Jff){_error_("Not implemented yet");};
		void  CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs);
		void  CreatePVector(Vector<IssmDouble>* pf);
		void  GetNodesSidList(int* sidlist);
		void  GetNodesLidList(int* lidlist);
		int   GetNumberOfNodes(void);
		bool  IsPenalty(void);
		void  PenaltyCreateJacobianMatrix(Matrix<IssmDouble>* Jff,IssmDouble kmax){_error_("Not implemented yet");};
		void  PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* kfs, IssmDouble kmax);
		void  PenaltyCreatePVector(Vector<IssmDouble>* pf, IssmDouble kmax);
		void  SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);
		void  SetwiseNodeConnectivity(int* d_nz,int* o_nz,Node* node,bool* flags,int* flagsindices,int* flagsindices_counter,int set1_enum,int set2_enum);
		void  ResetHooks();
		/*}}}*/

		ElementMatrix* CreateKMatrixHydrologyGlaDS(void);
		ElementVector* CreatePVectorHydrologyShakti(void);
		ElementVector* CreatePVectorHydrologyGlaDS(void);
		ElementVector* CreatePVectorHydrologyDCInefficient(void);
		ElementVector* CreatePVectorHydrologyDCEfficient(void);
};

#endif  /* _MOULIN_H_ */
