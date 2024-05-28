/*!\file Pengrid.h
 * \brief: header file for pengrid object */

#ifndef _PENGRID_H_
#define _PENGRID_H_

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

class Pengrid: public Load{

	private: 

		int id;

		/*Hooks*/
		Hook* hnode;  //hook to 1 node
		Hook* helement;  //hook to 1 element

		/*Corresponding fields*/
		Node    *node;
		Element *element;

		Parameters* parameters; //pointer to solution parameters

		/*internals: */
		int active;
		int zigzag_counter;

	public:

		/*Pengrid constructors, destructors {{{*/
		Pengrid();
		Pengrid(int id, int index, IoModel* iomodel);
		~Pengrid();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void  Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum();
		/*}}}*/
		/*Load virtual functions definitions: {{{*/
		void  Configure(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);
		void  CreateJacobianMatrix(Matrix<IssmDouble>* Jff){_error_("Not implemented yet");};
		void  CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs);
		void  CreatePVector(Vector<IssmDouble>* pf);
		void  GetNodesLidList(int* lidlist);
		void  GetNodesSidList(int* sidlist);
		int   GetNumberOfNodes(void);
		bool  IsPenalty(void);
		void  PenaltyCreateJacobianMatrix(Matrix<IssmDouble>* Jff,IssmDouble kmax){_error_("Not implemented yet");};
		void  PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* kfs, IssmDouble kmax);
		void  PenaltyCreatePVector(Vector<IssmDouble>* pf, IssmDouble kmax);
		void  ResetHooks();
		void  SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);
		void  SetwiseNodeConnectivity(int* d_nz,int* o_nz,Node* node,bool* flags,int* flagsindices,int* flagsindices_counter,int set1_enum,int set2_enum);
		/*}}}*/
		/*Pengrid management {{{*/
		void				ConstraintActivate(int* punstable);
		void           ConstraintActivateHydrologyDCInefficient(int* punstable);
		void           ConstraintActivateThermal(int* punstable);
		ElementMatrix* PenaltyCreateKMatrixHydrologyDCInefficient(IssmDouble kmax);
		ElementMatrix* PenaltyCreateKMatrixMelting(IssmDouble kmax);
		ElementMatrix* PenaltyCreateKMatrixThermal(IssmDouble kmax);
		ElementVector* PenaltyCreatePVectorHydrologyDCInefficient(IssmDouble kmax);
		ElementVector* PenaltyCreatePVectorMelting(IssmDouble kmax);
		ElementVector* PenaltyCreatePVectorThermal(IssmDouble kmax);
		void  ResetConstraint(void);
		void  ResetZigzagCounter(void);
		/*}}}*/

};

#endif  /* _PENGRID_H_ */
