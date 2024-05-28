/*!\file Riftfront.h
 * \brief: header file for riftfront object
 */

#ifndef _RIFTFRONT_H_
#define _RIFTFRONT_H_

/*Headers:*/
/*{{{*/
#include "./Load.h"
class Hook;
class Parameters;
class IoModel;
/*}}}*/

class Riftfront: public Load {

	public:
		int		id;

		/*properties*/
		int        type;
		int        fill;
		IssmDouble friction;
		IssmDouble fractionincrement;
		bool       shelf;

		/*hooks: */
		Hook* hnodes;
		Hook* hvertices;
		Hook* helements;

		/*Corresponding fields*/
		Node    **nodes;
		Vertex  **vertices;
		Element **elements;

		/*computational: */
		int         penalty_lock;
		bool        active;
		bool        frozen;
		int         counter;
		bool        prestable;
		bool        material_converged;
		IssmDouble  normal[2];
		IssmDouble  length;
		IssmDouble  fraction;
		int         state;

		Parameters *parameters;           //pointer to solution parameters

		/*Riftfrontconstructors,destructors: {{{*/
		Riftfront();
		Riftfront(int riftfront_id,int i, IoModel* iomodel);
		~Riftfront();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object*  copy();
		void     DeepEcho();
		void     Echo();
		int      Id(); 
		void		Marshall(MarshallHandle* marshallhandle);
		int      ObjectEnum();
		/*}}}*/
		/*Update virtual functions resolution: {{{*/
		void    InputUpdateFromConstant(IssmDouble constant, int name);
		void    InputUpdateFromConstant(int constant, int name){_error_("Not implemented yet!");}
		void    InputUpdateFromConstant(bool constant, int name);
		void    InputUpdateFromIoModel(int index, IoModel* iomodel){_error_("not implemented yet");};
		void    InputUpdateFromMatrixDakota(IssmDouble* matrix, int nrows,int ncols, int name, int type){_error_("Not implemented yet!");}
		void    InputUpdateFromVector(IssmDouble* vector, int name, int type);
		void    InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type){_error_("Not implemented yet!");}
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
		/*Riftfront specific routines: {{{*/
		int            Constrain(int* punstable);
		void           FreezeConstraints(void);
		bool           IsFrozen(void);
		ElementMatrix* PenaltyCreateKMatrixStressbalanceHoriz(IssmDouble kmax);
		ElementVector* PenaltyCreatePVectorStressbalanceHoriz(IssmDouble kmax);
		/*}}}*/
};
#endif  /* _RIFTFRONT_H_ */
