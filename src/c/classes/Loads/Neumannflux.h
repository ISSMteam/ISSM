/*!\file Neumannflux.h
 * \brief: header file for icefront object
 */

#ifndef _NEUMANNFLUX_H_
#define _NEUMANNFLUX_H_

/*Headers:*/
#include "./Load.h"
class Hook;
class Parameters;
class IoModel;
class Element;
class Vertex;
class ElementMatrix;
class ElementVector;

class Neumannflux: public Load {

	public: 
		int id;

		/*Hooks*/
		Hook *helement;
		Hook *hnodes;
		Hook *hvertices;

		/*Corresponding fields*/
		Element     *element;
		Vertex     **vertices;
		Node       **nodes;
		Parameters  *parameters;

		/*Neumannflux constructors,destructors {{{*/
		Neumannflux();
		Neumannflux(int numericalflux_id,int i,IoModel* iomodel,int* segments);
		~Neumannflux();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*Update virtual functions resolution: {{{*/
		void InputUpdateFromConstant(IssmDouble constant, int name){/*Do nothing*/};
		void InputUpdateFromConstant(int constant, int name){/*Do nothing*/};
		void InputUpdateFromConstant(bool constant, int name){/*Do nothing*/};
		void InputUpdateFromIoModel(int index, IoModel* iomodel){_error_("not implemented yet");};
		void InputUpdateFromMatrixDakota(IssmDouble* matrix, int nrows, int ncols, int name, int type){/*Do nothing*/}
		void InputUpdateFromVector(IssmDouble* vector, int name, int type){/*Do nothing*/}
		void InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type){/*Do nothing*/}
		/*}}}*/
		/*Load virtual functions definitions: {{{*/
		void Configure(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);
		void CreateJacobianMatrix(Matrix<IssmDouble>* Jff){_error_("Not implemented yet");};
		void CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs);
		void CreatePVector(Vector<IssmDouble>* pf);
		void GetNodesLidList(int* lidlist);
		void GetNodesSidList(int* sidlist);
		int  GetNumberOfNodes(void);
		bool IsPenalty(void);
		void PenaltyCreateJacobianMatrix(Matrix<IssmDouble>* Jff,IssmDouble kmax){_error_("Not implemented yet");};
		void PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* kfs, IssmDouble kmax);
		void PenaltyCreatePVector(Vector<IssmDouble>* pf, IssmDouble kmax);
		void ResetHooks();
		void SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);
		void SetwiseNodeConnectivity(int* d_nz,int* o_nz,Node* node,bool* flags,int* flagsindices,int* flagsindices_counter,int set1_enum,int set2_enum);
		/*}}}*/
		/*Neumannflux management:{{{*/
		ElementVector* CreatePVectorHydrologyShakti(void);
		ElementVector* CreatePVectorHydrologyGlaDS(void);
		/*}}}*/

};

#endif  /* _NEUMANNFLUX_H_ */
