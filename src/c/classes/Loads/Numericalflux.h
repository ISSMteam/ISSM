/*!\file Numericalflux.h
 * \brief: header file for icefront object
 */

#ifndef _NUMERICALFLUX_H_
#define _NUMERICALFLUX_H_

/*Headers:*/
#include "./Load.h"
class Hook;
class Parameters;
class IoModel;
class Element;
class Vertex;
class ElementMatrix;
class ElementVector;

class Numericalflux: public Load {

	public: 
		int id;
		int flux_type;
		int flux_degree;

		/*Hooks*/
		Hook *helement;
		Hook *hnodes;
		Hook *hvertices;

		/*Corresponding fields*/
		Element     *element;
		Vertex     **vertices;
		Node       **nodes;
		Parameters  *parameters;

		/*Numericalflux constructors,destructors {{{*/
		Numericalflux();
		Numericalflux(int numericalflux_id,int i,int index,IoModel* iomodel);
		~Numericalflux();
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
		void InputUpdateFromConstant(bool constant, int name){_error_("Not implemented yet!");}
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
		int  GetNumberOfNodesOneSide(void);
		bool IsPenalty(void);
		void PenaltyCreateJacobianMatrix(Matrix<IssmDouble>* Jff,IssmDouble kmax){_error_("Not implemented yet");};
		void PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* kfs, IssmDouble kmax);
		void PenaltyCreatePVector(Vector<IssmDouble>* pf, IssmDouble kmax);
		void ResetHooks();
		void SetCurrentConfiguration(Elements* elements,Loads* loads,Nodes* nodes,Vertices* vertices,Materials* materials,Parameters* parameters);
		void SetwiseNodeConnectivity(int* d_nz,int* o_nz,Node* node,bool* flags,int* flagsindices,int* flagsindices_counter,int set1_enum,int set2_enum);
		/*}}}*/
		/*Numericalflux management:{{{*/
		ElementMatrix* CreateKMatrixAdjointBalancethickness(void);
		ElementMatrix* CreateKMatrixAdjointBalancethicknessBoundary(void);
		ElementMatrix* CreateKMatrixAdjointBalancethicknessInternal(void);
		ElementMatrix* CreateKMatrixBalancethickness(void);
		ElementMatrix* CreateKMatrixBalancethicknessBoundary(void);
		ElementMatrix* CreateKMatrixBalancethicknessInternal(void);
		ElementMatrix* CreateKMatrixMasstransport(void);
		ElementMatrix* CreateKMatrixMasstransportBoundary(void);
		ElementMatrix* CreateKMatrixMasstransportInternal(void);
		ElementVector* CreatePVectorAdjointBalancethickness(void);
		ElementVector* CreatePVectorBalancethickness(void);
		ElementVector* CreatePVectorBalancethicknessBoundary(void);
		ElementVector* CreatePVectorBalancethicknessInternal(void);
		ElementVector* CreatePVectorMasstransport(void);
		ElementVector* CreatePVectorMasstransportBoundary(void);
		ElementVector* CreatePVectorMasstransportInternal(void);
		void           GetNormal(IssmDouble* normal,IssmDouble xyz_list[4][3]);
		/*}}}*/

};

#endif  /* _NUMERICALFLUX_H_ */
