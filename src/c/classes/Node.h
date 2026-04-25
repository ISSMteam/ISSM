/*!\file Node.h
 * \brief: header file for node object
 */

#ifndef _NODE_H_
#define _NODE_H_

/*Headers:*/
/*{{{*/
#include "../datastructures/datastructures.h"
#include "../shared/shared.h"
class  Inputs;
class  Hook;
class  IoModel;
class  DataSet;
class  Vertices;
template <class doubletype> class Vector;
template <class doubletype> class Matrix;
class ElementVector;
class ElementMatrix;
/*}}}*/

class Node: public Object{

	private:
		int  approximation; //For ice flow models, we need to know what ice flow approximation is employed on this node
		bool clone;  //this node is replicated from another one
		int  id;    // unique arbitrary id.
		int  sid;   // "serial" id (rank of this node if the dataset was serial on 1 cpu)
		int  lid;   // "local"  id (rank of this node in current partition)
		int  pid;   // parallel id (specific to this partition)

		/*Only this function can access these private fields*/
		friend class Nodes;
		friend class FemModel;

	public: 
		int        analysis_enum;
		bool       indexingupdate;
		bool       isrotated;
		IssmDouble coord_system[3][3];

		/*sizes: */
		int gsize;   //number of dofs for a node

		/*Activation*/
		bool active; //Is this node active or inactive (all dofs are constrained)
		bool freeze; //this is required for 2d solutions, we never activate nodes that are not on base

		/*boundary conditions sets: */
		bool       *f_set;     //is dof on f-set (on which we solve)
		bool       *s_set;     //is dof on s-set (on which boundary conditions -dirichlet- are applied)
		IssmDouble *svalues;   //list of constraint values. size g_size, for ease of use.

		/*types of dofs: */
		int  *doftype;   //approximation type of the dofs (used only for coupling), size g_size

		/*list of degrees of freedom: */
		int *gdoflist;
		int *fdoflist;
		int *sdoflist;
		int *gdoflist_local;
		int *fdoflist_local;
		int *sdoflist_local;

		/*Node constructors, destructors*/
		Node();
		Node(int node_id,int node_sid,int io_index,bool isclone,IoModel* iomodel,int analysis_enum,int approximation_in,bool isamr);
		~Node();

		/*Object virtual functions definitions:*/
		Object *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();

		/*Node numerical routines*/
		void  Activate(void);
		void  ApplyConstraint(int dof,IssmDouble value);
		void  CreateNodalConstraints(Vector<IssmDouble>* ys);
		void  Deactivate(void);
		void  DistributeLocalDofs(int* pdofcount,int setenum);
		void  DofInFSet(int dof);
		void  DofInSSet(int dof);
		void  FreezeDof(int dof);
		int   GetApproximation();
		void  GetCoordinateSystem(IssmDouble* coord_system_out);
		int   GetDof(int dofindex,int setenum);
		void  GetDofList(int* poutdoflist,int approximation_enum,int setenum,bool hideclones=0);
		void  GetDofListLocal(int* poutdoflist,int approximation_enum,int setenum);
		int   GetNumberOfDofs(int approximation_enum,int setenum);
		void  HardDeactivate(void);
		bool  IsActive(void);
		int   IsClone();
		int   Lid(void); 
		void  DistributeGlobalDofsMasters(int dofcount,int setenum);
		void  ReindexingDone(void);
		void  RelaxConstraint(int dof);
		bool  RequiresDofReindexing(void);
		void  SetCurrentConfiguration(DataSet* nodes,Vertices* vertices);
		void  ShowMasterDofs(int* truerows,int setenum);
		int   Sid(void); 
		int   Pid(void); 
		void  UpdateCloneDofs(int* alltruerows,int setenum);
		void  VecMerge(Vector<IssmDouble>* ug,IssmDouble* local_uf,int* indices_uf,IssmDouble* local_ys,int* indices_ys);
		void  VecReduce(Vector<IssmDouble>* uf, IssmDouble* local_ug,int* indices_ug);
		void  SetApproximation(int in_approximation);
		int   SSize();
		int   FSize();
};

/*Methods inherent to Node: */
int* GetGlobalDofList(Node** nodes,int numnodes,int setenum,int approximation);
int  GetNumberOfDofs(Node** nodes,int numnodes,int setenum,int approximation);

#endif  /* _NODE_H_ */
