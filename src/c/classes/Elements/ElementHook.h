/*!\file: ElementHook.h
 * \brief prototypes for ElementHook.h
 */ 

#ifndef _ELEMENTHOOK_H_
#define _ELEMENTHOOK_H_

class Hook;
class IoModel;

class ElementHook{

	public: 
		int    numanalyses;   //number of analysis types
		Hook **hnodes;        // set of nodes for each analysis type
		Hook  *hvertices;     // vertices
		Hook  *hmaterial;     // 1 ice material
		Hook  *hneighbors;    // 2 elements, first down, second up in 3d only

		/*constructors, destructors*/
		ElementHook();
		ElementHook(int in_numanalyses,int material_id,int numvertices,IoModel* iomodel);
		~ElementHook();
		void Marshall(MarshallHandle* marshallhandle);

		void DeepEcho();
		void Echo();
		void InitHookNeighbors(int* element_ids);               //3d only
		void SetHookNodes(int* node_ids,int numnodes,int analysis_counter);
		void SpawnSegHook(ElementHook* triahook,int ndex1,int index2); //2d only
		void SpawnTriaHook(ElementHook* triahook,int index1,int index2,int index3); //3d only
};

#endif //ifndef _ELEMENTHOOK_H_
