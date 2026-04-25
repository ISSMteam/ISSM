#ifndef _CONTAINER_DATASET_H_
#define _CONTAINER_DATASET_H_

#include <vector>
#include <cstring>

/*forward declarations */
class Object;
class MarshallHandle;

/*
 * Declaration of DataSet class.  A DataSet is a Container of Objects.
 */
class DataSet{

	public: 

		/*internals: */
		std::vector<Object*> objects;

		/*type of dataset: */
		int             enum_type;

		/*sorting: */
		int             sorted;
		int             presorted;
		int             numsorted;
		int*            sorted_ids;
		int*            id_offsets;

		/*constructors, destructors*/
		DataSet();
		DataSet(int enum_type);
		~DataSet();
		void  Marshall(MarshallHandle* marshallhandle);

		/*management*/
		int      GetEnum();
		int      GetEnum(int offset);
		void     Echo();
		void     DeepEcho();
		int      AddObject(Object *object);
		int      DeleteObject(int id);
		int      Size();
		void     clear();
		Object  *GetObjectByOffset(int  offset);
		Object  *GetObjectById(int *poffset,int eid);
		void     Presort();
		void     Sort();
		DataSet *Copy(void);
		int      DeleteObject(Object *object);

};

#endif
