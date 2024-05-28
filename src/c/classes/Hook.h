/*!\file Hook.h
 * \brief: header file for hook object.
 * A hook is a class  that can store the id, offset, and object corresponding to this id and offset into a dataset.
 * For example, an element has a hook to its nodes. A node has a hook to its vertex.  The hook abstracts the need for having
 * ids and offsets (necesarry for proper configuration of an object) in our objects. 
 */

#ifndef _HOOK_H_
#define _HOOK_H_

/*Headers:*/
/*{{{*/
#include "../datastructures/datastructures.h"
/*}}}*/

class Hook{

	private: 

		int     *ids;       //list of object ids, to go look for them in datasets.
		int      num;       //number of objects being hooked onto
		Object **objects;   //list of object pointers
		int     *offsets;   //list of object offsets into datasets, to speed up lookup.

	public:

		/*Hook constructors, destructors: {{{*/
		Hook();
		Hook(int* ids, int num);
		~Hook();
		/*}}}*/
		/*Object like functionality:{{{*/
		Object*    copy(void);
		void       DeepEcho(void);
		void       Echo(void);
		void       Marshall(MarshallHandle* marshallhandle);
		/*}}}*/
		/*Hook management: {{{*/
		void       configure(DataSet* dataset);
		Object**   deliverp(void); //deliver all objects
		Object*    delivers(void); //single object deliver
		int        GetNum(void);
		int*       Ids(void);
		void       reset(void);
		Hook*      Spawn(int* indices, int numindices);
		/*}}}*/
};

#endif  /* _HOOK_H_ */
