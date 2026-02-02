/*
 * This prototype describes the Object class. This is an abstract class, parent 
 * to any other objects (Quad, Tria, Node , etc ...), that can be included in a 
 * DataSet.
 */

#ifndef _OBJECT_H_
#define _OBJECT_H_

class Object {

	public: 
		virtual ~Object(){};
		virtual void    Echo()       = 0;
		virtual void    DeepEcho()   = 0;
		virtual int     Id()         = 0;
		virtual int     ObjectEnum() = 0;
		virtual Object *copy()       = 0;
		virtual void    Marshall(MarshallHandle *marshallhandle)=0;

};
#endif
