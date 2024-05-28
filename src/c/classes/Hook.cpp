/*!\file Hook.cpp
 * \brief: implementation of the Hook object: see Hook.h for more explanations.
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include <string.h>
#include "./classes.h"
#include "../shared/Enum/Enum.h"
#include "../shared/shared.h"

/*Constructor/Destructors*/
Hook::Hook(){/*{{{*/
	this->num     = 0;
	this->objects = NULL;
	this->ids     = NULL;
	this->offsets = NULL;
}
/*}}}*/
Hook::Hook(int* in_ids, int in_num){/*{{{*/

	/*Get number of objects to hook*/
	this->num=in_num;

	/*Get out if num=0*/
	if (this->num<=0){
		/*Empty hook*/
		this->ids     = NULL;
		this->objects = NULL;
		this->offsets = NULL;
		this->num = 0;
	}
	else{
		/*Allocate: */
		this->objects=xNew<Object*>(in_num);
		this->ids=xNew<int>(in_num);
		this->offsets=xNew<int>(in_num);

		/*Copy ids: */
		for(int i=0;i<this->num;i++){
			this->ids[i]     = in_ids[i];
			this->objects[i] = NULL;
			this->offsets[i] = 0;
		}
	}
}
/*}}}*/
Hook::~Hook(){/*{{{*/
	xDelete<Object*>(this->objects);
	xDelete<int>(this->ids);
	xDelete<int>(this->offsets);
}
/*}}}*/

/*Some of the Object functionality: */
Object* Hook::copy(void){/*{{{*/

	/*output: */
	Hook* output=NULL;

	/*initalize output: */
	output=new Hook(this->ids,this->num);

	for(int i=0;i<output->num;i++){
		output->objects[i] = this->objects[i];
		output->offsets[i] = this->offsets[i];
	}

	return (Object*)output;
}
/*}}}*/
void Hook::DeepEcho(void){/*{{{*/

	int i;
	if (num){
		_printf_("   Hook: \n");
		_printf_("      num=" << this->num << "\n");
		_printf_("      ids: ");
		for (i=0;i<this->num;i++) _printf_(this->ids[i] << " ");
		_printf_("\n");
		_printf_("      offsets: ");
		for (i=0;i<this->num;i++) _printf_(this->offsets[i] << " ");
		_printf_("\n");
		if (!objects) _printf_("      warning: object not hooked yet\n");
		else{
			_printf_("      objects:\n   ");
			for (i=0;i<this->num;i++){
				_printf_("         object " << i << "\n");
				if(objects[i]) objects[i]->DeepEcho();
				else           _printf_("            no object hooked yet (not configured)\n");
			}
		}
	}
	else{
		_printf_("   Hook: num=0 \n");
	}
}
/*}}}*/
void Hook::Echo(void){/*{{{*/
	_assert_(this);
	int i;
	if (num){
		_printf_("   Hook: \n");
		_printf_("      num=" << this->num << "\n");
		_printf_("      ids: ");
		for(i=0;i<this->num;i++) _printf_(this->ids[i] << " ");
		_printf_("\n");
		_printf_("      offsets: ");
		for (i=0;i<this->num;i++) _printf_(this->offsets[i] << " ");
		_printf_("\n");
	}
	else{
		_printf_("   Hook: num=0 \n");
	}
}
/*}}}*/
void Hook::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD) reset();

	int object_enum = HookEnum;
	marshallhandle->call(object_enum);
	marshallhandle->call(this->num);
	if (num<=0){
		/*Empty hook*/
		this->ids     = NULL;
		this->objects = NULL;
		this->offsets = NULL;
		this->num = 0;
	}
	else{
		marshallhandle->call(this->ids,num);
		marshallhandle->call(this->offsets,num);
		marshallhandle->call(this->objects,num);
	}

}
/*}}}*/

/*Hook management: */
void Hook::configure(DataSet* dataset){/*{{{*/

	/*intermediary: */
	Object* object=NULL;
	int i;

	/*Checks if debugging mode*/
	_assert_(this->num==0 || this->ids!=NULL);

	for(i=0;i<this->num;i++){

		/*is this object id -1? If so, drop this search, it was not requested: */
		if (this->ids[i]==-1) continue;

		/*Check whether existing this->objects are correct: */
		if(this->objects[i]){
			if(this->objects[i]->Id()==this->ids[i]) continue; //this node is good.
			else this->objects[i]=NULL; //this node was incorrect, reset it.
		}

		/*May be the object this->offsets into this->objects are valid?: */
		if(this->offsets[i]!=UNDEF){
			/* Look at the this->offsets[i]'th node in the nodes dataset. If it has the correct id, 
			 * we are good: */
			object=(Object*)dataset->GetObjectByOffset(this->offsets[i]);
			if (object->Id()==this->ids[i]){
				this->objects[i]=object;
				continue;
			}
			else this->offsets[i]=UNDEF; //object offset was wrong, reset it.
		}
		else this->offsets[i]=UNDEF;

		/*Now, for this->objects that did not get resolved, and for which we have no offset, chase them in the dataset, by id: */
		if(this->objects[i]==NULL){
			this->objects[i]=xDynamicCast<Object*>(dataset->GetObjectById(this->offsets+i,this->ids[i])); //remember the offset for later on.
			/*check the id is correct!: */
			if (this->objects[i]->Id()!=this->ids[i]) _error_("wrong id: " << this->objects[i]->Id() << " vs " << this->ids[i] << "  in resolved pointer!");
		}
	}
}
/*}}}*/
Object** Hook::deliverp(void){/*{{{*/
	return objects;
}
/*}}}*/
Object* Hook::delivers(void){/*{{{*/

	/*first, check that we only have one T object in our object list: */
	if (this->num!=1) _error_("trying to deliver a single hook object when hook holds " << this->num << " objects" << "\n");

	/*check NULL: */
	if (this->objects==NULL) _error_("hook is not pointing to any object, objects pointer is NULL");

	return *objects;
}

/*}}}*/
int Hook::GetNum(void){/*{{{*/
	return this->num;
}
/*}}}*/
int* Hook::Ids(void){/*{{{*/
	return this->ids;
}
/*}}}*/
void Hook::reset(){/*{{{*/

	/*intermediary: */
	Object* object=NULL;
	int i;

	for(i=0;i<this->num;i++){
			this->objects[i]=NULL; //reset this node.
	}
}
/*}}}*/
Hook* Hook::Spawn(int* indices, int numindices){/*{{{*/

	/*allocate output*/
	Hook* output=new Hook();

	/*If this Hook is empty, simply return*/
	if(this->num==0){
		output->num=0;
		return output;
	}

	/*Else, check that we are requesting a half of num*/
	if(numindices>this->num) _error_("Cannot spawn hook with " << numindices << " objects from a Hook of " << this->num << " objects");

	/*go pickup the correct objects, ids and offsets :*/
	output->num=numindices;
	if(output->num<1) _error_("Trying to spawn an empty ElementProperties!");

	output->objects = xNew<Object*>(output->num);
	output->ids     = xNew<int>(output->num);
	output->offsets = xNew<int>(output->num);

	for(int i=0;i<output->num;i++){
		output->objects[i] = this->objects[indices[i]];
		output->ids[i]     = this->ids[indices[i]];
		output->offsets[i] = this->offsets[indices[i]];
	}

	return output;
}
/*}}}*/
