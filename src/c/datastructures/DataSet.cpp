/*
 * \file DataSet.cpp
 * \brief: Implementation of DataSet class
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <cstring>
#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>

#include "../datastructures/datastructures.h"
#include "../shared/shared.h"
#include "../classes/classes.h"

using namespace std;
/*}}}*/

/*Constructors/Destructors*/
DataSet::DataSet(){/*{{{*/

	sorted=0;
	numsorted=0;
	presorted=0;
	enum_type=-1;
	sorted_ids=NULL;
	id_offsets=NULL;

}
/*}}}*/
DataSet::DataSet(int dataset_enum){/*{{{*/
	enum_type=dataset_enum;

	sorted=0;
	numsorted=0;
	presorted=0;
	sorted_ids=NULL;
	id_offsets=NULL;

}
/*}}}*/
DataSet* DataSet::Copy(void){/*{{{*/

	vector<Object*>::iterator obj;
	Object* object_copy=NULL;

	DataSet* copy=new DataSet(this->enum_type);

	copy->sorted=this->sorted;
	copy->numsorted=this->numsorted;
	copy->presorted=this->presorted;

	/*Now we need to deep copy the objects: */
	for ( obj=this->objects.begin() ; obj < this->objects.end(); obj++ ){
		/*Call copy on object: */
		object_copy = (*obj)->copy();
		copy->AddObject(object_copy);
	}

	/*Build id_offsets and sorted_ids*/
	int objsize = this->numsorted;
	if(this->sorted && objsize>0 && this->id_offsets){	
		/*Allocate new ids*/
		copy->id_offsets=xNew<int>(objsize);
		xMemCpy<int>(copy->id_offsets,this->id_offsets,objsize);
	}
	else copy->id_offsets=NULL;
	if(this->sorted && objsize>0 && this->sorted_ids){
		/*Allocate new ids*/
		copy->sorted_ids=xNew<int>(objsize);
		xMemCpy<int>(copy->sorted_ids,this->sorted_ids,objsize);
	}
	else copy->sorted_ids=NULL;

	return copy;
}
/*}}}*/
DataSet::~DataSet(){/*{{{*/
	clear();
	xDelete<int>(sorted_ids);
	xDelete<int>(id_offsets);
}
/*}}}*/

/*Specific methods*/
void  DataSet::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	vector<Object*>::iterator obj;
	int obj_enum=0;
	int i;

	int obj_size=0;
	if(marshallhandle->OperationNumber()!=MARSHALLING_LOAD){
		obj_size=objects.size();
	}
	else{
		/*FIXME: if the assert below does not go off, then remove else{}*/
		_assert_(this->Size()==0);
		//clear();
	}

	int object_enum = DataSetEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->sorted);
	marshallhandle->call(this->presorted);
	marshallhandle->call(this->numsorted);

	/*Now branch according to direction of marshalling: */
	if(marshallhandle->OperationNumber()!=MARSHALLING_LOAD){
		if(!(this->sorted && numsorted>0 && this->id_offsets)){
			this->sorted_ids=NULL;
			this->id_offsets=NULL;
		}
		marshallhandle->call(this->sorted_ids,numsorted);
		marshallhandle->call(this->id_offsets,numsorted);
		marshallhandle->call(obj_size);

		/*Go through our objects, and marshall them into the buffer: */
		for( obj=this->objects.begin() ; obj < this->objects.end(); obj++ ){
			obj_enum=(*obj)->ObjectEnum();
			marshallhandle->call(obj_enum);
			(*obj)->Marshall(marshallhandle);
		}
	}
	else{

		marshallhandle->call(this->sorted_ids,numsorted);
		marshallhandle->call(this->id_offsets,numsorted);
		if (!(this->sorted && numsorted>0)){
			sorted_ids=NULL;
			id_offsets=NULL;
		}
		marshallhandle->call(obj_size);

		/*This is the heart of the demashalling method. We have a buffer coming
		  in, and we are supposed to create a dataset out of it. No such thing
		  as class orientation for buffers, we need to key off the enum of each
		  object stored in the buffer. */
		for(i=0;i<obj_size;i++){

			/*Recover enum of object first: */
			marshallhandle->call(obj_enum); 

			/*Giant case statement to spin-up the right object, and demarshall into it the information 
			 *stored in the buffer: */
			if(obj_enum==NodeEnum){
				Node* node=new Node();
				node->Marshall(marshallhandle);
				this->AddObject(node);
			}
			else if(obj_enum==VertexEnum){
				Vertex* vertex=new Vertex();
				vertex->Marshall(marshallhandle);
				this->AddObject(vertex);
			}
			else if(obj_enum==MaticeEnum){
				Matice* matice=new Matice();
				matice->Marshall(marshallhandle);
				this->AddObject(matice);
			}
			else if(obj_enum==MatestarEnum){
				Matestar* matestar=new Matestar();
				matestar->Marshall(marshallhandle);
				this->AddObject(matestar);
			}
			else if(obj_enum==SpcStaticEnum){
				SpcStatic* spcstatic=new SpcStatic();
				spcstatic->Marshall(marshallhandle);
				this->AddObject(spcstatic);
			}
			else if(obj_enum==SpcDynamicEnum){
				SpcDynamic* spcdynamic=new SpcDynamic();
				spcdynamic->Marshall(marshallhandle);
				this->AddObject(spcdynamic);
			}
			else if(obj_enum==SpcTransientEnum){
				SpcTransient* spctransient=new SpcTransient();
				spctransient->Marshall(marshallhandle);
				this->AddObject(spctransient);
			}
			else if(obj_enum==TriaEnum){
				Tria* tria=new Tria();
				tria->Marshall(marshallhandle);
				this->AddObject(tria);
			}
			else if(obj_enum==PentaEnum){
				Penta* penta=new Penta();
				penta->Marshall(marshallhandle);
				this->AddObject(penta);
			}
			else if(obj_enum==TetraEnum){
				Tetra* tetra=new Tetra();
				tetra->Marshall(marshallhandle);
				this->AddObject(tetra);
			}
			else if(obj_enum==SegEnum){
				Seg* seg=new Seg();
				seg->Marshall(marshallhandle);
				this->AddObject(seg);
			}
			else if(obj_enum==RiftfrontEnum){
				Riftfront* rift=new Riftfront();
				rift->Marshall(marshallhandle);
				this->AddObject(rift);
			}
			else if(obj_enum==NumericalfluxEnum){
				Numericalflux* numflux=new Numericalflux();
				numflux->Marshall(marshallhandle);
				this->AddObject(numflux);
			}
			else if(obj_enum==PengridEnum){
				Pengrid* pengrid=new Pengrid();
				pengrid->Marshall(marshallhandle);
				this->AddObject(pengrid);
			}
			else if(obj_enum==PenpairEnum){
				Penpair* penpair=new Penpair();
				penpair->Marshall(marshallhandle);
				this->AddObject(penpair);
			}
			else if(obj_enum==DoubleExternalResultEnum){
				GenericExternalResult<double>* result=new GenericExternalResult<double>();
				result->Marshall(marshallhandle);
				this->AddObject(result);
			}
			else if(obj_enum==DependentObjectEnum){
				DependentObject* dep=new DependentObject();
				dep->Marshall(marshallhandle);
				this->AddObject(dep);
			}
			else if(obj_enum==GenericExternalResultEnum){
				_printf_("   WARNING: Could not load GenericExternalResult, need overload\n");
			}
			else if(obj_enum==IntExternalResultEnum){
				GenericExternalResult<int>* res=new GenericExternalResult<int>();
				res->Marshall(marshallhandle);
				this->AddObject(res);
			}
			else if(obj_enum==DoubleExternalResultEnum){
				GenericExternalResult<double>* res=new GenericExternalResult<double>();
				res->Marshall(marshallhandle);
				this->AddObject(res);
			}
			else if(obj_enum==CflevelsetmisfitEnum){
				Cflevelsetmisfit* Cflevelset=new Cflevelsetmisfit();
				Cflevelset->Marshall(marshallhandle);
				this->AddObject(Cflevelset);
			}
			else if(obj_enum==CfsurfacesquaretransientEnum){
				Cfsurfacesquaretransient* cfsurf=new Cfsurfacesquaretransient();
				cfsurf->Marshall(marshallhandle);
				this->AddObject(cfsurf);
			}
			else if(obj_enum==CfsurfacesquareEnum){
				Cfsurfacesquare* cfsurf=new Cfsurfacesquare();
				cfsurf->Marshall(marshallhandle);
				this->AddObject(cfsurf);
			}
			else if(obj_enum==CfsurfacelogvelEnum){
				Cfsurfacelogvel* cfsurf=new Cfsurfacelogvel();
				cfsurf->Marshall(marshallhandle);
				this->AddObject(cfsurf);
			}
			else if(obj_enum==CfdragcoeffabsgradEnum){
				Cfdragcoeffabsgrad* cfdragcoeff=new Cfdragcoeffabsgrad();
				cfdragcoeff->Marshall(marshallhandle);
				this->AddObject(cfdragcoeff);
			}
			else if(obj_enum==CfdragcoeffabsgradtransientEnum){
				Cfdragcoeffabsgradtransient* cfdragcoeff=new Cfdragcoeffabsgradtransient();
				cfdragcoeff->Marshall(marshallhandle);
				this->AddObject(cfdragcoeff);
			}
			else if(obj_enum==CfrheologybbarabsgradEnum){
				Cfrheologybbarabsgrad* cfrheologybbarabsgrad=new Cfrheologybbarabsgrad();
				cfrheologybbarabsgrad->Marshall(marshallhandle);
				this->AddObject(cfrheologybbarabsgrad);
			}
			else if(obj_enum==NodalvalueEnum){
				Nodalvalue* nodalvalue=new Nodalvalue();
				nodalvalue->Marshall(marshallhandle);
				this->AddObject(nodalvalue);
			}
			else if(obj_enum==MassfluxatgateEnum){
				Massfluxatgate<IssmDouble>* massfluxgate=new Massfluxatgate<IssmDouble>();
				massfluxgate->Marshall(marshallhandle);
				this->AddObject(massfluxgate);
			}
			else if(obj_enum==ChannelEnum){
				Channel* channel=new Channel();
				channel->Marshall(marshallhandle);
				this->AddObject(channel);
			}
			else if(obj_enum==MoulinEnum){
				Moulin* moulin=new Moulin();
				moulin->Marshall(marshallhandle);
				this->AddObject(moulin);
			}
			else if(obj_enum==NeumannfluxEnum){
				Neumannflux* neumannflux=new Neumannflux();
				neumannflux->Marshall(marshallhandle);
				this->AddObject(neumannflux);
			}
			else _error_("could not recognize enum type: " << obj_enum << ": " << EnumToStringx(obj_enum) ); 
		}
	}
}
/*}}}*/
int   DataSet::AddObject(Object* object){/*{{{*/

	_assert_(this);
	objects.push_back(object);

	return 1;
}
/*}}}*/
void  DataSet::clear(){/*{{{*/

/*  use reverse_iterator for efficiency in matlab memory manager
	(keeping old code in case it needs to revert back)  */

//	vector<Object*>::iterator object;
	vector<Object*>::reverse_iterator object;

//	for ( object=objects.begin() ; object < objects.end(); object++ ){
//		delete (*object);
//	}
	for ( object=objects.rbegin() ; object < objects.rend(); object++ ){
		delete (*object);
	}
	objects.clear();
}
/*}}}*/
int   DataSet::DeleteObject(Object* object){/*{{{*/

	vector<Object*>::iterator iterator;

	if(object){
		iterator = find(objects.begin(), objects.end(),object);
		delete *iterator;
		objects.erase(iterator);
	}

	return 1;

}
/*}}}*/
void  DataSet::DeepEcho(){/*{{{*/

	vector<Object*>::iterator object;

	_assert_(this);

	_printf0_("DataSet echo: " << objects.size() << " objects\n");

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		/*Call deep echo on object: */
		(*object)->DeepEcho();

	}
}
/*}}}*/
void  DataSet::Echo(){/*{{{*/

	vector<Object*>::iterator object;

	_assert_(this);

	_printf0_("DataSet echo: " << objects.size() << " objects\n");

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		/*Call echo on object: */
		(*object)->Echo();

	}
	return;
}
/*}}}*/
int   DataSet::GetEnum(){/*{{{*/
	return enum_type;
}
/*}}}*/
int   DataSet::GetEnum(int offset){/*{{{*/

	return objects[offset]->ObjectEnum();

}
/*}}}*/
Object* DataSet::GetObjectByOffset(int offset){/*{{{*/

	if(this->Size()<=offset) this->Echo();

	/*Check index in debugging mode*/
	_assert_(this!=NULL);
	_assert_(offset>=0);
	_assert_(offset<this->Size());

	return objects[offset];

}
/*}}}*/
Object* DataSet::GetObjectById(int* poffset,int eid){/*{{{*/

	int id_offset;
	int offset;

	_assert_(this);
	if(!sorted || objects.size()>numsorted)_error_("trying to binary search on a non-sorted dataset!");

	/*Carry out a binary search on the sorted_ids: */
	if(!binary_search(&id_offset,eid,sorted_ids,objects.size())){
		_error_("could not find object with id " << eid << " in DataSet " << EnumToStringx(enum_type));
	}

	/*Convert  the id offset into sorted offset: */
	offset=id_offsets[id_offset];

	/*Assign output pointers if requested:*/
	if(poffset)*poffset=offset;

	/*Return object at offset position in objects :*/
	return objects[offset];
}
/*}}}*/
void  DataSet::Presort(){/*{{{*/

	/*vector of objects is already sorted, just allocate the sorted ids and their
	 * offsets:*/
	if(objects.size()){

		/*Delete existing ids*/
		if(sorted_ids) xDelete<int>(sorted_ids);
		if(id_offsets) xDelete<int>(id_offsets);

		/*Allocate new ids*/
		sorted_ids=xNew<int>(objects.size());
		id_offsets=xNew<int>(objects.size());

		/*Build id_offsets and sorted_ids*/
		for(int i=0;i<objects.size();i++){
			id_offsets[i]=i;
			sorted_ids[i]=objects[i]->Id();

			/*In debugging mode, make sure Ids are ACTUALLY sorted...*/
			#ifdef _ISSM_DEBUG_
			//if(i>0) _assert_(sorted_ids[i]>sorted_ids[i-1]);
			#endif
		}
	}

	/*set sorted flag: */
	numsorted=objects.size();
	sorted=1;
}
/*}}}*/
int   DataSet::Size(void){/*{{{*/
	_assert_(this!=NULL);

	return objects.size();
}
/*}}}*/
void  DataSet::Sort(){/*{{{*/

	/*Only sort if we are not already sorted: */
	if(!sorted){
		_error_("not implemented yet!");
	}
}
/*}}}*/
