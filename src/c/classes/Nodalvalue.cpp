/*!\file Nodalvalue.h
 * \brief: header file for Nodalvalue object
 */
#ifdef HAVE_CONFIG_H
   #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Headers:*/
/*{{{*/
#include "./Elements/Element.h"
#include "./Elements/Elements.h"
#include "./FemModel.h"
#include "../modules/SurfaceAreax/SurfaceAreax.h"
#include "../classes/Params/Parameters.h"
#include "../classes/gauss/Gauss.h"
#include "./classes.h"
/*}}}*/

		/*Nodalvalue constructors, destructors :*/
Nodalvalue::Nodalvalue(){/*{{{*/

	this->definitionenum = -1;
	this->name = NULL;
	this->model_enum = UNDEF;
	this->node = -1;

}
/*}}}*/
Nodalvalue::Nodalvalue(char* in_name, int in_definitionenum, int in_model_enum, int in_node){/*{{{*/

	this->definitionenum=in_definitionenum;
	this->name   = xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->model_enum=in_model_enum;
	this->node=in_node;
}
/*}}}*/
Nodalvalue::~Nodalvalue(){/*{{{*/
	if(this->name)xDelete(this->name);
}
/*}}}*/
/*Object virtual function resolutoin: */
Object* Nodalvalue::copy() {/*{{{*/
	Nodalvalue* mf = new Nodalvalue(this->name,this->definitionenum, this->model_enum,this->node);
	return (Object*) mf;
}
/*}}}*/
void Nodalvalue::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void Nodalvalue::Echo(void){/*{{{*/
	_printf_(" Nodalvalue: " << name << " " << this->definitionenum << "\n");
	_printf_("    model_enum: " << model_enum << " " << EnumToStringx(model_enum) << "\n");
	_printf_("    node: " << node << "\n");
}
/*}}}*/
int Nodalvalue::Id(void){/*{{{*/
	return -1;
}
/*}}}*/
void Nodalvalue::Marshall(MarshallHandle* marshallhandle){/*{{{*/

	int object_enum=NodalvalueEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->definitionenum);
	marshallhandle->call(this->model_enum);
	marshallhandle->call(this->name);
	marshallhandle->call(this->node);
} 
/*}}}*/
int Nodalvalue::ObjectEnum(void){/*{{{*/
	return NodalvalueEnum;
}
/*}}}*/
/*Definition virtual function resolutoin: */
int Nodalvalue::DefinitionEnum(){/*{{{*/

	return this->definitionenum;
}
/*}}}*/
char* Nodalvalue::Name(){/*{{{*/

	char* name2=xNew<char>(strlen(this->name)+1);
	xMemCpy(name2,this->name,strlen(this->name)+1);

	return name2;
}
/*}}}*/
IssmDouble Nodalvalue::Response(FemModel* femmodel){/*{{{*/

	 /*output:*/
	 IssmDouble value;

	 /*set index, which will be used by the NodalValue module: */
	 femmodel->parameters->SetParam(node,IndexEnum);

	 /*call Nodalvalue:*/
	 NodalValuex(&value, model_enum, femmodel->elements, femmodel->nodes, femmodel->vertices, femmodel->loads, 
			 femmodel->materials, femmodel->parameters);

	 /*done:*/
	 return value;
 }
 /*}}}*/
