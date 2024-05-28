/*!\file SpcStatic.c
 * \brief: implementation of the SpcStatic object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "./Constraint.h"
#include "../../shared/shared.h"

/*SpcStatic constructors and destructor*/
SpcStatic::SpcStatic(){/*{{{*/
	return;
}
/*}}}*/
SpcStatic::SpcStatic(int spc_id,int spc_nodeid, int spc_dof,IssmDouble spc_value,int spc_analysis_type){/*{{{*/

	id           = spc_id;
	nodeid        = spc_nodeid;
	dof           = spc_dof;
	value         = spc_value;
	analysis_type = spc_analysis_type;
	penalty       = false;

	return;
}
/*}}}*/
SpcStatic::~SpcStatic(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* SpcStatic::copy() {/*{{{*/

	SpcStatic* spcstat = new SpcStatic(*this); 

	spcstat->id=this->id;
	spcstat->nodeid=this->nodeid;
	spcstat->dof=this->dof;
	spcstat->value=this->value;
	spcstat->analysis_type=this->analysis_type;

	return (Object*) spcstat;
}
/*}}}*/
void    SpcStatic::DeepEcho(void){/*{{{*/

	_printf_("SpcStatic:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   nodeid: " << nodeid << "\n");
	_printf_("   dof: " << dof << "\n");
	_printf_("   value: " << value << "\n");
	_printf_("   analysis_type: " << EnumToStringx(analysis_type) << "\n");
	return;
}		
/*}}}*/
void    SpcStatic::Echo(void){/*{{{*/

	_printf_("SpcStatic:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   nodeid: " << nodeid << "\n");
	_printf_("   dof: " << dof << "\n");
	_printf_("   value: " << value << "\n");
	_printf_("   analysis_type: " << EnumToStringx(analysis_type) << "\n");
	return;
}
/*}}}*/
int     SpcStatic::Id(void){ return id; }/*{{{*/
/*}}}*/
void    SpcStatic::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = SpcStaticEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->id);
	marshallhandle->call(this->nodeid);
	marshallhandle->call(this->dof);
	marshallhandle->call(this->value);
	marshallhandle->call(this->analysis_type);
	marshallhandle->call(this->penalty);

}
/*}}}*/
int     SpcStatic::ObjectEnum(void){/*{{{*/

	return SpcStaticEnum;

}
/*}}}*/

/*Constraint virtual functions definitions: */
void SpcStatic::ActivatePenaltyMethod(void){/*{{{*/
	   this->penalty = true;
}
/*}}}*/
void SpcStatic::ConstrainNode(Nodes* nodes,Parameters* parameters){/*{{{*/

	/*Chase through nodes and find the node to which this SpcStatic applys: */
	Node* node=(Node*)nodes->GetObjectById(NULL,nodeid);

	/*Apply constraint: */
	if(!this->penalty && node){ //in case the spc is dealing with a node on another cpu
		node->ApplyConstraint(dof,value);
	}
}
/*}}}*/
void SpcStatic::InputUpdateFromVectorDakota(IssmDouble* vector,Nodes* nodes,int name,int type){/*{{{*/

	/*Only update if this is a constraint parameter*/
	if(name != BalancethicknessSpcthicknessEnum) return;

	/*Chase through nodes and find the node to which this SpcStatic applies: */
	Node* node=(Node*)nodes->GetObjectById(NULL,nodeid);

	/*Apply constraint: */
	if(node){ //in case the spc is dealing with a node on another cpu
		int sid = node->Sid();
		this->value = vector[sid];
		_assert_(!xIsNan<IssmDouble>(this->value)); 
	}
}
/*}}}*/
void SpcStatic::PenaltyDofAndValue(int* pdof,IssmDouble* pvalue,Nodes* nodes,Parameters* parameters){/*{{{*/

	if(!this->penalty) _error_("cannot return dof and value for non penalty constraint");

	IssmDouble value_out = this->value;
	int gdof;

	/*Chase through nodes and find the node to which this SpcTransient applys: */
	Node* node=(Node*)nodes->GetObjectById(NULL,nodeid);

	if(node){ //in case the spc is dealing with a node on another cpu

		/*Get gdof */
		gdof = node->GetDof(dof,GsetEnum);
		if(xIsNan<IssmDouble>(value_out)) gdof = -1;
	}
	else{
		value_out = NAN;
		gdof = -1;
	}

	/*Assign output pointers*/
	*pdof   = gdof;
	*pvalue = value_out;
}
/*}}}*/

void SpcStatic::UpdateSpcThicknessAD(IssmDouble* vector,Nodes* nodes){/*{{{*/

	/*Chase through nodes and find the node to which this SpcStatic applies: */
	Node* node=(Node*)nodes->GetObjectById(NULL,nodeid);

	/*Apply constraint: */
	if(node){ //in case the spc is dealing with a node on another cpu
		int sid = node->Sid();
		this->value = vector[sid];
		_assert_(!xIsNan<IssmDouble>(this->value)); 
	}
}
/*}}}*/

/*SpcStatic functions*/
int        SpcStatic::GetDof(){/*{{{*/
	return dof;
}
/*}}}*/
int        SpcStatic::GetNodeId(){/*{{{*/

	return nodeid;
}
/*}}}*/
IssmDouble SpcStatic::GetValue(){/*{{{*/
	_assert_(!xIsNan<IssmDouble>(value));
	return value;
}
/*}}}*/
