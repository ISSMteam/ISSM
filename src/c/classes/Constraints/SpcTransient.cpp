/*!\file SpcTransient.c
 * \brief: implementation of the SpcTransient object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "./Constraint.h"
#include "shared/shared.h"

/*SpcTransient constructors and destructor*/
SpcTransient::SpcTransient(){/*{{{*/
	penalty       = false;
	id            = -1;
	nodeid        = -1;
	dof           = -1;
	values        = NULL;
	times         = NULL;
	nsteps        = -1;
	analysis_type = -1;
	return;
}
/*}}}*/
SpcTransient::SpcTransient(int spc_id,int spc_nodeid, int spc_dof,int spc_nsteps, IssmDouble* spc_times, IssmDouble* spc_values,int spc_analysis_type){/*{{{*/

	penalty = false;
	id     = spc_id;
	nodeid  = spc_nodeid;
	dof     = spc_dof;
	nsteps  = spc_nsteps;
	if(spc_nsteps){
		values = xNew<IssmDouble>(spc_nsteps);
		times  = xNew<IssmDouble>(spc_nsteps);
		xMemCpy<IssmDouble>(values,spc_values,nsteps);
		xMemCpy<IssmDouble>(times,spc_times,nsteps);
	}
	analysis_type=spc_analysis_type;
	return;
}
/*}}}*/
SpcTransient::~SpcTransient(){/*{{{*/
	xDelete<IssmDouble>(times);
	xDelete<IssmDouble>(values);
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* SpcTransient::copy() {/*{{{*/
	return new SpcTransient(id,nodeid,dof,nsteps,times,values,analysis_type);
}
/*}}}*/
void    SpcTransient::DeepEcho(void){/*{{{*/
	this->Echo();
}		
/*}}}*/
void    SpcTransient::Echo(void){/*{{{*/

	int i;
	_printf_("SpcTransient:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   nodeid: " << nodeid << "\n");
	_printf_("   dof: " << dof << "\n");
	_printf_("   nsteps: " << nsteps << "\n");
	_printf_("   analysis_type: " << EnumToStringx(analysis_type) << "\n");
	_printf_("   steps|times|values\n");
	for(i=0;i<nsteps;i++){
		_printf_(i << "-" << times[i] << ":" << values[i] << "\n");
	}
	return;
}
/*}}}*/
int     SpcTransient::Id(void){/*{{{*/
	return id;
}
/*}}}*/
void    SpcTransient::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = SpcTransientEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->id);
	marshallhandle->call(this->nodeid);
	marshallhandle->call(this->dof);
	marshallhandle->call(this->analysis_type);
	marshallhandle->call(this->penalty);
	marshallhandle->call(this->nsteps);
	if(nsteps){
		marshallhandle->call(this->values,nsteps);
		marshallhandle->call(this->times,nsteps);
	}
	else{
		this->values=NULL;
		this->times=NULL;
	}
}/*}}}*/
int     SpcTransient::ObjectEnum(void){/*{{{*/

	return SpcTransientEnum;

}
/*}}}*/

/*Constraint virtual functions definitions:*/
void SpcTransient::ActivatePenaltyMethod(void){/*{{{*/
	   this->penalty = true;
}
/*}}}*/
void SpcTransient::ConstrainNode(Nodes* nodes,Parameters* parameters){/*{{{*/

	Node       *node  = NULL;
	IssmDouble  time  = 0.;
	int         i;
	IssmDouble  alpha = -1.;
	IssmDouble  value;
	bool        found = false;

	/*Chase through nodes and find the node to which this SpcTransient applys: */
	node=(Node*)nodes->GetObjectById(NULL,nodeid);

	if(!this->penalty && node){ //in case the spc is dealing with a node on another cpu

		/*Retrieve time in parameters: */
		parameters->FindParam(&time,TimeEnum);

		/*Now, go fetch value for this time: */
		if (time<=times[0]){
			value=values[0];
			found=true;
		}
		else if (time>=times[nsteps-1]){
			value=values[nsteps-1];
			found=true;
		}
		else{
			for(i=0;i<nsteps-1;i++){
				if (times[i]<=time && time<times[i+1]){
					alpha=(time-times[i])/(times[i+1]-times[i]);
					value=(1-alpha)*values[i]+alpha*values[i+1];
					found=true;
					break;
				}
			}
		}

		if(!found)_error_("could not find time segment for constraint");

		/*Apply or relax constraint: */
		if(xIsNan<IssmDouble>(value)){
			node->RelaxConstraint(dof);
		}
		else node->ApplyConstraint(dof,value);
	}
}
/*}}}*/
void SpcTransient::PenaltyDofAndValue(int* pdof,IssmDouble* pvalue,Nodes* nodes,Parameters* parameters){/*{{{*/

	if(!this->penalty) _error_("cannot return dof and value for non penalty constraint");

	Node       *node  = NULL;
	IssmDouble  time  = 0.;
	int         i,gdof;
	IssmDouble  alpha = -1.;
	IssmDouble  value;
	bool        found = false;

	/*Chase through nodes and find the node to which this SpcTransient applys: */
	node=(Node*)nodes->GetObjectById(NULL,nodeid);

	if(node){ //in case the spc is dealing with a node on another cpu

		/*Retrieve time in parameters: */
		parameters->FindParam(&time,TimeEnum);

		/*Now, go fetch value for this time: */
		if (time<=times[0]){
			value=values[0];
			found=true;
		}
		else if (time>=times[nsteps-1]){
			value=values[nsteps-1];
			found=true;
		}
		else{
			for(i=0;i<nsteps-1;i++){
				if (times[i]<=time && time<times[i+1]){
					alpha=(time-times[i])/(times[i+1]-times[i]);
					value=(1-alpha)*values[i]+alpha*values[i+1];
					found=true;
					break;
				}
			}
		}
		if(!found)_error_("could not find time segment for constraint");

		/*Get gdof */
		gdof = node->GetDof(dof,GsetEnum);
		if(xIsNan<IssmDouble>(value)){
			gdof = -1;
		}
	}
	else{
		value = NAN;
		gdof = -1;
	}

	/*Assign output pointers*/
	*pdof   = gdof;
	*pvalue = value;
}
/*}}}*/

/*SpcTransient functions*/
int        SpcTransient::GetDof(){/*{{{*/
	return dof;
}
/*}}}*/
int        SpcTransient::GetNodeId(){/*{{{*/

	return nodeid;
}
/*}}}*/
IssmDouble SpcTransient::GetValue(){/*{{{*/
	return values[0];
}
/*}}}*/
