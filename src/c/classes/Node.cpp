/*!\file Node.c
 * \brief: implementation of the Node object
 */

/*Include files: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./classes.h"
#include "shared/shared.h"
#include "modules/ModelProcessorx/ModelProcessorx.h"
#include "../analyses/analyses.h"
/*}}}*/

/*Node constructors and destructors:*/
Node::Node(){/*{{{*/
	this->approximation  = 0;
	this->gsize          = -1;
	this->clone          = false;
	this->active         = true;
	this->freeze         = false;
	this->isrotated      = false;
	this->f_set          = NULL;
	this->s_set          = NULL;
	this->svalues        = NULL;
	this->doftype        = NULL;
	this->gdoflist       = NULL;
	this->fdoflist       = NULL;
	this->sdoflist       = NULL;
	this->gdoflist_local = NULL;
	this->fdoflist_local = NULL;
	this->sdoflist_local = NULL;
}
/*}}}*/
Node::Node(int node_id,int node_sid,int io_index,bool node_clone,IoModel* iomodel,int node_analysis,int in_approximation,bool isamr){/*{{{*/

	/*id: */
	this->id            = node_id;
	this->sid           = node_sid;
	this->lid           = -1; /*Assigned by Finalize*/
	this->pid           = -1; /*Assigned by Finalize*/
	this->analysis_enum = node_analysis;
	this->clone         = node_clone;
	this->active        = true;
	this->freeze        = false;

	/*Initialize coord_system: Identity matrix by default*/
	this->isrotated = false;
	for(int k=0;k<3;k++) for(int l=0;l<3;l++) this->coord_system[k][l]=0.0;
	for(int k=0;k<3;k++) this->coord_system[k][k]=1.0;

	this->approximation=0;
	if(analysis_enum==StressbalanceAnalysisEnum) this->approximation=in_approximation;

	/*indexing:*/
	this->indexingupdate = true;
	this->doftype        = NULL;
	Analysis *analysis = EnumToAnalysis(analysis_enum);
	this->gsize = analysis->DofsPerNode(&this->doftype,iomodel->domaintype,in_approximation);
	delete analysis;

	if(this->gsize>0){
		this->f_set          = xNew<bool>(this->gsize);
		this->s_set          = xNew<bool>(this->gsize);
		this->svalues        = xNew<IssmDouble>(this->gsize);
		this->gdoflist       = xNew<int>(this->gsize);
		this->gdoflist_local = xNew<int>(this->gsize);
		this->fdoflist       = xNew<int>(this->gsize);
		this->sdoflist       = xNew<int>(this->gsize);
		this->fdoflist_local = xNew<int>(this->gsize);
		this->sdoflist_local = xNew<int>(this->gsize);
	}
	else{
		this->f_set          = NULL;
		this->s_set          = NULL;
		this->svalues        = NULL;
		this->gdoflist       = NULL;
		this->gdoflist_local = NULL;
		this->fdoflist       = NULL;
		this->sdoflist       = NULL;
		this->fdoflist_local = NULL;
		this->sdoflist_local = NULL;
	}

	/*Assign values assuming no Dirichlet at this point*/
	for(int i=0;i<this->gsize;i++){
		this->f_set[i]    = true;
		this->s_set[i]    = false;
		this->svalues[i]  = 0.;
		this->gdoflist[i] = -1;
		this->fdoflist[i] = -1;
		this->sdoflist[i] = -1;
		this->gdoflist_local[i] = -1;
		this->fdoflist_local[i] = -1;
		this->sdoflist_local[i] = -1;
	}

	/*Stop here if AMR*/
	if(isamr) return;

	/*Stressbalance Horiz*/
	if(analysis_enum==StressbalanceAnalysisEnum){

		/*Coordinate system provided, convert to coord_system matrix*/
		_assert_(iomodel->Data("md.stressbalance.referential")); 
		this->isrotated = XZvectorsToCoordinateSystem(&this->coord_system[0][0],&iomodel->Data("md.stressbalance.referential")[io_index*6]);
		_assert_(sqrt( coord_system[0][0]*coord_system[0][0] + coord_system[1][0]*coord_system[1][0]) >1.e-4);

		if(iomodel->domaintype!=Domain2DhorizontalEnum && iomodel->domaintype!=Domain3DsurfaceEnum){
			/*We have a  3d mesh, we may have collapsed elements, hence dead nodes. Freeze them out: */
			_assert_(iomodel->Data("md.mesh.vertexonbase")); 
			_assert_(iomodel->Data("md.flowequation.vertex_equation"));
			if(in_approximation==SSAApproximationEnum && !reCast<int>(iomodel->Data("md.mesh.vertexonbase")[io_index])){
				this->HardDeactivate();
			}
			if(in_approximation==L1L2ApproximationEnum && !reCast<int>(iomodel->Data("md.mesh.vertexonbase")[io_index])){
				this->HardDeactivate();
			}
			if(in_approximation==MOLHOApproximationEnum && !reCast<int>(iomodel->Data("md.mesh.vertexonbase")[io_index])){
				this->HardDeactivate();
			}
			if(in_approximation==SSAHOApproximationEnum && reCast<int>(iomodel->Data("md.flowequation.borderSSA")[io_index])){
				if(!reCast<int>(iomodel->Data("md.mesh.vertexonbase")[io_index])){
					this->HardDeactivate();
				}
			}
			if(in_approximation==SSAFSApproximationEnum && reCast<int>(iomodel->Data("md.flowequation.borderSSA")[io_index])){
				if(!reCast<int>(iomodel->Data("md.mesh.vertexonbase")[io_index])){
					for(int k=0;k<=1;k++) this->FreezeDof(k);
				}
			}
		}
		/*spc all nodes on SIA*/
		if(in_approximation==SIAApproximationEnum){
			this->HardDeactivate();
		}
	}

	/*2d solutions in 3d, we need to constrain all the nodes that are not on base*/
	if(
				analysis_enum==FreeSurfaceBaseAnalysisEnum || 
				analysis_enum==MasstransportAnalysisEnum || 
				analysis_enum==MeltingAnalysisEnum || 
				analysis_enum==L2ProjectionBaseAnalysisEnum || 
				analysis_enum==BalancethicknessAnalysisEnum ||
				analysis_enum==HydrologyDCInefficientAnalysisEnum ||
				analysis_enum==HydrologyDCEfficientAnalysisEnum ||
				analysis_enum==HydrologyShaktiAnalysisEnum ||
				analysis_enum==HydrologyGlaDSAnalysisEnum ||
				analysis_enum==GLheightadvectionAnalysisEnum ||
				analysis_enum==LevelsetAnalysisEnum
	  ){
		if(iomodel->domaintype!=Domain2DhorizontalEnum & iomodel->domaintype!=Domain3DsurfaceEnum){
			/*On a 3d mesh, we may have collapsed elements, hence dead nodes. Freeze them out: */
			_assert_(iomodel->Data("md.mesh.vertexonbase"));
			if(!(reCast<bool>(iomodel->Data("md.mesh.vertexonbase")[io_index]))){
				this->HardDeactivate();
			}
		}
	}
	if(
			analysis_enum==FreeSurfaceTopAnalysisEnum ||
			analysis_enum==DebrisAnalysisEnum
	  ){
		if(iomodel->domaintype!=Domain2DhorizontalEnum){
			/*On a 3d mesh, we may have collapsed elements, hence dead nodes. Freeze them out: */
			_assert_(iomodel->Data("md.mesh.vertexonsurface"));
			if(!(reCast<bool>(iomodel->Data("md.mesh.vertexonsurface")[io_index]))){
				this->HardDeactivate();
			}
		}
	}

}
/*}}}*/
Node::~Node(){/*{{{*/

	if(this->f_set)          xDelete<bool>(f_set);
	if(this->s_set)          xDelete<bool>(s_set);
	if(this->svalues)        xDelete<IssmDouble>(svalues);
	if(this->doftype)        xDelete<int>(doftype);
	if(this->gdoflist)       xDelete<int>(gdoflist);
	if(this->fdoflist)       xDelete<int>(fdoflist);
	if(this->sdoflist)       xDelete<int>(sdoflist);
	if(this->gdoflist_local) xDelete<int>(gdoflist_local);
	if(this->fdoflist_local) xDelete<int>(fdoflist_local);
	if(this->sdoflist_local) xDelete<int>(sdoflist_local);
	return;
}
/*}}}*/
Object* Node::copy(void){/*{{{*/

	/*output: */
	Node* output=NULL;

	/*initalize output: */
	output=new Node();

	output->approximation = this->approximation;
	output->clone  = this->clone;
	output->id  = this->id;
	output->sid = this->sid;
	output->lid = this->lid;
	output->pid = this->pid;

	output->analysis_enum  = this->analysis_enum;
	output->indexingupdate = this->indexingupdate;
	output->isrotated      = this->isrotated;

	/*Initialize coord_system: */
	for(int k=0;k<3;k++) for(int l=0;l<3;l++) output->coord_system[k][l]=this->coord_system[k][l];

	/*indexing:*/
	output->gsize  = this->gsize;
	output->active = this->active;
	output->freeze = this->freeze;
	if(output->gsize>0){
		output->f_set=xNew<bool>(output->gsize);
		output->s_set=xNew<bool>(output->gsize);
		output->svalues=xNew<IssmDouble>(output->gsize);
		if(this->doftype) output->doftype=xNew<int>(output->gsize);
		output->gdoflist=xNew<int>(output->gsize);
		output->gdoflist_local=xNew<int>(output->gsize);
		output->fdoflist=xNew<int>(output->gsize);
		output->fdoflist_local=xNew<int>(output->gsize);
		output->sdoflist=xNew<int>(output->gsize);
		output->sdoflist_local=xNew<int>(output->gsize);
	}

	if(output->gsize>0){
		memcpy(output->f_set,this->f_set,output->gsize*sizeof(bool));
		memcpy(output->s_set,this->s_set,output->gsize*sizeof(bool));
		xMemCpy<IssmDouble>(output->svalues,this->svalues,output->gsize);
		if(output->doftype)memcpy(output->doftype,this->doftype,output->gsize*sizeof(int));
		memcpy(output->gdoflist,this->gdoflist,output->gsize*sizeof(int));
		memcpy(output->gdoflist_local,this->gdoflist_local,output->gsize*sizeof(int));
		memcpy(output->fdoflist,this->fdoflist,output->gsize*sizeof(int));
		memcpy(output->fdoflist_local,this->fdoflist_local,output->gsize*sizeof(int));
		memcpy(output->sdoflist,this->sdoflist,output->gsize*sizeof(int));
		memcpy(output->sdoflist_local,this->sdoflist_local,output->gsize*sizeof(int));
	}

	return (Object*)output; 
}
/*}}}*/
void Node::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = NodeEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->id);
	marshallhandle->call(this->sid);
	marshallhandle->call(this->lid);
	marshallhandle->call(this->pid);
	marshallhandle->call(this->indexingupdate);
	marshallhandle->call(this->analysis_enum);

	for(int k=0;k<3;k++) for(int l=0;l<3;l++) marshallhandle->call(this->coord_system[k][l]);

	marshallhandle->call(this->gsize);
	marshallhandle->call(this->clone);
	marshallhandle->call(this->active);
	marshallhandle->call(this->freeze);
	marshallhandle->call(this->f_set,gsize);
	marshallhandle->call(this->s_set,gsize);
	marshallhandle->call(this->svalues,gsize);
	marshallhandle->call(this->doftype,gsize);
	marshallhandle->call(this->gdoflist,gsize);
	marshallhandle->call(this->fdoflist,gsize);
	marshallhandle->call(this->sdoflist,gsize);
	marshallhandle->call(this->gdoflist_local,gsize);
	marshallhandle->call(this->fdoflist_local,gsize);
	marshallhandle->call(this->sdoflist_local,gsize);
} /*}}}*/

/*Object virtual functions definitions:*/
void Node::DeepEcho(void){/*{{{*/

	int i;
	_printf_("Node:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   sid: " << sid << "\n");
	_printf_("   analysis_enum: " << EnumToStringx(analysis_enum) << "\n");
	_printf_("   approximation: " << EnumToStringx(approximation) << "\n");
	_printf_("   indexingupdate: " << indexingupdate << "\n");
	_printf_("   gsize:  " << gsize << "\n");
	_printf_("   clone:  " << clone << "\n");
	_printf_("   active: " << active << "\n");
	_printf_("   freeze: " << freeze << "\n");
	_printf_("   f_set = [ ");
	for(i=0;i<gsize;i++) _printf_((f_set[i]?1:0)<< " ");
	_printf_("]\n");
	_printf_("   s_set = [ ");
	for(i=0;i<gsize;i++) _printf_((s_set[i]?1:0)<< " ");
	_printf_("]\n");
	_printf_("   svalues: |");
	for(i=0;i<this->gsize;i++){
		if(this->s_set[i])_printf_(" " << svalues[i] << " |");
	}
	_printf_("\n");
	if(doftype){
		_printf_("   doftype: |");
		for(i=0;i<gsize;i++){
			_printf_(" " << doftype[i] << " |");
		}
		_printf_("\n");
	}
	else _printf_("   doftype: NULL\n");

	_printf_("   g_doflist (" << this->gsize << "): |");
	for(i=0;i<this->gsize;i++) _printf_(" " << gdoflist[i] << " |");
	_printf_("\n");
	_printf_("   g_doflist_local (" << this->gsize << "): |");
	for(i=0;i<this->gsize;i++) _printf_(" " << gdoflist_local[i] << " |");
	_printf_("\n");

	_printf_("   f_doflist (" << this->gsize << "): |");
	for(i=0;i<this->gsize;i++) _printf_(" " << fdoflist[i] << " |");
	_printf_("\n");
	_printf_("   f_doflist_local (" << this->gsize << "): |");
	for(i=0;i<this->gsize;i++) _printf_(" " << fdoflist_local[i] << " |");
	_printf_("\n");

	_printf_("   s_doflist (" << this->gsize << "): |");
	for(i=0;i<this->gsize;i++) _printf_(" " << sdoflist[i] << " |");
	_printf_("\n");
	_printf_("   s_doflist_local (" << this->gsize << "): |");
	for(i=0;i<this->gsize;i++) _printf_(" " << sdoflist_local[i] << " |");
	_printf_("\n");

}
/*}}}*/
void Node::Echo(void){/*{{{*/

	_printf_("Node:\n");
	_printf_("   id : " << id << "\n");
	_printf_("   sid: " << sid << "\n");
	_printf_("   lid: " << lid << "\n");
	_printf_("   pid: " << pid << "\n");
	_printf_("   analysis_enum: " << EnumToStringx(analysis_enum) << "\n");
	_printf_("   approximation: " << EnumToStringx(approximation) << "\n");
	_printf_("   indexingupdate: " << indexingupdate << "\n");
	_printf_("   gsize:  " << gsize << "\n");
	_printf_("   clone:  " << clone << "\n");
	_printf_("   active: " << active << "\n");
	_printf_("   freeze: " << freeze << "\n");
}
/*}}}*/
int  Node::Id(void){ return id; }/*{{{*/
/*}}}*/
int  Node::ObjectEnum(void){/*{{{*/

	return NodeEnum;

}
/*}}}*/

/*Node management:*/
void Node::GetCoordinateSystem(IssmDouble* coord_system_out){/*{{{*/

	/*Copy coord_system*/
	for(int k=0;k<3;k++) for(int l=0;l<3;l++) coord_system_out[3*k+l]=this->coord_system[k][l];

}
/*}}}*/
int  Node::GetDof(int dofindex,int setenum){/*{{{*/

	_assert_(!this->indexingupdate);
	if(setenum==GsetEnum){
		_assert_(dofindex>=0 && dofindex<gsize);
		return gdoflist[dofindex];
	}
	else if(setenum==FsetEnum){
		_assert_(dofindex>=0 && dofindex<gsize);
		return fdoflist[dofindex];
	}
	else if(setenum==SsetEnum){
		_assert_(dofindex>=0 && dofindex<gsize);
		return sdoflist[dofindex];
	}
	else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");

} /*}}}*/
void Node::GetDofList(int* outdoflist,int approximation_enum,int setenum){/*{{{*/
	_assert_(!this->indexingupdate);
	int i;

	int* doflistpointer = NULL;
	if(setenum==GsetEnum) doflistpointer = gdoflist;
	else if(setenum==FsetEnum)for(i=0;i<this->gsize;i++) doflistpointer = fdoflist;
	else if(setenum==SsetEnum)for(i=0;i<this->gsize;i++) doflistpointer = sdoflist;
	else _error_("not supported");

	if(approximation_enum==NoneApproximationEnum){
		for(i=0;i<this->gsize;i++) outdoflist[i]=doflistpointer[i];
	}
	else{
		if(doftype){
			int count = 0;
			for(i=0;i<this->gsize;i++){
				if(doftype[i]==approximation_enum) outdoflist[count++]=doflistpointer[i];
			}
		}
		else for(i=0;i<this->gsize;i++) outdoflist[i]=doflistpointer[i];
	}
}/*}}}*/
void Node::GetDofListLocal(int* outdoflist,int approximation_enum,int setenum){/*{{{*/

	_assert_(!this->indexingupdate);
	int i;

	int* doflistpointer = NULL;
	if(setenum==GsetEnum) doflistpointer = gdoflist_local;
	else if(setenum==FsetEnum)for(i=0;i<this->gsize;i++) doflistpointer = fdoflist_local;
	else if(setenum==SsetEnum)for(i=0;i<this->gsize;i++) doflistpointer = sdoflist_local;
	else _error_("not supported");

	if(approximation_enum==NoneApproximationEnum){
		for(i=0;i<this->gsize;i++) outdoflist[i]=doflistpointer[i];
	}
	else{
		if(doftype){
			int count =0;
			for(i=0;i<this->gsize;i++){
				if(doftype[i]==approximation_enum) outdoflist[count++]=doflistpointer[i];
			}
		}
		else for(i=0;i<this->gsize;i++) outdoflist[i]=doflistpointer[i];
	}
}
/*}}}*/
int  Node::Lid(void){/*{{{*/
	return lid; 
}
/*}}}*/
int  Node::Sid(void){/*{{{*/
	return sid; 
}
/*}}}*/
int  Node::Pid(void){/*{{{*/
	return this->pid; 
}
/*}}}*/

/*Node numerics:*/
void Node::Activate(void){/*{{{*/

	if(!IsActive() && !this->freeze){
		this->indexingupdate = true;
		this->active = true;
		for(int i=0;i<this->gsize;i++){
			this->f_set[i]    = true;
			this->s_set[i]    = false;
			this->svalues[i]  = 0.; 
		}
	}

}
/*}}}*/
void Node::ApplyConstraint(int dof,IssmDouble value){/*{{{*/

	/*Dof should be added in the s set, describing which 
	 * dofs are constrained to a certain value (dirichlet boundary condition*/
	DofInSSet(dof);
	this->svalues[dof]=value;
}
/*}}}*/
void Node::CreateNodalConstraints(Vector<IssmDouble>* ys){/*{{{*/

	if(this->SSize()){
		/*Add values into constraint vector: */
		ys->SetValues(this->gsize,this->sdoflist,this->svalues,INS_VAL);
	}

}/*}}}*/
void Node::Deactivate(void){/*{{{*/

	if(IsActive() && !this->freeze){
		this->indexingupdate = true;
		this->active = false;
		/*Constrain to 0. at this point*/
		for(int i=0;i<this->gsize;i++){
			this->f_set[i]    = false;
			this->s_set[i]    = true;
			this->svalues[i]  = 0.; 
		}
	}
} /*}}}*/
void Node::DofInFSet(int dof){/*{{{*/

	/*Put dof for this node into the f set (ie, this dof will NOT be constrained 
	 * to a fixed value during computations. Only do this for active nodes. */
	_assert_(dof<this->gsize);
	_assert_(this->active);

	if(!this->f_set[dof]){
		if(this->freeze) _error_("Cannot change dof of frozen node");
		this->indexingupdate = true;
		this->f_set[dof]=true; 
		this->s_set[dof]=false;
	}
}
/*}}}*/
void Node::DofInSSet(int dof){/*{{{*/

	/*Put dof for this node into the s set (ie, this dof will be constrained 
	 * to a fixed value during computations. */
	_assert_(dof<this->gsize);

	if(this->f_set[dof]){
		//if(this->freeze) _error_("Cannot change dof of frozen node");
		this->indexingupdate = true;
		this->f_set[dof]=false; //n splits into f (for which we solve) and s (single point constraints)
		this->s_set[dof]=true;
	}
}
/*}}}*/
void Node::FreezeDof(int dof){/*{{{*/

	DofInSSet(dof); //with 0 displacement for this dof.
	//FIXME: for now we don't want this element to change so we use freeze
	this->freeze =true;

}
/*}}}*/
int  Node::GetApproximation(){/*{{{*/

	return approximation;
}
/*}}}*/
void Node::SetApproximation(int in_approximation){/*{{{*/

	this->approximation = in_approximation;
}
/*}}}*/
int  Node::GetNumberOfDofs(int approximation_enum,int setenum){/*{{{*/

	/*Get number of degrees of freedom in a node, for a certain set (g,f or s-set)
	 *and for a certain approximation type: */

	int i;
	int numdofs=0;

	if(approximation_enum==NoneApproximationEnum){
		if      (setenum==GsetEnum) numdofs=this->gsize;
		else if (setenum==FsetEnum) numdofs=this->FSize();
		else if (setenum==SsetEnum) numdofs=this->SSize();
		else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
	}
	else{
		if(setenum==GsetEnum){
			if(this->doftype){
				numdofs=0;
				for(i=0;i<this->gsize;i++){
					if(this->doftype[i]==approximation_enum) numdofs++;
				}
			}
			else numdofs=this->gsize;
		}
		else if (setenum==FsetEnum){
			if(this->doftype){
				numdofs=0;
				for(i=0;i<this->gsize;i++){
					if((this->doftype[i]==approximation_enum) && (this->f_set[i])) numdofs++;
				}
			}
			else numdofs=this->FSize();
		}
		else if (setenum==SsetEnum){
			if(this->doftype){
			numdofs=0;
				for(i=0;i<this->gsize;i++){
					if((this->doftype[i]==approximation_enum) && (this->s_set[i])) numdofs++;
				}
			}
			else numdofs=this->SSize();
		}
		else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
	}
	return numdofs;
}
/*}}}*/
void Node::HardDeactivate(void){/*{{{*/

	this->Deactivate();
	this->freeze =true;

}
/*}}}*/
bool Node::IsActive(void){/*{{{*/
	return active;
}/*}}}*/
int  Node::IsClone(){/*{{{*/
	return clone;
}/*}}}*/
void Node::ReindexingDone(void){/*{{{*/
	this->indexingupdate = false;
}/*}}}*/
void Node::RelaxConstraint(int dof){/*{{{*/

	/*Dof should be added to the f-set, and taken out of the s-set:*/
	if(this->active){
		DofInFSet(dof);
		this->svalues[dof]=0.;
	}
}
/*}}}*/
bool Node::RequiresDofReindexing(void){/*{{{*/

	return this->indexingupdate;

}
/*}}}*/
void Node::VecMerge(Vector<IssmDouble>* ug,IssmDouble* local_uf,int* indices_uf,IssmDouble* local_ys,int* indices_ys){/*{{{*/

	/*Only perform operation if not clone*/
	if(this->IsClone()) return;

	/*Get F size and S size*/
	int fsize = this->FSize();
	int ssize = this->SSize();

	if(fsize){
		int*        indices = xNew<int>(fsize);
		IssmDouble* values  = xNew<IssmDouble>(fsize);

		int count = 0;
		for(int i=0;i<this->gsize;i++){
			if(this->f_set[i]){
				_assert_(local_uf);
				_assert_(this->fdoflist[i]==indices_uf[this->fdoflist_local[i]]);

				values[count] =local_uf[this->fdoflist_local[i]];
				indices[count]=this->gdoflist[i];
				count++;
			}
		}
		ug->SetValues(fsize,indices,values,INS_VAL);

		xDelete<IssmDouble>(values);
		xDelete<int>(indices);
	}
	if(ssize){
		int*        indices = xNew<int>(ssize);
		IssmDouble* values  = xNew<IssmDouble>(ssize);

		int count = 0;
		for(int i=0;i<this->gsize;i++){
			if(this->s_set[i]){
				_assert_(local_ys);
				_assert_(this->sdoflist[i]==indices_ys[this->sdoflist_local[i]]);

				values[count] =local_ys[this->sdoflist_local[i]];
				indices[count]=this->gdoflist[i];
				count++;
			}
		}
		ug->SetValues(ssize,indices,values,INS_VAL);

		xDelete<IssmDouble>(values);
		xDelete<int>(indices);
	}
}
/*}}}*/
void Node::VecReduce(Vector<IssmDouble>* uf, IssmDouble* local_ug,int* indices_ug){/*{{{*/

	/*Only perform operation if not clone*/
	if(this->IsClone()) return;

	/*Get F size*/
	int fsize = this->FSize();

	if(fsize){
		int*        indices = xNew<int>(fsize);
		IssmDouble* values  = xNew<IssmDouble>(fsize);

		int count = 0;
		for(int i=0;i<this->gsize;i++){
			if(this->f_set[i]){
				_assert_(local_ug);
				_assert_(this->gdoflist[i]==indices_ug[this->gdoflist_local[i]]);
				_assert_(this->fdoflist[i]>=0);

				values[count] =local_ug[this->gdoflist_local[i]];
				indices[count]=this->fdoflist[i];
				count++;
			}
		}
		_assert_(count==fsize);
		uf->SetValues(fsize,indices,values,INS_VAL);

		xDelete<IssmDouble>(values);
		xDelete<int>(indices);
	}
}
/*}}}*/

/* indexing routines:*/
void Node::DistributeLocalDofs(int* pdofcount,int setenum){/*{{{*/

	/*Get current count*/
	int dofcount=*pdofcount;

	/*This node should distribute dofs for setenum set (eg, f_set or s_set), go ahead: */
	if(setenum==GsetEnum){
		_assert_(this->gsize==0 || this->gdoflist_local);
		for(int i=0;i<this->gsize;i++) gdoflist_local[i]=dofcount++;
	}
	else if(setenum==FsetEnum){
		for(int i=0;i<this->gsize;i++){
			if(this->f_set[i]) fdoflist_local[i]=dofcount++;
			else               fdoflist_local[i]=-1;
		}
	}
	else if(setenum==SsetEnum){
		for(int i=0;i<this->gsize;i++){
			if(this->s_set[i]) sdoflist_local[i]=dofcount++;
			else               sdoflist_local[i]=-1;
		}
	}
	else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");

	/*Assign output pointers: */
	*pdofcount=dofcount;
}/*}}}*/
void Node::DistributeGlobalDofsMasters(int dofcount,int setenum){/*{{{*/

	/*This node is a clone, don't offset the dofs!: */
	if(clone) return;

	/*This node should off_set the dofs, go ahead: */
	if(setenum==GsetEnum){
		_assert_(this->gsize==0 || this->gdoflist);
		for(int i=0;i<this->gsize;i++) this->gdoflist[i]=this->gdoflist_local[i]+dofcount;
	}
	else if(setenum==FsetEnum){
		for(int i=0;i<this->gsize;i++){
			if(this->f_set[i]) fdoflist[i]=this->fdoflist_local[i]+dofcount;
			else               fdoflist[i]=-1;
		}
	}
	else if(setenum==SsetEnum){
		for(int i=0;i<this->gsize;i++){
			if(this->s_set[i]) sdoflist[i]=this->sdoflist_local[i]+dofcount;
			else               sdoflist[i]=-1;
		}
	}
	else _error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
}
/*}}}*/
void Node::ShowMasterDofs(int* truedofs,int setenum){/*{{{*/

	_assert_(!this->clone);

	/*Ok, we are not a clone, just plug our dofs into truedofs: */
	switch(setenum){
		case GsetEnum:
			for(int j=0;j<this->gsize;j++) truedofs[j]=gdoflist[j];
			break;
		case FsetEnum:
			for(int j=0;j<this->gsize;j++) truedofs[j]=fdoflist[j];
			break;
		case SsetEnum:
			for(int j=0;j<this->gsize;j++) truedofs[j]=sdoflist[j];
			break;
		default:
			_error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
	}

}
/*}}}*/
void Node::UpdateCloneDofs(int* alltruedofs,int setenum){/*{{{*/

	_assert_(this->clone);

	/*Ok, we are a clone node, but we did not create the dofs for this node.
	 *Therefore, our doflist is garbage right now. Go pick it up in the alltruedofs: */
	switch(setenum){
		case GsetEnum:
			for(int j=0;j<this->gsize;j++) gdoflist[j]=alltruedofs[j];
			break;
		case FsetEnum:
			for(int j=0;j<this->gsize;j++) fdoflist[j]=alltruedofs[j];
			break;
		case SsetEnum:
			for(int j=0;j<this->gsize;j++) sdoflist[j]=alltruedofs[j];
			break;
		default:
			_error_("set of enum type " << EnumToStringx(setenum) << " not supported yet!");
	}
}
/*}}}*/
int  Node::FSize(void){/*{{{*/

	_assert_(this!=NULL && this->gdoflist);

	int fsize = 0;
	for(int i=0;i<this->gsize;i++) if(this->f_set[i]) fsize++;
	return fsize;
}
/*}}}*/
int  Node::SSize(void){/*{{{*/

	_assert_(this!=NULL && this->s_set);

	int ssize = 0;
	for(int i=0;i<this->gsize;i++) if(this->s_set[i]) ssize++;
	return ssize;
}
/*}}}*/

/*Methods inherent to Node: */
int* GetGlobalDofList(Node** nodes,int numnodes,int setenum,int approximation){/*{{{*/

	/*output*/
	int *doflist = NULL;

	if(numnodes){

		/*Allocate:*/
		int* ndof_list=xNew<int>(numnodes);

		/*First, figure out size of doflist: */
		int numdof=0;
		for(int i=0;i<numnodes;i++){
			ndof_list[i]=nodes[i]->GetNumberOfDofs(approximation,GsetEnum);
			numdof+=ndof_list[i];
		}

		if(numdof){
			/*Allocate: */
			doflist=xNew<int>(numdof);

			/*Populate: */
			int count=0;
			for(int i=0;i<numnodes;i++){
				nodes[i]->GetDofList(&doflist[count],approximation,setenum);
				count+=ndof_list[i];
			}
		}
		else doflist=NULL;

		/*Free resources:*/
		xDelete<int>(ndof_list);
	}

	return doflist;
}
/*}}}*/
int GetNumberOfDofs(Node** nodes,int numnodes,int setenum,int approximation){/*{{{*/

	/*output: */
	int numberofdofs=0;

	for(int i=0;i<numnodes;i++){
		numberofdofs+=nodes[i]->GetNumberOfDofs(approximation,setenum);
	}

	return numberofdofs;
}
/*}}}*/
