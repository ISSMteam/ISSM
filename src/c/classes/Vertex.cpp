/*!\file Vertex.c
 * \brief: implementation of the Vertex object
 */

/*Include files: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <string.h>
#include "classes.h"
#include "shared/shared.h"
/*}}}*/

/*Vertex constructors and destructor:*/
Vertex::Vertex(){/*{{{*/
	return;
}
/*}}}*/
Vertex::Vertex(int vertex_id, int vertex_sid,bool vertex_clone, IoModel* iomodel,bool isamr){/*{{{*/

	/*Checks in debugging mode*/
	_assert_(vertex_sid>=0 && vertex_sid<iomodel->numberofvertices);

	/*IDs*/
	this->id    = vertex_id;
	this->sid   = vertex_sid;
	this->pid   = -1; /*Assigned later*/
	this->lid   = -1; /*Assigned later*/
	this->clone = vertex_clone;

	/*Properties from iomodel*/
	_assert_(iomodel->numbernodetoelementconnectivity);
	this->connectivity = iomodel->numbernodetoelementconnectivity[vertex_sid];
	this->domaintype   = iomodel->domaintype;

	/*Coordinates, only if not AMR*/
	if(!isamr){
		_assert_(iomodel->Data("md.mesh.x") && iomodel->Data("md.mesh.y") && iomodel->Data("md.mesh.z"));
		this->x            = iomodel->Data("md.mesh.x")[vertex_sid];
		this->y            = iomodel->Data("md.mesh.y")[vertex_sid];
		this->z            = iomodel->Data("md.mesh.z")[vertex_sid];
		if(iomodel->Data("md.mesh.lat") && iomodel->Data("md.mesh.long")){
			this->latitute     = iomodel->Data("md.mesh.lat")[vertex_sid];
			this->longitude    = iomodel->Data("md.mesh.long")[vertex_sid];
		}

		switch(iomodel->domaintype){
			case Domain3DEnum:
				_assert_(iomodel->Data("md.geometry.base") && iomodel->Data("md.geometry.thickness"));
				this->sigma = (iomodel->Data("md.mesh.z")[vertex_sid]-iomodel->Data("md.geometry.base")[vertex_sid])/(iomodel->Data("md.geometry.thickness")[vertex_sid]);
				break;
			case Domain3DsurfaceEnum:
				_assert_(iomodel->Data("md.mesh.lat") && iomodel->Data("md.mesh.long") && iomodel->Data("md.mesh.r"));
				this->latitute     = iomodel->Data("md.mesh.lat")[vertex_sid];
				this->longitude    = iomodel->Data("md.mesh.long")[vertex_sid];
				this->R            = iomodel->Data("md.mesh.r")[vertex_sid];
				break;
			case Domain2DhorizontalEnum:
				this->sigma = 0.;
				break;
			case Domain2DverticalEnum:
				_assert_(iomodel->Data("md.geometry.base") && iomodel->Data("md.geometry.thickness"));
				this->sigma = (iomodel->Data("md.mesh.y")[vertex_sid]-iomodel->Data("md.geometry.base")[vertex_sid])/(iomodel->Data("md.geometry.thickness")[vertex_sid]);
				break;
		}
	}
	else{
		this->x         = 0.;
		this->y         = 0.;
		this->z         = 0.;
		this->latitute  = 0.;
		this->longitude = 0.;
		this->R         = 0.;
		this->sigma     = 0.;
	}

}/*}}}*/
Vertex::~Vertex(){/*{{{*/
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* Vertex::copy() {/*{{{*/

	return new Vertex(*this); 

}
/*}}}*/
void Vertex::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void Vertex::Echo(void){/*{{{*/

	_printf_("Vertex:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   sid: " << sid << "\n");
	_printf_("   pid: " << pid << "\n");
	_printf_("   lid: " << lid << "\n");
	_printf_("   x: " << x << "\n");
	_printf_("   y: " << y << "\n");
	_printf_("   z: " << z << "\n");
	_printf_("   sigma: " << sigma << "\n");
	_printf_("   connectivity: " << connectivity << "\n");
	_printf_("   clone: " << clone << "\n");

	return;
}
/*}}}*/
int Vertex::Id(void){ return id; }/*{{{*/
/*}}}*/
void Vertex::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = VertexEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->clone);
	marshallhandle->call(this->domaintype);
	marshallhandle->call(this->id);
	marshallhandle->call(this->sid);
	marshallhandle->call(this->pid);
	marshallhandle->call(this->lid);
	marshallhandle->call(this->x);
	marshallhandle->call(this->y);
	marshallhandle->call(this->z);
	marshallhandle->call(this->sigma);
	marshallhandle->call(this->connectivity);

}
/*}}}*/
int Vertex::ObjectEnum(void){/*{{{*/
	return VertexEnum;
}/*}}}*/

/*Vertex management: */
int        Vertex::Connectivity(void){return connectivity;}/*{{{*/
/*}}}*/
IssmDouble Vertex::GetLatitude(){/*{{{*/
	return this->latitute;
}
/*}}}*/
IssmDouble Vertex::GetLongitude(){/*{{{*/
	return this->longitude;
}
/*}}}*/
IssmDouble Vertex::GetRadius(){/*{{{*/
	return this->R;
}
/*}}}*/
IssmDouble Vertex::GetX(){/*{{{*/
	return this->x;
}
/*}}}*/
IssmDouble Vertex::GetY(){/*{{{*/
	return this->y;
}
/*}}}*/
IssmDouble Vertex::GetZ(){/*{{{*/
	return this->z;
}
/*}}}*/
int        Vertex::Pid(void){ return pid; }/*{{{*/
/*}}}*/
int        Vertex::Lid(void){ return lid; }/*{{{*/
/*}}}*/
int        Vertex::Sid(void){ return sid; }/*{{{*/
/*}}}*/
void       Vertex::UpdatePosition(Vector<IssmDouble>* vx,Vector<IssmDouble>* vy,Vector<IssmDouble>* vz,Parameters* parameters,IssmDouble* surface,IssmDouble* bed){/*{{{*/

	IssmDouble oldy,newy,vely;
	IssmDouble oldz,newz,velz;
	IssmDouble dt;

	/*Get time stepping*/
	parameters->FindParam(&dt,TimesteppingTimeStepEnum);

	/*sigma remains constant. z=bed+sigma*thickness*/
	switch(this->domaintype){
		case Domain2DhorizontalEnum:
		case Domain3DsurfaceEnum:
			/*Nothing*/
			return;
		case Domain2DverticalEnum:
			oldy = this->y;
			newy = bed[this->pid]+sigma*(surface[this->pid] - bed[this->pid]);
			vely = (newy-oldy)/dt;
			this->y = newy;
			vy->SetValue(this->pid,vely,INS_VAL);
			_assert_(!xIsNan<IssmDouble>(vely));
			return;
		case Domain3DEnum:
			oldz = this->z;
			newz = bed[this->pid]+sigma*(surface[this->pid] - bed[this->pid]);
			velz = (newz-oldz)/dt;
			this->z = newz;
			vz->SetValue(this->pid,velz,INS_VAL);
			_assert_(!xIsNan<IssmDouble>(velz));
			return;
		default:
			_error_("not implemented");
	}
}
/*}}}*/
void       Vertex::VertexCoordinates(Vector<IssmDouble>* vx,Vector<IssmDouble>* vy,Vector<IssmDouble>* vz,bool spherical){/*{{{*/

	if(this->clone==true) return;

	if(!spherical){
		vx->SetValue(this->sid,this->x,INS_VAL);
		vy->SetValue(this->sid,this->y,INS_VAL);
		vz->SetValue(this->sid,this->z,INS_VAL);
	}
	else{
		vx->SetValue(this->sid,this->latitute,INS_VAL);
		vy->SetValue(this->sid,this->longitude,INS_VAL);
		vz->SetValue(this->sid,this->R,INS_VAL);
	}
	return;
}
/*}}}*/

/*Methods relating to Vertex, but not internal methods: */
void GetVerticesCoordinates(IssmDouble* xyz,Vertex** vertices, int numvertices,bool spherical){ /*{{{*/

	_assert_(vertices);
	_assert_(xyz);

	if(!spherical){
		for(int i=0;i<numvertices;i++) {
			xyz[i*3+0]=vertices[i]->GetX();
			xyz[i*3+1]=vertices[i]->GetY();
			xyz[i*3+2]=vertices[i]->GetZ();
		}
	}
	else{
		for(int i=0;i<numvertices;i++) {
			xyz[i*3+0]=vertices[i]->GetLatitude();
			xyz[i*3+1]=vertices[i]->GetLongitude();
			xyz[i*3+2]=vertices[i]->GetRadius();
		}
	}
}/*}}}*/
