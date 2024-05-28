/*!\file IoModelToConstraintsx
 * \brief: used in ModelProcessor, to retrieve a nodal vector  and to create corresponding constraints using Constraints objects
 */

#include "./IoModelToConstraintsx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../ModelProcessorx/ModelProcessorx.h"

void IoModelToConstraintsx(Constraints* constraints,IoModel* iomodel,const char* spc_name,int analysis_type,int finite_element,int dof){

	/*intermediary: */
	int         code,vector_layout;
	IssmDouble *spcdata = NULL;
	int         M,N;

	/*First of, find the record for the enum, and get code  of data type: */
	iomodel->SetFilePointerToData(&code, &vector_layout,spc_name);
	if(code!=7)_error_("expecting a IssmDouble vector for constraints " << spc_name);
	if(vector_layout!=1)_error_("expecting a nodal vector for constraints " << spc_name);

	/*Fetch vector:*/
	iomodel->FetchData(&spcdata,&M,&N,spc_name);

	/*Call IoModelToConstraintsx*/
	IoModelToConstraintsx(constraints,iomodel,spcdata,M,N,analysis_type,finite_element,dof);

	/*Clean up*/
	xDelete<IssmDouble>(spcdata);
}
void IoModelToDynamicConstraintsx(Constraints* constraints,IoModel* iomodel,const char* spc_name,int analysis_type,int finite_element,int dof){

	/*intermediary: */
	int         code,vector_layout;
	IssmDouble *spcdata = NULL;
	int         M,N;

	/*First of, find the record for the enum, and get code  of data type: */
	iomodel->SetFilePointerToData(&code, &vector_layout,spc_name);
	if(code!=7)_error_("expecting a IssmDouble vector for constraints " << spc_name);
	if(vector_layout!=1)_error_("expecting a nodal vector for constraints " << spc_name);

	/*Fetch vector:*/
	iomodel->FetchData(&spcdata,&M,&N,spc_name);

	/*Call IoModelToConstraintsx*/
	IoModelToDynamicConstraintsx(constraints,iomodel,spcdata,M,N,analysis_type,finite_element,dof);

	/*Clean up*/
	xDelete<IssmDouble>(spcdata);
}

void IoModelToConstraintsx(Constraints* constraints,IoModel* iomodel,IssmDouble* spcdata,int M,int N,int analysis_type,int finite_element,int dof){/*{{{*/

	/*intermediary: */
	int         i,j,elementnbv,numfacevertices;
	IssmDouble  value;
	IssmDouble *times      = NULL;
	IssmDouble *values     = NULL;
	bool        spcpresent = false;

	/*Higher-order finite elements*/
	int   v1,v2;
	bool *boundaryedge = NULL;

	switch(finite_element){
		case P1Enum: 
			/*Nothing else to do*/
			break;
		case P1bubbleEnum:
			switch(iomodel->meshelementtype){
				case TriaEnum:  elementnbv = 3; break;
				case TetraEnum: elementnbv = 4; break;
				case PentaEnum: elementnbv = 6; break;
				default: _error_("mesh type not supported yet");
			}
			break;
		case P1bubblecondensedEnum: 
			/*Nothing else to do*/
			break;
		case P1xP2Enum:
		case P1xP3Enum:
		case P1xP4Enum:
			EdgesPartitioning(iomodel);
			break;
		case P2xP1Enum:
			EdgesPartitioning(iomodel);
			break;
		case P2Enum:
			EdgesPartitioning(iomodel);
	      if(iomodel->meshelementtype==PentaEnum) FacesPartitioning(iomodel);
			EdgeOnBoundaryFlags(&boundaryedge,iomodel);
			break;
		case P2bubbleEnum:
			EdgesPartitioning(iomodel);
			if(iomodel->meshelementtype==PentaEnum){
				FacesPartitioning(iomodel);
			}
			EdgeOnBoundaryFlags(&boundaryedge,iomodel);
			switch(iomodel->meshelementtype){
				case TriaEnum:  elementnbv = 3; break;
				case TetraEnum: elementnbv = 4; break;
				case PentaEnum: elementnbv = 6; break;
				default: _error_("mesh type not supported yet");
			}
			break;
		case P2xP4Enum:
			EdgesPartitioning(iomodel);
			FacesPartitioning(iomodel);
			break;
		default:
			_error_("Finite element "<<EnumToStringx(finite_element)<<" not supported yet");
	}

	int count =0;
	if(M==iomodel->numberofvertices){
		switch(finite_element){
			case P1Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){
						if (!xIsNan<IssmDouble>(spcdata[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,dof,spcdata[i],analysis_type));
							count++;
						}
					}
				}
				break;
			case P2Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){
						if (!xIsNan<IssmDouble>(spcdata[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,dof,spcdata[i],analysis_type));
							count++;
						}
					}
				}
				for(i=0;i<iomodel->numberofedges;i++){
					if(iomodel->my_edges[i] && boundaryedge[i]){
						v1 = iomodel->edges[3*i+0]-1;
						v2 = iomodel->edges[3*i+1]-1;
						if(!xIsNan<IssmDouble>(spcdata[v1]) && !xIsNan<IssmDouble>(spcdata[v2])){
							constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+i+1,
											dof,(spcdata[v1]+spcdata[v2])/2.,analysis_type));
							count++;
						}
					}
				}
				if(iomodel->meshelementtype==PentaEnum){
					for(i=0;i<iomodel->numberofverticalfaces;i++){
						if(iomodel->my_vfaces[i]){
							value=0.;
							for(j=0;j<4;j++) value += spcdata[iomodel->verticalfaces[i*4+j] -1]/4.;
							if(!xIsNan<IssmDouble>(value)){
								constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+iomodel->numberofedges+i+1,dof,value,analysis_type));
								count++;
							}
						}
					}
				}
				break;
			case P2bubbleEnum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){
						if (!xIsNan<IssmDouble>(spcdata[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,dof,spcdata[i],analysis_type));
							count++;
						}
					}
				}
				for(i=0;i<iomodel->numberofedges;i++){
					if(iomodel->my_edges[i] && boundaryedge[i]){
						v1 = iomodel->edges[3*i+0]-1;
						v2 = iomodel->edges[3*i+1]-1;
						if(!xIsNan<IssmDouble>(spcdata[v1]) && !xIsNan<IssmDouble>(spcdata[v2])){
							constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+i+1,
											dof,(spcdata[v1]+spcdata[v2])/2.,analysis_type));
							count++;
						}
					}
				}
				if(iomodel->meshelementtype==PentaEnum){
					for(i=0;i<iomodel->numberoffaces;i++){
						if(iomodel->faces[i*iomodel->facescols+2]==2){/*Vertical quads*/
							if(iomodel->my_faces[i]){
								numfacevertices = iomodel->faces[i*iomodel->facescols+3];
								value=0.;
								for(j=0;j<numfacevertices;j++){
									value += spcdata[iomodel->faces[i*iomodel->facescols+4+j] -1];
								}
								value = value/reCast<IssmDouble>(numfacevertices);
								if(!xIsNan<IssmDouble>(value)){
									constraints->AddObject(new SpcStatic(count+1,
													iomodel->numberofvertices+iomodel->numberofedges+i+1,
													dof,value,analysis_type));
									count++;
								}
							}
						}
					}
				}
				for(i=0;i<iomodel->numberofelements;i++){
					if(iomodel->my_elements[i]){
						value = spcdata[iomodel->elements[i*elementnbv+0]-1];
						for(j=1;j<elementnbv;j++) value += spcdata[iomodel->elements[i*elementnbv+j]-1];
						value = value/reCast<IssmDouble,int>(elementnbv+0);
						if(!xIsNan<IssmDouble>(value)){
							int nodeid = iomodel->numberofvertices+iomodel->numberofedges+i+1;
							if(iomodel->meshelementtype==PentaEnum){
								nodeid += iomodel->numberoffaces;
							}
							constraints->AddObject(new SpcStatic(count+1,nodeid,dof,value,analysis_type));
							count++;
						}
					}
				}
				break;
			case P2xP4Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){
						if (!xIsNan<IssmDouble>(spcdata[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,dof,spcdata[i],analysis_type));
							count++;
						}
					}
				}
				for(i=0;i<iomodel->numberofedges;i++){
					if(iomodel->my_edges[i]){
						v1 = iomodel->edges[3*i+0]-1;
						v2 = iomodel->edges[3*i+1]-1;
						if(!xIsNan<IssmDouble>(spcdata[v1]) && !xIsNan<IssmDouble>(spcdata[v2])){
							constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+i+1,
											dof,1./2.*spcdata[v1]+1./2.*spcdata[v2],analysis_type));
							count++;
						}
					}
				}
				for(i=0;i<iomodel->numberofverticaledges;i++){
					if(iomodel->my_vedges[i]){
						v1 = iomodel->verticaledges[2*i+0]-1;
						v2 = iomodel->verticaledges[2*i+1]-1;
						if(!xIsNan<IssmDouble>(spcdata[v1]) && !xIsNan<IssmDouble>(spcdata[v2])){
							/*FIXME: should not be 1/2 but 1/4 and 3/4!*/
							constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+iomodel->numberofedges+2*i+1,
											dof,1./2.*spcdata[v1]+1./2.*spcdata[v2],analysis_type));
							constraints->AddObject(new SpcStatic(count+2,iomodel->numberofvertices+iomodel->numberofedges+2*i+2,
											dof,1./2.*spcdata[v1]+1./2.*spcdata[v2],analysis_type));
							count=count+2;
						}
					}
				}
				for(i=0;i<iomodel->numberofverticalfaces;i++){
					if(iomodel->my_vfaces[i]){
						value=0.;
						for(j=0;j<4;j++) value += spcdata[iomodel->verticalfaces[i*4+j]-1]/4.;
						if(!xIsNan<IssmDouble>(value)){
							constraints->AddObject(new SpcStatic(count+1,
											iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*i+1,
											dof,value,analysis_type));
							constraints->AddObject(new SpcStatic(count+2,
											iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*i+2,
											dof,value,analysis_type));
							constraints->AddObject(new SpcStatic(count+3,
											iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*i+3,
											dof,value,analysis_type));
							count=count+3;
						}
					}
				}
				break;
			case P1bubbleEnum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){
						if (!xIsNan<IssmDouble>(spcdata[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,dof,spcdata[i],analysis_type));
							count++;
						}
					}
				}
				for(i=0;i<iomodel->numberofelements;i++){
					if(iomodel->my_elements[i]){
						value = spcdata[iomodel->elements[i*elementnbv+0]-1];
						for(j=1;j<elementnbv;j++) value += spcdata[iomodel->elements[i*elementnbv+j]-1];
						value = value/reCast<IssmDouble,int>(elementnbv+0);
						if(!xIsNan<IssmDouble>(value)){
							constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+i+1,
											dof,value,analysis_type));
							count++;
						}
					}
				}
				break;
			case P1bubblecondensedEnum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){
						if (!xIsNan<IssmDouble>(spcdata[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,dof,spcdata[i],analysis_type));
							count++;
						}
					}
				}
				break;
			case P1xP2Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){
						if (!xIsNan<IssmDouble>(spcdata[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,dof,spcdata[i],analysis_type));
							count++;
						}
					}
				}
				for(i=0;i<iomodel->numberofverticaledges;i++){
					if(iomodel->my_vedges[i]){
						v1 = iomodel->verticaledges[2*i+0]-1;
						v2 = iomodel->verticaledges[2*i+1]-1;
						if(!xIsNan<IssmDouble>(spcdata[v1]) && !xIsNan<IssmDouble>(spcdata[v2])){
							constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+i+1,dof,(spcdata[v1]+spcdata[v2])/2.,analysis_type));
							count++;
						}
					}
				}
				break;
			case P1xP3Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){
						if (!xIsNan<IssmDouble>(spcdata[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,dof,spcdata[i],analysis_type));
							count++;
						}
					}
				}
				for(i=0;i<iomodel->numberofverticaledges;i++){
					if(iomodel->my_vedges[i]){
						v1 = iomodel->verticaledges[2*i+0]-1;
						v2 = iomodel->verticaledges[2*i+1]-1;
						if(!xIsNan<IssmDouble>(spcdata[v1]) && !xIsNan<IssmDouble>(spcdata[v2])){
							constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+2*i+1,
											dof,2./3.*spcdata[v1]+1./3.*spcdata[v2],analysis_type));
							constraints->AddObject(new SpcStatic(count+2,iomodel->numberofvertices+2*i+2,
											dof,1./3.*spcdata[v1]+2./3.*spcdata[v2],analysis_type));
							count=count+2;
						}
					}
				}
				break;
			case P1xP4Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){
						if (!xIsNan<IssmDouble>(spcdata[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,dof,spcdata[i],analysis_type));
							count++;
						}
					}
				}
				for(i=0;i<iomodel->numberofverticaledges;i++){
					if(iomodel->my_vedges[i]){
						v1 = iomodel->verticaledges[2*i+0]-1;
						v2 = iomodel->verticaledges[2*i+1]-1;
						if(!xIsNan<IssmDouble>(spcdata[v1]) && !xIsNan<IssmDouble>(spcdata[v2])){
							constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+3*i+1,
											dof,1./2.*spcdata[v1]+1./2.*spcdata[v2],analysis_type));
							constraints->AddObject(new SpcStatic(count+2,iomodel->numberofvertices+3*i+2,
											dof,1./2.*spcdata[v1]+1./2.*spcdata[v2],analysis_type));
							constraints->AddObject(new SpcStatic(count+3,iomodel->numberofvertices+3*i+3,
											dof,1./2.*spcdata[v1]+1./2.*spcdata[v2],analysis_type));
							count=count+3;
						}
					}
				}
				break;
			case P2xP1Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){
						if (!xIsNan<IssmDouble>(spcdata[i])){
							constraints->AddObject(new SpcStatic(count+1,i+1,dof,spcdata[i],analysis_type));
							count++;
						}
					}
				}
				for(i=0;i<iomodel->numberofhorizontaledges;i++){
					if(iomodel->my_hedges[i]){
						v1 = iomodel->horizontaledges[2*i+0]-1;
						v2 = iomodel->horizontaledges[2*i+1]-1;
						if(!xIsNan<IssmDouble>(spcdata[v1]) && !xIsNan<IssmDouble>(spcdata[v2])){
							constraints->AddObject(new SpcStatic(count+1,iomodel->numberofvertices+i+1,dof,(spcdata[v1]+spcdata[v2])/2.,analysis_type));
							count++;
						}
					}
				}
				break;
			default:
				_error_("Finite element "<<EnumToStringx(finite_element)<<" not supported yet");
		}
	}
	else if (M==(iomodel->numberofvertices+1)){
		/*transient: create transient SpcTransient objects. Same logic, except we need to retrieve 
		 * various times and values to initialize an SpcTransient object: */

		/*figure out times: */
		times=xNew<IssmDouble>(N);
		for(j=0;j<N;j++) times[j]=spcdata[(M-1)*N+j];

		switch(finite_element){
			case P1Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){

						/*figure out times and values: */
						values=xNew<IssmDouble>(N);
						spcpresent=false;
						for(j=0;j<N;j++){
							values[j]=spcdata[i*N+j];
							if(!xIsNan<IssmDouble>(values[j]))spcpresent=true; //NaN means no spc by default
						}

						if(spcpresent){
							constraints->AddObject(new SpcTransient(count+1,i+1,dof,N,times,values,analysis_type));
							count++;
						}
						xDelete<IssmDouble>(values);
					}
				}
				break;
			case P2Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){

						/*figure out times and values: */
						values=xNew<IssmDouble>(N);
						spcpresent=false;
						for(j=0;j<N;j++){
							values[j]=spcdata[i*N+j];
							if(!xIsNan<IssmDouble>(values[j]))spcpresent=true; //NaN means no spc by default
						}
						if(iomodel->meshelementtype==PentaEnum){ _error_("need to add nodes on faces! (See Spcstatic above)");}

						if(spcpresent){
							constraints->AddObject(new SpcTransient(count+1,i+1,dof,N,times,values,analysis_type));
							count++;
						}
						xDelete<IssmDouble>(values);
					}
				}
				for(i=0;i<iomodel->numberofedges;i++){
					if(iomodel->my_edges[i]){
						v1 = iomodel->edges[3*i+0]-1;
						v2 = iomodel->edges[3*i+1]-1;
						values=xNew<IssmDouble>(N);
						spcpresent=false;
						for(j=0;j<N;j++){
							values[j]=(spcdata[v1*N+j]+spcdata[v2*N+j])/2.;
							if(!xIsNan<IssmDouble>(values[j])) spcpresent=true; //NaN means no spc by default
						}
						if(spcpresent){
							constraints->AddObject(new SpcTransient(count+1,iomodel->numberofvertices+i+1,dof,N,times,values,analysis_type));
							count++;
						}
						xDelete<IssmDouble>(values);
					}
				}
				break;
			case P1xP2Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){

						/*figure out times and values: */
						values=xNew<IssmDouble>(N);
						spcpresent=false;
						for(j=0;j<N;j++){
							values[j]=spcdata[i*N+j];
							if(!xIsNan<IssmDouble>(values[j]))spcpresent=true; //NaN means no spc by default
						}

						if(spcpresent){
							constraints->AddObject(new SpcTransient(count+1,i+1,dof,N,times,values,analysis_type));
							count++;
						}
						xDelete<IssmDouble>(values);
					}
				}
				for(i=0;i<iomodel->numberofverticaledges;i++){
					if(iomodel->my_vedges[i]){
						v1 = iomodel->verticaledges[2*i+0]-1;
						v2 = iomodel->verticaledges[2*i+1]-1;
						values=xNew<IssmDouble>(N);
						spcpresent=false;
						for(j=0;j<N;j++){
							values[j]=(spcdata[v1*N+j]+spcdata[v2*N+j])/2.;
							if(!xIsNan<IssmDouble>(values[j])) spcpresent=true; //NaN means no spc by default
						}
						if(spcpresent){
							constraints->AddObject(new SpcTransient(count+1,iomodel->numberofvertices+i+1,dof,N,times,values,analysis_type));
							count++;
						}
						xDelete<IssmDouble>(values);
					}
				}
				break;
			case P1xP3Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){

						/*figure out times and values: */
						values=xNew<IssmDouble>(N);
						spcpresent=false;
						for(j=0;j<N;j++){
							values[j]=spcdata[i*N+j];
							if(!xIsNan<IssmDouble>(values[j]))spcpresent=true; //NaN means no spc by default
						}

						if(spcpresent){
							constraints->AddObject(new SpcTransient(count+1,i+1,dof,N,times,values,analysis_type));
							count++;
						}
						xDelete<IssmDouble>(values);
					}
				}
				for(i=0;i<iomodel->numberofverticaledges;i++){
					if(iomodel->my_vedges[i]){
						v1 = iomodel->verticaledges[2*i+0]-1;
						v2 = iomodel->verticaledges[2*i+1]-1;
						values=xNew<IssmDouble>(N);
						spcpresent=false;
						for(j=0;j<N;j++){
							values[j]=2./3.*spcdata[v1*N+j]+1./3.*spcdata[v2*N+j];
							if(!xIsNan<IssmDouble>(values[j])) spcpresent=true; //NaN means no spc by default
						}
						if(spcpresent){
							constraints->AddObject(new SpcTransient(count+1,iomodel->numberofvertices+2*i+1,dof,N,times,values,analysis_type));
							count++;
						}
						spcpresent=false;
						for(j=0;j<N;j++){
							values[j]=1./3.*spcdata[v1*N+j]+2./3.*spcdata[v2*N+j];
							if(!xIsNan<IssmDouble>(values[j])) spcpresent=true; //NaN means no spc by default
						}
						if(spcpresent){
							constraints->AddObject(new SpcTransient(count+1,iomodel->numberofvertices+2*i+2,dof,N,times,values,analysis_type));
							count++;
						}
						xDelete<IssmDouble>(values);
					}
				}
				break;
			case P2xP1Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){

						/*figure out times and values: */
						values=xNew<IssmDouble>(N);
						spcpresent=false;
						for(j=0;j<N;j++){
							values[j]=spcdata[i*N+j];
							if(!xIsNan<IssmDouble>(values[j]))spcpresent=true; //NaN means no spc by default
						}

						if(spcpresent){
							constraints->AddObject(new SpcTransient(count+1,i+1,dof,N,times,values,analysis_type));
							count++;
						}
						xDelete<IssmDouble>(values);
					}
				}
				for(i=0;i<iomodel->numberofedges;i++){
					if(iomodel->edges[i*3+2]!=2){
						if(iomodel->my_edges[i]){
							v1 = iomodel->edges[3*i+0]-1;
							v2 = iomodel->edges[3*i+1]-1;
							values=xNew<IssmDouble>(N);
							spcpresent=false;
							for(j=0;j<N;j++){
								values[j]=(spcdata[v1*N+j]+spcdata[v2*N+j])/2.;
								if(!xIsNan<IssmDouble>(values[j])) spcpresent=true; //NaN means no spc by default
							}
							if(spcpresent){
								constraints->AddObject(new SpcTransient(count+1,iomodel->numberofvertices+i+1,dof,N,times,values,analysis_type));
								count++;
							}
							xDelete<IssmDouble>(values);
						}
					}
				}
				break;
			default:
				_error_("Finite element "<<EnumToStringx(finite_element)<<" not supported yet");
		}
	}
	else{
		_error_("Size of spc field not supported");
	}

	/*Free resources:*/
	xDelete<IssmDouble>(times);
	xDelete<IssmDouble>(values);
	xDelete<bool>(boundaryedge);
}/*}}}*/
void IoModelToDynamicConstraintsx(Constraints* constraints,IoModel* iomodel,IssmDouble* spcdata,int M,int N,int analysis_type,int finite_element,int dof){/*{{{*/

	/*intermediary: */
	int         i,j,elementnbv,numfacevertices;
	IssmDouble  value;
	IssmDouble *times            = NULL;
	IssmDouble *values           = NULL;
	bool        spcpresent       = false;

	/*Higher-order finite elements*/
	int   v1,v2;
	bool *boundaryedge = NULL;

	switch(finite_element){
		case P1Enum: 
			/*Nothing else to do*/
			break;
		case P1bubbleEnum:
			switch(iomodel->meshelementtype){
				case TriaEnum:  elementnbv = 3; break;
				case TetraEnum: elementnbv = 4; break;
				case PentaEnum: elementnbv = 6; break;
				default: _error_("mesh type not supported yet");
			}
			break;
		case P1bubblecondensedEnum: 
			/*Nothing else to do*/
			break;
		case P1xP2Enum:
			EdgesPartitioning(iomodel);
			break;
		case P1xP3Enum:
			EdgesPartitioning(iomodel);
			break;
		case P2xP1Enum:
			EdgesPartitioning(iomodel);
			break;
		case P2Enum:
			EdgesPartitioning(iomodel);
			if(iomodel->meshelementtype==PentaEnum){
				FacesPartitioning(iomodel);
			}
			EdgeOnBoundaryFlags(&boundaryedge,iomodel);
			break;
		case P2bubbleEnum:
			EdgesPartitioning(iomodel);
			if(iomodel->meshelementtype==PentaEnum){
				FacesPartitioning(iomodel);
			}
			EdgeOnBoundaryFlags(&boundaryedge,iomodel);
			switch(iomodel->meshelementtype){
				case TriaEnum:  elementnbv = 3; break;
				case TetraEnum: elementnbv = 4; break;
				case PentaEnum: elementnbv = 6; break;
				default: _error_("mesh type not supported yet");
			}
			break;
		case P2xP4Enum:
			EdgesPartitioning(iomodel);
			FacesPartitioning(iomodel);
			break;
		default:
			_error_("Finite element "<<EnumToStringx(finite_element)<<" not supported yet");
	}

	int count=0;
	if(M==iomodel->numberofvertices){
		switch(finite_element){
			case P1Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){
						if (!xIsNan<IssmDouble>(spcdata[i])){
							constraints->AddObject(new SpcDynamic(count+1,i+1,dof,spcdata[i],analysis_type));
							count++;
						}
					}
				}
				break;
			case P1xP2Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){
						if (!xIsNan<IssmDouble>(spcdata[i])){
							constraints->AddObject(new SpcDynamic(count+1,i+1,dof,spcdata[i],analysis_type));
							count++;
						}
					}
				}
				for(i=0;i<iomodel->numberofverticaledges;i++){
					if(iomodel->my_vedges[i]){
						v1 = iomodel->verticaledges[2*i+0]-1;
						v2 = iomodel->verticaledges[2*i+1]-1;
						if(!xIsNan<IssmDouble>(spcdata[v1]) && !xIsNan<IssmDouble>(spcdata[v2])){
							constraints->AddObject(new SpcDynamic(count+1,iomodel->numberofvertices+i+1,dof,(spcdata[v1]+spcdata[v2])/2.,analysis_type));
							count++;
						}
					}
				}
				break;
			case P1xP3Enum:
				for(i=0;i<iomodel->numberofvertices;i++){
					if((iomodel->my_vertices[i])){
						if (!xIsNan<IssmDouble>(spcdata[i])){
							constraints->AddObject(new SpcDynamic(count+1,i+1,dof,spcdata[i],analysis_type));
							count++;
						}
					}
				}
				for(i=0;i<iomodel->numberofverticaledges;i++){
					if(iomodel->my_vedges[i]){
						v1 = iomodel->verticaledges[2*i+0]-1;
						v2 = iomodel->verticaledges[2*i+1]-1;
						if(!xIsNan<IssmDouble>(spcdata[v1]) && !xIsNan<IssmDouble>(spcdata[v2])){
							constraints->AddObject(new SpcDynamic(count+1,iomodel->numberofvertices+2*i+1,
											dof,2./3.*spcdata[v1]+1./3.*spcdata[v2],analysis_type));
							constraints->AddObject(new SpcDynamic(count+2,iomodel->numberofvertices+2*i+2,
											dof,1./3.*spcdata[v1]+2./3.*spcdata[v2],analysis_type));
							count=count+2;
						}
					}
				}
				break;
			default:
				_error_("Finite element "<<EnumToStringx(finite_element)<<" not supported yet");
		}
	}
	else{
		_error_("Size of spc field not supported");
	}

	/*Free resources:*/
	xDelete<IssmDouble>(times);
	xDelete<IssmDouble>(values);
	xDelete<bool>(boundaryedge);
}/*}}}*/
