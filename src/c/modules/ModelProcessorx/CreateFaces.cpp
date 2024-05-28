/*!\file:  CreateFaces.cpp
 * \brief: create faces from 2d mesh
 */ 

#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "./ModelProcessorx.h"

void CreateFaces(IoModel* iomodel){/*{{{*/

	/*If faces are already present, exit*/
	if(iomodel->faces) return;

	/*Some checks*/
	if(iomodel->numberofvertices<3) _error_("not enough elements in mesh");
	_assert_(iomodel->elements);

	/*Check Iomodel properties*/
	if(iomodel->domaintype==Domain2DhorizontalEnum || iomodel->domaintype==Domain2DverticalEnum){
		/*Keep going*/
	}
	else if(iomodel->domaintype==Domain3DEnum){
		CreateFaces3d(iomodel);
		return;
	}
	else{
		_error_("mesh dimension not supported yet");
	}

	/*Intermediaries*/
	bool exist;
	int  i,j,v1,v2,v3;
	int  maxnbf,nbf;

	/*Maximum number of faces*/
	maxnbf = 3*iomodel->numberofelements;

	/*Initialize intermediaries*/
	int*  facestemp = xNew<int>(maxnbf*4);         /*format: [vertex1 vertex2 element1 element2]                */
	bool* exchange  = xNewZeroInit<bool>(maxnbf);  /*Faces are ordered, we need to keep track of vertex swapping*/
	for(i=0;i<maxnbf;i++) facestemp[i*4+3]=-1;     /*Initialize last column of faces as -1 (boundary edge)      */

	/*Initialize chain*/
	int* head_minv = xNew<int>(iomodel->numberofvertices);
	int* next_face = xNew<int>(maxnbf);
	for(i=0;i<iomodel->numberofvertices;i++) head_minv[i]=-1;

	/*Initialize number of faces*/
	nbf = 0;

	for(i=0;i<iomodel->numberofelements;i++){
		for(j=0;j<3;j++){

			/*Get the two indices of the edge number j of the ith triangle*/
			v1 = iomodel->elements[i*3+j];
			if(j==2)
			 v2 = iomodel->elements[i*3+0];
			else
			 v2 = iomodel->elements[i*3+j+1];

			/*v1 and v2 must be sorted*/
			if(v2<v1){
				v3=v2; v2=v1; v1=v3;
			}

			/*This edge a priori has not been processed yet*/
			exist = false;

			/*Go through all processed faces connected to v1 and check whether we have seen this edge yet*/
			_assert_(v1>=0 && v1<iomodel->numberofvertices);
			for(int e=head_minv[v1]; e!=-1; e=next_face[e]){
				if(facestemp[e*4+1]==v2){
					exist = true;
					facestemp[e*4+3]=i+1;
					break;
				}
			}

			/*If this edge is new, add it to the lists*/
			if(!exist){
				_assert_(nbf<maxnbf);

				/*Update faces*/
				facestemp[nbf*4+0] = v1;
				facestemp[nbf*4+1] = v2;
				facestemp[nbf*4+2] = i+1;
				if(v1!=iomodel->elements[i*3+j]) exchange[nbf]=true;

				/*Update chain*/
				next_face[nbf] = head_minv[v1];
				head_minv[v1]  = nbf;

				/*Increase number of faces*/
				nbf++;
			}
		}
	}

	/*Clean up*/
	xDelete<int>(head_minv);
	xDelete<int>(next_face);

	/*Create final faces*/
	int* faces = xNew<int>(nbf*4); /*vertex1 vertex2 element1 element2*/
	for(int i=0;i<nbf;i++){
		if(exchange[i]){
			faces[i*4+0]=facestemp[i*4+1];
			faces[i*4+1]=facestemp[i*4+0];
		}
		else{
			faces[i*4+0]=facestemp[i*4+0];
			faces[i*4+1]=facestemp[i*4+1];
		}
		faces[i*4+2]=facestemp[i*4+2];
		faces[i*4+3]=facestemp[i*4+3];
	}
	xDelete<int>(facestemp);
	xDelete<bool>(exchange);

	/*Assign output pointers*/
	iomodel->faces         = faces;
	iomodel->numberoffaces = nbf;
}/*}}}*/
void CreateFaces3d(IoModel* iomodel){/*{{{*/

	/*Intermediaries*/
	int  i,j,k,v0,cols,facemaxnbv;
	int  elementnbf,elementnbv,facenbv;
	int *elementfaces         = NULL;
	int *elementfaces_markers = NULL;

	/*Mesh specific face indexing per element*/
	switch(iomodel->meshelementtype){
		case PentaEnum:
			elementnbv = 6; /*Number of vertices per element*/
			elementnbf = 5; /*Number of faces per element*/
			facemaxnbv = 4; /*Maximum number of vertices per face*/
			cols       = facemaxnbv + 1;
			elementfaces         = xNew<int>(elementnbf*cols);
			elementfaces_markers = xNew<int>(elementnbf);
			/*2 triangles*/
			elementfaces_markers[0] = 1;
			elementfaces_markers[1] = 1;
			elementfaces[cols*0+0] = 3; elementfaces[cols*0+1] = 0;  elementfaces[cols*0+2] = 1; elementfaces[cols*0+3] = 2;
			elementfaces[cols*1+0] = 3; elementfaces[cols*1+1] = 3;  elementfaces[cols*1+2] = 4; elementfaces[cols*1+3] = 5;
			/*3 quads*/
			elementfaces_markers[2] = 2;
			elementfaces_markers[3] = 2;
			elementfaces_markers[4] = 2;
			elementfaces[cols*2+0] = 4; elementfaces[cols*2+1] = 1;  elementfaces[cols*2+2] = 2; elementfaces[cols*2+3] = 5;  elementfaces[cols*2+4] = 4;
			elementfaces[cols*3+0] = 4; elementfaces[cols*3+1] = 2;  elementfaces[cols*3+2] = 0; elementfaces[cols*3+3] = 3;  elementfaces[cols*3+4] = 5;
			elementfaces[cols*4+0] = 4; elementfaces[cols*4+1] = 0;  elementfaces[cols*4+2] = 1; elementfaces[cols*4+3] = 4;  elementfaces[cols*4+4] = 3;
			break;
		case TetraEnum:
			elementnbv = 4; /*Number of vertices per element*/
			elementnbf = 4; /*Number of faces per element*/
			facemaxnbv = 3; /*Maximum number of vertices per face*/
			cols       = facemaxnbv + 1;
			elementfaces         = xNew<int>(elementnbf*cols);
			elementfaces_markers = xNew<int>(elementnbf);
			/*4 triangles*/
			elementfaces_markers[0] = 1;
			elementfaces_markers[1] = 1;
			elementfaces_markers[2] = 1;
			elementfaces_markers[3] = 1;
			elementfaces[cols*0+0] = 3; elementfaces[cols*0+1] = 0;  elementfaces[cols*0+2] = 1; elementfaces[cols*0+3] = 2;
			elementfaces[cols*1+0] = 3; elementfaces[cols*1+1] = 0;  elementfaces[cols*1+2] = 3; elementfaces[cols*1+3] = 1;
			elementfaces[cols*2+0] = 3; elementfaces[cols*2+1] = 1;  elementfaces[cols*2+2] = 3; elementfaces[cols*2+3] = 2;
			elementfaces[cols*3+0] = 3; elementfaces[cols*3+1] = 0;  elementfaces[cols*3+2] = 2; elementfaces[cols*3+3] = 3;
			break;
		default:
		_error_("mesh "<< EnumToStringx(iomodel->meshelementtype) <<" not supported");
	}

	/*Allocate connectivity*/
	int *element_face_connectivity  = xNew<int>(iomodel->numberofelements*elementnbf); /*format: [face1 face2 ...] */
	int *element_vface_connectivity = NULL;
	if(iomodel->meshelementtype==PentaEnum){
		element_vface_connectivity = xNew<int>(iomodel->numberofelements*3); /*format: [face1 face2 face3] */
	}

	/*Maximum number of faces for initial allocation*/
	int maxnbf     = elementnbf*iomodel->numberofelements;
	int facescols  = 4+facemaxnbv; _assert_(facescols>6);

	/*Initialize intermediaries*/
	int* facestemp  = xNew<int>(maxnbf*facescols);        /*format: [element1 element2 marker nbv vertex1 vertex2 vertex3 ...]    */
	int* vfacestemp = xNew<int>(maxnbf*4);
	for(i=0;i<maxnbf;i++) facestemp[i*facescols+1]=-1;   /*Initialize second column of faces as -1 (boundary face)               */

	/*Initialize chain*/
	int* head_minv = xNew<int>(iomodel->numberofvertices);
	int* next_face = xNew<int>(maxnbf);
	for(i=0;i<iomodel->numberofvertices;i++) head_minv[i]=-1;

	/*Initialize number of faces and list of vertex indices*/
	int nbf = 0;
	int* v = xNew<int>(facemaxnbv);
	for(i=0;i<iomodel->numberofelements;i++){
		for(j=0;j<elementnbf;j++){

			/*Get indices of current face*/
			facenbv = elementfaces[cols*j+0];
			for(k=0;k<facenbv;k++){
				v[k] = iomodel->elements[i*elementnbv + elementfaces[cols*j+k+1]] - 1;
			}

			/*Sort list of vertices*/
			HeapSort(v,elementfaces[cols*j+0]);
			v0 = v[0]; _assert_(v0>=0 && v0<iomodel->numberofvertices);

			/*This face a priori has not been processed yet*/
			bool exist = false;

			/*Go through all processed faces connected to v0 and check whether we have seen this face yet*/
			for(int f=head_minv[v0]; f!=-1; f=next_face[f]){
				if(facestemp[f*facescols+5]==v[1]+1 && facestemp[f*facescols+6]==v[2]+1){
					exist = true;
					facestemp[f*facescols+1]=i+1;
					element_face_connectivity[i*elementnbf+j]=f;
					break;
				}
			}

			/*If this face is new, add it to the lists*/
			if(!exist){
				_assert_(nbf<maxnbf);

				/*Update faces*/
				facestemp[nbf*facescols+0] = i+1;
				facestemp[nbf*facescols+2] = elementfaces_markers[j];
				facestemp[nbf*facescols+3] = facenbv;
				for(k=0;k<facenbv;k++) facestemp[nbf*facescols+4+k] = v[k]+1;

				/*Update Connectivity*/
				element_face_connectivity[i*elementnbf+j]=nbf;

				/*Update chain*/
				next_face[nbf] = head_minv[v0];
				head_minv[v0]  = nbf;

				/*Increase number of faces*/
				nbf++;
			}
		}
	}

	/*Vertical faces*/
	int nbvf = 0;
	if(iomodel->meshelementtype==PentaEnum){
		for(i=0;i<iomodel->numberofvertices;i++) head_minv[i]=-1;
		for(i=0;i<iomodel->numberofelements;i++){
			for(j=2;j<5;j++){
				for(k=0;k<4;k++) v[k] = iomodel->elements[i*elementnbv + elementfaces[cols*j+k+1]] - 1;
				HeapSort(v,4);
				v0 = v[0]; _assert_(v0>=0 && v0<iomodel->numberofvertices);
				bool exist = false;
				for(int f=head_minv[v0]; f!=-1; f=next_face[f]){
					if(vfacestemp[f*4+1]==v[1]+1 && vfacestemp[f*4+2]==v[2]+1){
						exist = true;
						element_vface_connectivity[i*3+(j-2)]=f;
						break;
					}
				}
				if(!exist){ _assert_(nbvf<maxnbf);
					for(k=0;k<4;k++) vfacestemp[nbvf*4+k] = v[k]+1;
					element_vface_connectivity[i*3+(j-2)]=nbvf;
					next_face[nbvf] = head_minv[v0];
					head_minv[v0]   = nbvf;
					nbvf++;
				}
			}
		}
	}

	/*Clean up*/
	xDelete<int>(head_minv);
	xDelete<int>(next_face);
	xDelete<int>(v);
	xDelete<int>(elementfaces);
	xDelete<int>(elementfaces_markers);

	/*Create final faces (now that we have the correct size)*/
	int* faces = xNew<int>(nbf*facescols);
	xMemCpy<int>(faces,facestemp,nbf*facescols);
	xDelete<int>(facestemp);
	int* vfaces = xNew<int>(nbvf*4);
	xMemCpy<int>(vfaces,vfacestemp,nbvf*4);
	xDelete<int>(vfacestemp);

	/*Assign output pointers*/
	iomodel->faces                             = faces;
	iomodel->verticalfaces                     = vfaces;
	iomodel->numberoffaces                     = nbf;
	iomodel->numberofverticalfaces             = nbvf;
	iomodel->facescols                         = facescols;
	iomodel->elementtofaceconnectivity         = element_face_connectivity;
	iomodel->elementtoverticalfaceconnectivity = element_vface_connectivity;
}/*}}}*/
