/*\file WritePythonData.cpp:
 *\brief: general I/O interface to write data in matlab
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#define PY_ARRAY_UNIQUE_SYMBOL PythonIOSymbol
#define NO_IMPORT

#include "./pythonio.h"
#include "../../c/shared/shared.h"

/*Primitive data types*/
/*FUNCTION WriteData(PyObject* py_tuple,int index,int integer){{{*/
void WriteData(PyObject* py_tuple, int index, int integer){

	PyTuple_SetItem(py_tuple, index, PyLong_FromSsize_t((Py_ssize_t)integer));

}/*}}}*/
/*FUNCTION WriteData(PyObject* py_tuple,int index,char* string){{{*/
void WriteData(PyObject* py_tuple, int index, char* string){

	PyTuple_SetItem(py_tuple, index, PyUnicode_FromString(string));

}/*}}}*/
/*FUNCTION WriteData(PyObject* py_tuple,int index, double* matrix, int M, int N){{{*/
void WriteData(PyObject* tuple, int index, double* matrix, int M,int N){

	npy_intp dims[2]={0,0};
	PyObject* array=NULL;

	/*copy matrix: */
	double* matrix_python=xNew<double>(M*N);
	memcpy(matrix_python,matrix,M*N*sizeof(double));

	dims[0]=(npy_intp)M;
	dims[1]=(npy_intp)N;
	array=PyArray_SimpleNewFromData(2,dims,NPY_DOUBLE,matrix_python);

	PyTuple_SetItem(tuple, index, array);
}/*}}}*/
/*FUNCTION WriteData(PyObject* py_tuple,int index, int* matrix, int M, int N){{{*/
void WriteData(PyObject* tuple, int index, int* matrix, int M,int N){

	npy_intp dims[2]={0,0};
	PyObject* array=NULL;

	/*transform into long matrix: */
	long* lmatrix=xNew<long>(M*N);
	for(int i=0;i<M*N;i++)lmatrix[i]=(long)matrix[i];

	dims[0]=(npy_intp)M;
	dims[1]=(npy_intp)N;
	array=PyArray_SimpleNewFromData(2,dims,NPY_INT64,lmatrix);

	PyTuple_SetItem(tuple, index, array);
}/*}}}*/
/*FUNCTION WriteData(PyObject* py_tuple,int index, bool* matrix, int M, int N){{{*/
void WriteData(PyObject* tuple, int index, bool* matrix, int M,int N){

	npy_intp dims[2]={0,0};
	PyObject* array=NULL;

	/*copy matrix: */
	bool* matrix_python=xNew<bool>(M*N);
	memcpy(matrix_python,matrix,M*N*sizeof(bool));

	dims[0]=(npy_intp)M;
	dims[1]=(npy_intp)N;
	array=PyArray_SimpleNewFromData(2,dims,NPY_BOOL,matrix_python);

	PyTuple_SetItem(tuple, index, array);
}/*}}}*/
/*FUNCTION WriteData(PyObject* py_tuple,int index, double* vector, int M){{{*/
void WriteData(PyObject* py_tuple, int index, double* vector, int M){

	npy_intp  dim   = 10;
	PyObject *array = NULL;

	/*copy vector: */
	double* vector_python=xNew<double>(M);
	memcpy(vector_python,vector,M*sizeof(double));

	dim=(npy_intp)M;
	array=PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,vector_python);

	PyTuple_SetItem(py_tuple, index, array);

}/*}}}*/
/*FUNCTION WriteData(PyObject* py_tuple,int index, short* vector, int M){{{*/
void WriteData(PyObject* py_tuple, int index,short* vector, int M){

	long* lvector=NULL;
	npy_intp dim=10;
	PyObject* array=NULL;

	/*transform into long matrix: */
	lvector=xNew<long>(M);
	for(int i=0;i<M;i++)lvector[i]=(long)vector[i];

	dim=(npy_intp)M;
	array=PyArray_SimpleNewFromData(1,&dim,NPY_INT64,lvector);

	PyTuple_SetItem(py_tuple, index, array);

}/*}}}*/
/*FUNCTION WriteData(PyObject* py_tuple,int index, int* vector, int M){{{*/
void WriteData(PyObject* py_tuple, int index, int* vector, int M){

	long* lvector=NULL;
	npy_intp dim=10;
	PyObject* array=NULL;

	/*transform into long matrix: */
	lvector=xNew<long>(M);
	for(int i=0;i<M;i++)lvector[i]=(long)vector[i];

	dim=(npy_intp)M;
	array=PyArray_SimpleNewFromData(1,&dim,NPY_INT64,lvector);

	PyTuple_SetItem(py_tuple, index, array);

}/*}}}*/
/*FUNCTION WriteData(PyObject* py_tuple,int index){{{*/
void WriteData(PyObject* py_tuple, int index){

	PyTuple_SetItem(py_tuple, index, Py_None);

}/*}}}*/

/*ISSM objects*/
/*FUNCTION WriteData(PyObject* py_tuple,int index,BamgGeom* bamggeom){{{*/
void WriteData(PyObject* py_tuple,int index,BamgGeom* bamggeom){

	PyObject* dict=NULL;

	dict=PyDict_New();

	PyDict_SetItemString(dict,"Vertices",PyArrayFromCopiedData(bamggeom->VerticesSize,bamggeom->Vertices));
	PyDict_SetItemString(dict,"Edges",PyArrayFromCopiedData(bamggeom->EdgesSize,bamggeom->Edges));
	PyDict_SetItemString(dict,"TangentAtEdges",PyArrayFromCopiedData(bamggeom->TangentAtEdgesSize,bamggeom->TangentAtEdges));
	PyDict_SetItemString(dict,"Corners",PyArrayFromCopiedData(bamggeom->CornersSize,bamggeom->Corners));
	PyDict_SetItemString(dict,"RequiredVertices",PyArrayFromCopiedData(bamggeom->RequiredVerticesSize,bamggeom->RequiredVertices));
	PyDict_SetItemString(dict,"RequiredEdges",PyArrayFromCopiedData(bamggeom->RequiredEdgesSize,bamggeom->RequiredEdges));
	PyDict_SetItemString(dict,"CrackedEdges",PyArrayFromCopiedData(bamggeom->CrackedEdgesSize,bamggeom->CrackedEdges));
	PyDict_SetItemString(dict,"SubDomains",PyArrayFromCopiedData(bamggeom->SubDomainsSize,bamggeom->SubDomains));

	PyTuple_SetItem(py_tuple, index, dict);
}
/*}}}*/
/*FUNCTION WriteData(PyObject* py_tuple,int index,BamgMesh* bamgmesh){{{*/
void WriteData(PyObject* py_tuple,int index,BamgMesh* bamgmesh){

	PyObject* dict=NULL;

	dict=PyDict_New();

	PyDict_SetItemString(dict,"Vertices",PyArrayFromCopiedData(bamgmesh->VerticesSize,bamgmesh->Vertices));
	PyDict_SetItemString(dict,"Edges",PyArrayFromCopiedData(bamgmesh->EdgesSize,bamgmesh->Edges));
	PyDict_SetItemString(dict,"Triangles",PyArrayFromCopiedData(bamgmesh->TrianglesSize,bamgmesh->Triangles));
	PyDict_SetItemString(dict,"IssmEdges",PyArrayFromCopiedData(bamgmesh->IssmEdgesSize,bamgmesh->IssmEdges));
	PyDict_SetItemString(dict,"IssmSegments",PyArrayFromCopiedData(bamgmesh->IssmSegmentsSize,bamgmesh->IssmSegments));
	PyDict_SetItemString(dict,"VerticesOnGeomVertex",PyArrayFromCopiedData(bamgmesh->VerticesOnGeomVertexSize,bamgmesh->VerticesOnGeomVertex));
	PyDict_SetItemString(dict,"VerticesOnGeomEdge",PyArrayFromCopiedData(bamgmesh->VerticesOnGeomEdgeSize,bamgmesh->VerticesOnGeomEdge));
	PyDict_SetItemString(dict,"EdgesOnGeomEdge",PyArrayFromCopiedData(bamgmesh->EdgesOnGeomEdgeSize,bamgmesh->EdgesOnGeomEdge));
	PyDict_SetItemString(dict,"SubDomains",PyArrayFromCopiedData(bamgmesh->SubDomainsSize,bamgmesh->SubDomains));
	PyDict_SetItemString(dict,"SubDomainsFromGeom",PyArrayFromCopiedData(bamgmesh->SubDomainsFromGeomSize,bamgmesh->SubDomainsFromGeom));
	PyDict_SetItemString(dict,"ElementConnectivity",PyArrayFromCopiedData(bamgmesh->ElementConnectivitySize,bamgmesh->ElementConnectivity));
	PyDict_SetItemString(dict,"NodalConnectivity",PyArrayFromCopiedData(bamgmesh->NodalConnectivitySize,bamgmesh->NodalConnectivity));
	PyDict_SetItemString(dict,"NodalElementConnectivity",PyArrayFromCopiedData(bamgmesh->NodalElementConnectivitySize,bamgmesh->NodalElementConnectivity));
	PyDict_SetItemString(dict,"CrackedVertices",PyArrayFromCopiedData(bamgmesh->CrackedVerticesSize,bamgmesh->CrackedVertices));
	PyDict_SetItemString(dict,"CrackedEdges",PyArrayFromCopiedData(bamgmesh->CrackedEdgesSize,bamgmesh->CrackedEdges));

	PyTuple_SetItem(py_tuple, index, dict);
}
/*}}}*/
/*FUNCTION WriteData(PyObject* py_tuple,int index,IssmDenseMat<double>* matrix){{{*/
void WriteData(PyObject* py_tuple,int index,IssmDenseMat<double>* matrix){

	int M,N;
	double* buffer=NULL;
	npy_intp dims[2]={0,0};
	PyObject* array=NULL;

	matrix->GetSize(&M,&N);
	buffer=matrix->ToSerial();

	dims[0]=(npy_intp)M;
	dims[1]=(npy_intp)N;
	array=PyArray_SimpleNewFromData(2,dims,NPY_DOUBLE,buffer);

	PyTuple_SetItem(py_tuple, index, array);

}/*}}}*/
/*FUNCTION WriteData(PyObject* py_tuple,int index,IssmSeqVec<double>* vector){{{*/
void WriteData(PyObject* py_tuple,int index,IssmSeqVec<double>* vector){

	int M;
	double* buffer=NULL;
	npy_intp dim=10;
	PyObject* array=NULL;

	vector->GetSize(&M);
	buffer=vector->ToMPISerial();

	dim=(npy_intp)M;
	array=PyArray_SimpleNewFromData(1,&dim,NPY_DOUBLE,buffer);

	PyTuple_SetItem(py_tuple, index, array);
}
/*}}}*/
/*FUNCTION WriteData(PyObject* py_tuple,int index,RiftStruct* riftstruct){{{*/
void WriteData(PyObject* py_tuple,int index,RiftStruct* riftstruct){

	int i;
	PyObject* list=NULL;
	PyObject* dict=NULL;

	list=PyList_New((Py_ssize_t)0);

	for (i=0; i<riftstruct->numrifts; i++) {
		dict=PyDict_New();

		PyDict_SetItemString(dict,"numsegs"          ,PyLong_FromSsize_t((Py_ssize_t)riftstruct->riftsnumsegments[i]));
		PyDict_SetItemString(dict,"fill"             ,PyUnicode_FromString("Ice"));
		PyDict_SetItemString(dict,"friction"         ,PyLong_FromSsize_t((Py_ssize_t)0));

		PyDict_SetItemString(dict,"segments"         ,PyArrayFromCopiedData(riftstruct->riftsnumsegments[i]    ,3,riftstruct->riftssegments[i]));
		PyDict_SetItemString(dict,"pairs"            ,PyArrayFromCopiedData(riftstruct->riftsnumpairs[i]       ,2,riftstruct->riftspairs[i]));
		PyDict_SetItemString(dict,"tips"             ,PyArrayFromCopiedData(1                                  ,2,&riftstruct->riftstips[2*i]));
		PyDict_SetItemString(dict,"penaltypairs"     ,PyArrayFromCopiedData(riftstruct->riftsnumpenaltypairs[i],7,riftstruct->riftspenaltypairs[i]));
		PyDict_SetItemString(dict,"fraction"         ,PyFloat_FromDouble(0.));
		PyDict_SetItemString(dict,"fractionincrement",PyFloat_FromDouble(0.1));
		PyDict_SetItemString(dict,"state"            ,PyArrayFromCopiedData(riftstruct->riftsnumpenaltypairs[i],1,riftstruct->state[i]));

		PyList_Append(list, dict);
	}

	PyTuple_SetItem(py_tuple, index, list);
}
/*}}}*/

/*Utils*/
/*FUNCTION PyArrayFromCopiedData(int dims[2],double* data){{{*/
PyObject* PyArrayFromCopiedData(int dims[2],double* data){

	double* pydata;
	npy_intp pydims[2]={0,0};

	/*  note that PyArray_SimpleNewFromData does not copy the data, so that when the original
		 object (e.g. bamggeom,bamgmesh) is deleted, the data is gone.  */

	pydims[0]=(npy_intp)dims[0];
	pydims[1]=(npy_intp)dims[1];
	pydata=xNew<IssmDouble>(dims[0]*dims[1]);
	memcpy(pydata,data,dims[0]*dims[1]*sizeof(double));
	return PyArray_SimpleNewFromData(2,pydims,NPY_DOUBLE,pydata);
}
/*}}}*/
/*FUNCTION PyArrayFromCopiedData(int dimi,int dimj,double* data){{{*/
PyObject* PyArrayFromCopiedData(int dimi,int dimj,double* data){

	double* pydata;
	npy_intp pydims[2]={0,0};

	/*  note that PyArray_SimpleNewFromData does not copy the data, so that when the original
		 object (e.g. bamggeom,bamgmesh) is deleted, the data is gone.  */

	pydims[0]=(npy_intp)dimi;
	pydims[1]=(npy_intp)dimj;
	pydata=xNew<IssmDouble>(dimi*dimj);
	memcpy(pydata,data,dimi*dimj*sizeof(double));
	return PyArray_SimpleNewFromData(2,pydims,NPY_DOUBLE,pydata);
}
/*}}}*/
/*FUNCTION PyArrayFromCopiedData(int dimi,int dimj,int* data){{{*/
PyObject* PyArrayFromCopiedData(int dimi,int dimj,int* data){

	long* pydata;
	npy_intp pydims[2]={0,0};

	/*  note that PyArray_SimpleNewFromData does not copy the data, so that when the original
		 object (e.g. bamggeom,bamgmesh) is deleted, the data is gone.  */

	pydims[0]=(npy_intp)dimi;
	pydims[1]=(npy_intp)dimj;
	pydata=xNew<long>(dimi*dimj);
	for(int i=0;i<dimi*dimj;i++) pydata[i]=(long)data[i];
	return PyArray_SimpleNewFromData(2,pydims,NPY_INT64,pydata);
}
/*}}}*/
/*FUNCTION PyArrayFromCopiedData(int dimi,int dimj,bool* data){{{*/
PyObject* PyArrayFromCopiedData(int dimi,int dimj,bool* data){

	bool* pydata;
	npy_intp pydims[2]={0,0};

	/*  note that PyArray_SimpleNewFromData does not copy the data, so that when the original
		 object (e.g. bamggeom,bamgmesh) is deleted, the data is gone.  */

	pydims[0]=(npy_intp)dimi;
	pydims[1]=(npy_intp)dimj;
	pydata=xNew<bool>(dimi*dimj);
	memcpy(pydata,data,dimi*dimj*sizeof(bool));
	return PyArray_SimpleNewFromData(2,pydims,NPY_BOOL,pydata);
}
/*}}}*/
