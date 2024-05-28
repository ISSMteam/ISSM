/*!\file AmrBamg.cpp
 * \brief: implementation of the adaptive mesh refinement tool based on bamg
 */

#ifdef HAVE_CONFIG_H
    #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./AmrBamg.h"
#include "../bamg/bamgobjects.h"
#include "../modules/Bamgx/Bamgx.h"

using namespace bamg;
using namespace std;

/*Constructor, copy, clean up and destructor*/
AmrBamg::AmrBamg(){/*{{{*/

	/*These attributes MUST be setup by FemModel*/
	this->fieldenum								= -1;
	this->keepmetric								= -1;
	this->groundingline_resolution			= -1;
	this->groundingline_distance				= -1;
	this->icefront_resolution					= -1;
	this->icefront_distance						= -1;
	this->thicknesserror_resolution			= -1;
	this->thicknesserror_threshold			= -1;
	this->thicknesserror_groupthreshold 	= -1;
	this->thicknesserror_maximum				= -1;
	this->deviatoricerror_resolution			= -1;
	this->deviatoricerror_threshold			= -1;
	this->deviatoricerror_groupthreshold	= -1;
	this->deviatoricerror_maximum				= -1;

	/*Geometry and mesh as NULL*/
	this->geometry									= NULL;
	this->fathermesh								= NULL;
	this->previousmesh							= NULL;
	this->elementslist							= NULL;
	this->x											= NULL;
	this->y											= NULL;
	this->numberofvertices						= -1;
	this->numberofelements						= -1;

	/*Only initialize options for now (same as bamg.m)*/
	this->options									= new BamgOpts();
	this->options->anisomax						= 10.e30;
	this->options->cutoff						= 10.e-5;
	this->options->coeff							= 1;
	this->options->errg							= 0.1;
	this->options->gradation					= -1; //MUST be setup by the FemModel 
	this->options->Hessiantype					= 0;
	this->options->maxnbv						= 1e6;
	this->options->maxsubdiv					= 10;
	this->options->Metrictype					= 0;
	this->options->nbjacobi						= 1;
	this->options->nbsmooth						= 3;
	this->options->omega							= 1.8;
	this->options->power							= 1;
	this->options->verbose						= 0;
	this->options->Crack							= 0;
	this->options->KeepVertices				= 1; /*!!!!! VERY IMPORTANT !!!!! This avoid numerical errors when remeshing*/
	this->options->splitcorners				= 1;
	this->options->hmin							= -1;/*MUST be setup by the FemModel*/
	this->options->hmax							= -1;/*MUST be setup by the FemModel*/
	this->options->err							= xNew<IssmDouble>(1);
	this->options->err[0]						= -1;/*MUST be setup by the FemModel*/
	this->options->errSize[0]					= 1;
	this->options->errSize[1]					= 1;
}
/*}}}*/
AmrBamg::~AmrBamg(){/*{{{*/

	if(this->geometry) delete this->geometry;
	if(this->fathermesh) delete this->fathermesh;
	if(this->previousmesh) delete this->previousmesh;
	if(this->options) delete this->options;
	if(this->x) xDelete<IssmDouble>(this->x);
	if(this->y) xDelete<IssmDouble>(this->y);
	if(this->elementslist) xDelete<int>(this->elementslist);
}
/*}}}*/

/*Methods*/
void AmrBamg::SetMesh(int** elementslist_in,IssmDouble** x_in,IssmDouble** y_in,int* numberofvertices_in,int* numberofelements_in){/*{{{*/

	/*Delete previous mesh and keep the entire mesh*/
	if(this->elementslist) xDelete<int>(this->elementslist);
	if(this->x) xDelete<IssmDouble>(this->x);
	if(this->y) xDelete<IssmDouble>(this->y);

	this->elementslist		= *elementslist_in;
	this->x						= *x_in;
	this->y						= *y_in;
	this->numberofvertices	= *numberofvertices_in;
	this->numberofelements	= *numberofelements_in;
}/*}}}*/
void AmrBamg::GetMesh(int** elementslist_out,IssmDouble** x_out,IssmDouble** y_out,int* numberofvertices_out,int* numberofelements_out){/*{{{*/

	/*Get the entire mesh*/
	*elementslist_out		= this->elementslist;
	*x_out					= this->x;
	*y_out					= this->y;
	*numberofvertices_out= this->numberofvertices;
	*numberofelements_out= this->numberofelements;
}/*}}}*/
void AmrBamg::Initialize(){/*{{{*/

	/*Check options*/
	_assert_(this->options);
	this->options->Check();

	/*Read father mesh and create geometry*/
	Mesh* Th=new Mesh(this->elementslist,this->x,this->y,this->numberofvertices,this->numberofelements,this->options);

	/*Write geometry*/
	this->geometry = new BamgGeom();
	Th->Gh.WriteGeometry(this->geometry,this->options);

	/*Write father mesh*/
	this->fathermesh = new BamgMesh();
	Th->WriteMesh(this->fathermesh,this->options);

	/*Cleanup and return*/
	delete Th;
}/*}}}*/
void AmrBamg::ExecuteRefinementBamg(IssmDouble* field,IssmDouble* hmaxVertices,int** pdatalist,IssmDouble** pxylist,int** pelementslist){/*{{{*/

	/*Intermediaries*/
	BamgGeom* geomout=new BamgGeom();
	BamgMesh* meshout=new BamgMesh();

	/*Some checks*/
	_assert_(this->geometry);
	_assert_(this->options);
	_assert_(this->fathermesh);
	_assert_(field || hmaxVertices);//at least one is necessary

	/*Prepare field for metric*/
	this->options->field			 = field;
	this->options->hmaxVertices = hmaxVertices;

	/*remesh*/
	if(this->previousmesh){
		this->options->fieldSize[0]			= this->previousmesh->VerticesSize[0];
		this->options->fieldSize[1]			= 1;
		this->options->hmaxVerticesSize[0]	= this->previousmesh->VerticesSize[0];
		this->options->hmaxVerticesSize[1]	= 1;
		Bamgx(meshout,geomout,this->previousmesh,this->geometry,this->options);
	}
	else{
		this->options->fieldSize[0]			= this->fathermesh->VerticesSize[0];
		this->options->fieldSize[1]			= 1;
		this->options->hmaxVerticesSize[0]	= this->fathermesh->VerticesSize[0];
		this->options->hmaxVerticesSize[1]	= 1;
		Bamgx(meshout,geomout,this->fathermesh,this->geometry,this->options);
	}

	/*remove field and hmaxVertices for memory management (FemModel is taking care of deleting it)*/
	this->options->field = NULL;
	this->options->fieldSize[0] = 0;
	this->options->fieldSize[1] = 0;
	this->options->hmaxVertices = NULL;
	this->options->hmaxVerticesSize[0] = 0;
	this->options->hmaxVerticesSize[1] = 0;

	/*verify if the metric will be reseted or not*/
	if(this->keepmetric==0){
		if(this->options->metric) xDelete<IssmDouble>(this->options->metric);
		this->options->metricSize[0] = 0;
		this->options->metricSize[1] = 0;
	}

	/*Change previous mesh*/
	if(this->previousmesh) delete this->previousmesh;
	this->previousmesh = meshout;

	/*Prepare output*/
	int nbv				= meshout->VerticesSize[0];
	int nbt				= meshout->TrianglesSize[0];
	int *datalist		= xNew<int>(2);
	IssmDouble *xylist= xNew<IssmDouble>(nbv*2);
	int* elementslist = xNew<int>(nbt*3);

	datalist[0] = nbv;
	datalist[1] = nbt;

	for(int i=0;i<nbv;i++){
		xylist[2*i]		= meshout->Vertices[i*3+0];
		xylist[2*i+1]	= meshout->Vertices[i*3+1];
	}

	for(int i=0;i<nbt;i++){
		elementslist[3*i+0] = reCast<int>(meshout->Triangles[4*i+0]);
		elementslist[3*i+1] = reCast<int>(meshout->Triangles[4*i+1]);
		elementslist[3*i+2] = reCast<int>(meshout->Triangles[4*i+2]);
	}

	/*Assign pointers*/
	*pdatalist		= datalist;
	*pxylist			= xylist;
	*pelementslist = elementslist;

	/*Cleanup and return*/
	delete geomout;
}/*}}}*/
void AmrBamg::SetBamgOpts(IssmDouble hmin_in,IssmDouble hmax_in,IssmDouble err_in,IssmDouble gradation_in){/*{{{*/

	if(!this->options) _error_("AmrBamg->options is NULL!");

	this->options->hmin     = hmin_in; 
	this->options->hmax     = hmax_in; 
	this->options->err[0]	= err_in; 
	this->options->gradation= gradation_in; 
}/*}}}*/
