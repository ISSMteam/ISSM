/*!\file Bamgx
 * \brief: use Bamg capabilities.
 */
#include "./Bamgx.h"
#include "../../bamg/bamgobjects.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

using namespace bamg;
using namespace std;

int Bamgx(BamgMesh* bamgmesh_out,BamgGeom* bamggeom_out,BamgMesh* bamgmesh_in,BamgGeom* bamggeom_in,BamgOpts* bamgopts){

	/*Bamg options*/
	int    maxnbv;
	double coef;
	int    verbosity;
	int    nbsmooth;

	/*intermediary*/
	int i;
	int noerr=1;
	double hminaniso=1e-100;
	Mesh* Thr=NULL;
	Mesh* Thb=NULL;

	/*Bamg options*/
	nbsmooth =bamgopts->nbsmooth;
	coef     =bamgopts->coeff;
	maxnbv   =bamgopts->maxnbv;
	verbosity=bamgopts->verbose;

	// no metric -> no smoothing
	if (bamgopts->metric==NULL) nbsmooth=0;

	/*If no mesh in input, generate one*/
	if(bamgmesh_in->TrianglesSize[0]==0){
		/*Mesh generation {{{*/

		//Step1: generate geometry Gh
		if (verbosity>0) _printf_("Construction of a mesh from a given geometry\n");
		if (verbosity>1) _printf_("   Processing geometry...\n");
		Geometry Gh(bamggeom_in,bamgopts);

		//get hmin and hmax from geometry to generate the metric
		bamgopts->hmin = Max(bamgopts->hmin,Gh.MinimalHmin());
		bamgopts->hmax = Min(bamgopts->hmax,Gh.MaximalHmax());

		//build metric using geometry
		if (verbosity>1) _printf_("   Generating Metric...\n");
		for(i=0;i<Gh.nbv;i++){
			Metric M=Gh[i];
			EigenMetric Vp(M/coef);
			Vp.Maxh(bamgopts->hmax);
			Vp.Minh(bamgopts->hmin);
			Gh.vertices[i].m = Vp;
		}

		//generate mesh
		if (verbosity>1) _printf_("   Generating Mesh...\n");
		Mesh Th(maxnbv,Gh,bamgopts);

		//Split corners if requested
		if(bamgopts->splitcorners) Th.SplitInternalEdgeWithBorderVertices();

		//Renumbering
		Th.TrianglesRenumberBySubDomain();

		//Crack mesh if requested
		if(bamgopts->Crack) Th.CrackMesh(bamgopts);

		//Build output
		if (verbosity>1) _printf_("   Write Mesh...\n");
		Th.WriteMesh(bamgmesh_out,bamgopts);
		if (verbosity>1) _printf_("   Write Geometry...\n");
		Gh.WriteGeometry(bamggeom_out,bamgopts);

		//clean up
	//	delete &Th;
	//	delete &Gh;
		/*}}}*/
	}
	else{
		/*Anisotropic mesh adaptation {{{*/

		// read background mesh
		if (verbosity>0) _printf_("Anisotropic mesh adaptation\n");
		if (verbosity>1) _printf_("   Processing initial mesh and geometry...\n");
		Mesh BTh(bamggeom_in,bamgmesh_in,bamgopts);

		//Make Quadtree from background mesh
		BTh.MakeBamgQuadtree();

		//Bound hmin and hmax
		bamgopts->hmin=Max(bamgopts->hmin,BTh.MinimalHmin());
		bamgopts->hmax=Min(bamgopts->hmax,BTh.MaximalHmax());

		//Generate initial metric
		if (bamgopts->metric){
			if (verbosity>1) _printf_("   Processing Metric...\n");
			BTh.ReadMetric(bamgopts);
		}
		else {
			if (verbosity>1) _printf_("   Generating initial Metric...\n");
			Metric Mhmax(bamgopts->hmax);
			for (int iv=0;iv<BTh.nbv;iv++) BTh[iv].m = Mhmax;
		}

		//use present fields to generate metric if present
		if (bamgopts->field){
			if (verbosity>1) _printf_("   Merge metric with field provided...\n");
			BTh.AddMetric(bamgopts);
		}

		// change using hVertices if provided
		if(bamgopts->hVertices && bamgopts->hVerticesLength==BTh.nbv){
			if (verbosity>1) _printf_("   Merging Metric with hVertices...\n");
			for (i=0;i<BTh.nbv;i++){
				if (!xIsNan<IssmPDouble>(bamgopts->hVertices[i])){
					BTh[i].m=Metric((float)bamgopts->hVertices[i]);
				}
			}
		}

		// change using hminVertices if provided
		if (bamgopts->hminVertices){
			if (verbosity>1) _printf_("   Merging Metric with hminVertices...\n");
			for (i=0;i<BTh.nbv;i++){
				if (!xIsNan<IssmPDouble>(bamgopts->hminVertices[i])){
					Metric M=BTh.vertices[i].m;
					EigenMetric Vp(M/coef);
					Vp.Minh(bamgopts->hminVertices[i]);
					BTh.vertices[i].m=Vp;
				}
			}
		}

		// change using hmaxVertices if provided
		if (bamgopts->hmaxVertices){
			if (verbosity>1) _printf_("   Merging Metric with hmaxVertices...\n");
			for (i=0;i<BTh.nbv;i++){
				if (!xIsNan<IssmPDouble>(bamgopts->hmaxVertices[i])){
					Metric M=BTh.vertices[i].m;
					EigenMetric Vp(M/coef);
					Vp.Maxh(bamgopts->hmaxVertices[i]);
					BTh.vertices[i].m=Vp;
				}
			}
		}

		//Smoothe metric
		BTh.SmoothMetric(bamgopts,bamgopts->gradation);

		//Control element subdivision
		BTh.MaxSubDivision(bamgopts,bamgopts->maxsubdiv);

		//Bound anisotropy
		BTh.BoundAnisotropy(bamgopts,bamgopts->anisomax,hminaniso);

		//Build new mesh
		if (verbosity>1) _printf_("   Generating Mesh...\n");
		Thr=&BTh,Thb=0;
		Mesh & Th( *(0 ?  new Mesh(*Thr,&Thr->Gh,Thb,maxnbv) :  new Mesh(maxnbv,BTh,bamgopts,bamgopts->KeepVertices)));
		//if (Thr!=&BTh) delete Thr;

		//Split corners if requested
		if(bamgopts->splitcorners) Th.SplitInternalEdgeWithBorderVertices();

		//Renumber by subdomain
		Th.TrianglesRenumberBySubDomain();

		//Smooth vertices
		if(nbsmooth>0) Th.SmoothingVertex(bamgopts,nbsmooth,bamgopts->omega);

		//display info
		if(verbosity>0) {
			if (Th.nbt-Th.nbtout){
				_printf_("   new number of triangles = " << (Th.nbt-Th.nbtout) << "\n");
			}
		}

		//Build output
		if (verbosity>1) _printf_("   Write Mesh...\n");
		Th.WriteMesh(bamgmesh_out,bamgopts);
		if (verbosity>1) _printf_("   Write Geometry...\n");
		Th.Gh.WriteGeometry(bamggeom_out,bamgopts);
		if (verbosity>1) _printf_("   Write Metric...\n");
		BTh.WriteMetric(bamgopts);

		/*clean up*/
		delete &Th;
		//delete &BTh;
		/*}}}*/
	}

	/*No error return*/
	if (verbosity>1) _printf_("   Exiting Bamg.\n");
	return noerr;

}
