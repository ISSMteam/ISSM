/*!\file Trianglex
 * \brief: x code for Triangle mesher
 */

/*Header files: {{{*/
#include "./Trianglex.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
/*ANSI_DECLARATORS needed to call triangle library: */
#if defined(_HAVE_TRIANGLE_)
	#ifndef ANSI_DECLARATORS
	#define ANSI_DECLARATORS
	#endif
	#include "triangle.h"
	#undef ANSI_DECLARATORS
#endif
/*}}}*/

void Trianglex(int** pindex,IssmPDouble** px,IssmPDouble** py,int** psegments,int** psegmentmarkerlist,int* pnels,int* pnods, int* pnsegs,Contours* domain,Contours* rifts,double area){

#if !defined(_HAVE_TRIANGLE_)
	_error_("triangle has not been installed");
#else
	/*indexing: */
	int i,j;

	/*output: */
	int    *index             = NULL;
	double *x                 = NULL;
	double *y                 = NULL;
	int    *segments          = NULL;
	int    *segmentmarkerlist = NULL;

	/*intermediary: */
	int counter,counter2,backcounter;
	Contour<IssmPDouble> *contour = NULL;

	/* Triangle structures needed to call Triangle library routines: */
	struct triangulateio in,out;
	char   options[256];

	/*Create initial triangulation to call triangulate(). First number of points:*/
	in.numberofpoints=0;
	for (i=0;i<domain->Size();i++){
		contour=(Contour<IssmPDouble>*)domain->GetObjectByOffset(i);
		in.numberofpoints+=contour->nods-1;
	}
	for (i=0;i<rifts->Size();i++){
		contour=(Contour<IssmPDouble>*)rifts->GetObjectByOffset(i);
		in.numberofpoints+=contour->nods;
	}

	/*number of point attributes: */
	in.numberofpointattributes=1;

	/*fill in the point list: */
	in.pointlist = xNew<REAL>(in.numberofpoints*2);

	counter=0;
	for (i=0;i<domain->Size();i++){
		contour=(Contour<IssmPDouble>*)domain->GetObjectByOffset(i);
		for (j=0;j<contour->nods-1;j++){
			in.pointlist[2*counter+0]=contour->x[j];
			in.pointlist[2*counter+1]=contour->y[j];
			counter++;
		}
	}
	for (i=0;i<rifts->Size();i++){
		contour=(Contour<IssmPDouble>*)rifts->GetObjectByOffset(i);
		for (j=0;j<contour->nods;j++){
			in.pointlist[2*counter+0]=contour->x[j];
			in.pointlist[2*counter+1]=contour->y[j];
			counter++;
		}
	}

	/*fill in the point attribute list: */
	in.pointattributelist = xNew<REAL>(in.numberofpoints*in.numberofpointattributes);
	for (i=0;i<in.numberofpoints;i++) in.pointattributelist[i] = 0.0;

	/*fill in the point marker list: */
	in.pointmarkerlist = xNew<int>(in.numberofpoints);
	for(i=0;i<in.numberofpoints;i++) in.pointmarkerlist[i] = 0;

	/*Build segments. First figure out number of segments: holes and closed outlines have as many segments as vertices: */
	in.numberofsegments=0;
	for (i=0;i<domain->Size();i++){
		contour=(Contour<IssmPDouble>*)domain->GetObjectByOffset(i);
		in.numberofsegments+=contour->nods-1;
	}
	for(i=0;i<rifts->Size();i++){
		contour=(Contour<IssmPDouble>*)rifts->GetObjectByOffset(i);
		/*for rifts, we have one less segment as we have vertices*/
		in.numberofsegments+=contour->nods-1;
	}

	in.segmentlist = xNew<int>(in.numberofsegments*2);
	in.segmentmarkerlist = xNewZeroInit<int>(in.numberofsegments);
	counter=0;
	backcounter=0;
	for (i=0;i<domain->Size();i++){
		contour=(Contour<IssmPDouble>*)domain->GetObjectByOffset(i);
		for (j=0;j<contour->nods-2;j++){
			in.segmentlist[2*counter+0]=counter;
			in.segmentlist[2*counter+1]=counter+1;
			in.segmentmarkerlist[counter]=0;
			counter++;
		}
		/*Close this profile: */
		 in.segmentlist[2*counter+0]=counter;
		 in.segmentlist[2*counter+1]=backcounter;
		 in.segmentmarkerlist[counter]=0;
		 counter++;
		 backcounter=counter;
	}
	counter2=counter;
	for (i=0;i<rifts->Size();i++){
		contour=(Contour<IssmPDouble>*)rifts->GetObjectByOffset(i);
		for (j=0;j<(contour->nods-1);j++){
			in.segmentlist[2*counter2+0]=counter;
			in.segmentlist[2*counter2+1]=counter+1;
			in.segmentmarkerlist[counter2]=2+i;
			counter2++;
			counter++;
		}
		counter++;
	}

	/*Build regions: */
	in.numberofregions = 0;

	/*Build holes: */
	in.numberofholes = domain->Size()-1; /*everything is a hole, but for the first profile.*/
	if(in.numberofholes){
		in.holelist = xNew<REAL>(in.numberofholes*2);
		for (i=0;i<domain->Size()-1;i++){
			contour=(Contour<IssmPDouble>*)domain->GetObjectByOffset(i+1);
			GridInsideHole(&in.holelist[2*i+0],&in.holelist[2*i+1],contour->nods-1,contour->x,contour->y);
		}
	}

	/* Make necessary initializations so that Triangle can return a triangulation in `out': */
	out.pointlist             = (REAL*)NULL;
	out.pointattributelist    = (REAL*)NULL;
	out.pointmarkerlist       = (int *)NULL;
	out.trianglelist          = (int *)NULL;
	out.triangleattributelist = (REAL*)NULL;
	out.neighborlist          = (int *)NULL;
	out.segmentlist           = (int *)NULL;
	out.segmentmarkerlist     = (int *)NULL;
	out.edgelist              = (int *)NULL;
	out.edgemarkerlist        = (int *)NULL;

	/* Triangulate the points:.  Switches are chosen to read and write a  */
	/*   PSLG (p), preserve the convex hull (c), number everything from  */
	/*   zero (z), assign a regional attribute to each element (A), and  */
	/*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
	/*   neighbor list (n).                                              */
	snprintf(options, sizeof(options),"%s%lf","pQzDq30ia",area); /*replace V by Q to quiet down the logging*/
	triangulate(options, &in, &out, NULL);
	/*report(&out, 0, 1, 1, 1, 1, 0);*/

	/*Allocate index, x and y: */
	index=xNew<int>(3*out.numberoftriangles);
	x=xNew<double>(out.numberofpoints);
	y=xNew<double>(out.numberofpoints);
	segments=xNew<int>(3*out.numberofsegments);
	segmentmarkerlist=xNew<int>(out.numberofsegments);

	for (i = 0; i< out.numberoftriangles; i++) {
		for (j = 0; j < out.numberofcorners; j++) {
			index[3*i+j]=(int)out.trianglelist[i*out.numberofcorners+j]+1;
		}
	}
	for (i = 0; i< out.numberofpoints; i++){
		x[i]=(double)out.pointlist[i*2+0];
		y[i]=(double)out.pointlist[i*2+1];
	}
	for (i = 0; i<out.numberofsegments;i++){
		segments[3*i+0]=(int)out.segmentlist[i*2+0]+1;
		segments[3*i+1]=(int)out.segmentlist[i*2+1]+1;
		segmentmarkerlist[i]=(int)out.segmentmarkerlist[i];
	}

	/*Associate elements with segments: */
	AssociateSegmentToElement(&segments,out.numberofsegments,index,out.numberoftriangles);

	/*Order segments so that their normals point outside the domain: */
	OrderSegments(&segments,out.numberofsegments, index,out.numberoftriangles);

	/*Output : */
	*pindex=index;
	*px=x;
	*py=y;
	*psegments=segments;
	*psegmentmarkerlist=segmentmarkerlist;
	*pnels=out.numberoftriangles;
	*pnods=out.numberofpoints;
	*pnsegs=out.numberofsegments;
#endif
}
