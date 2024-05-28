/*! \file  MeshProfileIntersectionx.c
 */

#include "./MeshProfileIntersectionx.h"

void MeshProfileIntersectionx(double** psegments, int* pnumsegs, int* index, double* x, double* y, int nel, int nods,  Contour<IssmPDouble>** contours,int numcontours){/*{{{*/

	int i,j,k;

	/*Contour:*/
	Contour<IssmPDouble>* contouri=NULL;
	int      numnodes;
	double*  xc=NULL;
	double*  yc=NULL;

	/*output: */
	double* segments=NULL;
	int     numsegs;

	/*intermediary: */
	double** allsegments=NULL;
	double*  segmentsi=NULL;
	int*     allnumsegs=NULL;
	int      numsegsi;
	int      count;

	/*Allocate: */
	allsegments=xNew<double*>(numcontours);
	allnumsegs=xNew<int>(numcontours);

	/*Loop through all contours: */
	for (i=0;i<numcontours;i++){

		/*retrieve contour info: */
		contouri=*(contours+i);
		numnodes=contouri->nods;
		xc=contouri->x;
		yc=contouri->y;

		/*determine segmentsi and numsegsi for this contour and the mesh intersection: */
		MeshSegmentsIntersection(&segmentsi,&numsegsi,index,x,y,nel,nods,xc,yc,numnodes);

		/*save segmentsi: */
		allsegments[i]=segmentsi;
		allnumsegs[i]=numsegsi;
	}

	/*total number of segments? */
	numsegs=0;
	for(i=0;i<numcontours;i++)numsegs+=allnumsegs[i];

	/*Out of all segments, create one common array of segments: */
	segments=xNew<double>(5*numsegs);
	count=0;
	for(i=0;i<numcontours;i++){

		segmentsi=allsegments[i];
		numsegsi=allnumsegs[i];

		for(j=0;j<numsegsi;j++){
			for(k=0;k<5;k++){
				*(segments+count*5+k)=*(segmentsi+j*5+k);
			}
			count++;
		}
	}

	/*Assign output pointers:*/
	*psegments=segments;
	*pnumsegs=numsegs;
}/*}}}*/
void MeshSegmentsIntersection(double** psegments, int* pnumsegs,int* index, double* x, double* y, int nel, int nods, double* xc, double* yc, int numnodes){/*{{{*/

	int      i,j;

	/*output: */
	double* segments=NULL;
	int     numsegs;

	/*intermediary: */
	DataSet* segments_dataset=NULL;
	double   xnodes[3];
	double   ynodes[3];

	/*We don't know how many segments  we are going to get, so have a dynamic container: */
	segments_dataset=new DataSet();

	/*Go through elements, and call ElementSegmentsIntersection routine: */
	for(i=0;i<nel;i++){
		for(j=0;j<3;j++){
			xnodes[j]=x[*(index+3*i+j)];
			ynodes[j]=y[*(index+3*i+j)];
		}
		ElementSegmentsIntersection(segments_dataset,i,xnodes,ynodes,xc,yc,numnodes);
	}

	/*Using the segments_dataset dataset, create segments: */
	numsegs=segments_dataset->Size();
	segments=xNew<double>(5*numsegs);
	for(i=0;i<numsegs;i++){
		Segment<double>* segment=(Segment<double>*)segments_dataset->GetObjectByOffset(i);

		/*x1,y1,x2,y2 then element_id: */
		segments[5*i+0]=segment->x1;
		segments[5*i+1]=segment->y1;
		segments[5*i+2]=segment->x2;
		segments[5*i+3]=segment->y2;
		segments[5*i+4]=(double)segment->eid;
	}

	/*Free resources:*/
	delete segments_dataset;

	/*Assign output pointers:*/
	*psegments=segments;
	*pnumsegs=numsegs;
}/*}}}*/

/*Utilities*/
void ElementSegmentsIntersection(DataSet* segments_dataset,int el, double* xnodes,double* ynodes,double* xc,double* yc,int numnodes){/*{{{*/

	double xsegment[2];
	double ysegment[2];

	/*Loop through contour: */
	for(int i=0;i<numnodes-1;i++){
		xsegment[0]=xc[i];
		xsegment[1]=xc[i+1];
		ysegment[0]=yc[i];
		ysegment[1]=yc[i+1];
		/*if (el==318 && i==9){
			_printf_("contour: " << i << " " << xsegment[0] << " " << ysegment[0] << " " << xsegment[1] << " " << ysegment[1] 
				<< " " << xnodes[0] << " " << xnodes[1] << " " << xnodes[2] << " " << ynodes[0] << " " << ynodes[1] << " " << 
				ynodes[2] << "\n");
		}*/
		ElementSegment(segments_dataset,el, i, xnodes,ynodes,xsegment,ysegment);
	}
}/*}}}*/
void ElementSegment(DataSet* segments_dataset,int el, int contouri, double* xnodes,double* ynodes,double* xsegment,double* ysegment){/*{{{*/

	/*We have a tria element (xnodes,ynodes) and a segment (xsegment,ysegment). Find whether they intersect. 
	 * If they do, create a Segment object with the intersection, and add to segments_dataset dataset: */

	double alpha1,alpha2;
	double beta1,beta2;
	double gamma1,gamma2;

	int    edge1,edge2,edge3;

	double xel[2],yel[2];
	double coord1 = 0.;
	double coord2 = 0.;
	double xfinal[2],yfinal[2];

	/*edge 1: */
	xel[0]=xnodes[0];  yel[0]=ynodes[0]; xel[1]=xnodes[1];  yel[1]=ynodes[1];
	edge1=SegmentIntersect(&alpha1,&alpha2, xel,yel,xsegment,ysegment); //alpha1: segment coordinate of intersection. alpha2: same thing for second interesection if it exists (colinear edges)

	/*edge 2: */
	xel[0]=xnodes[1];  yel[0]=ynodes[1]; xel[1]=xnodes[2];  yel[1]=ynodes[2];
	edge2=SegmentIntersect(&beta1,&beta2, xel,yel,xsegment,ysegment);

	/*edge 3: */
	xel[0]=xnodes[2];  yel[0]=ynodes[2]; xel[1]=xnodes[0];  yel[1]=ynodes[0];
	edge3=SegmentIntersect(&gamma1,&gamma2, xel,yel,xsegment,ysegment);

	/*edge can be either IntersectEnum (one and only one intersection between the edge and the segment), ColinearEnum (edge and segment are collinear) and SeparateEnum (no intersection): */

	/*if (el==318 && contouri==9){
		_printf_(edge1 << " " << edge2 << " " << edge3 << " "  << alpha1 << " " << alpha2 << " " << beta1 << " " << beta2 << " " << gamma1 << " " << gamma2 << " " << xsegment[0] << " "  << xsegment[1] << " " << ysegment[0] << " " << ysegment[1] << " " << xnodes[0] << " " << xnodes[1] << " " << xnodes[2] << " " << ynodes[0] << " " << ynodes[1] << " " << ynodes[2]);

	_printf_("Bool" << (edge1==IntersectEnum) || (edge2==IntersectEnum) || (edge3==IntersectEnum));
	}*/

	if(    (edge1==IntersectEnum) && (edge2==IntersectEnum) && (edge3==IntersectEnum)   ){

		/*This can only be the case if the segment intersected through one vertex, meaning a pair from alpha1, beta1 or gamma1  is 0:*/
		if (alpha1!=0 && alpha1!=1){
			/*The vertex opposite edge 1 was intersected:*/
			xfinal[0]=xsegment[0]+alpha1*(xsegment[1]-xsegment[0]);
			yfinal[0]=ysegment[0]+alpha1*(ysegment[1]-ysegment[0]);
			xfinal[1]=xnodes[2];
			yfinal[1]=ynodes[2];
		}
		else if (beta1!=0 && beta1!=1){
			/*The vertex opposite edge 2 was intersected:*/
			xfinal[0]=xsegment[0]+beta1*(xsegment[1]-xsegment[0]);
			yfinal[0]=ysegment[0]+beta1*(ysegment[1]-ysegment[0]);
			xfinal[1]=xnodes[0];
			yfinal[1]=ynodes[0];
		}
		else if (gamma1!=0 && gamma1!=1){
			/*The vertex opposite edge 3 was intersected:*/
			xfinal[0]=xsegment[0]+gamma1*(xsegment[1]-xsegment[0]);
			yfinal[0]=ysegment[0]+gamma1*(ysegment[1]-ysegment[0]);
			xfinal[1]=xnodes[1];
			yfinal[1]=ynodes[1];
		}
		segments_dataset->AddObject(new  Segment<double>(el+1,xfinal[0],yfinal[0],xfinal[1],yfinal[1]));

		/*This case is impossible: not quite! */
		//_printf_(alpha1 << " " << alpha2 << " " << beta1 << " " << beta2 << " " << gamma1 << " " << gamma2 << " " << xsegment[0] << " "  << xsegment[1] << " " << ysegment[0] << " " << ysegment[1] << " " << xnodes[0] << " " << xnodes[1] << " " << xnodes[2] << " " << ynodes[0] << " " << ynodes[1] << " " << ynodes[2]);
		/* _error_("error: a line cannot go through 3 different vertices!");*/
	}
	else if(    ((edge1==IntersectEnum) && (edge2==IntersectEnum)) || ((edge2==IntersectEnum) && (edge3==IntersectEnum)) || ((edge3==IntersectEnum) && (edge1==IntersectEnum))   ){

		/*segment interscts 2 opposite edges of our triangle, at 2 segment coordinates, pick up the lowest (coord1) and highest (coord2): */
		if((edge1==IntersectEnum) && (edge2==IntersectEnum)) {coord1=min(alpha1,beta1); coord2=max(alpha1,beta1);}
		if((edge2==IntersectEnum) && (edge3==IntersectEnum)) {coord1=min(beta1,gamma1); coord2=max(beta1,gamma1);}
		if((edge3==IntersectEnum) && (edge1==IntersectEnum)) {coord1=min(gamma1,alpha1); coord2=max(gamma1,alpha1);}

		/*check this segment did not intersect at a vertex of the tria: */
		if(coord1!=coord2){

			xfinal[0]=xsegment[0]+coord1*(xsegment[1]-xsegment[0]);
			xfinal[1]=xsegment[0]+coord2*(xsegment[1]-xsegment[0]);
			yfinal[0]=ysegment[0]+coord1*(ysegment[1]-ysegment[0]);
			yfinal[1]=ysegment[0]+coord2*(ysegment[1]-ysegment[0]);

			segments_dataset->AddObject(new  Segment<double>(el+1,xfinal[0],yfinal[0],xfinal[1],yfinal[1]));
		}
		else{
			/*the segment intersected at the vertex, do not bother with this "0" length segment!:*/
		}
	}
	else if(  (edge1==IntersectEnum) || (edge2==IntersectEnum) || (edge3==IntersectEnum)   ){

		/*if (el==318 && contouri==9){
			_printf_("hello" <<  " NodeInElement 0 " << (NodeInElement(xnodes,ynodes,xsegment[0],ysegment[0])) <<  " NodeInElement 1 " << (NodeInElement(xnodes,ynodes,xsegment[1],ysegment[1])));
		}*/

		/*segment intersect only 1 edge. Figure out where the first point in the segment is, inside or outside the element, 
		 * this will decide the coordinate: */
		if (NodeInElement(xnodes,ynodes,xsegment[0],ysegment[0])){
			coord1=0;
			if(edge1==IntersectEnum){coord2=alpha1;}
			if(edge2==IntersectEnum){coord2=beta1;}
			if(edge3==IntersectEnum){coord2=gamma1;}
		}
		else if (NodeInElement(xnodes,ynodes,xsegment[1],ysegment[1])){
			if(edge1==IntersectEnum){coord1=alpha1;}
			if(edge2==IntersectEnum){coord1=beta1;}
			if(edge3==IntersectEnum){coord1=gamma1;}
			coord2=1.0;
		}
		else{
			double tolerance=1e-10;
			/*Ok, we have an issue here. Probably one of the segments' end is on a vertex, within a certain tolerance!*/
			if (IsIdenticalNode(xnodes[0],ynodes[0],xsegment[0],ysegment[0],tolerance) ||
				IsIdenticalNode(xnodes[1],ynodes[1],xsegment[0],ysegment[0],tolerance) ||
				IsIdenticalNode(xnodes[2],ynodes[2],xsegment[0],ysegment[0],tolerance)){

				/*ok, segments[0] is common to one of our vertices: */
				coord1=0;
				if(edge1==IntersectEnum){coord2=alpha1;}
				if(edge2==IntersectEnum){coord2=beta1;}
				if(edge3==IntersectEnum){coord2=gamma1;}
			}
			else if (IsIdenticalNode(xnodes[0],ynodes[0],xsegment[1],ysegment[1],tolerance) ||
				     IsIdenticalNode(xnodes[1],ynodes[1],xsegment[1],ysegment[1],tolerance) ||
				     IsIdenticalNode(xnodes[2],ynodes[2],xsegment[1],ysegment[1],tolerance)){

				/*ok, segments[1] is common to one of our vertices: */
				//if (el==318 && contouri==9){ _printf_("ok2" << "\n"); }
				if(edge1==IntersectEnum){coord1=alpha1;}
				if(edge2==IntersectEnum){coord1=beta1;}
				if(edge3==IntersectEnum){coord1=gamma1;}
				coord2=1.0;
			}
		}

		xfinal[0]=xsegment[0]+coord1*(xsegment[1]-xsegment[0]);
		xfinal[1]=xsegment[0]+coord2*(xsegment[1]-xsegment[0]);
		yfinal[0]=ysegment[0]+coord1*(ysegment[1]-ysegment[0]);
		yfinal[1]=ysegment[0]+coord2*(ysegment[1]-ysegment[0]);

		segments_dataset->AddObject(new  Segment<double>(el+1,xfinal[0],yfinal[0],xfinal[1],yfinal[1]));
	}
	else{
		/*No interesections, but the segment might be entirely inside this triangle!: */
		if ( (NodeInElement(xnodes,ynodes,xsegment[0],ysegment[0])) && (NodeInElement(xnodes,ynodes,xsegment[1],ysegment[1])) ){
			segments_dataset->AddObject(new  Segment<double>(el+1,xsegment[0],ysegment[0],xsegment[1],ysegment[1]));
		}
	}
}/*}}}*/
bool NodeInElement(double* xnodes, double* ynodes, double x, double y){/*{{{*/

	double x1,y1;
	double x2,y2;
	double x3,y3;
	double lambda1,lambda2,lambda3;
	double det;

	x1=xnodes[0];
	x2=xnodes[1];
	x3=xnodes[2];
	y1=ynodes[0];
	y2=ynodes[1];
	y3=ynodes[2];

	/*compute determinant: */
	det=x1*y2-x1*y3-x3*y2-x2*y1+x2*y3+x3*y1;

	/*area coordinates: */
	lambda1=((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/det;
	lambda2=((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/det;
	lambda3=1-lambda1-lambda2;

	if( ((lambda1<=1) && (lambda1>=0)) && ((lambda2<=1) && (lambda2>=0)) && ((lambda3<=1) && (lambda3>=0))  )return true;
	else return false;

}/*}}}*/
int SegmentIntersect(double* palpha, double* pbeta, double* x1, double* y1, double* x2, double* y2){/*{{{*/

	/*See ISSM_DIR/src/m/utils/Geometry/SegIntersect.m for matlab routine from which we take this routine: */

	/*output: */
	double alpha=-1;
	double beta=-1;

	double xA,xB,xC,xD,yA,yB,yC,yD;
	double O2A[2],O2B[2],O1C[2],O1D[2];
	double n1[2],n2[2];
	double test1, test2, test3, test4;
	double det;
	double O2O1[2];
	double pO1A,pO1B,pO1C,pO1D;

	xA=x1[0]; yA=y1[0];
	xB=x1[1]; yB=y1[1];
	xC=x2[0]; yC=y2[0];
	xD=x2[1]; yD=y2[1];

	O2A[0]=xA -(xD/2+xC/2); O2A[1]=yA -(yD/2+yC/2);
	O2B[0]=xB -(xD/2+xC/2); O2B[1]=yB -(yD/2+yC/2);
	O1C[0]=xC -(xA/2+xB/2); O1C[1]=yC -(yA/2+yB/2);
	O1D[0]=xD -(xA/2+xB/2); O1D[1]=yD -(yA/2+yB/2);

	n1[0]=yA-yB; n1[1]=xB-xA;  //normal vector to segA
	n2[0]=yC-yD; n2[1]=xD-xC;  //normal vector to segB

	test1=n2[0]*O2A[0]+n2[1]*O2A[1];
	test2=n2[0]*O2B[0]+n2[1]*O2B[1];

	if (test1*test2>0){
		return SeparateEnum;
	}

	test3=n1[0]*O1C[0]+n1[1]*O1C[1];
	test4=n1[0]*O1D[0]+n1[1]*O1D[1];

	if (test3*test4>0){
		return SeparateEnum;
	}

	/*If colinear: */
	det=n1[0]*n2[1]-n2[0]*n1[1];

	if(test1*test2==0 && test3*test4==0 && det==0){

		//projection on the axis O1O2
		O2O1[0]=(xA/2+xB/2)-(xD/2+xC/2);
		O2O1[1]=(yA/2+yB/2)-(yD/2+yC/2);

		pO1A=O2O1[0]*(O2A[0]-O2O1[0])+O2O1[1]*(O2A[1]-O2O1[1]);
		pO1B=O2O1[0]*(O2B[0]-O2O1[0])+O2O1[1]*(O2B[1]-O2O1[1]);
		pO1C=O2O1[0]*O1C[0]+O2O1[1]*O1C[1];
		pO1D=O2O1[0]*O1D[0]+O2O1[1]*O1D[1];

		//test if one point is included in the other segment (->intersects=true)
		if ((pO1C-pO1A)*(pO1D-pO1A)<0){
			alpha=0; beta=0;
			*palpha=alpha;*pbeta=beta;
			return ColinearEnum;
		}
		if ((pO1C-pO1B)*(pO1D-pO1B)<0){
			alpha=0; beta=0;
			*palpha=alpha;*pbeta=beta;
			return ColinearEnum;
		}
		if ((pO1A-pO1C)*(pO1B-pO1C)<0){
			alpha=0; beta=0;
			*palpha=alpha;*pbeta=beta;
			return ColinearEnum;
		}
		if ((pO1A-pO1D)*(pO1B-pO1D)<0){
			alpha=0; beta=0;
			*palpha=alpha;*pbeta=beta;
			return ColinearEnum;
		}

		//test if the 2 segments have the same middle (->intersects=true)
		if(O2O1[0]==0 && O2O1[1]){
			alpha=0; beta=0;
			*palpha=alpha;*pbeta=beta;
			return ColinearEnum;
		}

		//if we are here, both segments are colinear, but do not interset:
		alpha=-1; beta=-1;
		*palpha=alpha;*pbeta=beta;
		return SeparateEnum;
	}

	/*if we are here, both segments intersect. Determine where in the segment coordinate 
	 * system: */
	beta=-1;
	alpha=-(xA*yB-xC*yB+yC*xB-yC*xA+xC*yA-yA*xB)/(-xD*yB+xD*yA+xC*yB-xC*yA-yD*xA+yD*xB+yC*xA-yC*xB); //from intersect.m in formal calculus

	*palpha=alpha;*pbeta=beta;
	return IntersectEnum;
}/*}}}*/
bool IsIdenticalNode(double x1, double y1, double x2, double y2, double tolerance){ /*{{{*/

	if (sqrt(pow(x1-x2,2.0) + pow(y1-y2,2))<tolerance)return true;
	else return false;

}/*}}}*/
