/*\file ExpSimplify.c
 *\brief: exp to kml file conversion mex module.
 */
#include "./ExpSimplify.h"

void ExpSimplifyUsage(void){/*{{{*/
	_printf_("ExpSimplify - Simplify Exp contour\n");
	_printf_("\n");
	_printf_("   Recursive Douglas-Peucker Polygon Simplification\n");
	_printf_("\n");
	_printf_("   Usage:\n");
	_printf_("      ExpSimplify(expfile,tol);\n");
	_printf_("      - expfile: name of the exp file\n");
	_printf_("      - tol:  tolerance (maximal euclidean distance allowed between the new line and a vertex)\n");
	_printf_("      Additional options:\n");
	_printf_("      - 'min': minimum number of vertices to save contours in exp file (default is 3)\n");
	_printf_("\n");
	_printf_("   Example:\n");
	_printf_("      ExpSimplify('file.exp',100);\n");
	_printf_("      ExpSimplify('file.exp',100,'min',4);\n");
}/*}}}*/
void simplify(Contour<double>* contour,bool* flags,int ind0,int ind1,double tolerance){/*{{{*/

	bool    closed    = false;
	double  distance,beta,dx,dy;
	double  maxdistance;
	int     index  = -1;
	double *x      = contour->x;
	double *y      = contour->y;

	/*Some checks*/
	_assert_(ind0>=0 && ind0<contour->nods);
	_assert_(ind1>=0 && ind1<contour->nods);
	_assert_(ind1-ind0>=0);

	/*Check wether this portion is closed*/
	if(x[ind0]==x[ind1] && y[ind0]==y[ind1]) closed=true;

	if(closed){

		/*calculate the shortest distance of all vertices between ind0 and ind1*/
		for(int i=ind0;i<ind1-1;i++){
			distance = sqrt((x[ind0]-x[i+1])*(x[ind0]-x[i+1]) + (y[ind0]-y[i+1])*(y[ind0]-y[i+1]));
			if(i==ind0 || distance>maxdistance){
				maxdistance=distance;
				index = i + 1;
			}
		}
	}
	else{
		/*calculate shortest distance of all points to the line from ind0 to ind1
		 * subtract starting point from other locations
		 *
		 * d = || (x-x0) - beta (xend - x0) ||
		 * <x-x0,xend-x0>      = ||x-x0|| ||xend-x0|| cos(alpha)
		 * beta ||xend-x0|| = ||x-x0|| cos(alpha)
		 *
		 * So: beta = <x-x0,xend-x0>/<xend-x0,xend-x0>  */

		for(int i=ind0+1;i<=ind1;i++){
			beta = ((x[i]-x[ind0])*(x[ind1]-x[ind0]) + (y[i]-y[ind0])*(y[ind1]-y[ind0]))/((x[ind1]-x[ind0])*(x[ind1]-x[ind0])+(y[ind1]-y[ind0])*(y[ind1]-y[ind0]));
			dx   = x[i]-beta*x[ind1]+(beta-1.)*x[ind0];
			dy   = y[i]-beta*y[ind1]+(beta-1.)*y[ind0];
			distance = sqrt(dx*dx + dy*dy);
			if(i==ind0+1 || distance>maxdistance){
				maxdistance = distance;
				index       = i;
			}
		}

	}

	/*if the maximum distance is smaller than the tolerance remove vertices between ind0 and ind1*/
	if(maxdistance<tolerance){
		if(ind0!=ind1-1){
			for(int i=ind0+1;i<ind1;i++) flags[i]=false;
		}
	}
	else{
		/*if not, call simplifyrec for the segments between ind0 and index
		 * (index and ind1)*/
		_assert_(index!=-1);
		_assert_(index!=ind1);
		_assert_(index!=ind0);
		simplify(contour,flags,ind0 ,index,tolerance);
		simplify(contour,flags,index,ind1, tolerance);
	}

}/*}}}*/
WRAPPER(ExpSimplify_python){

	int i,verbose=1;

	/*input: */
	char*    expfile  = NULL;
	double   tolerance;
	Options *options      = NULL;
	double   minimumvertices_double;
	int      minimumvertices;

	/*output*/
	Contours* oldcontours = NULL;
	Contours* newcontours = NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments: */
	/*checks on arguments on the matlab side: */
	if (nrhs<NRHS || nlhs>NLHS){
		ExpSimplifyUsage(); _error_("ExpSimplify usage error");
	}

	/*Input datasets: */
	FetchData(&expfile,  EXPFILE);
	FetchData(&tolerance,TOLERANCE);
	FetchData(&options,NRHS,nrhs,prhs);

	/*some checks*/
	if(tolerance<0) _error_("tolerance must be a positivve scalar");

	/* Run core computations: */
	Contour<double>* contour = NULL;
	Contour<double>* newcontour = NULL;
	int     nods,newnods;
	bool*   flags = NULL;
	double* x = NULL;
	double* y = NULL;
	double distance;

	/*Process options*/
	options->Get(&minimumvertices_double,"min",3.);
	if(minimumvertices_double<1.) _error_("'min' (minimum number of verties) should be a positive integer");
	minimumvertices = int(minimumvertices_double);

	/*Read old contours and allocate new contours*/
	oldcontours=ExpRead<double>(expfile);
	newcontours=new Contours();
	for(int counter=0;counter<oldcontours->Size();counter++){

		/*Get single contour*/
		contour = (Contour<double>*)oldcontours->GetObjectByOffset(counter);
		nods    = contour->nods;
		x       = contour->x;
		y       = contour->y;
		_printf_("   Initial number of vertices in contour #"<<counter+1<<": "<<nods << "\n");

		/*Allocate flags (1=keep, 0=remove)*/
		if(nods>0) flags   = xNew<bool>(nods);

		if(nods==0){
			/*Don't do anything*/
		}
		else if(nods==1){
			flags[0] = true;
		}
		else if(nods==2){
			/*check if the distance between both is less than the tolerance
			 * If so, return the center*/
			distance = sqrt((x[1]-x[0])*(x[1]-x[0]) + (y[1]-y[0])*(y[1]-y[0]));
			if(distance<=tolerance){
				x[0]=(x[0]+x[1])/2.;
				y[0]=(y[0]+y[1])/2.;
				flags[0] = true;
				flags[1] = false;
			}
			else{
				flags[0] = true;
				flags[1] = true;
			}
		}
		else{
			/*Start recursive call to simplify*/
			for(int i=0;i<nods;i++) flags[i]=true;
			simplify(contour,flags,0,nods-1,tolerance);
		}

		/*Add new contour to newcontours*/
		newnods = 0;
		for(int i=0;i<nods;i++) if(flags[i]) newnods++;

		/*Do we save new profile?*/
		if(newnods>=minimumvertices){
			_printf_("   Final   number of vertices in contour #"<<counter+1<<": "<<newnods << "\n");
			newcontour       = new Contour<double>();
			newcontour->nods = newnods;
			newcontour->x    = xNew<double>(newnods);
			newcontour->y    = xNew<double>(newnods);
			newnods = 0;
			for(int i=0;i<nods;i++){
				if(flags[i]){
					newcontour->x[newnods] = contour->x[i];
					newcontour->y[newnods] = contour->y[i];
					newnods++;
				}
			}
			_assert_(newnods==newcontour->nods);

			/*Add to main dataset*/
			newcontours->AddObject(newcontour);
		}
		else{
			_printf_("   Final   number of vertices in contour #"<<counter+1<<": "<<newnods<<" (not saved)\n");
		}

		/*cleanup*/
		xDelete<bool>(flags);
	}
	_printf_("   Initial number of contours: "<<oldcontours->Size() << "\n");
	_printf_("   Final   number of contours: "<<newcontours->Size() << "\n");

	/*Write data: */
	ExpWrite(newcontours,expfile);

	/*Clean-up*/
	xDelete<char>(expfile);
	delete options;
	delete oldcontours;
	delete newcontours;

	/*end module: */
	MODULEEND();
}
