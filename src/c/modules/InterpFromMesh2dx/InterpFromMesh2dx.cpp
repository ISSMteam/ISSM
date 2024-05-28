/*!\file:  InterpFromMesh2dx.cpp
 * \brief  "c" core code for interpolating values from a structured grid.
 */ 

#include "./InterpFromMesh2dx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../ContourToNodesx/ContourToNodesx.h"

int InterpFromMesh2dx(IssmSeqVec<IssmPDouble>** pdata_prime,
			int* index_data, double* x_data, double* y_data, int nods_data,int nels_data, double* data, int data_length,
			double* x_prime, double* y_prime, int nods_prime,
			double* default_values,int num_default_values,Contour<IssmPDouble>** contours,int numcontours){

	/*Output*/
	IssmSeqVec<IssmPDouble>* data_prime=NULL;

	/*Intermediary*/
	int    i;
	int    interpolation_type;
	bool   debug;
	double xmin,xmax;
	double ymin,ymax;

	/*contours: */
	double *incontour     = NULL;

	/*some checks*/
	if (nels_data<1 || nods_data<3 || nods_prime==0){
		_error_("nothing to be done according to the mesh given in input");
	}

	/*Set debug to 1 if there are lots of elements*/
	debug=(bool)((double)nels_data*(double)nods_prime >= pow((double)10,(double)9));

	/*figure out what kind of interpolation is needed*/
	if (data_length==nods_data){
		interpolation_type=1;
	}
	else if (data_length==nels_data){
		interpolation_type=2;
	}
	else{
		_error_("length of vector data not supported yet. It should be of length (number of nodes) or (number of elements)!");
	}

	if((numcontours) && (interpolation_type==2)){
		_error_("element interpolation_type with contours not supported yet!");
	}

	/*Get prime mesh extrema coordinates*/
	xmin=x_prime[0]; xmax=x_prime[0];ymin=y_prime[0]; ymax=y_prime[0];
	for (i=1;i<nods_prime;i++){
		if (x_prime[i]<xmin) xmin=x_prime[i];
		if (x_prime[i]>xmax) xmax=x_prime[i];
		if (y_prime[i]<ymin) ymin=y_prime[i];
		if (y_prime[i]>ymax) ymax=y_prime[i];
	}

	/*Initialize output*/
	data_prime=new IssmSeqVec<IssmPDouble>(nods_prime);
	if(num_default_values){
		if(num_default_values==1)for (i=0;i<nods_prime;i++) data_prime->SetValue(i,default_values[0],INS_VAL);
		else for (i=0;i<nods_prime;i++) data_prime->SetValue(i,default_values[i],INS_VAL);
	}

	/*Build indices of contour: */
	if(numcontours){
		ContourToNodesx( &incontour,x_prime,y_prime,nods_prime,contours,numcontours,1);
	}
	else{
		 incontour=xNew<double>(nods_prime);
		 for (i=0;i<nods_prime;i++) incontour[i]=1.0;
	}

	/*initialize thread parameters: */
	InterpFromMesh2dxThreadStruct gate;
	gate.interpolation_type = interpolation_type;
	gate.debug              = debug;
	gate.nels_data          = nels_data;
	gate.index_data         = index_data;
	gate.x_data             = x_data;
	gate.y_data             = y_data;
	gate.data               = data;
	gate.xmin               = xmin;
	gate.xmax               = xmax;
	gate.ymin               = ymin;
	gate.ymax               = ymax;
	gate.nods_prime         = nods_prime;
	gate.data_prime         = data_prime;
	gate.x_prime            = x_prime;
	gate.y_prime            = y_prime;
	gate.default_values     = default_values;
	gate.num_default_values = num_default_values;
	gate.incontour          = incontour;

	/*launch the thread manager with InterpFromGridToMeshxt as a core: */
	LaunchThread(InterpFromMesh2dxt,(void*)&gate,_NUMTHREADS_);

	/*Assign output pointers:*/
	 xDelete<double>(incontour);
	*pdata_prime=data_prime;
	return 1;
}
