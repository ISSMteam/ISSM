/*!\file GaussPenta.c
 * \brief: implementation of the GaussPenta object
 */

#include "./GaussPenta.h"
#include "./GaussTria.h"
#include "../../shared/io/Print/Print.h"
#include "../../shared/Exceptions/exceptions.h"
#include "../../shared/MemOps/MemOps.h"
#include "../../shared/Numerics/recast.h"
#include "../../shared/Enum/Enum.h"
#include "../../shared/Numerics/GaussPoints.h"
#include "../../shared/Numerics/constants.h"

/*GaussPenta constructors and destructors:*/
GaussPenta::GaussPenta(){/*{{{*/

	ig = -1;
	numgauss=-1;

	weights=NULL;
	coords1=NULL;
	coords2=NULL;
	coords3=NULL;
	coords4=NULL;

	weight=UNDEF;
	coord1=UNDEF;
	coord2=UNDEF;
	coord3=UNDEF;
	coord4=UNDEF;
}
/*}}}*/
GaussPenta::GaussPenta(int order_horiz,int order_vert){/*{{{*/

	/*Intermediaries*/
	int         numgauss_horiz;
	int         numgauss_vert;
	IssmDouble *coords1_horiz  = NULL;
	IssmDouble *coords2_horiz  = NULL;
	IssmDouble *coords3_horiz  = NULL;
	IssmDouble *weights_horiz  = NULL;
	double     *coords_vert    = NULL;
	double     *weights_vert   = NULL;

	/*Get gauss points*/
	GaussLegendreTria(&numgauss_horiz,&coords1_horiz,&coords2_horiz,&coords3_horiz,&weights_horiz,order_horiz);
	GaussLegendreLinear(&coords_vert,&weights_vert,order_vert);
	numgauss_vert=order_vert;

	/*Allocate GaussPenta fields*/
   ig      = -1;
	numgauss=numgauss_horiz*numgauss_vert;
	coords1=xNew<IssmDouble>(numgauss);
	coords2=xNew<IssmDouble>(numgauss);
	coords3=xNew<IssmDouble>(numgauss);
	coords4=xNew<IssmDouble>(numgauss);
	weights=xNew<IssmDouble>(numgauss);

	/*Combine Horizontal and vertical points*/
	for(int ighoriz=0; ighoriz<numgauss_horiz; ighoriz++){
		for(int igvert=0; igvert<numgauss_vert; igvert++){
			coords1[numgauss_vert*ighoriz+igvert]=coords1_horiz[ighoriz];
			coords2[numgauss_vert*ighoriz+igvert]=coords2_horiz[ighoriz];
			coords3[numgauss_vert*ighoriz+igvert]=coords3_horiz[ighoriz];
			coords4[numgauss_vert*ighoriz+igvert]=coords_vert[igvert];
			weights[numgauss_vert*ighoriz+igvert]=weights_horiz[ighoriz]*weights_vert[igvert];
		}
	}

	/*Initialize static fields as undefinite*/
	weight=UNDEF;
	coord1=UNDEF;
	coord2=UNDEF;
	coord3=UNDEF;
	coord4=UNDEF;

	/*Clean up*/
	xDelete<IssmDouble>(coords1_horiz);
	xDelete<IssmDouble>(coords2_horiz);
	xDelete<IssmDouble>(coords3_horiz);
	xDelete<double>(coords_vert);
	xDelete<IssmDouble>(weights_horiz);
	xDelete<double>(weights_vert);
}
/*}}}*/
GaussPenta::GaussPenta(int index1, int index2,int order){/*{{{*/

	/*Intermediaties*/
	double *seg_coords  = NULL;
	double *seg_weights = NULL;
	int     i;

	/*Get Segment gauss points*/
   ig      = -1;
	numgauss=order;
	GaussLegendreLinear(&seg_coords,&seg_weights,numgauss);

	/*Allocate GaussPenta fields*/
	coords1=xNew<IssmDouble>(numgauss);
	coords2=xNew<IssmDouble>(numgauss);
	coords3=xNew<IssmDouble>(numgauss);
	coords4=xNew<IssmDouble>(numgauss);
	weights=xNew<IssmDouble>(numgauss);

	if(index1==0 && index2==3){
		for(i=0;i<numgauss;i++) coords1[i]=1.0;
		for(i=0;i<numgauss;i++) coords2[i]=0.0;
		for(i=0;i<numgauss;i++) coords3[i]=0.0;
		for(i=0;i<numgauss;i++) coords4[i]=seg_coords[i];
		for(i=0;i<numgauss;i++) weights[i]=seg_weights[i];
	}
	else if (index1==1 && index2==4){
		for(i=0;i<numgauss;i++) coords1[i]=0.0;
		for(i=0;i<numgauss;i++) coords2[i]=1.0;
		for(i=0;i<numgauss;i++) coords3[i]=0.0;
		for(i=0;i<numgauss;i++) coords4[i]=seg_coords[i];
		for(i=0;i<numgauss;i++) weights[i]=seg_weights[i];
	}
	else if (index1==2 && index2==5){
		for(i=0;i<numgauss;i++) coords1[i]=0.0;
		for(i=0;i<numgauss;i++) coords2[i]=0.0;
		for(i=0;i<numgauss;i++) coords3[i]=1.0;
		for(i=0;i<numgauss;i++) coords4[i]=seg_coords[i];
		for(i=0;i<numgauss;i++) weights[i]=seg_weights[i];
	}
	else{
		_error_("Penta not supported yet");
	}

	/*Initialize static fields as undefined*/
	weight=UNDEF;
	coord1=UNDEF;
	coord2=UNDEF;
	coord3=UNDEF;
	coord4=UNDEF;

	/*clean up*/
	xDelete<double>(seg_coords);
	xDelete<double>(seg_weights);

}
/*}}}*/
GaussPenta::GaussPenta(int index1, int index2, int index3, int order){/*{{{*/

	/*Basal Tria*/
	if(index1==0 && index2==1 && index3==2){

		/*Get GaussTria*/
		GaussLegendreTria(&numgauss,&coords1,&coords2,&coords3,&weights,order);

		/*compute z coordinate*/
		coords4=xNew<IssmDouble>(numgauss);
		for(int i=0;i<numgauss;i++) coords4[i]=-1.0;
	}
	/*Upper surface Tria*/
	else if(index1==3 && index2==4 && index3==5){

		/*Get GaussTria*/
		GaussLegendreTria(&numgauss,&coords1,&coords2,&coords3,&weights,order);

		/*compute z coordinate*/
		coords4=xNew<IssmDouble>(numgauss);
		for(int i=0;i<numgauss;i++) coords4[i]=1.0;
	}
	else{
		_error_("Tria not supported yet");
	}

   this->ig = -1;

}
/*}}}*/
GaussPenta::GaussPenta(int index1, int index2, int index3, int index4,int order_horiz,int order_vert){/*{{{*/

	/*Intermediaties*/
	double *seg_horiz_coords  = NULL;
	double *seg_horiz_weights = NULL;
	double *seg_vert_coords   = NULL;
	double *seg_vert_weights  = NULL;
	int     i,j;

	/*get the gauss points using the product of two line rules*/
	GaussLegendreLinear(&seg_horiz_coords,&seg_horiz_weights,order_horiz);
	GaussLegendreLinear(&seg_vert_coords, &seg_vert_weights, order_vert);

	/*Allocate GaussPenta fields*/
   this->ig = -1;
	numgauss=order_horiz*order_vert;
	coords1=xNew<IssmDouble>(numgauss);
	coords2=xNew<IssmDouble>(numgauss);
	coords3=xNew<IssmDouble>(numgauss);
	coords4=xNew<IssmDouble>(numgauss);
	weights=xNew<IssmDouble>(numgauss);

	/*Quads: get the gauss points using the product of two line rules  */
	if(index1==0 && index2==1 && index3==4 && index4==3){
		for(i=0;i<order_horiz;i++){
			for(j=0;j<order_vert;j++){
				coords1[i*order_vert+j]=  0.5*(1-seg_horiz_coords[i]);
				coords2[i*order_vert+j]=1-0.5*(1-seg_horiz_coords[i]);
				coords3[i*order_vert+j]=0.0;
				coords4[i*order_vert+j]=seg_vert_coords[j];
				weights[i*order_vert+j]=seg_horiz_weights[i]*seg_vert_weights[j];
			}
		}
	}
	else if(index1==1 && index2==2 && index3==5 && index4==4){
		for(i=0;i<order_horiz;i++){
			for(j=0;j<order_vert;j++){
				coords1[i*order_vert+j]=0.0;
				coords2[i*order_vert+j]=  0.5*(1-seg_horiz_coords[i]);
				coords3[i*order_vert+j]=1-0.5*(1-seg_horiz_coords[i]);
				coords4[i*order_vert+j]=seg_vert_coords[j];
				weights[i*order_vert+j]=seg_horiz_weights[i]*seg_vert_weights[j];
			}
		}
	}
	else if(index1==2 && index2==0 && index3==3 && index4==5){
		for(i=0;i<order_horiz;i++){
			for(j=0;j<order_vert;j++){
				coords1[i*order_vert+j]=1-0.5*(1-seg_horiz_coords[i]);
				coords2[i*order_vert+j]=0.0;
				coords3[i*order_vert+j]=  0.5*(1-seg_horiz_coords[i]);
				coords4[i*order_vert+j]=seg_vert_coords[j];
				weights[i*order_vert+j]=seg_horiz_weights[i]*seg_vert_weights[j];
			}
		}
	}
	else{
		_error_("Tria not supported yet (user provided indices " << index1 << " " << index2 << " " << index3 << " " << index4 << ")");
	}

	/*clean-up*/
	xDelete<double>(seg_horiz_coords);
	xDelete<double>(seg_horiz_weights);
	xDelete<double>(seg_vert_coords);
	xDelete<double>(seg_vert_weights);
}
/*}}}*/
GaussPenta::GaussPenta(int index,IssmDouble r1,IssmDouble r2,bool mainlyfloating,int order){/*{{{*/

	IssmDouble x,y;
	IssmDouble xy_list[3][2];

	if(mainlyfloating){
		/*Get gauss points*/
		GaussLegendreTria(&this->numgauss,&this->coords1,&this->coords2,&this->coords3,&this->weights,order);

		xy_list[0][0]=0;  xy_list[0][1]=0; 
		xy_list[1][0]=r1; xy_list[1][1]=0; 
		xy_list[2][0]=0;  xy_list[2][1]=r2; 

		for(int ii=0;ii<this->numgauss;ii++){
			x = this->coords1[ii]*xy_list[0][0] + this->coords2[ii]*xy_list[1][0] + this->coords3[ii]*xy_list[2][0];
			y = this->coords1[ii]*xy_list[0][1] + this->coords2[ii]*xy_list[1][1] + this->coords3[ii]*xy_list[2][1];

			switch(index){
				case 0:
					this->coords1[ii] = 1.-x-y;
					this->coords2[ii] = x;
					this->coords3[ii] = y;
					break;
				case 1:
					this->coords1[ii] = y;
					this->coords2[ii] = 1.-x-y;
					this->coords3[ii] = x;
					break;
				case 2:
					this->coords1[ii] = x;
					this->coords2[ii] = y;
					this->coords3[ii] = 1.-x-y;
					break;
				default:
					_error_("index "<<index<<" not supported yet");
			}
			this->weights[ii] = this->weights[ii]*r1*r2;
		}
		this->coords4=xNew<IssmDouble>(numgauss);
		for(int ii=0;ii<numgauss;ii++) this->coords4[ii]=-1.0;
	}
	else{
		/*Double number of gauss points*/
		GaussPenta *gauss1    = NULL;
		GaussPenta *gauss2    = NULL;
		gauss1=new GaussPenta(0,1,2,order);
		gauss2=new GaussPenta(0,1,2,order);

		xy_list[0][0]=r1; xy_list[0][1]=0; 
		xy_list[1][0]=0;  xy_list[1][1]=1.; 
		xy_list[2][0]=0;  xy_list[2][1]=r2; 

			//gauss1->Echo();
		for(int ii=0;ii<gauss1->numgauss;ii++){
			x = gauss1->coords1[ii]*xy_list[0][0] + gauss1->coords2[ii]*xy_list[1][0] + gauss1->coords3[ii]*xy_list[2][0];
			y = gauss1->coords1[ii]*xy_list[0][1] + gauss1->coords2[ii]*xy_list[1][1] + gauss1->coords3[ii]*xy_list[2][1];

			switch(index){
				case 0:
					gauss1->coords1[ii] = 1.-x-y;
					gauss1->coords2[ii] = x;
					gauss1->coords3[ii] = y;
					break;
				case 1:
					gauss1->coords1[ii] = y;
					gauss1->coords2[ii] = 1.-x-y;
					gauss1->coords3[ii] = x;
					break;
				case 2:
					gauss1->coords1[ii] = x;
					gauss1->coords2[ii] = y;
					gauss1->coords3[ii] = 1.-x-y;
					break;
				default:
					_error_("index "<<index<<" not supported yet");
			}
			gauss1->weights[ii] = gauss1->weights[ii]*r1*(1-r2);
		}
			//gauss1->Echo();
		xy_list[0][0]=r1; xy_list[0][1]=0; 
		xy_list[1][0]=1.; xy_list[1][1]=0; 
		xy_list[2][0]=0;  xy_list[2][1]=1.; 

			//gauss2->Echo();
		for(int ii=0;ii<gauss2->numgauss;ii++){
			x = gauss2->coords1[ii]*xy_list[0][0] + gauss2->coords2[ii]*xy_list[1][0] + gauss2->coords3[ii]*xy_list[2][0];
			y = gauss2->coords1[ii]*xy_list[0][1] + gauss2->coords2[ii]*xy_list[1][1] + gauss2->coords3[ii]*xy_list[2][1];

			switch(index){
				case 0:
					gauss2->coords1[ii] = 1.-x-y;
					gauss2->coords2[ii] = x;
					gauss2->coords3[ii] = y;
					break;
				case 1:
					gauss2->coords1[ii] = y;
					gauss2->coords2[ii] = 1.-x-y;
					gauss2->coords3[ii] = x;
					break;
				case 2:
					gauss2->coords1[ii] = x;
					gauss2->coords2[ii] = y;
					gauss2->coords3[ii] = 1.-x-y;
					break;
				default:
					_error_("index "<<index<<" not supported yet");
			}
			gauss2->weights[ii] = gauss2->weights[ii]*(1-r1);
		}

		this->numgauss = gauss1->numgauss + gauss2->numgauss;
		this->coords1=xNew<IssmDouble>(this->numgauss);
		this->coords2=xNew<IssmDouble>(this->numgauss);
		this->coords3=xNew<IssmDouble>(this->numgauss);
		this->coords4=xNew<IssmDouble>(this->numgauss);
		this->weights=xNew<IssmDouble>(this->numgauss);

		for(int ii=0;ii<gauss1->numgauss;ii++){ // Add the first triangle gauss points
			this->coords1[ii]=gauss1->coords1[ii];
			this->coords2[ii]=gauss1->coords2[ii];
			this->coords3[ii]=gauss1->coords3[ii];
			this->coords4[ii]=gauss1->coords4[ii];
			this->weights[ii]=gauss1->weights[ii];
		}
		for(int ii=0;ii<gauss2->numgauss;ii++){ // Add the second triangle gauss points
			this->coords1[gauss1->numgauss+ii]=gauss2->coords1[ii];
			this->coords2[gauss1->numgauss+ii]=gauss2->coords2[ii];
			this->coords3[gauss1->numgauss+ii]=gauss2->coords3[ii];
			this->coords4[gauss1->numgauss+ii]=gauss2->coords4[ii];
			this->weights[gauss1->numgauss+ii]=gauss2->weights[ii];
		}

		/*Delete gauss points*/
		delete gauss1;
		delete gauss2;
	}

	/*Initialize static fields as undefined*/
	weight=UNDEF;
	coord1=UNDEF;
	coord2=UNDEF;
	coord3=UNDEF;
	coord4=UNDEF;
   ig    = -1;
}
/*}}}*/
GaussPenta::GaussPenta(IssmDouble area_coordinates[4][3],int order_horiz,int order_vert){/*{{{*/

	/*Intermediaties*/
	IssmPDouble *seg_horiz_coords  = NULL;
	IssmPDouble *seg_horiz_weights = NULL;
	IssmPDouble *seg_vert_coords   = NULL;
	IssmPDouble *seg_vert_weights  = NULL;

	/*get the gauss points using the product of two line rules*/
	GaussLegendreLinear(&seg_horiz_coords,&seg_horiz_weights,order_horiz);
	GaussLegendreLinear(&seg_vert_coords, &seg_vert_weights, order_vert);

	/*Allocate GaussPenta fields*/
   ig      = -1;
	numgauss=order_horiz*order_vert;
	coords1=xNew<IssmDouble>(numgauss);
	coords2=xNew<IssmDouble>(numgauss);
	coords3=xNew<IssmDouble>(numgauss);
	coords4=xNew<IssmDouble>(numgauss);
	weights=xNew<IssmDouble>(numgauss);

	/*Quads: get the gauss points using the product of two line rules  */
	for(int i=0;i<order_horiz;i++){
		for(int j=0;j<order_vert;j++){
			coords1[i*order_vert+j]=0.5*(area_coordinates[0][0]+area_coordinates[1][0]) + 0.5*seg_horiz_coords[i]*(area_coordinates[1][0]-area_coordinates[0][0]);
			coords2[i*order_vert+j]=0.5*(area_coordinates[0][1]+area_coordinates[1][1]) + 0.5*seg_horiz_coords[i]*(area_coordinates[1][1]-area_coordinates[0][1]);
			coords3[i*order_vert+j]=0.5*(area_coordinates[0][2]+area_coordinates[1][2]) + 0.5*seg_horiz_coords[i]*(area_coordinates[1][2]-area_coordinates[0][2]);
			coords4[i*order_vert+j]=seg_vert_coords[j];
			weights[i*order_vert+j]=seg_horiz_weights[i]*seg_vert_weights[j];
		}
	}

	/*clean-up*/
	xDelete<IssmPDouble>(seg_horiz_coords);
	xDelete<IssmPDouble>(seg_horiz_weights);
	xDelete<IssmPDouble>(seg_vert_coords);
	xDelete<IssmPDouble>(seg_vert_weights);
}
/*}}}*/
GaussPenta::GaussPenta(IssmDouble area_coordinates[2][3],int order_horiz){/*{{{*/

	/*Intermediaties*/
	IssmPDouble *seg_horiz_coords  = NULL;
	IssmPDouble *seg_horiz_weights = NULL;

	/*get the gauss points using the product of two line rules*/
	GaussLegendreLinear(&seg_horiz_coords,&seg_horiz_weights,order_horiz);

	/*Allocate GaussPenta fields*/
   ig = -1;
	numgauss=order_horiz;
	coords1=xNew<IssmDouble>(numgauss);
	coords2=xNew<IssmDouble>(numgauss);
	coords3=xNew<IssmDouble>(numgauss);
	coords4=xNew<IssmDouble>(numgauss);
	weights=xNew<IssmDouble>(numgauss);

	/*Quads: get the gauss points using the product of two line rules  */
	for(int i=0;i<order_horiz;i++){
		coords1[i]=0.5*(area_coordinates[0][0]+area_coordinates[1][0]) + 0.5*seg_horiz_coords[i]*(area_coordinates[1][0]-area_coordinates[0][0]);
		coords2[i]=0.5*(area_coordinates[0][1]+area_coordinates[1][1]) + 0.5*seg_horiz_coords[i]*(area_coordinates[1][1]-area_coordinates[0][1]);
		coords3[i]=0.5*(area_coordinates[0][2]+area_coordinates[1][2]) + 0.5*seg_horiz_coords[i]*(area_coordinates[1][2]-area_coordinates[0][2]);
		coords4[i]=0.;
		weights[i]=seg_horiz_weights[i];
	}

	/*clean-up*/
	xDelete<IssmPDouble>(seg_horiz_coords);
	xDelete<IssmPDouble>(seg_horiz_weights);
}
/*}}}*/
GaussPenta::~GaussPenta(){/*{{{*/
	xDelete<IssmDouble>(weights);
	xDelete<IssmDouble>(coords1);
	xDelete<IssmDouble>(coords2);
	xDelete<IssmDouble>(coords3);
	xDelete<IssmDouble>(coords4);
}
/*}}}*/

/*Methods*/
void GaussPenta::Echo(void){/*{{{*/

	_printf_("GaussPenta:\n");
	_printf_("   numgauss: " << numgauss << "\n");

	if (weights){
	 _printf_("   weights = ["); 
	 for(int i=0;i<numgauss;i++) _printf_(" " << weights[i] << "\n");
	 _printf_("]\n");
	}
	else _printf_("weights = NULL\n");
	if (coords1){
	 _printf_("   coords1 = ["); 
	 for(int i=0;i<numgauss;i++) _printf_(" " << coords1[i] << "\n");
	 _printf_("]\n");
	}
	else _printf_("coords1 = NULL\n");
	if (coords2){
	 _printf_("   coords2 = ["); 
	 for(int i=0;i<numgauss;i++) _printf_(" " << coords2[i] << "\n");
	 _printf_("]\n");
	}
	else _printf_("coords2 = NULL\n");
	if (coords3){
	 _printf_("   coords3 = ["); 
	 for(int i=0;i<numgauss;i++) _printf_(" " << coords3[i] << "\n");
	 _printf_("]\n");
	}
	else _printf_("coords3 = NULL\n");
	if (coords4){
		_printf_("   coords4 = ["); 
		for(int i=0;i<numgauss;i++) _printf_(" " << coords4[i] << "\n");
		_printf_("]\n");
	}
	else _printf_("coords4 = NULL\n");

	_printf_("   weight = " << weight << "\n");
	_printf_("   coord1 = " << coord1 << "\n");
	_printf_("   coord2 = " << coord2 << "\n");
	_printf_("   coord3 = " << coord3 << "\n");
	_printf_("   coord4 = " << coord4 << "\n");

}
/*}}}*/
int GaussPenta::Enum(void){/*{{{*/
	return GaussPentaEnum;
}
/*}}}*/
void GaussPenta::GaussFaceTria(int index1, int index2, int index3, int order){/*{{{*/

	/*in debugging mode: check that the default constructor has been called*/
	_assert_(numgauss==-1);

	/*Basal Tria*/
	if(index1==0 && index2==1 && index3==2){
		GaussLegendreTria(&numgauss,&coords1,&coords2,&coords3,&weights,order);
		coords4=xNew<IssmDouble>(numgauss);
		for(int i=0;i<numgauss;i++) coords4[i]=-1.0;
	}
	else{
		_error_("Tria not supported yet");
	}

}
/*}}}*/
void GaussPenta::GaussPoint(int ig){/*{{{*/

	/*Check input in debugging mode*/
	 _assert_(ig>=0 && ig< numgauss);

	 /*update static arrays*/
	 weight=weights[ig];
	 coord1=coords1[ig];
	 coord2=coords2[ig];
	 coord3=coords3[ig];
	 coord4=coords4[ig];

}
/*}}}*/
bool GaussPenta::next(void){/*{{{*/

	/*Increment Gauss point*/
	this->ig++;

	/*Have we reached the end?*/
	if(this->ig==this->numgauss) return false;

	/*If not let's go to the next point*/
	_assert_(this->ig>=0 && this->ig< numgauss);
	weight=weights[ig];
	coord1=coords1[ig];
	coord2=coords2[ig];
	coord3=coords3[ig];
	coord4=coords4[ig];

	return true;
}/*}}}*/
void GaussPenta::GaussNode(int finiteelement,int iv){/*{{{*/

	/*in debugging mode: check that the default constructor has been called*/
	_assert_(numgauss==-1);

	/*update static arrays*/
	switch(finiteelement){
		case P1Enum: case P1DGEnum:
		case P1P1GLSEnum: case P1P1Enum: /* added to allow P1-P1 GLS */
			switch(iv){
				case 0: coord1=1.; coord2=0.; coord3=0.; coord4=-1.; break;
				case 1: coord1=0.; coord2=1.; coord3=0.; coord4=-1.; break;
				case 2: coord1=0.; coord2=0.; coord3=1.; coord4=-1.; break;
				case 3: coord1=1.; coord2=0.; coord3=0.; coord4=+1.; break;
				case 4: coord1=0.; coord2=1.; coord3=0.; coord4=+1.; break;
				case 5: coord1=0.; coord2=0.; coord3=1.; coord4=+1.; break;
				default: _error_("node index should be in [0 5]");
			}
			break;
		case P1xP2Enum: 
			switch(iv){
				case 0: coord1=1.; coord2=0.; coord3=0.; coord4=-1.; break;
				case 1: coord1=0.; coord2=1.; coord3=0.; coord4=-1.; break;
				case 2: coord1=0.; coord2=0.; coord3=1.; coord4=-1.; break;
				case 3: coord1=1.; coord2=0.; coord3=0.; coord4=+1.; break;
				case 4: coord1=0.; coord2=1.; coord3=0.; coord4=+1.; break;
				case 5: coord1=0.; coord2=0.; coord3=1.; coord4=+1.; break;

				case 6: coord1=1.; coord2=0.; coord3=0.; coord4=0.; break;
				case 7: coord1=0.; coord2=1.; coord3=0.; coord4=0.; break;
				case 8: coord1=0.; coord2=0.; coord3=1.; coord4=0.; break;
				default: _error_("node index should be in [0 8]");
			}
			break;
		case P1xP3Enum: 
			switch(iv){
				case 0 : coord1=1.; coord2=0.; coord3=0.; coord4=-1.; break;
				case 1 : coord1=0.; coord2=1.; coord3=0.; coord4=-1.; break;
				case 2 : coord1=0.; coord2=0.; coord3=1.; coord4=-1.; break;
				case 3 : coord1=1.; coord2=0.; coord3=0.; coord4=+1.; break;
				case 4 : coord1=0.; coord2=1.; coord3=0.; coord4=+1.; break;
				case 5 : coord1=0.; coord2=0.; coord3=1.; coord4=+1.; break;

				case 6 : coord1=1.; coord2=0.; coord3=0.; coord4=-1./3.; break;
				case 7 : coord1=0.; coord2=1.; coord3=0.; coord4=-1./3.; break;
				case 8 : coord1=0.; coord2=0.; coord3=1.; coord4=-1./3.; break;
				case 9 : coord1=1.; coord2=0.; coord3=0.; coord4=+1./3.; break;
				case 10: coord1=0.; coord2=1.; coord3=0.; coord4=+1./3.; break;
				case 11: coord1=0.; coord2=0.; coord3=1.; coord4=+1./3.; break;
				default: _error_("node index should be in [0 11]");
			}
			break;
		case P1xP4Enum: 
			switch(iv){
				case 0 : coord1=1.; coord2=0.; coord3=0.; coord4=-1.; break;
				case 1 : coord1=0.; coord2=1.; coord3=0.; coord4=-1.; break;
				case 2 : coord1=0.; coord2=0.; coord3=1.; coord4=-1.; break;
				case 3 : coord1=1.; coord2=0.; coord3=0.; coord4=+1.; break;
				case 4 : coord1=0.; coord2=1.; coord3=0.; coord4=+1.; break;
				case 5 : coord1=0.; coord2=0.; coord3=1.; coord4=+1.; break;

				case 6 : coord1=1.; coord2=0.; coord3=0.; coord4=0.; break;
				case 7 : coord1=0.; coord2=1.; coord3=0.; coord4=0.; break;
				case 8 : coord1=0.; coord2=0.; coord3=1.; coord4=0.; break;

				case 9 : coord1=1.; coord2=0.; coord3=0.; coord4=-0.5; break;
				case 10: coord1=0.; coord2=1.; coord3=0.; coord4=-0.5; break;
				case 11: coord1=0.; coord2=0.; coord3=1.; coord4=-0.5; break;

				case 12: coord1=1.; coord2=0.; coord3=0.; coord4=+0.5; break;
				case 13: coord1=0.; coord2=1.; coord3=0.; coord4=+0.5; break;
				case 14: coord1=0.; coord2=0.; coord3=1.; coord4=+0.5; break;
				default: _error_("node index should be in [0 14]");
			}
			break;
		case P2xP1Enum: 
			switch(iv){
				case 0: coord1=1.; coord2=0.; coord3=0.; coord4=-1.; break;
				case 1: coord1=0.; coord2=1.; coord3=0.; coord4=-1.; break;
				case 2: coord1=0.; coord2=0.; coord3=1.; coord4=-1.; break;
				case 3: coord1=1.; coord2=0.; coord3=0.; coord4=+1.; break;
				case 4: coord1=0.; coord2=1.; coord3=0.; coord4=+1.; break;
				case 5: coord1=0.; coord2=0.; coord3=1.; coord4=+1.; break;

				case  6: coord1=0.; coord2=.5; coord3=.5; coord4=-1.;break;
				case  7: coord1=.5; coord2=0.; coord3=.5; coord4=-1.;break;
				case  8: coord1=.5; coord2=.5; coord3=0.; coord4=-1.;break;
				case  9: coord1=0.; coord2=.5; coord3=.5; coord4=+1.;break;
				case 10: coord1=.5; coord2=0.; coord3=.5; coord4=+1.;break;
				case 11: coord1=.5; coord2=.5; coord3=0.; coord4=+1.;break;
				default: _error_("node index should be in [0 11]");
			}
			break;
		case P1bubbleEnum:  case P1bubblecondensedEnum:
			switch(iv){
				case 0: coord1=1.;    coord2=0.;    coord3=0.;    coord4=-1.; break;
				case 1: coord1=0.;    coord2=1.;    coord3=0.;    coord4=-1.; break;
				case 2: coord1=0.;    coord2=0.;    coord3=1.;    coord4=-1.; break;
				case 3: coord1=1.;    coord2=0.;    coord3=0.;    coord4=+1.; break;
				case 4: coord1=0.;    coord2=1.;    coord3=0.;    coord4=+1.; break;
				case 5: coord1=0.;    coord2=0.;    coord3=1.;    coord4=+1.; break;
				case 6: coord1=1./3.; coord2=1./3.; coord3=1./3.; coord4=0.;  break;
				default: _error_("node index should be in [0 6]");
			}
			break;
		case P2Enum:
			switch(iv){
				case 0: coord1=1.; coord2=0.; coord3=0.; coord4=-1.; break;
				case 1: coord1=0.; coord2=1.; coord3=0.; coord4=-1.; break;
				case 2: coord1=0.; coord2=0.; coord3=1.; coord4=-1.; break;
				case 3: coord1=1.; coord2=0.; coord3=0.; coord4=+1.; break;
				case 4: coord1=0.; coord2=1.; coord3=0.; coord4=+1.; break;
				case 5: coord1=0.; coord2=0.; coord3=1.; coord4=+1.; break;

				case 6: coord1=1.; coord2=0.; coord3=0.; coord4=0.; break;
				case 7: coord1=0.; coord2=1.; coord3=0.; coord4=0.; break;
				case 8: coord1=0.; coord2=0.; coord3=1.; coord4=0.; break;

				case  9: coord1=0.; coord2=.5; coord3=.5; coord4=-1.;break;
				case 10: coord1=.5; coord2=0.; coord3=.5; coord4=-1.;break;
				case 11: coord1=.5; coord2=.5; coord3=0.; coord4=-1.;break;
				case 12: coord1=0.; coord2=.5; coord3=.5; coord4=+1.;break;
				case 13: coord1=.5; coord2=0.; coord3=.5; coord4=+1.;break;
				case 14: coord1=.5; coord2=.5; coord3=0.; coord4=+1.;break;

				case 15: coord1=0.; coord2=.5; coord3=.5; coord4=0.;break;
				case 16: coord1=.5; coord2=0.; coord3=.5; coord4=0.;break;
				case 17: coord1=.5; coord2=.5; coord3=0.; coord4=0.;break;
				default: _error_("node index should be in [0 17]");
			}
			break;
		case P2bubbleEnum: case P2bubblecondensedEnum:
			switch(iv){
				case 0: coord1=1.; coord2=0.; coord3=0.; coord4=-1.; break;
				case 1: coord1=0.; coord2=1.; coord3=0.; coord4=-1.; break;
				case 2: coord1=0.; coord2=0.; coord3=1.; coord4=-1.; break;
				case 3: coord1=1.; coord2=0.; coord3=0.; coord4=+1.; break;
				case 4: coord1=0.; coord2=1.; coord3=0.; coord4=+1.; break;
				case 5: coord1=0.; coord2=0.; coord3=1.; coord4=+1.; break;

				case 6: coord1=1.; coord2=0.; coord3=0.; coord4=0.; break;
				case 7: coord1=0.; coord2=1.; coord3=0.; coord4=0.; break;
				case 8: coord1=0.; coord2=0.; coord3=1.; coord4=0.; break;

				case  9: coord1=0.; coord2=.5; coord3=.5; coord4=-1.;break;
				case 10: coord1=.5; coord2=0.; coord3=.5; coord4=-1.;break;
				case 11: coord1=.5; coord2=.5; coord3=0.; coord4=-1.;break;
				case 12: coord1=0.; coord2=.5; coord3=.5; coord4=+1.;break;
				case 13: coord1=.5; coord2=0.; coord3=.5; coord4=+1.;break;
				case 14: coord1=.5; coord2=.5; coord3=0.; coord4=+1.;break;

				case 15: coord1=0.; coord2=.5; coord3=.5; coord4=0.;break;
				case 16: coord1=.5; coord2=0.; coord3=.5; coord4=0.;break;
				case 17: coord1=.5; coord2=.5; coord3=0.; coord4=0.;break;

				case 18: coord1=1./3.; coord2=1./3.; coord3=1./3.; coord4=0.;  break;
				default: _error_("node index should be in [0 18]");
			}
			break;
		case P2xP4Enum:
			switch(iv){
				case 0: coord1=1.; coord2=0.; coord3=0.; coord4=-1.; break;
				case 1: coord1=0.; coord2=1.; coord3=0.; coord4=-1.; break;
				case 2: coord1=0.; coord2=0.; coord3=1.; coord4=-1.; break;
				case 3: coord1=1.; coord2=0.; coord3=0.; coord4=+1.; break;
				case 4: coord1=0.; coord2=1.; coord3=0.; coord4=+1.; break;
				case 5: coord1=0.; coord2=0.; coord3=1.; coord4=+1.; break;

				case 6: coord1=1.; coord2=0.; coord3=0.; coord4=0.; break;
				case 7: coord1=0.; coord2=1.; coord3=0.; coord4=0.; break;
				case 8: coord1=0.; coord2=0.; coord3=1.; coord4=0.; break;

				case  9: coord1=0.; coord2=.5; coord3=.5; coord4=-1.;break;
				case 10: coord1=.5; coord2=0.; coord3=.5; coord4=-1.;break;
				case 11: coord1=.5; coord2=.5; coord3=0.; coord4=-1.;break;
				case 12: coord1=0.; coord2=.5; coord3=.5; coord4=+1.;break;
				case 13: coord1=.5; coord2=0.; coord3=.5; coord4=+1.;break;
				case 14: coord1=.5; coord2=.5; coord3=0.; coord4=+1.;break;

				case 15: coord1=1.; coord2=0.; coord3=0.; coord4=-.5; break;
				case 16: coord1=0.; coord2=1.; coord3=0.; coord4=-.5; break;
				case 17: coord1=0.; coord2=0.; coord3=1.; coord4=-.5; break;
				case 18: coord1=1.; coord2=0.; coord3=0.; coord4=+.5; break;
				case 19: coord1=0.; coord2=1.; coord3=0.; coord4=+.5; break;
				case 20: coord1=0.; coord2=0.; coord3=1.; coord4=+.5; break;

				case 21: coord1=0.; coord2=.5; coord3=.5; coord4=-.5;break;
				case 22: coord1=.5; coord2=0.; coord3=.5; coord4=-.5;break;
				case 23: coord1=.5; coord2=.5; coord3=0.; coord4=-.5;break;
				case 24: coord1=0.; coord2=.5; coord3=.5; coord4=0.;break;
				case 25: coord1=.5; coord2=0.; coord3=.5; coord4=0.;break;
				case 26: coord1=.5; coord2=.5; coord3=0.; coord4=0.;break;
				case 27: coord1=0.; coord2=.5; coord3=.5; coord4=+.5;break;
				case 28: coord1=.5; coord2=0.; coord3=.5; coord4=+.5;break;
				case 29: coord1=.5; coord2=.5; coord3=0.; coord4=+.5;break;
				default: _error_("node index should be in [0 29]");
			}
			break;
		default: _error_("Finite element "<<EnumToStringx(finiteelement)<<" not supported");
	}

}
/*}}}*/
void GaussPenta::GaussVertex(int iv){/*{{{*/

	/*in debugging mode: check that the default constructor has been called*/
	_assert_(numgauss==-1);

	/*update static arrays*/
	switch(iv){
		case 0: coord1=1.; coord2=0.; coord3=0.; coord4= -1.; break;
		case 1: coord1=0.; coord2=1.; coord3=0.; coord4= -1.; break;
		case 2: coord1=0.; coord2=0.; coord3=1.; coord4= -1.; break;
		case 3: coord1=1.; coord2=0.; coord3=0.; coord4= +1.; break;
		case 4: coord1=0.; coord2=1.; coord3=0.; coord4= +1.; break;
		case 5: coord1=0.; coord2=0.; coord3=1.; coord4= +1.; break;
		default: _error_("vertex index should be in [0 5]");

	}

}
/*}}}*/
void GaussPenta::Reset(void){/*{{{*/

	/*Check that this has been initialized*/
	_assert_(numgauss>0);

	/*Reset counter*/
	this->ig=-1;
} /*}}}*/
void GaussPenta::SynchronizeGaussBase(Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	GaussTria* gauss_tria = xDynamicCast<GaussTria*>(gauss);

	gauss_tria->coord1=this->coord1;
	gauss_tria->coord2=this->coord2;
	gauss_tria->coord3=this->coord3;
	gauss_tria->weight=UNDEF;
}
/*}}}*/
