/*!\file GaussTetra.c
 * \brief: implementation of the GaussTetra object
 */

#include <math.h>
#include "./GaussTetra.h"
#include "../../shared/io/Print/Print.h"
#include "../../shared/Exceptions/exceptions.h"
#include "../../shared/MemOps/MemOps.h"
#include "../../shared/Enum/Enum.h"
#include "../../shared/Numerics/GaussPoints.h"
#include "../../shared/Numerics/constants.h"

/*GaussTetra constructors and destructors:*/
GaussTetra::GaussTetra(){/*{{{*/

	ig     = -1;
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
GaussTetra::GaussTetra(int order){/*{{{*/

	/*Get gauss points*/
	GaussLegendreTetra(&numgauss,&coords1,&coords2,&coords3,&coords4,&weights,order);

	/*Rescale weights if necessary*/
	IssmDouble sumweights = 0.;
	for(int i=0;i<numgauss;i++) sumweights += this->weights[i];
	if(sumweights==1.){
		for(int i=0;i<numgauss;i++) this->weights[i] = this->weights[i]/6.; /*rescale volume to 1/6*/
	}

	/*Check final weights in debugging mode*/
	#ifdef _ISSM_DEBUG_
	sumweights = 0.; for(int i=0;i<numgauss;i++) sumweights += this->weights[i];
	_assert_(sumweights>1./6.-1e-10);
	_assert_(sumweights<1./6.+1e-10);
	#endif

	/*Initialize static fields as undefinite*/
	weight=UNDEF;
	coord1=UNDEF;
	coord2=UNDEF;
	coord3=UNDEF;
	ig    = -1;
}
/*}}}*/
GaussTetra::GaussTetra(int index1,int index2,int index3,int order){/*{{{*/

	/*Basal Tria*/
	if(index1==0 && index2==1 && index3==2){
		GaussLegendreTria(&numgauss,&coords1,&coords2,&coords3,&weights,order);
		coords4=xNew<IssmDouble>(numgauss);
		for(int i=0;i<numgauss;i++) coords4[i]=0.;
	}
	else if(index1==0 && index2==3 && index3==1){
		GaussLegendreTria(&numgauss,&coords1,&coords2,&coords4,&weights,order);
		coords3=xNew<IssmDouble>(numgauss);
		for(int i=0;i<numgauss;i++) coords3[i]=0.;
	}
	else if(index1==1 && index2==3 && index3==2){
		GaussLegendreTria(&numgauss,&coords2,&coords3,&coords4,&weights,order);
		coords1=xNew<IssmDouble>(numgauss);
		for(int i=0;i<numgauss;i++) coords1[i]=0.;
	}
	else if(index1==0 && index2==2 && index3==3){
		GaussLegendreTria(&numgauss,&coords1,&coords3,&coords4,&weights,order);
		coords2=xNew<IssmDouble>(numgauss);
		for(int i=0;i<numgauss;i++) coords2[i]=0.;
	}
	else{
		_error_(index1 <<" "<<index2 <<" "<<index3 <<" Not supported yet");
	}
	this->ig = -1;
}
/*}}}*/
GaussTetra::~GaussTetra(){/*{{{*/
	xDelete<IssmDouble>(weights);
	xDelete<IssmDouble>(coords1);
	xDelete<IssmDouble>(coords2);
	xDelete<IssmDouble>(coords3);
	xDelete<IssmDouble>(coords4);
}
/*}}}*/

/*Methods*/
void GaussTetra::Echo(void){/*{{{*/

	_printf_("GaussTetra:\n");
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
int GaussTetra::Enum(void){/*{{{*/
	return GaussTetraEnum;
}
/*}}}*/
void GaussTetra::GaussPoint(int ig){/*{{{*/

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
bool GaussTetra::next(void){/*{{{*/

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
void GaussTetra::GaussNode(int finiteelement,int iv){/*{{{*/

	/*in debugging mode: check that the default constructor has been called*/
	_assert_(numgauss==-1);

	/*update static arrays*/
	switch(finiteelement){
		case P1Enum: case P1DGEnum:
			switch(iv){
				case 0: coord1=1.; coord2=0.; coord3=0.; coord4=0.; break;
				case 1: coord1=0.; coord2=1.; coord3=0.; coord4=0.; break;
				case 2: coord1=0.; coord2=0.; coord3=1.; coord4=0.; break;
				case 3: coord1=0.; coord2=0.; coord3=0.; coord4=1.; break;
				default: _error_("node index should be in [0 3]");
			}
			break;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			switch(iv){
				case 0: coord1=1.; coord2=0.; coord3=0.; coord4=0.; break;
				case 1: coord1=0.; coord2=1.; coord3=0.; coord4=0.; break;
				case 2: coord1=0.; coord2=0.; coord3=1.; coord4=0.; break;
				case 3: coord1=0.; coord2=0.; coord3=0.; coord4=1.; break;
				case 4: coord1=1./4.; coord2=1./4.; coord3=1./4.; coord4=1./4.; break;
				default: _error_("node index should be in [0 4]");
			}
			break;
		case P2Enum:
			switch(iv){
				case 0: coord1=1.; coord2=0.; coord3=0.; coord4=0.; break;
				case 1: coord1=0.; coord2=1.; coord3=0.; coord4=0.; break;
				case 2: coord1=0.; coord2=0.; coord3=1.; coord4=0.; break;
				case 3: coord1=0.; coord2=0.; coord3=0.; coord4=1.; break;

				case 4: coord1=0.; coord2=.5; coord3=.5; coord4=0.; break;
				case 5: coord1=.5; coord2=0.; coord3=.5; coord4=0.; break;
				case 6: coord1=.5; coord2=.5; coord3=0.; coord4=0.; break;
				case 7: coord1=.5; coord2=0.; coord3=0.; coord4=.5; break;
				case 8: coord1=0.; coord2=.5; coord3=0.; coord4=.5; break;
				case 9: coord1=0.; coord2=0.; coord3=.5; coord4=.5; break;
				default: _error_("node index should be in [0 9]");
			}
			break;
		default: _error_("Finite element "<<EnumToStringx(finiteelement)<<" not supported");
	}

}
/*}}}*/
void GaussTetra::GaussVertex(int iv){/*{{{*/

	/*in debugging mode: check that the default constructor has been called*/
	_assert_(numgauss==-1);

	/*update static arrays*/
	switch(iv){
		case 0: coord1=1.; coord2=0.; coord3=0.; coord4=0.; break;
		case 1: coord1=0.; coord2=1.; coord3=0.; coord4=0.; break;
		case 2: coord1=0.; coord2=0.; coord3=1.; coord4=0.; break;
		case 3: coord1=0.; coord2=0.; coord3=0.; coord4=1.; break;
		default: _error_("vertex index should be in [0 3]");

	}

}
/*}}}*/
void GaussTetra::Reset(void){/*{{{*/

	/*Check that this has been initialized*/
	_assert_(numgauss>0);

	/*Reset counter*/
	this->ig=-1;
} /*}}}*/
void GaussTetra::SynchronizeGaussBase(Gauss* gauss){/*{{{*/

	_error_("not supported");
}
/*}}}*/
