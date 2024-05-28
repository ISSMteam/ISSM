#include "./bamgobjects.h"
#include "../shared/shared.h"

/*Constructors/Destructors*/
BamgOpts::BamgOpts(){/*{{{*/

	this->anisomax          = 0;
	this->cutoff            = 0;
	this->coeff             = 0;
	this->errg              = 0;
	this->gradation         = 0;
	this->Hessiantype       = 0;
	this->maxnbv            = 0;
	this->maxsubdiv         = 0;
	this->Metrictype        = 0;
	this->nbjacobi          = 0;
	this->nbsmooth          = 0;
	this->omega             = 0;
	this->power             = 0;
	this->verbose           = 0;

	this->Crack             = 0;
	this->KeepVertices      = 0;
	this->splitcorners      = 0;

	this->hmin              = 0;
	this->hmax              = 0;
	this->hminVertices=NULL; this->hminVerticesSize[0]=this->hminVerticesSize[1]=0;
	this->hmaxVertices=NULL; this->hmaxVerticesSize[0]=this->hmaxVerticesSize[1]=0;
	this->hVertices=NULL;    this->hVerticesLength=0;
	this->metric=NULL;       this->metricSize[0]=this->metricSize[1]=0;
	this->field=NULL;        this->fieldSize[0]=this->fieldSize[1]=0;
	this->err=NULL;          this->errSize[0]=this->errSize[1]=0;

}
/*}}}*/
BamgOpts::~BamgOpts(){/*{{{*/

	xDelete<double>(this->hminVertices);
	xDelete<double>(this->hmaxVertices);
	xDelete<double>(this->hVertices);
	xDelete<double>(this->metric);
	xDelete<double>(this->field);
	xDelete<double>(this->err);

}
/*}}}*/

/*Methods*/
void BamgOpts::Check(void){/*{{{*/

	if(this->anisomax<1) _error_("'anisomax' option should be >=1");
	if(this->coeff==0) _error_("'coeff' should be positive");
	if(this->errg<0) _error_("'errg' option should be >0");
	if(this->gradation<1) _error_("'gradation' option should be >=1");
	if(this->Hessiantype!=0  && this->Hessiantype!=1) _error_("'Hessiantype' supported options are 0 and 1");
	if(this->maxnbv<3) _error_("'maxnbv' option should be >3");
	if(this->maxsubdiv<=1) _error_("'maxsubdiv' should be >1");
	if(this->Metrictype!=0   && this->Metrictype!=1 && this->Metrictype!=2) _error_("'Metrictype' supported options are 0, 1 and 2");
	if(this->nbjacobi<=0) _error_("'nbjacobi' option should be >0");
	if(this->nbsmooth<=0) _error_("'nbsmooth' option should be >0");

	if(this->Crack!=0  && this->Crack!=1) _error_("'Crack' supported options are 0 and 1");
	if(this->KeepVertices!=0 && this->KeepVertices!=1) _error_("'KeepVertices' supported options are 0 and 1");

	if(this->hmin<=0) _error_("'hmin' option should be >0");
	if(this->hmax<=0 || this->hmax<this->hmin) _error_("'hmax' option should be between 0 and hmin=" << this->hmin);
	if(this->hminVertices && this->hminVerticesSize[1]!=1) _error_("'hminVertices' should be a column");
	if(this->hmaxVertices && this->hmaxVerticesSize[1]!=1) _error_("'hmaxVertices' should be a column");
	if(this->hVertices && this->hVerticesLength<2) _error_("'hVertices' should be a vector");
	if(this->metric && (this->metricSize[1]!=1 && this->metricSize[1]!=3)) _error_("'metric' should have either 1 (iso) or 3 (aniso) columns.");
	if(this->field){
		if(this->errSize[1]!=this->fieldSize[1]) _error_("'err' should have "<< this->fieldSize[1]<<" columns");
		for(int i=0;i<this->fieldSize[1];i++){if (this->err[i]<=0) _error_("'err' option should be >0");};
	}

}
/*}}}*/
void BamgOpts::Echo(){/*{{{*/

	_printf_("anisomax: "    << this->anisomax<<"\n");
	_printf_("cutoff: "      << this->cutoff<<"\n");
	_printf_("coeff: "       << this->coeff<<"\n");
	_printf_("errg: "        << this->errg<<"\n");
	_printf_("gradation: "   << this->gradation<<"\n");
	_printf_("Hessiantype: " << this->Hessiantype<<"\n");
	_printf_("maxnbv: "      << this->maxnbv<<"\n");
	_printf_("maxsubdiv: "   << this->maxsubdiv<<"\n");
	_printf_("Metrictype: "  << this->Metrictype<<"\n");
	_printf_("nbjacobi: "    << this->nbjacobi<<"\n");
	_printf_("nbsmooth: "    << this->nbsmooth<<"\n");
	_printf_("omega: "       << this->omega<<"\n");
	_printf_("power: "       << this->power<<"\n");
	_printf_("verbose: "     << this->verbose<<"\n");
   _printf_("\n"<<"\n");
	_printf_("Crack: "        << this->Crack<<"\n");
	_printf_("KeepVertices: " << this->KeepVertices<<"\n");
	_printf_("splitcorners: " << this->splitcorners<<"\n");

	_printf_("hmin: "<<this->hmin<<"\n");
	_printf_("hmax: "<<this->hmax<<"\n");
   _printf_("hVerticesLength: "  << this->hVerticesLength     << "\n");
   _printf_("hminVerticesSize: " << this->hminVerticesSize[0] << "x"    << this->hminVerticesSize[1] << "\n");
   _printf_("hmaxVerticesSize: " << this->hmaxVerticesSize[0] << "x"    << this->hmaxVerticesSize[1] << "\n");
   _printf_("metricSize: "       << this->metricSize[0]       << "x"    << this->metricSize[1]       << "\n");
   _printf_("fieldSize: "        << this->fieldSize[0]        << "x"    << this->fieldSize[1]        << "\n");
   _printf_("errSize: "          << this->errSize[0]          << "x"    << this->errSize[1]          << "\n");

}
/*}}}*/
