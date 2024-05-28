/*\file ShpRead.c
 *\brief: shp to exp file conversion mex module.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./ShpRead.h"
#ifdef _HAVE_SHAPELIB_ //only works if Shapelib library has been compiled in.
#include "shapefil.h"
#endif

void ShpReadUsage(void){/*{{{*/
	_printf0_("ShpRead - Read shapefile\n");
	_printf0_("\n");
	_printf0_("   This module reads shapefiles and converts them to matlab/python structures\n");
	_printf0_("\n");
	_printf0_("   Usage:\n");
	_printf0_("      ShpRead(filename);\n");
	_printf0_("      filexp      file name of exp file to be written\n");
	_printf0_("\n");
	_printf0_("   Examples:\n");
	_printf0_("      ShpRead('file.shp');\n");
}/*}}}*/
WRAPPER(ShpRead_python){

	/*input: */
	char *filename= NULL;

	/*Boot module: */
	MODULEBOOT();

	#ifndef _HAVE_SHAPELIB_ //only works if shapelib library has been compiled in.
	_error_("Shapelib not available! Cannot carry out shp file translation!");
	#else

	/*checks on arguments on the matlab side: */
	if(nlhs != NLHS){ShpReadUsage(); _error_("ShpRead usage error");}
	if(nrhs != NRHS){ShpReadUsage(); _error_("ShpRead usage error");}

	/*Input datasets: */
	FetchData(&filename,SHP_IN);

	/*Open shapefile*/
	SHPHandle hSHP = SHPOpen( filename, "rb" );
	if(!hSHP) _error_("Error opening shp/shx files.");

	/*read header and print out file bounds*/
	int         nShapeType,nEntities;
	IssmPDouble adfMinBound[4], adfMaxBound[4];
	SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );
	_printf_("Shapefile Type: "<<SHPTypeName(nShapeType)<<"   number of Shapes: "<< nEntities<<"\n\n");

	/*Initialize output*/
	Contours* contours = new Contours();

	/*Read all objects*/
	for(int i=0; i<nEntities;i++ ){

		SHPObject* psShape = SHPReadObject(hSHP,i);
		_printf_( "Shape #"<<i<<" ("<<SHPTypeName(psShape->nSHPType)<<") nVertices="<<psShape->nVertices<<", nParts="<<psShape->nParts<<"\n");

		Contour<double> *contour = NULL;

		switch(psShape->nSHPType){
			case SHPT_POINTZ:
				contour=new Contour<double>(0,psShape->nVertices,psShape->padfX,psShape->padfY,false);
				break;
			case SHPT_ARC:
				contour=new Contour<double>(0,psShape->nVertices,psShape->padfX,psShape->padfY,false);
				break;
			default:
				_printf_("Shape type "<<SHPTypeName(psShape->nSHPType)<<" not supported yet, skipping...\n");
		}

		/*Add to contours and clean up*/
		if(contour) contours->AddObject(contour);
		SHPDestroyObject(psShape);
	}

	/*Write output*/
	WriteData(SHP_OUT,contours);

	/*Clean-up*/
	delete contours;
	xDelete<char>(filename);

	#endif
	/*end module: */
	MODULEEND();
}
