/*!\file Shp2Kmlx
 * \brief shp to kml conversion routines.
 */

#include "./Shp2Kmlx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../kml/kmlobjects.h"

int Shp2Kmlx(char* filshp,char* filkml,int sgn){

	#ifdef _HAVE_SHAPELIB_ //only works if Shapelib library has been compiled in.

	double  cm,sp;
	Xy2lldef(&cm,&sp,sgn);

	return(Shp2Kmlx(filshp,filkml,sgn,cm,sp));

	#else //ifdef _HAVE_SHAPELIB_
	return 0;
	#endif
}

int Shp2Kmlx(char* filshp,char* filkml,int sgn,double cm,double sp){

	#ifdef _HAVE_SHAPELIB_ //only works if Shapelib library has been compiled in.

	int     i,j,k,iret=0;
	int     lwidth=1;
	double  popac=0.50;
	int     nshape,ncoord;
	double  cpsum;
	int     *pstype = NULL, *pnpart=NULL,**ppstrt=NULL,**pptype=NULL,*pnvert=NULL;
	double **pshapx = NULL,**pshapy=NULL,**pshapz=NULL,**pshapm=NULL;
	double  *lat    = NULL, *lon=NULL;

	SHPHandle   hSHP;
	int     nShapeType, nEntities, iPart, bValidate = 0,nInvalidCount=0;
	const char  *pszPlus;
	double  adfMinBound[4], adfMaxBound[4];

	char    indent[81]="";
	KML_File          *kfile  = NULL;
	KML_Document      *kdoc   = NULL;
	KML_Style         *kstyle = NULL;
	KML_LineStyle     *klsty  = NULL;
	KML_PolyStyle     *kpsty  = NULL;
	KML_Folder        *kfold  = NULL;
	KML_Placemark     *kplace = NULL;
	KML_MultiGeometry *kmulti = NULL;
	KML_Polygon       *kpoly  = NULL;
	KML_LinearRing    *kring  = NULL;
	KML_LineString    *kline  = NULL;
	KML_Point         *kpoint = NULL;
	FILE              *fid    = NULL;

	clock_t clock0,clock1;
	time_t  time0, time1;

	clock0=clock();
	time0 =time(NULL);
	_printf0_("\nShp2Kmlx Module -- " << ctime(&time0));

/*  note that much of the following code is taken from shpdump.c in shapelib.  */

/*  open shp/shx files  */

	hSHP = SHPOpen( filshp, "rb" );
	if (!hSHP) _error_("Error opening shp/shx files.");

/*  read header and print out file bounds  */

	SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

	printf( "Shapefile Type: %s   # of Shapes: %d\n\n",
			SHPTypeName( nShapeType ), nEntities );

	printf( "File Bounds: (%12.3f,%12.3f,%g,%g)\n"
			"         to  (%12.3f,%12.3f,%g,%g)\n",
			adfMinBound[0],
			adfMinBound[1],
			adfMinBound[2],
			adfMinBound[3],
			adfMaxBound[0],
			adfMaxBound[1],
			adfMaxBound[2],
			adfMaxBound[3] );

	nshape=nEntities;
	pstype=xNew<int>(nshape);
	pnpart=xNew<int>(nshape);
	ppstrt=xNew<int*>(nshape);
	pptype=xNew<int*>(nshape);
	pnvert=xNew<int>(nshape);
	pshapx=xNew<double*>(nshape);
	pshapy=xNew<double*>(nshape);
	pshapz=xNew<double*>(nshape);
	pshapm=xNew<double*>(nshape);

	/* loop over the list of shapes  */
	for(i=0;i<nEntities;i++ ){
		SHPObject   *psShape;

	psShape = SHPReadObject( hSHP, i );

	printf( "\nShape:%d (%s)  nVertices=%d, nParts=%d\n"
				"  Bounds:(%12.3f,%12.3f, %g, %g)\n"
				"      to (%12.3f,%12.3f, %g, %g)\n",
			i, SHPTypeName(psShape->nSHPType),
				psShape->nVertices, psShape->nParts,
				psShape->dfXMin, psShape->dfYMin,
				psShape->dfZMin, psShape->dfMMin,
				psShape->dfXMax, psShape->dfYMax,
				psShape->dfZMax, psShape->dfMMax );

	pstype[i]=psShape->nSHPType;
	pnpart[i]=psShape->nParts;
	if (pnpart[i]) {
		ppstrt[i]=xNew<int>(pnpart[i]);
		pptype[i]=xNew<int>(pnpart[i]);
	}
	else {
		ppstrt[i]=NULL;
		pptype[i]=NULL;
	}
	pnvert[i]=psShape->nVertices;
	if (pnvert[i]) {
		pshapx[i]=xNew<double>(pnvert[i]);
		pshapy[i]=xNew<double>(pnvert[i]);
		pshapz[i]=xNew<double>(pnvert[i]);
		pshapm[i]=xNew<double>(pnvert[i]);
	}
	else {
		pshapx[i]=NULL;
		pshapy[i]=NULL;
		pshapz[i]=NULL;
		pshapm[i]=NULL;
	}

	for( j = 0, iPart = 1; j < psShape->nVertices; j++ )
	{
			const char  *pszPartType = "";

			if( j == 0 && psShape->nParts > 0 )
			{
				pszPartType = SHPPartTypeName( psShape->panPartType[0] );
				ppstrt[i][0]=psShape->panPartStart[0];
				pptype[i][0]=psShape->panPartType[0];
			}

		if( iPart < psShape->nParts
				&& psShape->panPartStart[iPart] == j )
		{
				pszPartType = SHPPartTypeName( psShape->panPartType[iPart] );
				ppstrt[i][iPart]=psShape->panPartStart[iPart];
				pptype[i][iPart]=psShape->panPartType[iPart];
		iPart++;
		pszPlus = "+";
		}
		else
			pszPlus = " ";

//		printf("   %s (%12.3f,%12.3f, %g, %g) %s \n",
//				   pszPlus,
//				   psShape->padfX[j],
//				   psShape->padfY[j],
//				   psShape->padfZ[j],
//				   psShape->padfM[j],
//				   pszPartType );

		pshapx[i][j]=psShape->padfX[j];
		pshapy[i][j]=psShape->padfY[j];
		pshapz[i][j]=psShape->padfZ[j];
		pshapm[i][j]=psShape->padfM[j];
	}

		if( bValidate )
		{
			int nAltered = SHPRewindObject( hSHP, psShape );

			if( nAltered > 0 )
			{
				printf( "  %d rings wound in the wrong direction.\n",
						nAltered );
				nInvalidCount++;
			}
		}

		SHPDestroyObject( psShape );
	}

/*  close shp/shx files  */

	SHPClose( hSHP );

/*  construct kml file  */

	kfile =new KML_File();
	kfile->AddAttrib("xmlns","http://www.opengis.net/kml/2.2");

/*  construct kml document  */

	kdoc  =new KML_Document();
	sprintf(kdoc->name      ,"Shp2Kmlx Module -- %s",ctime(&time0));
	kdoc->open      =1;

/*  construct style templates for defaults  */

	klsty =new KML_LineStyle();
	sprintf(klsty->color     ,"ff000000");
	sprintf(klsty->colormode ,"normal");
	klsty->width     =lwidth;
	kpsty =new KML_PolyStyle();
	sprintf(kpsty->color     ,"%02xffffff",(int)floor(popac*255+0.5));
	sprintf(kpsty->colormode ,"random");
	kstyle=new KML_Style();
	kstyle->AddAttrib("id","BlackLineRandomPoly");
	kstyle->line      =klsty;
	kstyle->poly      =kpsty;
	(kdoc->style     )->AddObject((Object*)kstyle);

	klsty =new KML_LineStyle();
	sprintf(klsty->color     ,"ff000000");
	sprintf(klsty->colormode ,"normal");
	klsty->width     =lwidth;
	kpsty =new KML_PolyStyle();
	sprintf(kpsty->color     ,"00ffffff");
	sprintf(kpsty->colormode ,"random");
	kstyle=new KML_Style();
	kstyle->AddAttrib("id","BlackLineEmptyPoly");
	kstyle->line      =klsty;
	kstyle->poly      =kpsty;
	(kdoc->style     )->AddObject((Object*)kstyle);

	klsty =new KML_LineStyle();
	sprintf(klsty->color     ,"%02xffffff",(int)floor(popac*255+0.5));
	sprintf(klsty->colormode ,"random");
	klsty->width     =lwidth*2;
	kpsty =new KML_PolyStyle();
	sprintf(kpsty->color     ,"00ffffff");
	sprintf(kpsty->colormode ,"random");
	kstyle=new KML_Style();
	kstyle->AddAttrib("id","RandomLineEmptyPoly");
	kstyle->line      =klsty;
	kstyle->poly      =kpsty;
	(kdoc->style     )->AddObject((Object*)kstyle);

/*  construct kml folder for shapes  */

	kfold =new KML_Folder();
	sprintf(kfold->name      ,"Shapefile: %s  Type: %s  nShapes: %d",
			filshp, SHPTypeName( nShapeType ), nEntities );
	kfold->open      =1;

/*  loop over the list of shapes  */

	for (i=0; i<nshape; i++) {

/*  null type  */

		if      (pstype[i] == SHPT_NULL) {
			;
		}

/*  point types  */

		else if (pstype[i] == SHPT_POINT ||
			  	 pstype[i] == SHPT_POINTZ ||
			 	 pstype[i] == SHPT_POINTM) {
			kplace=new KML_Placemark();

			sprintf(kplace->name      ,"Shape:%d (%s)  nVertices=%d, nParts=%d",
					i,SHPTypeName(pstype[i]),pnvert[i],pnpart[i]);
			kplace->visibility=true;
			sprintf(kplace->styleurl  ,"#RandomLineEmptyPoly");

			if (pnpart[i] > 0)
				_printf_("Warning -- Shape "<< i << " of type \"" << SHPTypeName( pstype[i] ) << "\" should not have " << pnpart[i] << " > 0 parts.\n\n");
			if (pnvert[i] > 1)
				_printf_("Warning -- Shape " << i << " of type \"" << SHPTypeName( pstype[i] ) << "\" should not have " << pnpart[i] << " > 1 vertices.\n\n");

			kpoint=new KML_Point();

			lat=xNew<double>(pnvert[i]);
			lon=xNew<double>(pnvert[i]);
			if (sgn) {
				Xy2llx(lat,lon,pshapx[i],pshapy[i],pnvert[i],sgn,cm,sp);
			}
			else  {
				memcpy(lon,pshapx[i],pnvert[i]*sizeof(double));
				memcpy(lat,pshapy[i],pnvert[i]*sizeof(double));
			}

			kpoint->coords[0]=lon      [0];
			kpoint->coords[1]=lat      [0];
			kpoint->coords[2]=pshapz[i][0];

			xDelete<double>(lon);
			xDelete<double>(lat);

			(kplace->geometry  )->AddObject((Object*)kpoint);
			kpoint=NULL;
			(kfold ->feature   )->AddObject((Object*)kplace);
			kplace=NULL;
		}

/*  polyline types  */

		else if (pstype[i] == SHPT_ARC ||
				 pstype[i] == SHPT_ARCZ ||
				 pstype[i] == SHPT_ARCM) {
			kplace=new KML_Placemark();

			sprintf(kplace->name      ,"Shape:%d (%s)  nVertices=%d, nParts=%d",
					i,SHPTypeName(pstype[i]),pnvert[i],pnpart[i]);
			kplace->visibility=true;
			sprintf(kplace->styleurl  ,"#RandomLineEmptyPoly");

/*  create a multigeometry to hold all the lines  */

			kmulti=new KML_MultiGeometry();

/*  convert to lat/lon, if necessary  */

			lat=xNew<double>(pnvert[i]);
			lon=xNew<double>(pnvert[i]);
			if (sgn) {
				Xy2llx(lat,lon,pshapx[i],pshapy[i],pnvert[i],sgn,cm,sp);
			}
			else  {
				memcpy(lon,pshapx[i],pnvert[i]*sizeof(double));
				memcpy(lat,pshapy[i],pnvert[i]*sizeof(double));
			}

/*  loop over the lines  */

			for (j=0; j<pnpart[i]; j++) {
				kline =new KML_LineString();

				kline->ncoord    =(j<pnpart[i]-1 ? ppstrt[i][j+1]-ppstrt[i][j] : pnvert[i]-ppstrt[i][j]);
				kline->coords    =xNew<double>(kline->ncoord*3);
				for (k=0; k<kline->ncoord; k++) {
					kline->coords[3*k+0]=lon      [ppstrt[i][j]+k];
					kline->coords[3*k+1]=lat      [ppstrt[i][j]+k];
					kline->coords[3*k+2]=pshapz[i][ppstrt[i][j]+k];
				}
				(kmulti->geometry  )->AddObject((Object*)kline);
				kline = NULL;
			}

			xDelete<double>(lon);
			xDelete<double>(lat);

			(kplace->geometry)->AddObject((Object*)kmulti);
			kmulti=NULL;
			(kfold ->feature )->AddObject((Object*)kplace);
			kplace=NULL;
		}

/*  polygon types  */

		else if (pstype[i] == SHPT_POLYGON ||
				 pstype[i] == SHPT_POLYGONZ ||
				 pstype[i] == SHPT_POLYGONM) {

/*  the shp format specifies that outer rings are cw, while inner rings are ccw.  there
	may be multiple outer rings and inner rings in any order.  the kml format specifies
	all rings are ccw (right-hand rule).  there may be only one outer ring with multiple
	inner rings, and rings are differentiated by keyword.

	at least for now, assume that each cw ring forms a new kml polygon, and each ccw
	ring forms an inner ring for the most recent outer ring.  a more elaborate solution
	would be for each inner ring to search in which outer ring it occurs.  */

			kplace=new KML_Placemark();

			sprintf(kplace->name      ,"Shape:%d (%s)  nVertices=%d, nParts=%d",
					i,SHPTypeName(pstype[i]),pnvert[i],pnpart[i]);
			kplace->visibility=true;
			sprintf(kplace->styleurl  ,"#BlackLineRandomPoly");

/*  create a multigeometry to hold all the polygons  */

			kmulti=new KML_MultiGeometry();

/*  convert to lat/lon, if necessary  */

			lat=xNew<double>(pnvert[i]);
			lon=xNew<double>(pnvert[i]);
			if (sgn) {
				Xy2llx(lat,lon,pshapx[i],pshapy[i],pnvert[i],sgn,cm,sp);
			}
			else  {
				memcpy(lon,pshapx[i],pnvert[i]*sizeof(double));
				memcpy(lat,pshapy[i],pnvert[i]*sizeof(double));
			}

/*  loop over the polygons  */

			for (j=0; j<pnpart[i]; j++) {

/*  check if polygon is ccw or cw by computing sum of cross products (twice the area)  */

				ncoord=(j<pnpart[i]-1 ? ppstrt[i][j+1]-ppstrt[i][j] : pnvert[i]-ppstrt[i][j]);
				cpsum =0.;

				for (k=ppstrt[i][j]; k<ppstrt[i][j]+ncoord-1; k++)
					cpsum +=pshapx[i][k]*pshapy[i][k+1         ]-pshapy[i][k]*pshapx[i][k+1         ];
				cpsum +=pshapx[i][k]*pshapy[i][ppstrt[i][j]]-pshapy[i][k]*pshapx[i][ppstrt[i][j]];

/*  outer ring (cw) (allow exception for single-part shapes)  */

				if (cpsum < 0 || pnpart[i] == 1) {
					if (kpoly) {
						(kmulti->geometry  )->AddObject((Object*)kpoly);
						kpoly =NULL;
					}

/*  create a new polygon from the outer ring (reversing cw to ccw)  */

					kpoly =new KML_Polygon();
					kring =new KML_LinearRing();

					kring->ncoord    =(j<pnpart[i]-1 ? ppstrt[i][j+1]-ppstrt[i][j] : pnvert[i]-ppstrt[i][j]);
					kring->coords    =xNew<double>(kring->ncoord*3);
					if (cpsum < 0)
						for (k=0; k<kring->ncoord; k++) {
							kring->coords[3*(kring->ncoord-1-k)+0]=lon      [ppstrt[i][j]+k];
							kring->coords[3*(kring->ncoord-1-k)+1]=lat      [ppstrt[i][j]+k];
							kring->coords[3*(kring->ncoord-1-k)+2]=pshapz[i][ppstrt[i][j]+k];
						}
					else
						for (k=0; k<kring->ncoord; k++) {
							kring->coords[3*k+0]=lon      [ppstrt[i][j]+k];
							kring->coords[3*k+1]=lat      [ppstrt[i][j]+k];
							kring->coords[3*k+2]=pshapz[i][ppstrt[i][j]+k];
						}

					(kpoly ->outer     )->AddObject((Object*)kring);
					kring =NULL;
				}

/*  inner ring (ccw)  */

				else {
					if (!kpoly) {
						_printf_("Warning -- Shape " << i << " of type \"" << SHPTypeName( pstype[i] ) << "\", part " << j << ", expected to be outer loop (cw).\n\n");
						continue;
					}

/*  add the inner ring to the current polygon  */

					kring =new KML_LinearRing();

					kring->ncoord    =(j<pnpart[i]-1 ? ppstrt[i][j+1]-ppstrt[i][j] : pnvert[i]-ppstrt[i][j]);
					kring->coords    =xNew<double>(kring->ncoord*3);
					for (k=0; k<kring->ncoord; k++) {
						kring->coords[3*k+0]=lon      [ppstrt[i][j]+k];
						kring->coords[3*k+1]=lat      [ppstrt[i][j]+k];
						kring->coords[3*k+2]=pshapz[i][ppstrt[i][j]+k];
					}

					(kpoly ->inner     )->AddObject((Object*)kring);
					kring =NULL;
				}
			}

			if (kpoly) {
				(kmulti->geometry  )->AddObject((Object*)kpoly);
				kpoly =NULL;
			}

			xDelete<double>(lon);
			xDelete<double>(lat);

			(kplace->geometry  )->AddObject((Object*)kmulti);
			kmulti=NULL;
			(kfold ->feature   )->AddObject((Object*)kplace);
			kplace=NULL;
		}

/*  multipoint types  */

		else if (pstype[i] == SHPT_MULTIPOINT ||
				 pstype[i] == SHPT_MULTIPOINTZ ||
				 pstype[i] == SHPT_MULTIPOINTM) {
			kplace=new KML_Placemark();

			sprintf(kplace->name      ,"Shape:%d (%s)  nVertices=%d, nParts=%d",
					i,SHPTypeName(pstype[i]),pnvert[i],pnpart[i]);
			kplace->visibility=true;
			sprintf(kplace->styleurl  ,"#RandomLineEmptyPoly");

			if (pnpart[i] > 0)
				_printf_("Warning -- Shape " << i << " of type \"" << SHPTypeName( pstype[i] ) << "\" should not have " << pnpart[i] << " > 0 parts.\n\n");

/*  create a multigeometry to hold all the points  */

			kmulti=new KML_MultiGeometry();

/*  convert to lat/lon, if necessary  */

			lat=xNew<double>(pnvert[i]);
			lon=xNew<double>(pnvert[i]);
			if (sgn) {
				Xy2llx(lat,lon,pshapx[i],pshapy[i],pnvert[i],sgn,cm,sp);
			}
			else  {
				memcpy(lon,pshapx[i],pnvert[i]*sizeof(double));
				memcpy(lat,pshapy[i],pnvert[i]*sizeof(double));
			}

/*  loop over the points  */

			for (j=0; j<pnvert[i]; j++) {
				kpoint=new KML_Point();

				kpoint->coords[0]=lon      [j];
				kpoint->coords[1]=lat      [j];
				kpoint->coords[2]=pshapz[i][j];

				(kmulti->geometry  )->AddObject((Object*)kpoint);
				kpoint=NULL;
			}

			xDelete<double>(lon);
			xDelete<double>(lat);

			(kplace->geometry  )->AddObject((Object*)kmulti);
			kmulti=NULL;
			(kfold ->feature   )->AddObject((Object*)kplace);
			kplace=NULL;
		}

/*  multipatch types  */

		else if (pstype[i] == SHPT_MULTIPATCH) {
			_printf_("Warning -- Shape " << i << " of type \"" << SHPTypeName( pstype[i] ) << "\" will be ignored.\n\n");
			continue;
		}

/*  unknown type  */

		else {
			_printf_("Warning -- Shape " << i << " of type \"" << SHPTypeName( pstype[i] ) << "\" will be ignored.\n\n");
		}
	}

/*  assemble the rest of the kml hierarchy  */

	(kdoc ->feature   )->AddObject((Object*)kfold);
	kfold=NULL;
	(kfile->kmlobj    )->AddObject((Object*)kdoc);
	kdoc =NULL;

/*  write kml file  */

	_printf0_("Writing kml document to file.\n");
	fid=fopen(filkml,"w");
	fprintf(fid,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
	kfile->Write(fid,indent);
	fclose(fid);

	delete kfile;
	for (i=nshape-1; i>=0; i--) {
		xDelete<double>((pshapm[i]));
		xDelete<double>((pshapz[i]));
		xDelete<double>((pshapy[i]));
		xDelete<double>((pshapx[i]));
	}
	xDelete<double*>(pshapm);
	xDelete<double*>(pshapz);
	xDelete<double*>(pshapy);
	xDelete<double*>(pshapx);
	xDelete<int>(pnvert);
	for (i=nshape-1; i>=0; i--) {
		xDelete<int>((pptype[i]));
		xDelete<int>((ppstrt[i]));
	}
	xDelete<int*>(pptype);
	xDelete<int*>(ppstrt);
	xDelete<int>(pnpart);
	xDelete<int>(pstype);

	clock1=clock();
	time1 =time(NULL);
	_printf_("Shp2Kmlx Module -- " << ((double)(clock1-clock0))/CLOCKS_PER_SEC << " CPU seconds; " << difftime(time1,time0) << " elapsed seconds.\n\n\n");

	return(iret);

	#else //ifdef _HAVE_SHAPELIB_
	return 0;
	#endif
}
