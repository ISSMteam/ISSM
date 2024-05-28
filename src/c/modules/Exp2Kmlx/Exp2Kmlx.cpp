/*!\file Exp2Kmlx
 * \brief exp to kml conversion routines.
 */

#include "./Exp2Kmlx.h"
#include "../../shared/shared.h"
#include "../../kml/kmlobjects.h"
#include "../../toolkits/toolkits.h"

int Exp2Kmlx(char* filexp,char* filkml,int sgn,bool holes){

	double  cm,sp;

	Xy2lldef(&cm,&sp,sgn);
	return(Exp2Kmlx(filexp,filkml,sgn,cm,sp,holes));
}

int Exp2Kmlx(char* filexp,char* filkml,int sgn,double cm,double sp,bool holes){

	int      i        ,j,iret=0;
	int      lwidth = 1;
	double   popac  = 0.50;
	int      nprof;
	int     *pnvert = NULL;
	double **pprofx = NULL,**pprofy=NULL;
	bool    *closed = NULL;
	double  *lat    = NULL, *lon=NULL;

	char    indent[81]="";
	KML_File*          kfile =NULL;
	KML_Document*      kdoc  =NULL;
	KML_Style*         kstyle=NULL;
	KML_LineStyle*     klsty =NULL;
	KML_PolyStyle*     kpsty =NULL;
	KML_Folder*        kfold =NULL;
	KML_Placemark*     kplace=NULL;
	KML_Polygon*       kpoly =NULL;
	KML_LinearRing*    kring =NULL;
	KML_LineString*    kline =NULL;
	KML_Point*         kpoint=NULL;

	FILE*   fid=NULL;

	clock_t clock0,clock1;
	time_t  time0, time1;

	clock0=clock();
	time0 =time(NULL);
	_printf0_("\nExp2Kmlx Module -- " << ctime(&time0));

	/*read exp file  */

	if (!ExpRead(&nprof,&pnvert,&pprofx,&pprofy,&closed,filexp))
		_error_("Error reading exp file.");
	_printf0_("Exp2Kmlx -- Reading " << nprof << " exp profiles from file \"" << filexp << "\".\n");
//	for (i=0; i<nprof; i++)
//		_printf_("i=" << i << "; nvert=" << pnvert[i] << ", closed=" << closed[i] << "\n");

/*  construct kml file  */

	kfile =new KML_File();
	kfile->AddAttrib("xmlns","http://www.opengis.net/kml/2.2");

/*  construct kml document  */

	kdoc  =new KML_Document();
	sprintf(kdoc->name      ,"Exp2Kmlx Module -- %s",ctime(&time0));
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

/*  construct kml folder for polygons  */

	kfold =new KML_Folder();
	sprintf(kfold->name      ,"Profiles translated from file \"%s\".",filexp);
	kfold->open      =1;

/*  polygon with multiple holes  */

	if (holes && nprof && (pnvert[0] <= 1 || pprofx[0][pnvert[0]-1] != pprofx[0][0] || pprofy[0][pnvert[0]-1] != pprofy[0][0])) {
		_printf0_("Warning -- Outer profile is not closed, so \"holes\" option will be ignored.\n");
		holes=false;
	}

	if (holes) {
		i=0;
		kplace=new KML_Placemark();
		sprintf(kplace->name      ,"Polygon with Holes");
		kplace->visibility=true;
		sprintf(kplace->styleurl  ,"#BlackLineRandomPoly");

		kpoly =new KML_Polygon();
		kring =new KML_LinearRing();

		kring->ncoord    =pnvert[i]-1;
		lat=xNew<double>(kring->ncoord);
		lon=xNew<double>(kring->ncoord);
		Xy2llx(lat,lon,pprofx[i],pprofy[i],kring->ncoord,sgn,cm,sp);
		kring->coords=xNew<double>(kring->ncoord*3);
		for (j=0; j<kring->ncoord; j++) {
			kring->coords[3*j+0]=lon[j];
			kring->coords[3*j+1]=lat[j];
			kring->coords[3*j+2]=0.;
		}
		xDelete<double>(lon);
		xDelete<double>(lat);

		(kpoly ->outer     )->AddObject((Object*)kring);
		kring =NULL;

		for (i=1; i<nprof; i++) {
			if (pnvert[i] <= 1 || pprofx[i][pnvert[i]-1] != pprofx[i][0] || pprofy[i][pnvert[i]-1] != pprofy[i][0]) {
				_printf0_("Warning -- Inner profile " << i+1 << " is not closed with \"holes\" specified, so it will be ignored.\n");
				continue;
			}

			kring =new KML_LinearRing();

			kring->ncoord    =pnvert[i]-1;
			lat=xNew<double>(kring->ncoord);
			lon=xNew<double>(kring->ncoord);
			Xy2llx(lat,lon,pprofx[i],pprofy[i],kring->ncoord,sgn,cm,sp);
			kring->coords    =xNew<double>(kring->ncoord*3);
			for (j=0; j<kring->ncoord; j++) {
				kring->coords[3*j+0]=lon[j];
				kring->coords[3*j+1]=lat[j];
				kring->coords[3*j+2]=0.;
			}
			xDelete<double>(lon);
			xDelete<double>(lat);

			(kpoly ->inner     )->AddObject((Object*)kring);
			kring =NULL;
		}

		(kplace->geometry  )->AddObject((Object*)kpoly);
		kpoly =NULL;
		(kfold ->feature   )->AddObject((Object*)kplace);
		kplace=NULL;
	}

/*  multiple polygons or linestrings  */

	else {
		for (i=0; i<nprof; i++) {
			kplace=new KML_Placemark();

			if     (pnvert[i] > 1 && pprofx[i][pnvert[i]-1] == pprofx[i][0] && pprofy[i][pnvert[i]-1] == pprofy[i][0]) {
				sprintf(kplace->name      ,"Polygon %d",i+1);
				kplace->visibility=true;
				sprintf(kplace->styleurl  ,"#BlackLineRandomPoly");

				kpoly =new KML_Polygon();
				kring =new KML_LinearRing();

				kring->ncoord    =pnvert[i]-1;
				lat=xNew<double>(kring->ncoord);
				lon=xNew<double>(kring->ncoord);
				Xy2llx(lat,lon,pprofx[i],pprofy[i],kring->ncoord,sgn,cm,sp);
				kring->coords    =xNew<double>(kring->ncoord*3);
				for (j=0; j<kring->ncoord; j++) {
					kring->coords[3*j+0]=lon[j];
					kring->coords[3*j+1]=lat[j];
					kring->coords[3*j+2]=0.;
				}
				xDelete<double>(lon);
				xDelete<double>(lat);

				(kpoly ->outer     )->AddObject((Object*)kring);
				kring =NULL;

				(kplace->geometry  )->AddObject((Object*)kpoly);
				kpoly =NULL;
			}

			else if (pnvert[i] > 1) {
				sprintf(kplace->name      ,"LineString %d",i+1);
				kplace->visibility=true;
				sprintf(kplace->styleurl  ,"#RandomLineEmptyPoly");

				kline =new KML_LineString();

				kline->ncoord    =pnvert[i];
				lat=xNew<double>(kline->ncoord);
				lon=xNew<double>(kline->ncoord);
				Xy2llx(lat,lon,pprofx[i],pprofy[i],kline->ncoord,sgn,cm,sp);
				kline->coords    =xNew<double>(kline->ncoord*3);
				for (j=0; j<kline->ncoord; j++) {
					kline->coords[3*j+0]=lon[j];
					kline->coords[3*j+1]=lat[j];
					kline->coords[3*j+2]=0.;
				}
				xDelete<double>(lon);
				xDelete<double>(lat);

				(kplace->geometry  )->AddObject((Object*)kline);
				kline =NULL;
			}

			else if (pnvert[i]) {
				sprintf(kplace->name      ,"Point %d",i+1);
				kplace->visibility=true;
				sprintf(kplace->styleurl  ,"#RandomLineEmptyPoly");
				int one=1;

				kpoint=new KML_Point();

				lat=xNew<double>(one);
				lon=xNew<double>(one);
				Xy2llx(lat,lon,pprofx[i],pprofy[i],1,sgn,cm,sp);
				kpoint->coords[0]=lon[0];
				kpoint->coords[1]=lat[0];
				kpoint->coords[2]=0.;
				xDelete<double>(lon);
				xDelete<double>(lat);

				(kplace->geometry  )->AddObject((Object*)kpoint);
				kpoint =NULL;
			}

			(kfold ->feature   )->AddObject((Object*)kplace);
			kplace=NULL;
		}
	}

/*  assemble the rest of the kml hierarchy  */

	(kdoc ->feature   )->AddObject((Object*)kfold);
	kfold=NULL;
	(kfile->kmlobj    )->AddObject((Object*)kdoc);
	kdoc =NULL;

/*  write kml file  */

	_printf0_("Exp2Kmlx -- Writing kml document to file \"" << filkml << "\".\n");
	fid=fopen(filkml,"w");
	fprintf(fid,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
	kfile->Write(fid,indent);
	fclose(fid);

	delete kfile;
	for (i=nprof-1; i>=0; i--) {
		xDelete<double>(pprofy[i]);
		xDelete<double>(pprofx[i]);
	}
	xDelete<int>(pnvert);

	clock1=clock();
	time1 =time(NULL);
	_printf_("Exp2Kmlx Module -- " <<((double)(clock1-clock0))/CLOCKS_PER_SEC << " CPU seconds; " <<difftime(time1,time0) << " elapsed seconds.\n\n\n");

	return(iret);
}
