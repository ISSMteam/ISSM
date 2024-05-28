/*!\file KMLMeshWritex
 */

#include "./KMLMeshWritex.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void KMLMeshWritex(int* ierror,char* name,char* notes,int* elem,int melem,int nelem,int* nodecon,int mncon,int nncon,double* lat, double* lng,int* part,double* data, int mdata, int ndata,double* cmap, int mcmap, int ncmap,FILE* fid){

	int                 i,j,k,ipt=0,jpt=0,nnodes;
	int                 mxepg      = 25;
	int                 lwidth     = 1;
	double              popac      = 0.50;
	char                indent[81] = " ";
	char                cstr[81];
	double             *edata = NULL;
	bool ncfree=false, edfree=false;
	KML_Document       *kdoc = NULL;
	KML_Style          *kstyle;
	KML_LineStyle      *klsty;
	KML_PolyStyle      *kpsty;

	clock_t clock0,clock1,clock0a,clock0b,clock0c;
	time_t  time0, time1, time0a, time0b, time0c;

	clock0=clock();
	time0 =time(NULL);
	_printf0_("\nKMLMeshWritex Module -- " << ctime(&time0));

/*  construct kml document  */

	kdoc=new KML_Document();
	sprintf(kdoc->name      ,"ISSM Mesh: %s",name);
	kdoc->open      =1;
	sprintf(kdoc->descript  ,"%s",notes);

/*  write style templates for defaults and for each color of the matlab
	colormap (note that matlab colormap format is rgb, where each varies
	from 0 to 1, whereas the kml color format is aabbggrr, where each
	varies from 00 to ff.)  */

	klsty=new KML_LineStyle();
	sprintf(klsty->color     ,"ff000000");
	sprintf(klsty->colormode ,"normal");
	klsty->width     =lwidth;
	kpsty=new KML_PolyStyle();
	sprintf(kpsty->color     ,"%02xffffff",(int)floor(popac*255+0.5));
	sprintf(kpsty->colormode ,"random");
	kstyle=new KML_Style();
	kstyle->AddAttrib("id","BlackLineRandomPoly");
	kstyle->line      =klsty;
	kstyle->poly      =kpsty;
	(kdoc->style     )->AddObject((Object*)kstyle);

	klsty=new KML_LineStyle();
	sprintf(klsty->color     ,"ff000000");
	sprintf(klsty->colormode ,"normal");
	klsty->width     =lwidth;
	kpsty=new KML_PolyStyle();
	sprintf(kpsty->color     ,"00ffffff");
	sprintf(kpsty->colormode ,"random");
	kstyle=new KML_Style();
	kstyle->AddAttrib("id","BlackLineEmptyPoly");
	kstyle->line      =klsty;
	kstyle->poly      =kpsty;
	(kdoc->style     )->AddObject((Object*)kstyle);

	klsty=new KML_LineStyle();
	sprintf(klsty->color     ,"ff0000ff");
	sprintf(klsty->colormode ,"normal");
	klsty->width     =lwidth;
	kpsty=new KML_PolyStyle();
	sprintf(kpsty->color     ,"%02x0000ff",(int)floor(popac*255+0.5));
	sprintf(kpsty->colormode ,"random");
	kstyle=new KML_Style();
	kstyle->AddAttrib("id","RedLineRedPoly");
	kstyle->line      =klsty;
	kstyle->poly      =kpsty;
	(kdoc->style     )->AddObject((Object*)kstyle);

	if (cmap) {
		_printf0_("Writing " << mcmap << " Matlab colors as KML style templates.\n");
		ipt=0;
		for (i=0; i<mcmap; i++) {
			klsty=new KML_LineStyle();
//			sprintf(klsty->color     ,"ff000000");
			sprintf(klsty->color     ,"%02x%02x%02x%02x",
					(int)255,
					(int)floor(cmap[ipt+2]*255+0.5),
					(int)floor(cmap[ipt+1]*255+0.5),
					(int)floor(cmap[ipt  ]*255+0.5));
			sprintf(klsty->colormode ,"normal");
			klsty->width     =lwidth;
			kpsty=new KML_PolyStyle();
			sprintf(kpsty->color     ,"%02x%02x%02x%02x",
					(int)floor(popac*255+0.5),
					(int)floor(cmap[ipt+2]*255+0.5),
					(int)floor(cmap[ipt+1]*255+0.5),
					(int)floor(cmap[ipt  ]*255+0.5));
			sprintf(kpsty->colormode ,"normal");
			kstyle=new KML_Style();
			sprintf(cstr,"MatlabColor%d",i+1);
			kstyle->AddAttrib("id",cstr);
			kstyle->line      =klsty;
			kstyle->poly      =kpsty;
			(kdoc->style     )->AddObject((Object*)kstyle);
			ipt+=ncmap;
		}
	}
//	kdoc->DeepEcho();

/*  create the node connectivity table, if necessary
	(noting that rows do not need to be sorted, since the elements
	are consecutively numbered)  */

	if (!nodecon) {
		_printf0_("Creating the node connectivity table.\n");
		nncon=mxepg+1;
		nodecon=xNewZeroInit<int>(mncon*nncon);
		ncfree=true;

		jpt=0;
		for (i=0; i<melem; i++) {
			for (j=0; j<nelem; j++) {
				if (elem[jpt]) {
					ipt=(elem[jpt]-1)*nncon;
					if (nodecon[ipt+(nncon-1)] < mxepg) {
						nodecon[ipt+nodecon[ipt+(nncon-1)]]=i+1;
						nodecon[ipt+(nncon-1)]++;
					}
					else
						_error_("Nodal connectivity table needs more than specified " << mxepg << " columns.\n");
				}
				jpt++;
			}
		}
	}

/*  average nodal data to element data, if necessary
	(noting that multiple columns of data are handled here, but not
	yet below)  */

	if (data) {
		if      (mdata == melem)
			edata=data;

		else if (mdata == mncon) {
			_printf0_("Averaging nodal data to element data.\n");
			edata=xNewZeroInit<double>(melem*ndata);
			edfree=true;

			ipt=0;
			jpt=0;
			for (i=0; i<melem; i++) {
				nnodes=0;
				for (j=0; j<nelem; j++) {
					if (elem[jpt]) {
						for (k=0; k<ndata; k++)
							edata[ipt+k]+=data[(elem[jpt]-1)*ndata+k];
						nnodes++;
					}
					jpt++;
				}
				if (nnodes)
					for (k=0; k<ndata; k++)
						edata[ipt+k]/=(double)nnodes;
				ipt+=ndata;
			}
		}

		else
			_error_("Data matrix has incorrect number of " << mdata << " rows.\n");
	}

/*  write folder for mesh  */

	(kdoc ->feature   )->AddObject((Object*)KMLMeshElem(elem,melem,nelem,
														nodecon,mncon,nncon,
														lat,lng,
														edata,
														cmap,mcmap,ncmap));

	if(edfree) xDelete<double>(edata);
	if(ncfree) xDelete<int>(nodecon);
	clock0a=clock();
	time0a =time(NULL);
	_printf_("  Constructed kml document -- " << ((double)(clock0a-clock0))/CLOCKS_PER_SEC << " CPU seconds; " << difftime(time0a,time0) << " elapsed seconds.\n\n\n");

/*  write kml file  */

	_printf0_("Writing kml document to file.\n");
	fprintf(fid,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
	fprintf(fid,"<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n");
	kdoc->Write(fid,indent);
	fprintf(fid,"</kml>\n");
	clock0b=clock();
	time0b =time(NULL);
	_printf_("  Wrote kml file -- " << ((double)(clock0b-clock0a))/CLOCKS_PER_SEC << " CPU seconds; " << difftime(time0b,time0a) << " elapsed seconds.\n\n\n");

	_printf0_("Deleting kml document.\n");
	delete kdoc;
	clock0c=clock();
	time0c =time(NULL);
	_printf_("  Deleted kml document -- " << ((double)(clock0c-clock0b))/CLOCKS_PER_SEC << " CPU seconds; " << difftime(time0c,time0b) << " elapsed seconds.\n\n\n");

	clock1=clock();
	time1 =time(NULL);
	_printf_("KMLMeshWritex Module -- " << ((double)(clock1-clock0))/CLOCKS_PER_SEC << " CPU seconds; " << difftime(time1,time0) << " elapsed seconds.\n\n\n");

	return;
}

KML_Folder* KMLMeshElem(int* elem,int melem,int nelem,
						int* nodecon,int mncon,int nncon,
						double* lat, double* lng,
						double* edata,
						double* cmap, int mcmap, int ncmap){

	int     i,j,ipt=0;
	double  alt=0;
	double  cmin= DBL_MAX,
			cmax=-DBL_MAX;
	int     imap;
	KML_Folder*     kfold =NULL;
	KML_Placemark*  kplace=NULL;
	KML_Polygon*    kpoly =NULL;
	KML_LinearRing* kring =NULL;

/*  write folder for mesh  */

	kfold=new KML_Folder();
//	sprintf(kfold->name      ,"Mesh");
	sprintf(kfold->name      ,"ISSM Targets");
	kfold->visibility=1;
//	sprintf(kfold->descript  ,"Elements=%d, Nodes=%d",melem,mncon);
	sprintf(kfold->descript  ,"campaign{\n");
	strcat(kfold->descript  ,"  evaluator ClaspTargetEvaluator;\n");
	strcat(kfold->descript  ,"  solver IssmSolver;\n");
	strcat(kfold->descript  ,"  spacecraft airplane ClaspSpacecraft(\n");
	strcat(kfold->descript  ,"    dutyCycleDuration=0,\n");
	strcat(kfold->descript  ,"    dutyCycleOnDuration=0,\n");
	strcat(kfold->descript  ,"    memoryInit=15000,\n");
	strcat(kfold->descript  ,"    memoryLimit=40000000,\n");
	strcat(kfold->descript  ,"    maxDataCollectionRate=27,\n");
	strcat(kfold->descript  ,"    maxDataDownlinkRate=10);\n");
	strcat(kfold->descript  ,"\n");
	strcat(kfold->descript  ,"  //sensor names\n");
	strcat(kfold->descript  ,"  sensor qqp_swath = 2,102,1002,1102;\n");
	strcat(kfold->descript  ,"\n");
	strcat(kfold->descript  ,"  //sensor ids to modes\n");
	strcat(kfold->descript  ,"  low_bandwidth_single_pol = 2,102,1002,1102;\n");
	strcat(kfold->descript  ,"  single_pol = 2,102,1002,1102;\n");
	strcat(kfold->descript  ,"  dual_pol = 2,102,1002,1102;\n");
	strcat(kfold->descript  ,"  quad_pol = 2,102,1002,1102;\n");
	strcat(kfold->descript  ,"\n");
	strcat(kfold->descript  ,"  //LRAD\n");
	strcat(kfold->descript  ,"  //Note all targets are \"ascending right\"-- i.e. mode=2\n");
	strcat(kfold->descript  ,"  left = 1002,1102;\n");
	strcat(kfold->descript  ,"  right = 2,102;\n");
	strcat(kfold->descript  ,"  ascending = 2,1002;\n");
	strcat(kfold->descript  ,"  descending = 102,1102;\n");
	strcat(kfold->descript  ,"\n");
	strcat(kfold->descript  ,"  //data rates\n");
	strcat(kfold->descript  ,"  low_bandwidth_single_pol datarate = 0.896;\n");
	strcat(kfold->descript  ,"  single_pol datarate = 4.214;\n");
	strcat(kfold->descript  ,"  dual_pol datarate = 8.428;\n");
	strcat(kfold->descript  ,"  quad_pol datarate = 16.856;\n");
	strcat(kfold->descript  ,"\n");
	strcat(kfold->descript  ,"  //mode domination relationships\n");
	strcat(kfold->descript  ,"  quad_pol dominates low_bandwidth_single_pol;\n");
	strcat(kfold->descript  ,"  quad_pol dominates single_pol;\n");
	strcat(kfold->descript  ,"  quad_pol dominates dual_pol;\n");
	strcat(kfold->descript  ,"  dual_pol dominates low_bandwidth_single_pol;\n");
	strcat(kfold->descript  ,"  dual_pol dominates single_pol;\n");
	strcat(kfold->descript  ,"  single_pol dominates low_bandwidth_single_pol;\n");
	strcat(kfold->descript  ,"\n");
	strcat(kfold->descript  ,"  //sensor styles\n");
	strcat(kfold->descript  ,"  2 0xff00ffff 0xff000000;\n");
	strcat(kfold->descript  ,"  102 0x7f00ffff 0xff00ffff;\n");
	strcat(kfold->descript  ,"  1002 0xffffff00 0xffffff00;\n");
	strcat(kfold->descript  ,"  1102 0x7fffff00 0xffffff00;\n");
	strcat(kfold->descript  ,"\n");
	strcat(kfold->descript  ,"  //discipline styles\n");
	strcat(kfold->descript  ,"  deformation 0xff006090 0xff006090 0xff0000ff 0xff1010ff;\n");
	strcat(kfold->descript  ,"  vegetation  0xff00ff00 0xff00ff00 0xff0000ff 0xff0020ff;\n");
	strcat(kfold->descript  ,"  ice         0xffff0000 0xffff0000 0xff0000ff 0xff2000ff;\n");
	strcat(kfold->descript  ,"}");

	if (edata)
		for (i=0; i<melem; i++) {
			if (edata[i] < cmin)
				cmin=edata[i];
			if (edata[i] > cmax)
				cmax=edata[i];
		}

/*  write each element as a polygon placemark  */

	_printf0_("Writing " << melem << " tria elements as KML polygons.\n");

	for (i=0; i<melem; i++) {
		kplace=new KML_Placemark();
		sprintf(kplace->name      ,"Element %d",(i+1));
		kplace->visibility=1;
		if (edata) {
//			sprintf(kplace->descript  ,"Element data: %g",edata[i]);
			sprintf(kplace->descript  ,"campaign{\n  deformation 1 %g quad_pol ascending right asap;\n}",edata[i]);
			imap = (int)floor((edata[i]-cmin)/(cmax-cmin)*mcmap+0.5)+1;
			if      ((imap >= 1) && (imap <= mcmap))
				sprintf(kplace->styleurl  ,"#MatlabColor%d",imap);
			else if (edata[i] == cmax)
				sprintf(kplace->styleurl  ,"#MatlabColor%d",mcmap);
			else
				sprintf(kplace->styleurl  ,"#BlackLineEmptyPoly");
		}
		else {
			sprintf(kplace->descript  ,"");
			sprintf(kplace->styleurl  ,"#BlackLineRandomPoly");
		}
//		kplace->DeepEcho();

		kpoly=new KML_Polygon();
		kpoly->extrude   =1;
		sprintf(kpoly->altmode   ,"clampToGround");
//		kpoly->DeepEcho();

		kring=new KML_LinearRing();
		kring->ncoord    =nelem+1;
		kring->coords =xNew<double>((nelem+1)*3);

/*  write the nodal coordinates as a linear ring  */

		for (j=0; j<nelem; j++) {
			kring->coords[3*j+0]=lng[elem[ipt]-1];
			kring->coords[3*j+1]=lat[elem[ipt]-1];
			kring->coords[3*j+2]=alt;
			ipt++;
		}
		kring->coords[3*nelem+0]=kring->coords[3*0+0];
		kring->coords[3*nelem+1]=kring->coords[3*0+1];
		kring->coords[3*nelem+2]=kring->coords[3*0+2];
//		kring->DeepEcho();

/*  assemble the linear ring into polygon into placemark into folder  */

		(kpoly ->outer   )->AddObject((Object*)kring);
		(kplace->geometry)->AddObject((Object*)kpoly);
		(kfold ->feature )->AddObject((Object*)kplace);

//		if (!(int)fmod((double)(i+1),1000))
//			_printf0_("  " << (i+1) << " tria elements written.\n");
	}
	_printf0_("  " << melem << " tria elements written.\n");

	return(kfold);
}
