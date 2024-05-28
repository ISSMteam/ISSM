/*!\file KMLOverlayx
 */

#include "./KMLOverlayx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../kml/kmlobjects.h"

void KMLOverlayx(int* ierror,
				 double* lataxis, double* longaxis,
				 int nimages, char** pimages,
				 FILE* fid){

	int     i;
	char    indent[81]="";
	KML_File*          kfile=NULL;
	KML_Document*      kdoc=NULL;
	KML_Folder*        kfold=NULL;
	KML_GroundOverlay* kgover=NULL;
	KML_Icon*          kicon=NULL;
	KML_LatLonBox*     kllbox=NULL;

	clock_t clock0,clock1;
	time_t  time0, time1;

	clock0=clock();
	time0 =time(NULL);
	_printf0_("\nKMLOverlayx Module -- " << ctime(&time0));

/*  construct kml file  */

	kfile=new KML_File();
	kfile->AddAttrib("xmlns","http://www.opengis.net/kml/2.2");

/*  construct kml document  */

	kdoc=new KML_Document();
	sprintf(kdoc->name      ,"Ground Overlays from ISSM");
	kdoc->open      =1;

/*  construct kml folder for overlays  */

	kfold=new KML_Folder();
	sprintf(kfold->name      ,"Ground Overlays");
	kfold->open      =1;

/*  construct ground overlay, icon, and lat/long box for each image  */

	for (i=0; i<nimages; i++) {
		kgover=new KML_GroundOverlay();
		sprintf(kgover->name      ,"%s",pimages[i]);
		kgover->visibility=0;

		kicon=new KML_Icon();
		sprintf(kicon->href      ,"%s",pimages[i]);
		kgover->icon      =kicon;
		kicon=NULL;

		kllbox=new KML_LatLonBox();
		kllbox->north     =lataxis[1];
		kllbox->south     =lataxis[0];
		kllbox->east      =longaxis[1];
		kllbox->west      =longaxis[0];
		kllbox->rotation  = 0.;
		kgover->llbox     =kllbox;
		kllbox=NULL;

		(kfold->feature   )->AddObject((Object*)kgover);
		kgover=NULL;
	}

/*  assemble the rest of the kml hierarchy  */

	(kdoc->feature   )->AddObject((Object*)kfold);
	kfold=NULL;
	(kfile->kmlobj    )->AddObject((Object*)kdoc);
	kdoc=NULL;

/*  write kml file  */

	_printf0_("Writing kml document to file.\n");
	fprintf(fid,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
	kfile->Write(fid,indent);

	delete kfile;

	clock1=clock();
	time1 =time(NULL);
	_printf_("KMLOverlayx Module -- " << ((double)(clock1-clock0))/CLOCKS_PER_SEC << " CPU seconds; " << difftime(time1,time0) << " elapsed seconds.\n\n\n");

	return;
}
