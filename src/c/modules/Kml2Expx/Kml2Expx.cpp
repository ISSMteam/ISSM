/*!\file Kml2Expx
 * \brief kml to exp conversion routines.
 */

#include "./Kml2Expx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../KMLFileReadx/KMLFileReadx.h"

int Kml2Expx(char* filkml,char* filexp,int sgn){

	double  cm,sp;
	Ll2xydef(&cm,&sp,sgn);

	return(Kml2Expx(filkml,filexp,sgn,cm,sp));
}

int Kml2Expx(char* filkml,char* filexp,int sgn,double cm,double sp){

	int         iret   = 0;
	KML_Object *kobj   = NULL;
	FILE       *fidi   = NULL;
	FILE       *fido   = NULL;
	clock_t     clock0,clock1;
	time_t      time0 ,time1;

	clock0=clock();
	time0 =time(NULL);
	_printf0_("\nKml2Expx Module -- " << ctime(&time0));

	/*read kml file*/
	fidi=fopen(filkml,"r");
	if (!(kobj=KMLFileReadx(fidi)))
	 _error_("Error reading kml file.");
	fclose(fidi);

	/*open exp file*/
	_printf0_("Writing exp profiles to file.\n");
	fido=fopen(filexp,"w");

	/*write the polygons and linestrings  */
	kobj->WriteExp(fido,"",sgn,cm,sp);

	/*close exp file  */
	fclose(fido);
	delete kobj;

	clock1=clock();
	time1 =time(NULL);
	_printf_("Kml2Expx Module -- " << ((double)(clock1-clock0))/CLOCKS_PER_SEC << " CPU seconds; " << difftime(time1,time0) << " elapsed seconds.\n\n\n");

	return(iret);
}
