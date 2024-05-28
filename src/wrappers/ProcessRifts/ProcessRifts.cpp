/*!\file:  ProcessRifts.cpp
 * \brief split a mesh where a rift (or fault) is present
 */ 

#include "./ProcessRifts.h"

void ProcessRiftsUsage(void){/*{{{*/
	_printf_("\n");
	_printf_("   usage: [index2,x2,y2,segments2,segmentmarkers2,rifts2]=ProcessRifts(index1,x1,y1,segments1,segmentmarkers1) \n");
	_printf_("      where: (index1,x1,y1,segments1,segmentmarkers1) is an initial triangulation.\n");
	_printf_("      index2,x2,y2,segments2,segmentmarkers2,rifts2 is the resulting triangulation where rifts have been processed.\n");
}/*}}}*/
WRAPPER(ProcessRifts_python){

	/* returned quantities: */
	RiftStruct *riftstruct = NULL;

	/* input: */
	int     nel,nods;
	int    *index          = NULL;
	double *x              = NULL;
	double *y              = NULL;
	int    *segments       = NULL;
	int    *segmentmarkers = NULL;
	int     num_seg;

	/*Boot module*/
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	CHECKARGUMENTS(NLHS,NRHS,&ProcessRiftsUsage);

	/*Fetch data */
	FetchData(&index,&nel,NULL,INDEXIN);
	FetchData(&x,&nods,XIN);
	FetchData(&y,NULL,YIN);
	FetchData(&segments,&num_seg,NULL,SEGMENTSIN);
	FetchData(&segmentmarkers,NULL,SEGMENTMARKERSIN);

	/*call x layer*/
	ProcessRiftsx(&index,&nel,&x,&y,&nods,&segments,&segmentmarkers,&num_seg,&riftstruct);

	/*Output : */
	WriteData(INDEXOUT,index,nel,3);
	WriteData(XOUT,x,nods,1);
	WriteData(YOUT,y,nods,1);
	WriteData(SEGMENTSOUT,segments,num_seg,3);
	WriteData(SEGMENTMARKERSOUT,segmentmarkers,num_seg,1);
	WriteData(RIFTSTRUCT,riftstruct);

	/*end module: */
	delete riftstruct;
	xDelete<int>(index);
	xDelete<double>(x);
	xDelete<double>(y);
	xDelete<int>(segments);
	xDelete<int>(segmentmarkers );
	MODULEEND();
}
